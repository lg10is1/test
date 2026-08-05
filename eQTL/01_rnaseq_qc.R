#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Usage: Rscript 01_rnaseq_qc.R COUNT.csv SAMPLE_NAMES.txt",
    "GENE_COORDINATES.tsv BASE_COVARIATES.tsv PCA.sscore OUTPUT_DIR"
  ))
}

suppressPackageStartupMessages({
  library(data.table)
  library(peer)
})

count_file <- args[1]
sample_file <- args[2]
coordinate_file <- args[3]
base_covariate_file <- args[4]
pca_file <- args[5]
output_dir <- args[6]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- fread(sample_file, header = FALSE)[[1]]
counts <- fread(count_file)
setnames(counts, 2:ncol(counts), sample_ids)

gene_parts <- strsplit(as.character(counts[[1]]), "\\|")
gene_ids <- vapply(
  gene_parts,
  function(value) if (length(value) >= 2) value[2] else NA_character_,
  character(1)
)

valid_gene <- !is.na(gene_ids)
count_matrix <- as.matrix(counts[valid_gene, 2:ncol(counts), with = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- gene_ids[valid_gene]

keep_gene <- rowSums(count_matrix > 6) > ncol(count_matrix) * 0.8
count_matrix <- count_matrix[keep_gene, , drop = FALSE]
sample_by_gene <- t(count_matrix)

count_qc_file <- file.path(output_dir, "rna_count_qc.csv")
count_qc <- as.data.table(sample_by_gene, keep.rownames = "sample_id")
fwrite(count_qc, count_qc_file)

tmm_file <- file.path(output_dir, "rna_tmm.csv")
status <- system2("rnanorm", c("tmm", count_qc_file), stdout = tmm_file)
if (status != 0) stop("rnanorm tmm failed")

tmm <- fread(tmm_file)
tmm_samples <- tmm[[1]]
tmm_matrix <- as.matrix(tmm[, 2:ncol(tmm), with = FALSE])
storage.mode(tmm_matrix) <- "numeric"
rownames(tmm_matrix) <- tmm_samples

inverse_normal <- function(value) {
  qnorm((rank(value, na.last = "keep") - 0.5) / sum(!is.na(value)))
}

inverse_sample_by_gene <- apply(tmm_matrix, 2, inverse_normal)
inverse_sample_by_gene <- as.matrix(inverse_sample_by_gene)
rownames(inverse_sample_by_gene) <- tmm_samples
inverse_gene_by_sample <- t(inverse_sample_by_gene)

inverse_file <- file.path(output_dir, "rna_invt.tsv")
inverse_table <- as.data.table(inverse_gene_by_sample, keep.rownames = "gene")
fwrite(inverse_table, inverse_file, sep = "\t")

coordinates <- fread(coordinate_file)
phenotypes <- merge(coordinates, inverse_table, by = "gene")
chromosome <- ifelse(
  startsWith(as.character(phenotypes$chr), "chr"),
  as.character(phenotypes$chr),
  paste0("chr", phenotypes$chr)
)

bed <- data.table(
  "#Chr" = chromosome,
  start = phenotypes$start,
  end = phenotypes$end,
  pid = phenotypes$gene,
  gid = phenotypes$gene,
  strand = phenotypes$strand
)
bed <- cbind(bed, phenotypes[, ..tmm_samples])
setorderv(bed, c("#Chr", "start", "end"))

bed_file <- file.path(output_dir, "rna_tmm_inv.bed")
fwrite(bed, bed_file, sep = "\t")
system2("bgzip", c("-f", bed_file))
system2("tabix", c("-f", "-p", "bed", paste0(bed_file, ".gz")))

peer_model <- PEER()
PEER_setPhenoMean(peer_model, inverse_sample_by_gene)
PEER_setNk(peer_model, 15)
PEER_setNmax_iterations(peer_model, 1000)
PEER_setAdd_mean(peer_model, TRUE)
PEER_update(peer_model)

peer_factors <- as.data.table(PEER_getX(peer_model))
setnames(peer_factors, paste0("PEER", seq_len(ncol(peer_factors))))
peer_factors[, IID := tmm_samples]
setcolorder(peer_factors, c("IID", setdiff(names(peer_factors), "IID")))

base_covariates <- fread(base_covariate_file)
if (!"IID" %in% names(base_covariates)) setnames(base_covariates, 1, "IID")

pca <- fread(pca_file)
iid_column <- grep("^#?IID$", names(pca), ignore.case = TRUE, value = TRUE)[1]
setnames(pca, iid_column, "IID")
pc_columns <- grep("^PC[0-9]+", names(pca), value = TRUE)
pca <- pca[, c("IID", pc_columns), with = FALSE]

covariates <- merge(peer_factors, base_covariates, by = "IID")
covariates <- merge(covariates, pca, by = "IID")
covariates <- covariates[match(tmm_samples, IID)]

covariate_columns <- setdiff(names(covariates), "IID")
covariate_matrix <- t(as.matrix(covariates[, ..covariate_columns]))
covariate_output <- as.data.table(covariate_matrix, keep.rownames = "id")
setnames(covariate_output, c("id", as.character(covariates$IID)))
fwrite(covariate_output, file.path(output_dir, "covariates_qtltools.tsv"), sep = "\t")
