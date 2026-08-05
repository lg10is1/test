options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
})

analysis_date <- "2026-07-16"
args_cli <- commandArgs(trailingOnly = TRUE)
get_cli_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args_cli[startsWith(args_cli, prefix)]
  if (length(hit) > 0) return(sub(prefix, "", hit[1], fixed = TRUE))
  default
}
project_root <- normalizePath(
  get_cli_arg("project-root", Sys.getenv("EOSCZ_PROJECT_ROOT", getwd())),
  winslash = "/",
  mustWork = FALSE
)
cnv_root <- get_cli_arg("cnv-root", file.path(project_root, "cnv_analysis/protein_coding_genes"))
burden_root <- get_cli_arg("burden-root", file.path(project_root, "cnv_analysis/pathway_burden"))
out_dir <- file.path(burden_root, "CNV_PCA_3cohorts_case_cohort_26-7-16")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PUBLIC_BINARY_SAMPLE_CN <- "sample_CN_binary.xlsx"
LEGACY_BINARY_SAMPLE_CN <- "sample_copy_number_presence.xlsx"

sample_binary_path <- function(cohort_dir) {
  public <- file.path(cohort_dir, PUBLIC_BINARY_SAMPLE_CN)
  if (file.exists(public)) return(public)
  legacy <- file.path(cohort_dir, LEGACY_BINARY_SAMPLE_CN)
  if (file.exists(legacy)) return(legacy)
  public
}

inputs <- list(
  sample = list(
    `case_cohort` = sample_binary_path(file.path(cnv_root, "slurm_scripts_case_cohort")),
    comparison_cohort = sample_binary_path(file.path(cnv_root, "slurm_scripts_comparison_site")),
    public_reference = sample_binary_path(file.path(cnv_root, "slurm_scripts_public_reference"))
  ),
  haplotype = list(
    `case_cohort` = file.path(cnv_root, "slurm_scripts_case_cohort/haplotype_CN.xlsx"),
    comparison_cohort = file.path(cnv_root, "slurm_scripts_comparison_site/haplotype_CN.xlsx"),
    public_reference = file.path(cnv_root, "slurm_scripts_public_reference/haplotype_CN.xlsx")
  )
)

clean_gene <- function(x) toupper(trimws(as.character(x)))
empty_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == "" | toupper(trimws(as.character(x))) %in% c("NA", "N/A", "UNKNOWN")
}

read_binary_matrix <- function(path) {
  if (!file.exists(path)) stop("Input not found: ", path)
  sheet <- readxl::excel_sheets(path)[1]
  raw <- readxl::read_excel(path, sheet = sheet, .name_repair = "minimal")
  ids <- trimws(as.character(raw[[1]]))
  mat <- as.data.frame(raw[, -1, drop = FALSE], check.names = FALSE)
  names(mat) <- clean_gene(names(mat))
  mat[] <- lapply(mat, function(x) {
    y <- suppressWarnings(as.numeric(x))
    y[is.na(y)] <- 0
    as.integer(y > 0)
  })
  keep <- !empty_value(ids)
  mat <- mat[keep, , drop = FALSE]
  ids <- ids[keep]
  if (anyDuplicated(names(mat))) {
    split_cols <- split(seq_along(names(mat)), names(mat))
    mat <- as.data.frame(
      lapply(split_cols, function(idx) as.integer(rowSums(mat[, idx, drop = FALSE]) > 0)),
      check.names = FALSE
    )
  }
  rownames(mat) <- ids
  if (anyDuplicated(rownames(mat))) {
    current_genes <- colnames(mat)
    split_rows <- split(seq_len(nrow(mat)), rownames(mat))
    mat <- as.data.frame(
      t(vapply(split_rows, function(idx) as.integer(colSums(mat[idx, , drop = FALSE]) > 0), integer(ncol(mat)))),
      check.names = FALSE
    )
    colnames(mat) <- current_genes
  }
  mat
}

align_matrices <- function(mats) {
  all_genes <- sort(unique(unlist(lapply(mats, colnames))))
  lapply(mats, function(mat) {
    missing <- setdiff(all_genes, colnames(mat))
    if (length(missing) > 0) {
      zeros <- as.data.frame(matrix(0L, nrow = nrow(mat), ncol = length(missing)), check.names = FALSE)
      colnames(zeros) <- missing
      mat <- cbind(mat, zeros)
    }
    mat[, all_genes, drop = FALSE]
  })
}

run_pca <- function(level) {
  mats <- lapply(inputs[[level]], read_binary_matrix)
  mats <- align_matrices(mats)
  cohort <- rep(names(mats), vapply(mats, nrow, integer(1)))
  record_id <- unlist(lapply(mats, rownames), use.names = FALSE)
  x <- as.matrix(do.call(rbind, unname(mats)))
  storage.mode(x) <- "numeric"
  vars <- apply(x, 2, var)
  keep <- is.finite(vars) & vars > 0
  x_keep <- x[, keep, drop = FALSE]
  pca <- prcomp(x_keep, center = TRUE, scale. = FALSE)
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  scores <- data.frame(
    level = level,
    record_id = record_id,
    cohort = cohort,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    total_CNV_gene_count = rowSums(x),
    stringsAsFactors = FALSE
  )
  meta <- data.frame(
    level = level,
    n_records = nrow(x),
    n_cohorts = length(mats),
    n_genes_union = ncol(x),
    n_genes_used_nonzero_variance = ncol(x_keep),
    PC1_variance = var_exp[1],
    PC2_variance = var_exp[2],
    stringsAsFactors = FALSE
  )
  list(scores = scores, meta = meta)
}

sample_res <- run_pca("sample")
hap_res <- run_pca("haplotype")
all_scores <- rbind(sample_res$scores, hap_res$scores)
all_meta <- rbind(sample_res$meta, hap_res$meta)

write.table(all_scores, file.path(out_dir, "CNV_PCA_scores_sample_haplotype_3cohorts_26-7-16.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(all_meta, file.path(out_dir, "CNV_PCA_meta_sample_haplotype_3cohorts_26-7-16.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

plot_one <- function(scores, meta, main_title) {
  cols <- c(`case_cohort` = "#D55E00", comparison_cohort = "#0072B2", public_reference = "#009E73")
  pch_map <- c(`case_cohort` = 16, comparison_cohort = 17, public_reference = 15)
  xlab <- sprintf("PC1 (%.2f%%)", 100 * meta$PC1_variance[1])
  ylab <- sprintf("PC2 (%.2f%%)", 100 * meta$PC2_variance[1])
  plot(scores$PC1, scores$PC2,
       col = cols[scores$cohort], pch = pch_map[scores$cohort],
       xlab = xlab, ylab = ylab, main = main_title,
       cex = 0.75, las = 1)
  grid(col = "#DDDDDD")
  legend("topright", legend = names(cols), col = cols, pch = pch_map, bty = "n", cex = 0.8)
  mtext(sprintf("n=%s; genes used=%s", meta$n_records[1], meta$n_genes_used_nonzero_variance[1]), side = 3, line = 0.2, cex = 0.75, col = "#444444")
}

pdf_file <- file.path(out_dir, "CNV_PCA_3cohorts_sample_and_haplotype_26-7-16.pdf")
png_file <- file.path(out_dir, "CNV_PCA_3cohorts_sample_and_haplotype_26-7-16.png")

pdf(pdf_file, width = 12, height = 5.8, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(4.3, 4.3, 3.4, 1.2), oma = c(0, 0, 2, 0))
plot_one(sample_res$scores, sample_res$meta, "Sample-level CNV PCA")
plot_one(hap_res$scores, hap_res$meta, "Haplotype-level CNV PCA")
mtext("CNV PCA across case_cohort, comparison_cohort and public_reference", outer = TRUE, cex = 1.1, font = 2)
dev.off()

png(png_file, width = 3600, height = 1740, res = 300)
par(mfrow = c(1, 2), mar = c(4.3, 4.3, 3.4, 1.2), oma = c(0, 0, 2, 0))
plot_one(sample_res$scores, sample_res$meta, "Sample-level CNV PCA")
plot_one(hap_res$scores, hap_res$meta, "Haplotype-level CNV PCA")
mtext("CNV PCA across case_cohort, comparison_cohort and public_reference", outer = TRUE, cex = 1.1, font = 2)
dev.off()

cat("DONE\n")
cat("Output directory:", out_dir, "\n")
print(all_meta)


