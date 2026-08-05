#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || identical(x, "")) y else x
parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) { out[[sub("^--", "", key)]] <- TRUE; i <- i + 1L }
    else { out[[sub("^--", "", key)]] <- x[[i + 1L]]; i <- i + 2L }
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))

source_name <- args[["source"]] %||% stop("Missing --source")
component <- args[["component"]] %||% stop("Missing --component")
method <- args[["method"]] %||% stop("Missing --method")
scores_file <- args[["scores"]] %||% stop("Missing --scores")
eigenval_file <- args[["eigenval"]] %||% stop("Missing --eigenval")
unrelated_file <- args[["unrelated"]] %||% ""
pheno_file <- args[["pheno"]] %||% stop("Missing --pheno")
batch_file <- args[["batch"]] %||% stop("Missing --batch")
out_prefix <- args[["out-prefix"]] %||% stop("Missing --out-prefix")
if (!method %in% c("unrelated_projected", "full_sample_projected_fallback", "full_grm_fallback")) stop("Invalid PCA method: ", method)
for (f in c(scores_file, eigenval_file, pheno_file, batch_file)) if (!file.exists(f)) stop("Missing file: ", f)

read_scores <- function(file, method) {
  if (method == "full_grm_fallback") {
    x <- fread(file, header = FALSE)
    if (ncol(x) < 22L) stop("GCTA eigenvec must contain FID, IID and at least 20 PCs: ", file)
    setnames(x, 1:2, c("FID", "IID"))
    setnames(x, 3:ncol(x), paste0("PC", seq_len(ncol(x) - 2L)))
  } else {
    x <- fread(file, header = TRUE, check.names = FALSE)
    id1 <- intersect(c("#FID", "FID"), names(x))[1]
    id2 <- intersect(c("IID", "#IID"), names(x))[1]
    if (is.na(id1) || is.na(id2)) stop("PLINK2 .sscore lacks FID/IID: ", file)
    setnames(x, c(id1, id2), c("FID", "IID"))
    candidates <- grep("^PC[0-9]+(_AVG)?$", names(x), value = TRUE)
    candidates <- candidates[order(as.integer(sub("^PC([0-9]+).*$", "\\1", candidates)))]
    if (length(candidates) < 20L) stop("Projected score file has fewer than 20 PC columns: ", file)
    x <- x[, c("FID", "IID", candidates), with = FALSE]
    setnames(x, candidates, paste0("PC", seq_along(candidates)))
  }
  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  if (anyDuplicated(x[, .(FID, IID)])) stop("Duplicated FID/IID in scores: ", file)
  pc <- paste0("PC", 1:20)
  for (nm in pc) set(x, j = nm, value = as.numeric(x[[nm]]))
  if (any(!is.finite(as.matrix(x[, ..pc])))) stop("Non-finite projected PC value in: ", file)
  x[, c("FID", "IID", pc), with = FALSE]
}

read_three_col <- function(file, value_name) {
  x <- fread(file, header = "auto")
  id_names <- toupper(names(x)[1:2])
  if (length(id_names) < 2L || !all(id_names %in% c("FID", "IID", "#FID", "#IID"))) x <- fread(file, header = FALSE)
  if (ncol(x) < 3L) stop("Expected at least three columns: ", file)
  x <- x[, 1:3]
  setnames(x, c("FID", "IID", value_name))
  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  if (anyDuplicated(x[, .(FID, IID)])) stop("Duplicated IDs in: ", file)
  x
}

scores <- read_scores(scores_file, method)
pheno <- read_three_col(pheno_file, "phenotype")
batch <- read_three_col(batch_file, "batch")

if (nzchar(unrelated_file)) {
  if (!file.exists(unrelated_file)) stop("Missing unrelated keep: ", unrelated_file)
  unrelated <- fread(unrelated_file, header = FALSE)
  if (ncol(unrelated) < 2L) stop("Unrelated keep has fewer than two columns: ", unrelated_file)
  unrelated <- unrelated[, 1:2]; setnames(unrelated, c("FID", "IID"))
  unrelated[, `:=`(FID = as.character(FID), IID = as.character(IID), pca_reference = "unrelated_reference")]
  scores <- merge(scores, unrelated, by = c("FID", "IID"), all.x = TRUE)
  scores[is.na(pca_reference), pca_reference := "projected_other"]
} else {
  scores[, pca_reference := if (method == "full_grm_fallback") "full_sample_grm_pca" else "full_sample_genotype_pca"]
}
scores <- merge(scores, pheno, by = c("FID", "IID"), all.x = TRUE)
scores <- merge(scores, batch, by = c("FID", "IID"), all.x = TRUE)
scores[, `:=`(source = source_name, component = component, pca_method = method)]

eigen <- fread(eigenval_file, header = FALSE)
if (!nrow(eigen)) stop("Empty eigenvalue file: ", eigenval_file)
vals <- suppressWarnings(as.numeric(eigen[[1]]))
if (any(!is.finite(vals))) stop("Non-numeric eigenvalue in: ", eigenval_file)
eigen_dt <- data.table(
  source = source_name, component = component, pca_method = method,
  PC = paste0("PC", seq_along(vals)), eigenvalue = vals,
  percent = if (sum(vals) > 0) 100 * vals / sum(vals) else NA_real_
)

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
pc10 <- scores[, c("FID", "IID", paste0("PC", 1:10)), with = FALSE]
pc20 <- scores[, c("FID", "IID", paste0("PC", 1:20)), with = FALSE]
fwrite(pc10, paste0(out_prefix, ".pc10.qcovar.tsv"), sep = "\t")
fwrite(pc20, paste0(out_prefix, ".pc20.qcovar.tsv"), sep = "\t")
fwrite(scores, paste0(out_prefix, ".plot.tsv"), sep = "\t", na = "NA")
fwrite(eigen_dt, paste0(out_prefix, ".eigenval.tsv"), sep = "\t", na = "NA")

meta <- data.table(
  source = source_name, component = component, pca_method = method,
  n_scores = nrow(scores), n_unrelated_reference = sum(scores$pca_reference == "unrelated_reference"),
  score_file = scores_file, eigenval_file = eigenval_file
)
fwrite(meta, paste0(out_prefix, ".pca_summary.tsv"), sep = "\t")
message("[DONE] ", source_name, " ", component, " PCA outputs | method=", method, " | N=", nrow(scores))
