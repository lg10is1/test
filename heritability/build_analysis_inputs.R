#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || identical(x, "")) y else x
parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Invalid argument near: ", x[[i]])
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("source", "pheno", "batch", "snv-pc10", "snv-pc20", "sv-pc10", "sv-pc20",
  "tr-pc10", "tr-pc20", "snv-grm-id", "sv-grm-id", "tr-grm-id", "out-dir")
missing_args <- required[!required %in% names(args)]
if (length(missing_args)) stop("Missing argument(s): ", paste(paste0("--", missing_args), collapse = ", "))
source_name <- args[["source"]]; out_dir <- args[["out-dir"]]
for (nm in setdiff(required, c("source", "out-dir"))) if (!file.exists(args[[nm]])) stop("Missing file: ", args[[nm]])

read_ids <- function(file, label, min_cols = 2L) {
  x <- fread(file, header = FALSE)
  if (ncol(x) < min_cols) stop(label, " has fewer than ", min_cols, " columns: ", file)
  x <- x[, 1:2]; setnames(x, c("FID", "IID"))
  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  if (anyDuplicated(x)) stop("Duplicated FID/IID in ", label, ": ", file)
  unique(x)
}

read_ids_auto <- function(file, label) {
  x <- fread(file, header = "auto")
  ids_are_headers <- ncol(x) >= 2L && all(toupper(names(x)[1:2]) %in% c("FID", "IID", "#FID", "#IID"))
  if (!ids_are_headers) x <- fread(file, header = FALSE)
  if (ncol(x) < 2L) stop(label, " has fewer than two columns: ", file)
  x <- x[, 1:2]; setnames(x, c("FID", "IID"))
  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  if (anyDuplicated(x)) stop("Duplicated FID/IID in ", label, ": ", file)
  unique(x)
}

read_pc <- function(file, prefix, expected_n) {
  x <- fread(file)
  if (ncol(x) != expected_n + 2L) stop("Expected ", expected_n, " PCs in: ", file)
  setnames(x, 1:2, c("FID", "IID"))
  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  if (anyDuplicated(x[, .(FID, IID)])) stop("Duplicated FID/IID in PC file: ", file)
  pc_cols <- names(x)[-(1:2)]
  setnames(x, pc_cols, paste0(prefix, "_PC", seq_along(pc_cols)))
  for (nm in names(x)[-(1:2)]) set(x, j = nm, value = as.numeric(x[[nm]]))
  if (any(!is.finite(as.matrix(x[, -(1:2)])))) stop("Non-finite PC value in: ", file)
  x
}

id_inputs <- list(
  phenotype = read_ids_auto(args[["pheno"]], "phenotype"),
  batch = read_ids_auto(args[["batch"]], "batch"),
  snv_grm = read_ids(args[["snv-grm-id"]], "SNV GRM"),
  sv_grm = read_ids(args[["sv-grm-id"]], "SV GRM"),
  tr_grm = read_ids(args[["tr-grm-id"]], "TR GRM")
)

pc <- list(
  pc10 = list(
    SNV = read_pc(args[["snv-pc10"]], "SNV", 10L),
    SV = read_pc(args[["sv-pc10"]], "SV", 10L),
    TR = read_pc(args[["tr-pc10"]], "TR", 10L)
  ),
  pc20 = list(
    SNV = read_pc(args[["snv-pc20"]], "SNV", 20L),
    SV = read_pc(args[["sv-pc20"]], "SV", 20L),
    TR = read_pc(args[["tr-pc20"]], "TR", 20L)
  )
)

all_id_tables <- c(id_inputs, lapply(pc$pc20, function(x) x[, .(FID, IID)]))
common <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID"), all = FALSE), all_id_tables)
n_common_before_removal <- nrow(common)
n_removed <- 0L
remove_file <- args[["remove-file"]] %||% ""
if (nzchar(remove_file)) {
  if (!file.exists(remove_file)) stop("Missing remove file: ", remove_file)
  remove <- fread(remove_file, header = FALSE, fill = TRUE)
  if (ncol(remove) < 1L) stop("Empty remove file: ", remove_file)
  if (ncol(remove) >= 2L) {
    remove <- remove[, 1:2]; setnames(remove, c("FID", "IID"))
    remove[, `:=`(FID = as.character(FID), IID = as.character(IID))]
    before <- nrow(common)
    common <- common[!remove, on = .(FID, IID)]
    n_removed <- before - nrow(common)
  } else {
    remove_iid <- unique(as.character(remove[[1]]))
    before <- nrow(common)
    common <- common[!IID %chin% remove_iid]
    n_removed <- before - nrow(common)
  }
}
setorder(common, FID, IID)
if (nrow(common) < 50L) stop("Unified analysis intersection is unexpectedly small (N=", nrow(common), ")")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
keep_file <- file.path(out_dir, paste0(source_name, ".analysis.keep"))
fwrite(common, keep_file, sep = "\t", col.names = FALSE)

audit <- rbindlist(c(
  lapply(names(id_inputs), function(nm) data.table(source = source_name, input = nm, n_unique = nrow(id_inputs[[nm]]))),
  lapply(names(pc$pc20), function(nm) data.table(source = source_name, input = paste0(tolower(nm), "_pc"), n_unique = nrow(pc$pc20[[nm]])))
))
audit[, `:=`(n_common = nrow(common), n_excluded = n_unique - nrow(common))]
fwrite(audit, file.path(out_dir, paste0(source_name, ".sample_intersection.tsv")), sep = "\t")
fwrite(audit, file.path(out_dir, paste0(source_name, ".sample_intersection.csv")))
removal_summary <- data.table(source = source_name, remove_file = remove_file,
  n_common_before_removal = n_common_before_removal, n_removed = n_removed,
  n_common_after_removal = nrow(common))
fwrite(removal_summary, file.path(out_dir, paste0(source_name, ".sample_removal_summary.tsv")), sep = "\t")
fwrite(removal_summary, file.path(out_dir, paste0(source_name, ".sample_removal_summary.csv")))

strategy_components <- list(
  snv_pc = "SNV", sv_pc = "SV", tr_pc = "TR",
  snv_sv_pc = c("SNV", "SV"), snv_tr_pc = c("SNV", "TR"),
  sv_tr_pc = c("SV", "TR"), snv_sv_tr_pc = c("SNV", "SV", "TR")
)

manifest <- list()
for (dim_name in names(pc)) {
  n_pc <- as.integer(sub("pc", "", dim_name))
  restricted <- lapply(pc[[dim_name]], function(x) merge(common, x, by = c("FID", "IID"), all.x = TRUE, sort = FALSE))
  for (strategy in names(strategy_components)) {
    comps <- strategy_components[[strategy]]
    q <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID"), all = FALSE, sort = FALSE), restricted[comps])
    q <- merge(common, q, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
    numeric_cols <- setdiff(names(q), c("FID", "IID"))
    mat <- as.matrix(q[, ..numeric_cols])
    finite <- all(is.finite(mat))
    rank <- if (finite) qr(scale(mat, center = TRUE, scale = FALSE))$rank else NA_integer_
    zero_var <- numeric_cols[vapply(q[, ..numeric_cols], function(z) is.na(var(z)) || var(z) == 0, logical(1))]
    expected_rank <- length(numeric_cols)
    valid <- finite && !length(zero_var) && rank == expected_rank
    path <- file.path(out_dir, sprintf("%s.pc%d.%s.qcovar", source_name, n_pc, strategy))
    fwrite(q, path, sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")
    manifest[[length(manifest) + 1L]] <- data.table(
      source = source_name, pc_n = n_pc, strategy = strategy, qcovar_file = path,
      n_samples = nrow(q), n_covariates = expected_rank, matrix_rank = rank,
      zero_variance_columns = paste(zero_var, collapse = ","), valid = valid,
      validation_message = if (valid) "OK" else "non-finite, zero-variance, or rank-deficient qcovar"
    )
  }
}
manifest_dt <- rbindlist(manifest)
fwrite(manifest_dt, file.path(out_dir, paste0(source_name, ".qcovar_manifest.tsv")), sep = "\t")
fwrite(manifest_dt, file.path(out_dir, paste0(source_name, ".qcovar_manifest.csv")))
message("[DONE] Analysis inputs: ", source_name, " | common N=", nrow(common), " | qcovars=", nrow(manifest_dt))
