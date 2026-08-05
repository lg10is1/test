#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (is.character(x) && length(x) == 1L && (is.na(x) || x == ""))) y else x
}

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[i], "--")) stop("Unexpected argument: ", x[i])
    key <- sub("^--", "", x[i])
    if (i == length(x) || startsWith(x[i + 1L], "--")) {
      out[[key]] <- TRUE; i <- i + 1L
    } else {
      out[[key]] <- x[i + 1L]; i <- i + 2L
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
base <- "/path/to/EOSCZ_PROJECT"
matches_file <- args[["matches"]] %||% file.path(
  base, "figure_analysis", "02.meQTL", "public", "SV", "matching",
  "sig_sv_to_pav_sv_len50.best.tsv"
)
pav_vcf <- args[["pav-sv-vcf"]] %||% file.path(
  base, "TGS_SV_merge_SCZ", "truvari_single_sample",
  "truvari_merged_sort_pP0.5.sv_len_gt50.sorted.vcf.gz"
)
sample_info_file <- args[["sample-info"]] %||% file.path(base, "Sample_info", "sample_case.txt")
genomic_pc_file <- args[["genomic-pc"]] %||% file.path(
  base, "TGS_callset", "glnexus_merge", "autosome_pca", "unrelated_projected_pca",
  "all_samples.projected_pca.PC1_PC2.projected_scores.tsv"
)
epigenetic_pc_file <- args[["epigenetic-pc"]] %||%
  "/path/to/SCZ/methy_pbcpg/pca_600k/methylation_uniform_600k_pca.scores.tsv"
bcftools <- path.expand(args[["bcftools"]] %||% Sys.which("bcftools"))
out_prefix <- path.expand(args[["out-prefix"]] %||% file.path(
  base, "figure_analysis", "02.meQTL", "public", "SV", "tables",
  "sv_meqtl.genotype_debug"
))
min_n <- as.integer(args[["min-n"]] %||% 50L)
match_types <- trimws(strsplit(
  args[["match-types"]] %||% "exact_ref_alt,fuzzy_ref_same", ",", fixed = TRUE
)[[1L]])

check_file <- function(x, label) {
  if (!file.exists(path.expand(x)) || dir.exists(path.expand(x)) || file.info(path.expand(x))$size == 0) {
    stop("Missing or empty ", label, ": ", x)
  }
}
for (z in list(
  c(matches_file, "matches"), c(pav_vcf, "TGS SV VCF"),
  c(sample_info_file, "sample info"), c(genomic_pc_file, "genomic PCs"),
  c(epigenetic_pc_file, "epigenetic PCs")
)) check_file(z[1L], z[2L])
if (bcftools == "" || !file.exists(bcftools)) stop("bcftools not found; pass --bcftools")
if (is.na(min_n) || min_n < 1L) stop("--min-n must be positive")

vcf_has_index <- function(vcf) any(file.exists(paste0(path.expand(vcf), c(".tbi", ".csi"))))
view_args <- function(vcf, region) {
  if (vcf_has_index(vcf)) c("view", "-H", "-r", region, path.expand(vcf))
  else c("view", "-H", "-t", region, path.expand(vcf))
}

.sample_cache <- new.env(parent = emptyenv())
get_samples <- function(vcf) {
  key <- path.expand(vcf)
  if (!exists(key, envir = .sample_cache, inherits = FALSE)) {
    x <- system2(bcftools, c("query", "-l", key), stdout = TRUE)
    if (!length(x)) stop("No VCF samples: ", key)
    assign(key, x, envir = .sample_cache)
  }
  get(key, envir = .sample_cache, inherits = FALSE)
}

read_record <- function(vcf, chr, pos, id, samples) {
  region <- paste0(chr, ":", pos, "-", pos)
  lines <- system2(bcftools, view_args(vcf, region), stdout = TRUE, stderr = FALSE)
  if (!length(lines)) stop("No VCF record at ", region)
  x <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t", data.table = FALSE)
  hit <- x[as.character(x[[3L]]) == as.character(id), , drop = FALSE]
  if (nrow(hit) != 1L) stop("Expected one record ID=", id, "; observed ", nrow(hit))
  if (ncol(hit) - 9L != length(samples)) stop("VCF sample/header column count mismatch")
  fmt <- strsplit(as.character(hit[[9L]]), ":", fixed = TRUE)[[1L]]
  gt_index <- match("GT", fmt)
  if (is.na(gt_index)) stop("FORMAT lacks GT for ", id)
  fields <- as.character(hit[1L, 10:ncol(hit)])
  gt <- vapply(
    strsplit(fields, ":", fixed = TRUE),
    function(z) if (length(z) >= gt_index) z[[gt_index]] else NA_character_,
    character(1L)
  )
  names(gt) <- samples
  list(ref = as.character(hit[[4L]]), alt = as.character(hit[[5L]]), gt = gt)
}

parse_dosage <- function(gt, target, impute_missing_ref = TRUE) {
  gt <- sub(":.*$", "", as.character(gt))
  out <- rep(NA_real_, length(gt))
  ok <- !is.na(gt) & gt != ""
  alleles <- strsplit(gt[ok], "[/|]")
  out[ok] <- vapply(alleles, function(z) {
    if (any(z == ".") && !impute_missing_ref) return(NA_real_)
    if (impute_missing_ref) z[z == "."] <- "0"
    if (identical(as.character(target), "ANY_ALT")) return(sum(z != "0"))
    sum(z == as.character(target))
  }, numeric(1L))
  out
}

count_missing_alleles <- function(gt) {
  gt <- sub(":.*$", "", as.character(gt))
  out <- rep(NA_integer_, length(gt))
  ok <- !is.na(gt) & gt != ""
  out[ok] <- vapply(strsplit(gt[ok], "[/|]"), function(z) sum(z == "."), integer(1L))
  out
}

impute_gt_display <- function(gt) {
  gt <- as.character(gt)
  ok <- !is.na(gt) & gt != ""
  gt[ok] <- gsub(".", "0", gt[ok], fixed = TRUE)
  gt
}

is_true <- function(x) length(x) == 1L && !is.na(x) && tolower(as.character(x)) %in% c("true", "t", "1")
choose_target <- function(row, ref, alt) {
  if (identical(as.character(row$best_match_type), "exact_ref_alt")) {
    target <- toupper(as.character(row$lead_a2))
    ref2 <- toupper(ref); alts <- strsplit(toupper(alt), ",", fixed = TRUE)[[1L]]
    if (identical(target, ref2)) return(list(index = 0L, rule = "exact_lead_a2_is_ref"))
    idx <- match(target, alts)
    if (is.na(idx)) stop("Exact match target allele absent from VCF REF/ALT")
    return(list(index = as.integer(idx), rule = "exact_lead_a2_is_alt"))
  }
  if (is_true(row$best_ref_same_as_a2)) {
    return(list(index = 0L, rule = "fuzzy_shared_ref_is_lead_a2"))
  }
  if (is_true(row$best_ref_same_as_a1)) {
    return(list(index = "ANY_ALT", rule = "fuzzy_shared_ref_is_lead_a1_any_alt_proxy_for_a2"))
  }
  stop("Fuzzy match lacks shared-reference orientation")
}

count_string <- function(x) {
  x <- as.character(x); x[is.na(x) | x == ""] <- "<missing>"
  tab <- table(x)
  paste(paste0(names(tab), "=", as.integer(tab)), collapse = ";")
}

# Recreate the covariate intersection used by the main SV meQTL script.
sample_info <- fread(sample_info_file, showProgress = FALSE)
if (!all(c("SampleID", "Sex", "Age") %in% names(sample_info))) stop("sample_info lacks SampleID/Sex/Age")
sample_info <- sample_info[, .(sample = as.character(SampleID))]
gpc <- fread(genomic_pc_file, select = "IID", showProgress = FALSE)
setnames(gpc, "IID", "sample"); gpc[, sample := as.character(sample)]
epc <- fread(epigenetic_pc_file, select = "sample", showProgress = FALSE)
epc[, sample := as.character(sample)]
cov_samples <- unique(merge(merge(sample_info, gpc, by = "sample"), epc, by = "sample")$sample)

matches <- fread(matches_file, showProgress = FALSE)
matches <- matches[best_match_type %in% match_types]
matches[, lead_chr_filter := sub("^chr", "", as.character(lead_chr), ignore.case = TRUE)]
matches <- matches[lead_chr_filter %in% as.character(1:22)]
matches[, lead_chr_filter := NULL]
if (!nrow(matches)) stop("No exact/fuzzy matches selected")

summary_rows <- list(); sample_rows <- list()
for (i in seq_len(nrow(matches))) {
  row <- matches[i]
  vcf <- pav_vcf
  if ("best_vcf_file" %in% names(matches) && !is.na(row$best_vcf_file) && row$best_vcf_file != "") {
    vcf <- as.character(row$best_vcf_file)
  }
  message("[", i, "/", nrow(matches), "] ", row$lead_id, " -> ", row$best_vcf_id)
  result <- tryCatch({
    samples <- get_samples(vcf)
    rec <- read_record(vcf, row$best_vcf_chr, row$best_vcf_pos, row$best_vcf_id, samples)
    target <- choose_target(row, rec$ref, rec$alt)
    dosage_strict <- parse_dosage(rec$gt, target$index, impute_missing_ref = FALSE)
    dosage_imputed <- parse_dosage(rec$gt, target$index, impute_missing_ref = TRUE)
    missing_alleles <- count_missing_alleles(rec$gt)
    detail <- data.table(
      source_set = row$source_set,
      lead_id = row$lead_id,
      match_type = row$best_match_type,
      vcf_id = row$best_vcf_id,
      sample = samples,
      in_covariate_set = samples %in% cov_samples,
      gt_original = as.character(rec$gt),
      gt_after_reference_imputation = impute_gt_display(rec$gt),
      n_missing_alleles_imputed_as_ref = missing_alleles,
      dosage_strict_missing = dosage_strict,
      dosage_after_reference_imputation = dosage_imputed
    )
    analysis <- detail[in_covariate_set == TRUE]
    strict_nonmissing <- analysis[!is.na(dosage_strict_missing)]
    imputed_nonmissing <- analysis[!is.na(dosage_after_reference_imputation)]
    n_overlap <- nrow(analysis)
    n_nonmissing_strict <- nrow(strict_nonmissing)
    n_nonmissing_imputed <- nrow(imputed_nonmissing)
    n_unique_strict <- uniqueN(strict_nonmissing$dosage_strict_missing)
    n_unique_imputed <- uniqueN(imputed_nonmissing$dosage_after_reference_imputation)
    classify <- function(n_nonmissing, n_unique) if (n_overlap < min_n) {
      "vcf_covariate_overlap_below_min_n"
    } else if (n_nonmissing < min_n) {
      "nonmissing_genotype_below_min_n"
    } else if (n_unique < 2L) {
      "monomorphic_in_analysis_samples"
    } else {
      "PASS"
    }
    reason_strict <- classify(n_nonmissing_strict, n_unique_strict)
    reason_imputed <- classify(n_nonmissing_imputed, n_unique_imputed)
    af_strict <- if (n_nonmissing_strict) sum(strict_nonmissing$dosage_strict_missing) / (2 * n_nonmissing_strict) else NA_real_
    af_imputed <- if (n_nonmissing_imputed) sum(imputed_nonmissing$dosage_after_reference_imputation) / (2 * n_nonmissing_imputed) else NA_real_
    summary <- data.table(
      status_before_imputation = if (reason_strict == "PASS") "PASS" else "FAIL",
      failure_reason_before_imputation = reason_strict,
      status_after_imputation = if (reason_imputed == "PASS") "PASS" else "FAIL",
      failure_reason_after_imputation = reason_imputed,
      rescued_by_reference_imputation = reason_strict != "PASS" & reason_imputed == "PASS",
      source_set = row$source_set,
      lead_id = row$lead_id,
      match_type = row$best_match_type,
      vcf_file = path.expand(vcf),
      vcf_id = row$best_vcf_id,
      vcf_chr = row$best_vcf_chr,
      vcf_pos = row$best_vcf_pos,
      vcf_ref = rec$ref,
      vcf_alt = rec$alt,
      dosage_target = as.character(target$index),
      dosage_allele_rule = target$rule,
      min_n = min_n,
      n_vcf_samples = length(samples),
      n_covariate_samples = length(cov_samples),
      n_vcf_covariate_overlap = n_overlap,
      n_samples_with_missing_gt_allele = sum(analysis$n_missing_alleles_imputed_as_ref > 0L, na.rm = TRUE),
      n_reference_alleles_imputed = sum(analysis$n_missing_alleles_imputed_as_ref, na.rm = TRUE),
      n_nonmissing_before_imputation = n_nonmissing_strict,
      n_nonmissing_after_imputation = n_nonmissing_imputed,
      n_unique_dosage_before_imputation = n_unique_strict,
      n_unique_dosage_after_imputation = n_unique_imputed,
      n_dosage_0_before = sum(strict_nonmissing$dosage_strict_missing == 0),
      n_dosage_1_before = sum(strict_nonmissing$dosage_strict_missing == 1),
      n_dosage_2_before = sum(strict_nonmissing$dosage_strict_missing == 2),
      n_dosage_0_after = sum(imputed_nonmissing$dosage_after_reference_imputation == 0),
      n_dosage_1_after = sum(imputed_nonmissing$dosage_after_reference_imputation == 1),
      n_dosage_2_after = sum(imputed_nonmissing$dosage_after_reference_imputation == 2),
      target_allele_frequency_before = af_strict,
      target_allele_frequency_after = af_imputed,
      gt_counts_all_vcf_samples = count_string(rec$gt),
      gt_counts_covariate_samples_before = count_string(analysis$gt_original),
      gt_counts_covariate_samples_after = count_string(analysis$gt_after_reference_imputation),
      dosage_counts_before = count_string(analysis$dosage_strict_missing),
      dosage_counts_after = count_string(analysis$dosage_after_reference_imputation),
      genotype_missing_allele_policy = "dot_to_reference: ./.=0/0; .|1=0|1; 1|.=1|0"
    )
    list(summary = summary, detail = detail)
  }, error = function(e) {
    list(
      summary = data.table(
        status_before_imputation = "ERROR", failure_reason_before_imputation = "record_or_dosage_error",
        status_after_imputation = "ERROR", failure_reason_after_imputation = "record_or_dosage_error",
        source_set = row$source_set, lead_id = row$lead_id,
        match_type = row$best_match_type, vcf_file = path.expand(vcf),
        vcf_id = row$best_vcf_id, error = conditionMessage(e)
      ),
      detail = data.table()
    )
  })
  summary_rows[[length(summary_rows) + 1L]] <- result$summary
  if (nrow(result$detail)) sample_rows[[length(sample_rows) + 1L]] <- result$detail
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
detail_dt <- rbindlist(sample_rows, fill = TRUE)
reason_dt <- summary_dt[, .(n_sv = .N), by = .(
  status_before_imputation, failure_reason_before_imputation,
  status_after_imputation, failure_reason_after_imputation,
  rescued_by_reference_imputation
)][order(status_after_imputation, failure_reason_after_imputation)]
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(summary_dt, paste0(out_prefix, ".tsv"), sep = "\t", na = "NA")
fwrite(detail_dt, paste0(out_prefix, ".samples.tsv.gz"), sep = "\t", na = "NA")
fwrite(reason_dt, paste0(out_prefix, ".reason_summary.tsv"), sep = "\t", na = "NA")

message("[DONE] Per-SV debug: ", paste0(out_prefix, ".tsv"))
message("[DONE] Per-sample debug: ", paste0(out_prefix, ".samples.tsv.gz"))
message("[DONE] Reason summary: ", paste0(out_prefix, ".reason_summary.tsv"))
print(reason_dt)
