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
base <- "/path/to/EOSCZ_PROJECT/figure_analysis/02.meQTL/public/SV/tables"
input_file <- path.expand(args[["input"]] %||% file.path(base, "meqtl_cis.all_results.tsv.gz"))
out_dir <- path.expand(args[["out-dir"]] %||% base)
mac_thresholds <- as.integer(strsplit(args[["mac-thresholds"]] %||% "2,5", ",", fixed = TRUE)[[1L]])
fdr_threshold <- as.numeric(args[["fdr-threshold"]] %||% 0.05)
p_threshold <- as.numeric(args[["p-threshold"]] %||% 1e-8)

if (!file.exists(input_file) || dir.exists(input_file) || file.info(input_file)$size == 0) {
  stop("Missing or empty input: ", input_file)
}
if (anyNA(mac_thresholds) || any(mac_thresholds < 0L)) stop("Invalid --mac-thresholds")
mac_thresholds <- sort(unique(mac_thresholds))
if (!is.finite(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold > 1) stop("Invalid --fdr-threshold")
if (!is.finite(p_threshold) || p_threshold <= 0 || p_threshold > 1) stop("Invalid --p-threshold")

dt <- fread(input_file, showProgress = FALSE)
required <- c("source_set", "lead_id", "p", "genotype_nonmissing", "dosage_mean")
missing <- setdiff(required, names(dt))
if (length(missing)) stop("Input lacks required columns: ", paste(missing, collapse = ", "))
dt[, p := as.numeric(p)]
dt[, genotype_nonmissing := as.integer(genotype_nonmissing)]
dt[, dosage_mean := as.numeric(dosage_mean)]
dt <- dt[is.finite(p) & p >= 0 & p <= 1]
if (!nrow(dt)) stop("No valid association P values in input")

optional_first <- function(x) if (length(x)) as.character(x[1L]) else NA_character_
per_sv <- dt[, .(
  n_associations_original = .N,
  n_distinct_genotype_n = uniqueN(genotype_nonmissing),
  n_distinct_dosage_mean = uniqueN(round(dosage_mean, 12)),
  genotype_n = genotype_nonmissing[1L],
  dosage_mean = dosage_mean[1L],
  match_type = if ("match_type" %in% names(.SD)) optional_first(match_type) else NA_character_,
  vcf_id = if ("vcf_id" %in% names(.SD)) optional_first(vcf_id) else NA_character_
), by = .(source_set, lead_id)]

bad <- per_sv[
  n_distinct_genotype_n != 1L | n_distinct_dosage_mean != 1L |
    is.na(genotype_n) | genotype_n < 1L | !is.finite(dosage_mean)
]
if (nrow(bad)) {
  stop("Inconsistent or missing genotype summaries for ", nrow(bad), " SV(s); cannot calculate auditable MAC")
}

# Dosage is the count of the modeled target allele. Convert its mean back to
# an integer allele count, then take the smaller of target and other allele
# counts. This definition remains valid whether the modeled target is REF or ALT.
per_sv[, target_allele_count := as.integer(round(dosage_mean * genotype_n))]
per_sv[, total_allele_count := 2L * genotype_n]
per_sv[, other_allele_count := total_allele_count - target_allele_count]
if (per_sv[target_allele_count < 0L | target_allele_count > total_allele_count, .N]) {
  stop("Calculated target allele count falls outside [0, 2N]")
}
per_sv[, minor_allele_count := pmin(target_allele_count, other_allele_count)]
per_sv[, minor_allele_frequency := minor_allele_count / total_allele_count]
setorder(per_sv, source_set, lead_id)

dt <- merge(
  dt,
  per_sv[, .(
    source_set, lead_id, genotype_n, target_allele_count,
    other_allele_count, minor_allele_count, minor_allele_frequency
  )],
  by = c("source_set", "lead_id"), all.x = TRUE, sort = FALSE
)

count_unique_cpg <- function(x) {
  cols <- intersect(c("cpg_chr", "cpg_start", "cpg_end"), names(x))
  if (!length(cols)) return(NA_integer_)
  uniqueN(x, by = cols)
}

summary_rows <- list()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(per_sv, file.path(out_dir, "meqtl_cis.per_sv_mac.tsv"), sep = "\t", na = "NA")

for (mac_threshold in mac_thresholds) {
  label <- paste0("mac_ge", mac_threshold)
  kept <- copy(dt[minor_allele_count >= mac_threshold])
  if (nrow(kept)) {
    kept[, fdr_bh_within_mac_filter := p.adjust(p, method = "BH")]
  } else {
    kept[, fdr_bh_within_mac_filter := numeric()]
  }
  fdr_hits <- kept[is.finite(fdr_bh_within_mac_filter) & fdr_bh_within_mac_filter < fdr_threshold]
  p_hits <- kept[p < p_threshold]
  both_hits <- kept[
    p < p_threshold & is.finite(fdr_bh_within_mac_filter) &
      fdr_bh_within_mac_filter < fdr_threshold
  ]

  fwrite(
    kept,
    file.path(out_dir, paste0("meqtl_cis.", label, ".all_results.tsv.gz")),
    sep = "\t", na = "NA"
  )
  fwrite(
    fdr_hits,
    file.path(out_dir, paste0("meqtl_cis.", label, ".fdr_lt_0.05.tsv.gz")),
    sep = "\t", na = "NA"
  )
  fwrite(
    p_hits,
    file.path(out_dir, paste0("meqtl_cis.", label, ".p_lt_1e-8.tsv.gz")),
    sep = "\t", na = "NA"
  )

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    mac_threshold = mac_threshold,
    mac_definition = "min(sum(target_allele_dosage), 2N-sum(target_allele_dosage))",
    fdr_method = "Benjamini-Hochberg",
    fdr_scope = paste0("all SV-CpG associations retained after MAC >= ", mac_threshold),
    fdr_threshold = fdr_threshold,
    p_threshold = p_threshold,
    n_sv_before_mac_filter = uniqueN(dt, by = c("source_set", "lead_id")),
    n_sv_after_mac_filter = uniqueN(kept, by = c("source_set", "lead_id")),
    n_sv_excluded_by_mac = uniqueN(dt, by = c("source_set", "lead_id")) - uniqueN(kept, by = c("source_set", "lead_id")),
    n_associations_tested = nrow(kept),
    n_unique_cpg_tested = count_unique_cpg(kept),
    n_associations_fdr_lt_0.05 = nrow(fdr_hits),
    n_sv_with_fdr_lt_0.05 = uniqueN(fdr_hits, by = c("source_set", "lead_id")),
    n_unique_cpg_fdr_lt_0.05 = count_unique_cpg(fdr_hits),
    n_associations_p_lt_1e_8 = nrow(p_hits),
    n_sv_with_p_lt_1e_8 = uniqueN(p_hits, by = c("source_set", "lead_id")),
    n_unique_cpg_p_lt_1e_8 = count_unique_cpg(p_hits),
    n_associations_both_fdr_and_p = nrow(both_hits),
    n_sv_with_both_fdr_and_p = uniqueN(both_hits, by = c("source_set", "lead_id")),
    minimum_p = if (nrow(kept)) min(kept$p) else NA_real_
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(out_dir, "meqtl_cis.mac_filter_summary.tsv"), sep = "\t", na = "NA")

message("[DONE] Original meQTL table preserved: ", input_file)
message("[DONE] Per-SV MAC audit: ", file.path(out_dir, "meqtl_cis.per_sv_mac.tsv"))
message("[DONE] MAC filter summary: ", file.path(out_dir, "meqtl_cis.mac_filter_summary.tsv"))
print(summary_dt)
