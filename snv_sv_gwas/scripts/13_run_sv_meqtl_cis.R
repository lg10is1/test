#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (is.character(x) && length(x) == 1L && (is.na(x) || identical(x, ""))) return(y)
  x
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[sub("^--", "", key)]] <- TRUE
      i <- i + 1L
    } else {
      out[[sub("^--", "", key)]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

norm_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  paste0("chr", x)
}

path_exists_file <- function(x) {
  !is.na(x) && x != "" && file.exists(path.expand(x)) && !dir.exists(path.expand(x))
}

safe_id <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

parse_gt_dosage <- function(gt, target_index) {
  sample_names <- names(gt)
  gt <- sub(":.*$", "", as.character(gt))
  out <- rep(NA_real_, length(gt))
  names(out) <- sample_names
  ok <- !is.na(gt) & gt != ""
  if (!any(ok)) return(out)
  alleles <- strsplit(gt[ok], "[/|]")
  out[ok] <- vapply(
    alleles,
    function(z) {
      # Project policy: a missing allele is imputed as the reference allele.
      # Examples: ./. -> 0/0, .|1 -> 0|1, 1|. -> 1|0.
      z[z == "."] <- "0"
      if (identical(as.character(target_index), "ANY_ALT")) return(sum(z != "0"))
      sum(z == as.character(target_index))
    },
    numeric(1L)
  )
  out
}

count_missing_gt_alleles <- function(gt) {
  gt <- sub(":.*$", "", as.character(gt))
  out <- rep(NA_integer_, length(gt))
  ok <- !is.na(gt) & gt != ""
  out[ok] <- vapply(strsplit(gt[ok], "[/|]"), function(z) sum(z == "."), integer(1L))
  out
}

vcf_has_index <- function(vcf_file) {
  vcf_file <- path.expand(vcf_file)
  any(file.exists(paste0(vcf_file, c(".tbi", ".csi"))))
}

bcftools_view_region_args <- function(vcf_file, region) {
  if (vcf_has_index(vcf_file)) {
    c("view", "-H", "-r", region, path.expand(vcf_file))
  } else {
    c("view", "-H", "-t", region, path.expand(vcf_file))
  }
}

read_vcf_record_gt <- function(vcf, bcftools, chr, pos, vcf_id, sample_ids) {
  region <- paste0(chr, ":", pos, "-", pos)
  lines <- tryCatch(
    system2(
      path.expand(bcftools),
      bcftools_view_region_args(vcf, region),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
  if (!length(lines)) stop("No VCF record returned for ", vcf_id, " at ", region)

  dt <- fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    sep = "\t",
    showProgress = FALSE,
    data.table = FALSE
  )
  hit <- dt[as.character(dt[[3L]]) == vcf_id, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop("Expected exactly one VCF record with ID=", vcf_id, " at ", region, "; found ", nrow(hit))
  }
  if (ncol(hit) < 10L) stop("VCF record has no sample genotype columns: ", vcf_id)

  fmt <- strsplit(as.character(hit[[9L]]), ":", fixed = TRUE)[[1L]]
  gt_idx <- match("GT", fmt)
  if (is.na(gt_idx)) stop("FORMAT has no GT field for VCF record: ", vcf_id)

  sample_fields <- as.character(hit[1L, 10:ncol(hit)])
  names(sample_fields) <- sample_ids
  gt <- vapply(
    strsplit(sample_fields, ":", fixed = TRUE),
    function(z) if (length(z) >= gt_idx) z[[gt_idx]] else NA_character_,
    character(1L)
  )

  list(
    chrom = as.character(hit[[1L]]),
    pos = as.integer(hit[[2L]]),
    id = as.character(hit[[3L]]),
    ref = as.character(hit[[4L]]),
    alt = as.character(hit[[5L]]),
    gt = gt
  )
}

.vcf_sample_cache <- new.env(parent = emptyenv())

get_vcf_samples <- function(vcf, bcftools) {
  key <- path.expand(vcf)
  if (!exists(key, envir = .vcf_sample_cache, inherits = FALSE)) {
    samples <- system2(path.expand(bcftools), c("query", "-l", key), stdout = TRUE)
    if (!length(samples)) stop("No samples found in VCF: ", vcf)
    assign(key, samples, envir = .vcf_sample_cache)
  }
  get(key, envir = .vcf_sample_cache, inherits = FALSE)
}

target_allele_index <- function(vcf_ref, vcf_alt, target_allele) {
  target <- toupper(target_allele)
  ref <- toupper(vcf_ref)
  alts <- strsplit(toupper(vcf_alt), ",", fixed = TRUE)[[1L]]
  if (identical(target, ref)) return(0L)
  idx <- match(target, alts)
  if (is.na(idx)) {
    stop("Target allele not found in VCF REF/ALT. target=", target_allele,
         "; REF=", vcf_ref, "; ALT=", vcf_alt)
  }
  as.integer(idx)
}

is_true_value <- function(x) {
  length(x) == 1L && !is.na(x) && tolower(as.character(x)) %in% c("true", "t", "1")
}

select_dosage_target <- function(sv, vcf_ref, vcf_alt) {
  match_type <- as.character(sv$best_match_type)
  if (identical(match_type, "exact_ref_alt")) {
    return(list(
      index = target_allele_index(vcf_ref, vcf_alt, sv$lead_a2),
      rule = "exact_lead_a2_in_vcf_ref_or_alt"
    ))
  }
  if (!identical(match_type, "fuzzy_ref_same")) {
    stop("Unsupported meQTL match type: ", match_type)
  }
  if (is_true_value(sv$best_ref_same_as_a2)) {
    return(list(index = 0L, rule = "fuzzy_shared_ref_is_lead_a2"))
  }
  if (is_true_value(sv$best_ref_same_as_a1)) {
    return(list(index = "ANY_ALT", rule = "fuzzy_shared_ref_is_lead_a1_any_alt_proxy_for_a2"))
  }
  stop("Fuzzy match lacks a recorded shared-reference orientation for ", sv$lead_id)
}

extract_methyl_rows <- function(value_file, row_indices) {
  if (!length(row_indices)) return(data.table())
  keep_file <- tempfile("meqtl_rows_")
  on.exit(unlink(keep_file), add = TRUE)
  writeLines(as.character(as.integer(row_indices)), keep_file)
  cmd <- paste(
    "awk",
    shQuote('NR==FNR {keep[$1+1]=1; next} FNR==1 {print; next} (FNR in keep)'),
    shQuote(keep_file),
    shQuote(value_file)
  )
  fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
}

fit_one_cpg <- function(y, x_full, dosage_col, min_n) {
  if (length(y) != nrow(x_full)) {
    stop("Internal error: methylation vector length does not match design matrix rows.")
  }
  ok <- !is.na(y)
  n <- sum(ok)
  p <- ncol(x_full)
  if (n < max(min_n, p + 2L)) return(NULL)

  y <- as.numeric(y[ok])
  if (!is.finite(var(y)) || var(y) == 0) return(NULL)
  x <- x_full[ok, , drop = FALSE]
  qrx <- qr(x)
  if (qrx$rank < p) return(NULL)

  fit <- lm.fit(x, y)
  df <- n - p
  if (df <= 0L) return(NULL)
  rss <- sum(fit$residuals * fit$residuals)
  sigma2 <- rss / df
  xtx_inv <- chol2inv(qr.R(qrx))
  beta <- fit$coefficients[[dosage_col]]
  se <- sqrt(sigma2 * xtx_inv[dosage_col, dosage_col])
  if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)
  t_stat <- beta / se
  p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)

  list(n = n, beta = beta, se = se, t = t_stat, p = p_value, df = df)
}

plot_local_manhattan <- function(sv_res, lead_id, sv_chr, sv_pos, pdf_file, png_file) {
  if (!nrow(sv_res)) return(FALSE)
  dt <- copy(sv_res[is.finite(p) & p > 0 & !is.na(cpg_start)])
  if (!nrow(dt)) return(FALSE)
  dt[, log10_p := -log10(p)]
  setorder(dt, cpg_start)

  ymax <- max(dt$log10_p, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  xlab <- paste0(sv_chr, " CpG position")
  ylab <- "-log10(P value)"
  main <- paste0(lead_id, " cis meQTL")

  draw_plot <- function() {
    plot(
      dt$cpg_start,
      dt$log10_p,
      pch = 20,
      cex = 0.55,
      col = "grey25",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ylim = c(0, ymax * 1.08)
    )
    abline(v = sv_pos, col = "#D55E00", lty = 2, lwd = 1.2)
    best <- dt[which.min(p)]
    points(best$cpg_start, best$log10_p, pch = 20, cex = 0.9, col = "#0072B2")
    legend(
      "topright",
      legend = c("CpG", "SV position", "min P"),
      col = c("grey25", "#D55E00", "#0072B2"),
      pch = c(20, NA, 20),
      lty = c(NA, 2, NA),
      lwd = c(NA, 1.2, NA),
      bty = "n",
      cex = 0.8
    )
  }

  pdf(pdf_file, width = 7.2, height = 4.8)
  draw_plot()
  dev.off()

  png(png_file, width = 1800, height = 1200, res = 220)
  draw_plot()
  dev.off()

  TRUE
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

match_file <- args[["matches"]] %||%
  "/path/to/EOSCZ_PROJECT/figure_analysis/02.meQTL/public/SV/matching/sig_sv_to_pav_sv_len50.best.tsv"
pav_sv_vcf <- args[["pav-sv-vcf"]] %||% file.path(
  "/path/to/EOSCZ_PROJECT/TGS_SV_merge_SCZ/truvari_single_sample",
  paste0("tru", "vari_merged_sort_pP0.5.sv_len_gt50.sorted.vcf.gz")
)
methyl_dir <- args[["methyl-dir"]] %||%
  "/path/to/SCZ/methy_pbcpg/qc_results"
sample_info_file <- args[["sample-info"]] %||%
  "/path/to/EOSCZ_PROJECT/Sample_info/sample_case.txt"
genomic_pc_file <- args[["genomic-pc"]] %||%
  "/path/to/EOSCZ_PROJECT/TGS_callset/glnexus_merge/autosome_pca/unrelated_projected_pca/all_samples.projected_pca.PC1_PC2.projected_scores.tsv"
epigenetic_pc_file <- args[["epigenetic-pc"]] %||%
  "/path/to/SCZ/methy_pbcpg/pca_600k/methylation_uniform_600k_pca.scores.tsv"
bcftools_bin <- args[["bcftools"]] %||% Sys.which("bcftools")
out_dir <- args[["out-dir"]] %||%
  "/path/to/EOSCZ_PROJECT/figure_analysis/02.meQTL/public/SV"
window_bp <- as.integer(args[["window-bp"]] %||% 1000000L)
min_n <- as.integer(args[["min-n"]] %||% 50L)
genomic_pc_n <- as.integer(args[["genomic-pc-n"]] %||% 20L)
epigenetic_pc_n <- as.integer(args[["epigenetic-pc-n"]] %||% 20L)
match_types <- trimws(strsplit(
  args[["match-types"]] %||% "exact_ref_alt,fuzzy_ref_same", ",", fixed = TRUE
)[[1L]])

if (!path_exists_file(match_file)) stop("Missing --matches: ", match_file)
if (!path_exists_file(pav_sv_vcf)) stop("Missing --pav-sv-vcf: ", pav_sv_vcf)
if (!path_exists_file(sample_info_file)) stop("Missing --sample-info: ", sample_info_file)
if (!path_exists_file(genomic_pc_file)) stop("Missing --genomic-pc: ", genomic_pc_file)
if (!path_exists_file(epigenetic_pc_file)) stop("Missing --epigenetic-pc: ", epigenetic_pc_file)
if (bcftools_bin == "" || !file.exists(path.expand(bcftools_bin))) stop("bcftools not found; pass --bcftools")
if (is.na(window_bp) || window_bp < 1L) stop("--window-bp must be positive")
if (is.na(min_n) || min_n < 1L) stop("--min-n must be positive")

dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)
per_sv_dir <- file.path(path.expand(out_dir), "per_lead")
dir.create(per_sv_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(path.expand(out_dir), "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(path.expand(out_dir), "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] Loading matches: ", match_file)
matches <- fread(path.expand(match_file), showProgress = FALSE)
matches <- matches[best_match_type %in% match_types]
if (!"lead_chr" %in% names(matches)) stop("Match table lacks lead_chr required for chr1-22 filtering")
matches[, lead_chr_filter := sub("^chr", "", as.character(lead_chr), ignore.case = TRUE)]
matches <- matches[lead_chr_filter %in% as.character(1:22)]
matches[, lead_chr_filter := NULL]
if (!nrow(matches)) stop("No rows selected from match file for match_types=", paste(match_types, collapse = ","))

message("[INFO] Selected SVs: ", nrow(matches), " | match_types=", paste(match_types, collapse = ","))
fwrite(
  matches,
  file.path(table_dir, "meqtl_selected_matches.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  data.table(
    parameter = c("genotype_source", "missing_allele_policy", "examples", "match_types"),
    value = c(
      path.expand(pav_sv_vcf),
      "Each missing GT allele (.) is imputed as reference allele 0 before dosage calculation",
      "./. -> 0/0; .|1 -> 0|1; 1|. -> 1|0",
      paste(match_types, collapse = ",")
    )
  ),
  file.path(table_dir, "meqtl_cis.genotype_policy.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fallback_vcf_samples <- get_vcf_samples(pav_sv_vcf, bcftools_bin)
message("[INFO] Fallback VCF samples: ", length(fallback_vcf_samples), " | ", pav_sv_vcf)

sample_info <- fread(path.expand(sample_info_file), showProgress = FALSE)
need_cov <- c("SampleID", "Sex", "Age")
missing_cov <- setdiff(need_cov, names(sample_info))
if (length(missing_cov)) stop("Sample info missing columns: ", paste(missing_cov, collapse = ", "))
sample_info <- sample_info[, .(sample = as.character(SampleID), Sex = as.factor(Sex), Age = as.numeric(Age))]

gpc <- fread(path.expand(genomic_pc_file), showProgress = FALSE)
if (!"IID" %in% names(gpc)) stop("Genomic PC file must contain IID")
gpc_cols <- paste0("PC", seq_len(genomic_pc_n))
missing_gpc <- setdiff(gpc_cols, names(gpc))
if (length(missing_gpc)) stop("Genomic PC file missing: ", paste(missing_gpc, collapse = ", "))
gpc <- gpc[, c("IID", gpc_cols), with = FALSE]
setnames(gpc, c("IID", gpc_cols), c("sample", paste0("gPC", seq_len(genomic_pc_n))))

epc <- fread(path.expand(epigenetic_pc_file), showProgress = FALSE)
if (!"sample" %in% names(epc)) stop("Epigenetic PC file must contain sample")
epc_cols <- paste0("PC", seq_len(epigenetic_pc_n))
missing_epc <- setdiff(epc_cols, names(epc))
if (length(missing_epc)) stop("Epigenetic PC file missing: ", paste(missing_epc, collapse = ", "))
epc <- epc[, c("sample", epc_cols), with = FALSE]
setnames(epc, c("sample", epc_cols), c("sample", paste0("ePC", seq_len(epigenetic_pc_n))))

cov <- merge(sample_info, gpc, by = "sample", all = FALSE)
cov <- merge(cov, epc, by = "sample", all = FALSE)
if (!nrow(cov)) stop("No overlapping sample_info/genomic_pc/epigenetic_pc samples")
message("[INFO] Samples with covariates and PCs: ", nrow(cov))

all_results <- list()
summary_rows <- list()

skip_summary <- function(sv, lead_id, chr, pos, reason) {
  data.table(
    status = "SKIP",
    variant_class = "SV",
    source_set = sv$source_set,
    lead_id = lead_id,
    sv_chr = chr,
    sv_pos = pos,
    match_type = sv$best_match_type,
    error = reason
  )
}

for (sv_i in seq_len(nrow(matches))) {
  sv <- matches[sv_i]
  lead_id <- sv$lead_id
  chr <- norm_chr(sv$lead_chr)
  chr_no <- sub("^chr", "", chr)
  sv_pos <- as.integer(sv$lead_pos)
  value_file <- file.path(path.expand(methyl_dir), paste0(chr, "_col9_val.tsv"))
  pos_file <- file.path(path.expand(methyl_dir), paste0(chr, "_positions.tsv"))
  if (!file.exists(value_file) && chr_no != chr) {
    value_file <- file.path(path.expand(methyl_dir), paste0("chr", chr_no, "_col9_val.tsv"))
  }
  if (!file.exists(pos_file) && chr_no != chr) {
    pos_file <- file.path(path.expand(methyl_dir), paste0("chr", chr_no, "_positions.tsv"))
  }
  if (!file.exists(value_file) || !file.exists(pos_file)) {
    warning("Skipping ", lead_id, ": missing methylation file or position file for ", chr)
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, "missing methylation value or position file")
    next
  }

  sv_vcf <- pav_sv_vcf
  sv_vcf_source <- "fallback_pav_sv_vcf"
  if ("best_vcf_file" %in% names(matches) &&
      !is.na(sv$best_vcf_file) &&
      as.character(sv$best_vcf_file) != "") {
    sv_vcf <- as.character(sv$best_vcf_file)
    sv_vcf_source <- if ("best_vcf_source" %in% names(matches)) as.character(sv$best_vcf_source) else "best_vcf_file"
  }
  if (!path_exists_file(sv_vcf)) {
    warning("Skipping ", lead_id, ": matched VCF file is missing: ", sv_vcf)
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, paste0("matched VCF missing: ", sv_vcf))
    next
  }
  vcf_samples <- get_vcf_samples(sv_vcf, bcftools_bin)

  message("[SV ", sv_i, "/", nrow(matches), "] ", lead_id, " | ", chr, ":", sv_pos,
          " | VCF=", sv$best_vcf_id, " | source=", sv_vcf_source)

  rec <- read_vcf_record_gt(
    sv_vcf,
    bcftools_bin,
    sv$best_vcf_chr,
    sv$best_vcf_pos,
    sv$best_vcf_id,
    vcf_samples
  )
  dosage_target <- select_dosage_target(sv, rec$ref, rec$alt)
  target_idx <- dosage_target$index
  dosage <- parse_gt_dosage(rec$gt, target_idx)
  missing_gt_alleles <- count_missing_gt_alleles(rec$gt)
  if (is.null(names(dosage)) || any(names(dosage) == "")) {
    stop("Internal error: genotype dosage vector has missing sample names for ", lead_id)
  }
  dosage_dt <- data.table(
    sample = names(dosage),
    dosage = as.numeric(dosage),
    gt_missing_alleles_imputed_as_ref = missing_gt_alleles
  )

  sv_cov <- merge(cov[sample %in% vcf_samples], dosage_dt, by = "sample", all = FALSE)
  sv_cov <- sv_cov[!is.na(dosage)]
  if (uniqueN(sv_cov$dosage) < 2L || nrow(sv_cov) < min_n) {
    reason <- if (nrow(sv_cov) < min_n) "genotype sample size below min_n after reference imputation" else "monomorphic dosage after reference imputation"
    warning("Skipping ", lead_id, ": ", reason)
    skipped <- skip_summary(sv, lead_id, chr, sv_pos, reason)
    skipped[, `:=`(
      vcf_source = sv_vcf_source,
      vcf_id = sv$best_vcf_id,
      vcf_file = path.expand(sv_vcf),
      n_genotype_nonmissing = nrow(sv_cov),
      n_genotype_samples_with_imputed_ref = sum(sv_cov$gt_missing_alleles_imputed_as_ref > 0L, na.rm = TRUE),
      n_reference_alleles_imputed = sum(sv_cov$gt_missing_alleles_imputed_as_ref, na.rm = TRUE),
      dosage_values = paste(sort(unique(sv_cov$dosage)), collapse = ","),
      genotype_missing_allele_policy = "dot_to_reference: ./.=0/0; .|1=0|1; 1|.=1|0"
    )]
    summary_rows[[length(summary_rows) + 1L]] <- skipped
    next
  }

  positions <- fread(pos_file, showProgress = FALSE)
  required_pos_cols <- c("chr", "start", "end")
  if (!all(required_pos_cols %in% names(positions))) {
    stop("Position file missing chr/start/end columns: ", pos_file)
  }
  positions[, row_index := .I]
  cis <- positions[
    start >= (sv_pos - window_bp) &
      start <= (sv_pos + window_bp)
  ]
  if (!nrow(cis)) {
    warning("Skipping ", lead_id, ": no CpGs in cis window")
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, "no CpGs in cis window")
    next
  }

  meth <- extract_methyl_rows(value_file, cis$row_index)
  if (!nrow(meth)) {
    warning("Skipping ", lead_id, ": methylation extraction returned no rows")
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, "methylation extraction returned no rows")
    next
  }

  meth_samples <- intersect(names(meth), sv_cov$sample)
  if (length(meth_samples) < min_n) {
    warning("Skipping ", lead_id, ": too few methylation samples overlap covariates")
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, "too few methylation/covariate overlaps")
    next
  }
  sv_cov <- sv_cov[match(meth_samples, sample)]
  covariate_cols <- c(
    "dosage",
    "Sex",
    "Age",
    paste0("gPC", seq_len(genomic_pc_n)),
    paste0("ePC", seq_len(epigenetic_pc_n))
  )
  complete_design <- complete.cases(sv_cov[, covariate_cols, with = FALSE])
  if (!all(complete_design)) {
    sv_cov <- sv_cov[complete_design]
    meth_samples <- sv_cov$sample
  }
  if (length(meth_samples) < min_n) {
    warning("Skipping ", lead_id, ": too few samples after complete covariate filtering")
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, "too few complete covariate samples")
    next
  }
  if (uniqueN(sv_cov$dosage) < 2L) {
    warning("Skipping ", lead_id, ": insufficient genotype variation after complete covariate filtering")
    summary_rows[[length(summary_rows) + 1L]] <- skip_summary(sv, lead_id, chr, sv_pos, "insufficient genotype variation after filtering")
    next
  }
  n_imputed_samples_final <- sum(sv_cov$gt_missing_alleles_imputed_as_ref > 0L, na.rm = TRUE)
  n_imputed_alleles_final <- sum(sv_cov$gt_missing_alleles_imputed_as_ref, na.rm = TRUE)
  design_formula <- as.formula(
    paste(
      "~ dosage + Sex + Age +",
      paste(c(paste0("gPC", seq_len(genomic_pc_n)), paste0("ePC", seq_len(epigenetic_pc_n))), collapse = " + ")
    )
  )
  x_full <- model.matrix(design_formula, data = sv_cov)
  dosage_col <- match("dosage", colnames(x_full))
  if (is.na(dosage_col)) stop("Design matrix does not contain dosage column")

  meth_mat <- as.matrix(meth[, meth_samples, with = FALSE])
  storage.mode(meth_mat) <- "double"
  if (nrow(meth_mat) != nrow(cis)) {
    stop("Extracted methylation rows do not match cis position rows for ", lead_id)
  }
  if (ncol(meth_mat) != nrow(x_full)) {
    stop("Internal error: methylation sample columns do not match design matrix rows for ", lead_id)
  }

  res_list <- vector("list", nrow(meth_mat))
  for (j in seq_len(nrow(meth_mat))) {
    fit <- fit_one_cpg(meth_mat[j, ], x_full, dosage_col, min_n)
    if (is.null(fit)) next
    res_list[[j]] <- data.table(
      variant_class = "SV",
      lead_id = lead_id,
      source_set = sv$source_set,
      match_type = sv$best_match_type,
      allele_source_orientation = sv$best_allele_source_orientation,
      vcf_source = sv_vcf_source,
      vcf_file = path.expand(sv_vcf),
      vcf_id = sv$best_vcf_id,
      sv_chr = chr,
      sv_pos = sv_pos,
      vcf_chr = sv$best_vcf_chr,
      vcf_pos = sv$best_vcf_pos,
      vcf_ref = sv$best_vcf_ref,
      vcf_alt = sv$best_vcf_alt,
      target_allele = sv$lead_a2,
      target_allele_index = target_idx,
      dosage_allele_rule = dosage_target$rule,
      genotype_missing_allele_policy = "dot_to_reference",
      genotype_samples_with_imputed_ref = n_imputed_samples_final,
      reference_alleles_imputed = n_imputed_alleles_final,
      cpg_chr = cis$chr[j],
      cpg_start = cis$start[j],
      cpg_end = cis$end[j],
      distance_bp = cis$start[j] - sv_pos,
      n = fit$n,
      beta = fit$beta,
      se = fit$se,
      t = fit$t,
      p = fit$p,
      df = fit$df,
      methyl_missing = sum(is.na(meth_mat[j, ])),
      genotype_nonmissing = nrow(sv_cov),
      dosage_mean = mean(sv_cov$dosage),
      dosage_sd = sd(sv_cov$dosage)
    )
  }

  sv_res <- rbindlist(res_list, fill = TRUE)
  source_per_dir <- file.path(per_sv_dir, safe_id(sv$source_set))
  source_plot_dir <- file.path(plot_dir, safe_id(sv$source_set))
  dir.create(source_per_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_plot_dir, recursive = TRUE, showWarnings = FALSE)
  out_sv <- file.path(source_per_dir, paste0(safe_id(lead_id), ".cis_meqtl.tsv.gz"))
  out_plot_pdf <- file.path(source_plot_dir, paste0(safe_id(lead_id), ".local_manhattan.pdf"))
  out_plot_png <- file.path(source_plot_dir, paste0(safe_id(lead_id), ".local_manhattan.png"))
  if (nrow(sv_res)) {
    sv_res[, p_bonferroni_sv := pmin(p * .N, 1)]
    sv_res[, fdr_bh_sv := p.adjust(p, method = "BH")]
    fwrite(sv_res, out_sv, sep = "\t", quote = FALSE, na = "NA")
    plot_local_manhattan(sv_res, lead_id, chr, sv_pos, out_plot_pdf, out_plot_png)
    all_results[[length(all_results) + 1L]] <- sv_res
  } else {
    fwrite(data.table(), out_sv, sep = "\t")
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    status = "PASS",
    variant_class = "SV",
    lead_id = lead_id,
    source_set = sv$source_set,
    match_type = sv$best_match_type,
    vcf_source = sv_vcf_source,
    vcf_id = sv$best_vcf_id,
    vcf_file = path.expand(sv_vcf),
    sv_chr = chr,
    sv_pos = sv_pos,
    n_cpg_in_window = nrow(cis),
    n_cpg_tested = nrow(sv_res),
    n_genotype_nonmissing = nrow(sv_cov),
    dosage_values = paste(sort(unique(sv_cov$dosage)), collapse = ","),
    dosage_allele_rule = dosage_target$rule,
    genotype_missing_allele_policy = "dot_to_reference: ./.=0/0; .|1=0|1; 1|.=1|0",
    n_genotype_samples_with_imputed_ref = n_imputed_samples_final,
    n_reference_alleles_imputed = n_imputed_alleles_final,
    out_file = out_sv,
    plot_pdf = if (nrow(sv_res)) out_plot_pdf else NA_character_,
    plot_png = if (nrow(sv_res)) out_plot_png else NA_character_,
    error = NA_character_
  )
}

combined <- rbindlist(all_results, fill = TRUE)
summary_dt <- rbindlist(summary_rows, fill = TRUE)

out_all <- file.path(table_dir, "meqtl_cis.all_results.tsv.gz")
out_summary <- file.path(table_dir, "meqtl_cis.summary.tsv")

if (nrow(combined)) {
  combined[, p_bonferroni_all := pmin(p * .N, 1)]
  combined[, fdr_bh_all := p.adjust(p, method = "BH")]
  setorder(combined, p)
}

fwrite(combined, out_all, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary_dt, out_summary, sep = "\t", quote = FALSE, na = "NA")
if (nrow(combined)) {
  top_per_lead <- combined[, .SD[which.min(p)], by = .(source_set, lead_id)]
  significant <- combined[fdr_bh_all < 0.05]
  fwrite(top_per_lead, file.path(table_dir, "meqtl_cis.top_per_lead.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(significant, file.path(table_dir, "meqtl_cis.global_fdr05.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
}

message("[DONE] Combined results: ", out_all)
message("[DONE] Summary: ", out_summary)
message("[DONE] Per-SV directory: ", per_sv_dir)
message("[DONE] Plot directory: ", plot_dir)
