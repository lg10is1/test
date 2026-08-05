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

parse_info_value <- function(info, key) {
  info <- as.character(info)
  out <- rep(NA_character_, length(info))
  pattern <- paste0("(^|;)", key, "=")
  hit <- grepl(pattern, info)
  if (any(hit)) {
    out[hit] <- sub(
      paste0("^.*(^|;)", key, "=([^;]*).*$"),
      "\\2",
      info[hit]
    )
  }
  out
}

path_exists_file <- function(x) {
  !is.na(x) && x != "" && file.exists(path.expand(x)) && !dir.exists(path.expand(x))
}

split_alt <- function(x) {
  strsplit(x, ",", fixed = TRUE)
}

get_vcf_chroms <- function(vcf_file, bcftools_bin) {
  idx <- tryCatch(
    system2(path.expand(bcftools_bin), c("index", "-s", path.expand(vcf_file)), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(idx)) {
    chroms <- strsplit(idx, "\t", fixed = TRUE)
    chroms <- vapply(chroms, `[`, character(1L), 1L)
    chroms <- chroms[chroms != ""]
    if (length(chroms)) return(unique(chroms))
  }

  hdr <- tryCatch(
    system2(path.expand(bcftools_bin), c("view", "-h", path.expand(vcf_file)), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  contig <- grep("^##contig=<ID=", hdr, value = TRUE)
  chroms <- sub("^##contig=<ID=([^,>]+).*$", "\\1", contig)
  chroms <- chroms[chroms != ""]
  unique(chroms)
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

make_vcf_chr_mapper <- function(vcf_chroms) {
  if (!length(vcf_chroms)) {
    warning("Could not detect VCF contigs; using chr-prefixed query names.")
    return(function(x) norm_chr(x))
  }
  has_chr <- any(grepl("^chr", vcf_chroms, ignore.case = TRUE))
  no_chr <- any(!grepl("^chr", vcf_chroms, ignore.case = TRUE))

  function(x) {
    x_norm <- norm_chr(x)
    if (has_chr && x_norm %in% vcf_chroms) return(x_norm)

    x_nochr <- sub("^chr", "", x_norm, ignore.case = TRUE)
    if (no_chr && x_nochr %in% vcf_chroms) return(x_nochr)

    if (has_chr) x_norm else x_nochr
  }
}

parse_bim_map <- function(x) {
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  parts <- parts[parts != ""]
  out <- character()
  for (part in parts) {
    kv <- strsplit(part, "=", fixed = TRUE)[[1L]]
    if (length(kv) != 2L) stop("Invalid --bim-map entry: ", part)
    out[kv[[1L]]] <- kv[[2L]]
  }
  out
}

read_bim_one <- function(file, label) {
  if (!path_exists_file(file)) {
    warning("BIM not found for ", label, ": ", file)
    return(data.table())
  }
  x <- fread(
    path.expand(file),
    header = FALSE,
    col.names = c("bim_chr", "lead_id", "cm", "bim_bp", "bim_a1", "bim_a2"),
    showProgress = FALSE
  )
  x[, source_set := label]
  x[, bim_chr := norm_chr(bim_chr)]
  x[, bim_bp := as.integer(bim_bp)]
  x[, bim_a1 := as.character(bim_a1)]
  x[, bim_a2 := as.character(bim_a2)]
  x[]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

sig_file <- args[["sig-file"]] %||%
  "/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public/tables/lead_sig_from_gwas.canonical_1000kb.tsv"
pav_sv_vcf <- args[["pav-sv-vcf"]] %||% file.path(
  "/path/to/EOSCZ_PROJECT/TGS_SV_merge_SCZ/truvari_single_sample",
  paste0("tru", "vari_merged_sort_pP0.5.sv_len_gt50.sorted.vcf.gz")
)
paragraph_vcf <- args[["paragraph-vcf"]] %||% NA_character_
bcftools_bin <- args[["bcftools"]] %||% Sys.which("bcftools")
out_dir <- args[["out-dir"]] %||%
  "/path/to/EOSCZ_PROJECT/figure_analysis/02.meQTL/public/SV/matching"
window_bp <- as.integer(args[["window-bp"]] %||% 1000L)
allow_table_fallback <- !isTRUE(args[["no-table-fallback"]])

if (is.na(window_bp) || window_bp < 0L) stop("--window-bp must be a non-negative integer.")
if (!path_exists_file(sig_file)) stop("Missing --sig-file: ", sig_file)
if (!path_exists_file(pav_sv_vcf)) stop("Missing --pav-sv-vcf: ", pav_sv_vcf)
has_paragraph_vcf <- !is.na(paragraph_vcf) && paragraph_vcf != ""
if (has_paragraph_vcf && !path_exists_file(paragraph_vcf)) {
  stop("Missing --paragraph-vcf: ", paragraph_vcf)
}
if (bcftools_bin == "" || !file.exists(path.expand(bcftools_bin))) {
  stop("bcftools not found. Pass --bcftools /path/to/bcftools")
}

default_bims <- c(
  ## Public example: set00 only. Add set01/set02 here if available.
  set00 = "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/NGS.QCsite.QCind.bim",
  paragraph = "/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph/chr_all2.strict_step2_genimi.common_samples.merged.bim"
)
bim_map <- c(default_bims, parse_bim_map(args[["bim-map"]] %||% NA_character_))
bim_map <- bim_map[!duplicated(names(bim_map), fromLast = TRUE)]

dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)
out_prefix <- file.path(path.expand(out_dir), "sig_sv_to_pav_sv_len50")
out_bim <- paste0(out_prefix, ".from_bim.tsv")
out_candidates <- paste0(out_prefix, ".candidates.tsv")
out_best <- paste0(out_prefix, ".best.tsv")
out_summary <- paste0(out_prefix, ".summary.tsv")

message("[INFO] Significant SV table: ", sig_file)
message("[INFO] PAV SV VCF: ", pav_sv_vcf)
if (has_paragraph_vcf) message("[INFO] Paragraph VCF: ", paragraph_vcf)
message("[INFO] Window: +/-", window_bp, " bp")

vcf_chroms <- get_vcf_chroms(pav_sv_vcf, bcftools_bin)
to_vcf_chr <- make_vcf_chr_mapper(vcf_chroms)
message("[INFO] Detected VCF contigs: ", paste(head(vcf_chroms, 10), collapse = ", "),
        if (length(vcf_chroms) > 10L) ", ..." else "")
paragraph_chroms <- character()
to_paragraph_vcf_chr <- to_vcf_chr
if (has_paragraph_vcf) {
  paragraph_chroms <- get_vcf_chroms(paragraph_vcf, bcftools_bin)
  to_paragraph_vcf_chr <- make_vcf_chr_mapper(paragraph_chroms)
  message("[INFO] Detected paragraph VCF contigs: ", paste(head(paragraph_chroms, 10), collapse = ", "),
          if (length(paragraph_chroms) > 10L) ", ..." else "")
}

sig <- fread(path.expand(sig_file), showProgress = FALSE)
required <- c("lead_id", "variant_type", "lead_chr")
missing <- setdiff(required, names(sig))
if (length(missing)) stop("Significant table missing columns: ", paste(missing, collapse = ", "))
sig[, lead_chr_filter := sub("^chr", "", as.character(lead_chr), ignore.case = TRUE)]
sig_sv <- sig[tolower(variant_type) == "sv" & lead_chr_filter %in% as.character(1:22)]
sig_sv[, lead_chr_filter := NULL]
if (!nrow(sig_sv)) stop("No variant_type == sv rows found in: ", sig_file)
sig_sv[, lead_id := as.character(lead_id)]
if (!"source_set" %in% names(sig_sv)) sig_sv[, source_set := NA_character_]

bim <- rbindlist(
  lapply(names(bim_map), function(label) read_bim_one(bim_map[[label]], label)),
  fill = TRUE
)

if (!nrow(bim) && !allow_table_fallback) {
  stop("No BIM rows were loaded and --no-table-fallback was set.")
}

setkey(bim, source_set, lead_id)
from_bim <- merge(
  sig_sv,
  bim,
  by = c("source_set", "lead_id"),
  all.x = TRUE,
  allow.cartesian = TRUE
)

if (!"source_info" %in% names(from_bim)) {
  from_bim[, source_info := NA_character_]
}
from_bim[!is.na(bim_bp), source_info := "bim"]

if (allow_table_fallback) {
  missing_bim <- is.na(from_bim$bim_bp)
  if (any(missing_bim)) {
    first_present <- function(candidates) {
      hit <- candidates[candidates %in% names(from_bim)]
      if (length(hit)) hit[1L] else NA_character_
    }
    chr_col <- first_present(c("lead_chr", "CHR", "CHROM", "#CHROM"))
    pos_col <- first_present(c("lead_pos", "START", "POS", "BP"))
    a1_col <- first_present(c("A1", "original_REF_or_A2", "REF"))
    a2_col <- first_present(c("A2", "original_ALT_or_A1", "ALT"))

    if (is.na(chr_col) || is.na(pos_col)) {
      warning(
        "Table fallback is needed for ", sum(missing_bim),
        " row(s), but no compatible chromosome/position columns are available."
      )
    } else {
      from_bim[missing_bim, bim_chr := norm_chr(get(chr_col))]
      from_bim[missing_bim, bim_bp := suppressWarnings(as.integer(get(pos_col)))]
      if (!is.na(a1_col)) from_bim[missing_bim, bim_a1 := as.character(get(a1_col))]
      if (!is.na(a2_col)) from_bim[missing_bim, bim_a2 := as.character(get(a2_col))]
      from_bim[missing_bim, source_info := "sig_table_fallback"]
    }
  }
}

from_bim[is.na(source_info), source_info := "bim"]
from_bim[, bim_found := !is.na(bim_bp) & !is.na(bim_a1) & !is.na(bim_a2)]

fwrite(
  from_bim[, .(
    source_set, lead_id, source_info, bim_found,
    bim_chr, bim_bp, bim_a1, bim_a2,
    variant_type
  )],
  out_bim,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

query_vcf_window <- function(chr, pos, vcf_file, chr_mapper) {
  region_start <- max(1L, as.integer(pos) - window_bp)
  region_end <- max(region_start, as.integer(pos) + window_bp)
  query_chr <- chr_mapper(chr)
  region <- paste0(query_chr, ":", region_start, "-", region_end)
  lines <- tryCatch(
    system2(
      path.expand(bcftools_bin),
      bcftools_view_region_args(vcf_file, region),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
  if (!length(lines)) return(data.table())
  out <- fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    sep = "\t",
    select = 1:8,
    col.names = c("vcf_chr", "vcf_pos", "vcf_id", "vcf_ref", "vcf_alt", "vcf_qual", "vcf_filter", "vcf_info"),
    showProgress = FALSE
  )
  out[, vcf_chr := norm_chr(vcf_chr)]
  out[, vcf_pos := as.integer(vcf_pos)]
  out[, vcf_end := parse_info_value(vcf_info, "END")]
  out[, vcf_svtype := parse_info_value(vcf_info, "SVTYPE")]
  out[, vcf_svlen := parse_info_value(vcf_info, "SVLEN")]
  out[, c("vcf_qual", "vcf_filter", "vcf_info") := NULL]
  out[]
}

get_query_vcf <- function(source_set) {
  if (has_paragraph_vcf && identical(as.character(source_set), "paragraph")) {
    return(list(
      source = "paragraph_vcf",
      file = paragraph_vcf,
      mapper = to_paragraph_vcf_chr
    ))
  }
  list(
    source = "pav_sv_vcf",
    file = pav_sv_vcf,
    mapper = to_vcf_chr
  )
}

candidate_list <- list()
query_rows <- from_bim[bim_found == TRUE]
if (!nrow(query_rows)) {
  stop("No SV leads have usable chromosome, position, and allele information after BIM/fallback mapping.")
}
for (i in seq_len(nrow(query_rows))) {
  lead <- query_rows[i]
  query_vcf <- get_query_vcf(lead$source_set)
  cand <- query_vcf_window(lead$bim_chr, lead$bim_bp, query_vcf$file, query_vcf$mapper)
  if (!nrow(cand)) {
    candidate_list[[length(candidate_list) + 1L]] <- data.table(
      source_set = lead$source_set,
      lead_id = lead$lead_id,
      source_info = lead$source_info,
      vcf_source = query_vcf$source,
      vcf_file = path.expand(query_vcf$file),
      lead_chr = lead$bim_chr,
      lead_pos = lead$bim_bp,
      lead_a1 = lead$bim_a1,
      lead_a2 = lead$bim_a2,
      vcf_chr = NA_character_,
      vcf_pos = NA_integer_,
      vcf_id = NA_character_,
      vcf_ref = NA_character_,
      vcf_alt = NA_character_,
      vcf_end = NA_character_,
      vcf_svtype = NA_character_,
      vcf_svlen = NA_character_,
      distance_bp = NA_integer_,
      id_same = FALSE,
      ref_same_as_a1 = FALSE,
      ref_same_as_a2 = FALSE,
      alt_contains_a1 = FALSE,
      alt_contains_a2 = FALSE,
      exact_ref_alt = FALSE,
      fuzzy_ref_same = FALSE,
      match_type = "none",
      dosage_flip_for_a2 = NA,
      allele_source_orientation = "none"
    )
    next
  }

  cand[, source_set := lead$source_set]
  cand[, lead_id := lead$lead_id]
  cand[, source_info := lead$source_info]
  cand[, vcf_source := query_vcf$source]
  cand[, vcf_file := path.expand(query_vcf$file)]
  cand[, lead_chr := lead$bim_chr]
  cand[, lead_pos := lead$bim_bp]
  cand[, lead_a1 := lead$bim_a1]
  cand[, lead_a2 := lead$bim_a2]
  cand[, distance_bp := abs(vcf_pos - lead$bim_bp)]
  cand[, id_same := as.character(vcf_id) == as.character(lead$lead_id)]
  cand[, vcf_ref_cmp := toupper(vcf_ref)]
  cand[, vcf_alt_cmp := toupper(vcf_alt)]

  eval_orientation <- function(dt, a1, a2, label) {
    a1_cmp <- toupper(a1)
    a2_cmp <- toupper(a2)
    out <- copy(dt)
    out[, ref_same_as_a1 := vcf_ref_cmp == a1_cmp]
    out[, ref_same_as_a2 := vcf_ref_cmp == a2_cmp]
    out[, alt_contains_a1 := vapply(split_alt(vcf_alt_cmp), function(z) a1_cmp %in% z, logical(1L))]
    out[, alt_contains_a2 := vapply(split_alt(vcf_alt_cmp), function(z) a2_cmp %in% z, logical(1L))]

    # For a BIM allele pair, one allele must be VCF REF. The other allele is then
    # represented as VCF ALT. If the target is lead_a2, ref_same_as_a1 means VCF
    # ALT dosage already points to lead_a2; ref_same_as_a2 means dosage must be
    # flipped to keep lead_a2-oriented betas.
    out[, exact_ref_alt := distance_bp == 0L & (
      (ref_same_as_a1 & alt_contains_a2) |
        (ref_same_as_a2 & alt_contains_a1)
    )]
    out[, fuzzy_ref_same := distance_bp <= window_bp & (ref_same_as_a1 | ref_same_as_a2)]
    out[, match_type := fifelse(exact_ref_alt, "exact_ref_alt",
                                fifelse(fuzzy_ref_same, "fuzzy_ref_same", "none"))]
    out[, dosage_flip_for_a2 := fifelse(
      ref_same_as_a1 & alt_contains_a2,
      FALSE,
      fifelse(ref_same_as_a2 & alt_contains_a1, TRUE, NA)
    )]
    out[, allele_source_orientation := label]
    out
  }

  orientation_original <- eval_orientation(cand, lead$bim_a1, lead$bim_a2, "as_is")
  try_swapped <- identical(as.character(lead$source_set), "paragraph") &&
    identical(as.character(lead$source_info), "sig_table_fallback")
  if (try_swapped) {
    orientation_swapped <- eval_orientation(cand, lead$bim_a2, lead$bim_a1, "paragraph_a1_a2_swapped")
    cand <- rbindlist(list(orientation_original, orientation_swapped), fill = TRUE)
  } else {
    cand <- orientation_original
  }

  candidate_list[[length(candidate_list) + 1L]] <- cand[
    ,
    .(
      source_set, lead_id, source_info, vcf_source, vcf_file,
      lead_chr, lead_pos, lead_a1, lead_a2,
      vcf_chr, vcf_pos, vcf_id, vcf_ref, vcf_alt, vcf_end, vcf_svtype, vcf_svlen,
      distance_bp, id_same,
      ref_same_as_a1, ref_same_as_a2,
      alt_contains_a1, alt_contains_a2,
      exact_ref_alt, fuzzy_ref_same, match_type,
      dosage_flip_for_a2, allele_source_orientation
    )
  ]
}

candidates <- rbindlist(candidate_list, fill = TRUE)
setorder(candidates, source_set, lead_id, -exact_ref_alt, -fuzzy_ref_same, -id_same, distance_bp)
fwrite(candidates, out_candidates, sep = "\t", quote = FALSE, na = "NA")

best <- candidates[
  ,
  {
    cand <- .SD[!is.na(vcf_pos)]
    if (!nrow(cand)) {
      .(
        n_candidates = 0L,
        n_exact_ref_alt = 0L,
        n_fuzzy_ref_same = 0L,
        best_match_type = "none",
        best_vcf_source = NA_character_,
        best_vcf_file = NA_character_,
        best_vcf_chr = NA_character_,
        best_vcf_pos = NA_integer_,
        best_vcf_id = NA_character_,
        best_vcf_ref = NA_character_,
        best_vcf_alt = NA_character_,
        best_distance_bp = NA_integer_,
        best_id_same = FALSE,
        best_ref_same_as_a1 = FALSE,
        best_ref_same_as_a2 = FALSE,
        best_alt_contains_a1 = FALSE,
        best_alt_contains_a2 = FALSE,
        best_dosage_flip_for_a2 = NA,
        best_allele_source_orientation = "none"
      )
    } else {
      setorder(cand, -exact_ref_alt, -fuzzy_ref_same, -id_same, distance_bp)
      b <- cand[1L]
      .(
        n_candidates = nrow(cand),
        n_exact_ref_alt = sum(cand$exact_ref_alt == TRUE, na.rm = TRUE),
        n_fuzzy_ref_same = sum(cand$fuzzy_ref_same == TRUE, na.rm = TRUE),
        best_match_type = b$match_type,
        best_vcf_source = b$vcf_source,
        best_vcf_file = b$vcf_file,
        best_vcf_chr = b$vcf_chr,
        best_vcf_pos = b$vcf_pos,
        best_vcf_id = b$vcf_id,
        best_vcf_ref = b$vcf_ref,
        best_vcf_alt = b$vcf_alt,
        best_distance_bp = b$distance_bp,
        best_id_same = b$id_same,
        best_ref_same_as_a1 = b$ref_same_as_a1,
        best_ref_same_as_a2 = b$ref_same_as_a2,
        best_alt_contains_a1 = b$alt_contains_a1,
        best_alt_contains_a2 = b$alt_contains_a2,
        best_dosage_flip_for_a2 = b$dosage_flip_for_a2,
        best_allele_source_orientation = b$allele_source_orientation
      )
    }
  },
  by = .(source_set, lead_id, source_info, lead_chr, lead_pos, lead_a1, lead_a2)
]
setorder(best, source_set, lead_chr, lead_pos)
fwrite(best, out_best, sep = "\t", quote = FALSE, na = "NA")

summary_dt <- rbindlist(list(
  data.table(metric = "sig_sv_rows", value = nrow(sig_sv)),
  data.table(metric = "sig_sv_with_bim_or_fallback", value = sum(from_bim$bim_found == TRUE, na.rm = TRUE)),
  data.table(metric = "best_exact_ref_alt", value = sum(best$best_match_type == "exact_ref_alt", na.rm = TRUE)),
  data.table(metric = "best_fuzzy_ref_same", value = sum(best$best_match_type == "fuzzy_ref_same", na.rm = TRUE)),
  data.table(metric = "best_no_match", value = sum(best$best_match_type == "none", na.rm = TRUE)),
  data.table(metric = "window_bp", value = window_bp),
  data.table(metric = "vcf_contigs_detected_head", value = paste(head(vcf_chroms, 20), collapse = ",")),
  data.table(metric = "paragraph_vcf_contigs_detected_head", value = if (has_paragraph_vcf) paste(head(paragraph_chroms, 20), collapse = ",") else NA_character_),
  data.table(metric = "pav_sv_vcf", value = path.expand(pav_sv_vcf)),
  data.table(metric = "paragraph_vcf", value = if (has_paragraph_vcf) path.expand(paragraph_vcf) else NA_character_),
  data.table(metric = "out_bim", value = out_bim),
  data.table(metric = "out_candidates", value = out_candidates),
  data.table(metric = "out_best", value = out_best)
), fill = TRUE)
fwrite(summary_dt, out_summary, sep = "\t", quote = FALSE, na = "NA")

message("[DONE] BIM-derived lead info: ", out_bim)
message("[DONE] Candidate matches: ", out_candidates)
message("[DONE] Best matches: ", out_best)
message("[DONE] Summary: ", out_summary)
print(summary_dt)
