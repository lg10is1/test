#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(qqman)
})

# ============================================================
# Config
# ============================================================

workdir <- "/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD"
deepvariant_paragraph_dir <- "/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph"

outdir <- file.path(workdir, "LD_decay_public")
ids_dir <- file.path(outdir, "ids")
cmd_dir <- file.path(outdir, "cmd")
ld_dir  <- file.path(outdir, "ld")
maf_dir <- file.path(outdir, "maf01")
rdata_dir <- file.path(outdir, "rdata")
fig_dir <- file.path(outdir, "figures")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(ids_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cmd_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rdata_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
manual_exclusion_file <- file.path(script_dir, "manual_excluded_leads.tsv")

parse_seed_arg <- function(args, default = 1L) {
  seed <- as.integer(default)
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--seed") {
      i <- i + 1L
      if (i > length(args)) stop("--seed requires a value")
      seed <- suppressWarnings(as.integer(args[[i]]))
    } else {
      stop("Unknown argument: ", args[[i]])
    }
    i <- i + 1L
  }
  if (is.na(seed)) stop("--seed must be an integer")
  seed
}

plink_bin <- "plink"
max_plink_jobs <- 2L

set_names <- c("set00")
pangenie_maf <- 0.01
pangenie_maf_lists <- setNames(
  file.path(maf_dir, paste0(set_names, ".maf0.01.snplist")),
  set_names
)
deepvariant_paragraph_null_maf <- 0.02
deepvariant_paragraph_null_maf_list <- file.path(
  outdir,
  "maf02",
  "deepvariant_paragraph.maf0.02.snplist"
)

pangenie_base <- "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas"

clump_dir <- "/path/to/EOSCZ_PROJECT/figure_analysis/01.GWAS_figure.public/clumping_by_set_subtype"
canonical_lead_file <- file.path(
  outdir, "tables", "lead_sig_from_gwas.canonical_1000kb.tsv"
)
canonical_mapping_file <- file.path(
  outdir, "tables", "lead_sig_from_gwas.canonical_1000kb.mapping.tsv"
)

paragraph_bfile <- file.path(
  deepvariant_paragraph_dir,
  "chr_all2.strict_step2_genimi.common_samples.merged"
)

deepvariant_snv_clump_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/deepvar/04_clumping/SCZ.deepvar.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped"
)

paragraph_sv_clump_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/paragraph_test/04_clumping/SCZ.paragraph_test.pcsrc_deepvar_pc20_grm_deepvar_with_batch.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped"
)

pipeline_seed <- parse_seed_arg(commandArgs(trailingOnly = TRUE), default = 1L)
set.seed(pipeline_seed)

null_ratio <- 100
sig_p <- 5e-6

# LD command parameters
ld_window_kb <- 1000
ld_window_variant_count <- 999999
ld_window_r2 <- 0

# ============================================================
# Helper functions
# ============================================================

check_file <- function(x) {
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }
  invisible(TRUE)
}

check_bfile <- function(prefix) {
  required <- paste0(prefix, c(".bed", ".bim", ".fam"))
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop("Missing bfile components:\n", paste(missing, collapse = "\n"))
  }
  invisible(TRUE)
}

read_bim <- function(bfile_prefix) {
  check_bfile(bfile_prefix)
  bim_file <- paste0(bfile_prefix, ".bim")
  bim <- fread(bim_file, header = FALSE)
  if (ncol(bim) < 6) {
    stop("BIM file has fewer than 6 columns: ", bim_file)
  }
  colnames(bim)[1:6] <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
  bim[, CHR := as.character(CHR)]
  bim[, CHR := sub("^chr", "", CHR, ignore.case = TRUE)]
  bim <- bim[CHR %in% as.character(1:22)]
  bim[, SNP := as.character(SNP)]
  bim[, A1 := as.character(A1)]
  bim[, A2 := as.character(A2)]

  # Pangenie IDs do not encode SV length, so use allele length.
  # This is also valid for paragraph/deepvariant bims.
  bim[, allele_len_diff := abs(nchar(A1) - nchar(A2))]
  bim[, max_allele_len := pmax(nchar(A1), nchar(A2))]

  # For paragraph/deepvariant, use ID type only as extra information.
  bim[, id_has_sv_tag := grepl("(^|[-_])(INS|DEL|DUP|INV)([-_]|$)", SNP)]

  bim
}

read_clumped_ids <- function(file, id_col = "SNP") {
  check_file(file)
  x <- fread(file)
  if (!id_col %in% colnames(x)) {
    stop("Column ", id_col, " not found in clumped file: ", file)
  }
  ids <- unique(as.character(x[[id_col]]))
  ids <- ids[!is.na(ids) & ids != ""]
  ids
}

read_tsv_ids <- function(file, id_col = "ID") {
  check_file(file)
  x <- fread(file)
  if (!id_col %in% colnames(x)) {
    stop("Column ", id_col, " not found in file: ", file)
  }
  ids <- unique(as.character(x[[id_col]]))
  ids <- ids[!is.na(ids) & ids != ""]
  ids
}

read_canonical_ids <- function(canonical, source_name, variant_group) {
  variant_group <- match.arg(variant_group, c("SV", "SNV_INDEL"))
  rows <- canonical[
    source_set == source_name & canonical_variant_group == tolower(variant_group)
  ]
  ids <- unique(as.character(rows$lead_id))
  ids[!is.na(ids) & ids != ""]
}

read_manual_excluded_ids <- function(path, source_names, variant_groups = NULL) {
  if (!file.exists(path)) {
    message("[MANUAL EXCLUSION] No manual exclusion file found for LD null exclusion: ", path)
    return(character())
  }
  x <- fread(path, showProgress = FALSE)
  required <- c("source_set", "variant_type", "lead_id")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Manual exclusion file lacks: ", paste(missing, collapse = ", "))
  }
  x[, source_set := tolower(trimws(as.character(source_set)))]
  x[, variant_type := tolower(trimws(as.character(variant_type)))]
  x[, lead_id := trimws(as.character(lead_id))]
  source_names <- tolower(as.character(source_names))
  keep <- x$source_set %in% source_names
  if (!is.null(variant_groups)) {
    keep <- keep & x$variant_type %in% tolower(as.character(variant_groups))
  }
  ids <- unique(x[keep, lead_id])
  ids[!is.na(ids) & ids != ""]
}

read_gwas_sig_ids <- function(file, sig_p = 5e-6) {
  check_file(file)
  x <- fread(file)

  if (!"P" %in% colnames(x)) {
    stop("Column P not found in GWAS file: ", file)
  }

  candidate_cols <- c(
    "SNP",
    "ID",
    "MarkerName",
    "MARKERNAME",
    "variant_id",
    "VARIANT_ID",
    "rsid",
    "RSID"
  )

  id_cols <- intersect(candidate_cols, colnames(x))

  if (length(id_cols) == 0) {
    stop(
      "No candidate ID column found in GWAS file: ", file,
      "\nColumns are: ", paste(colnames(x), collapse = ", ")
    )
  }

  id_col <- id_cols[1]
  x[, P := as.numeric(P)]
  sig <- x[!is.na(P) & P < sig_p]

  ids <- unique(as.character(sig[[id_col]]))
  ids <- ids[!is.na(ids) & ids != ""]

  message("[GWAS SIG] ", basename(file), " | P<", sig_p,
          " | id_col=", id_col, " | n_ids=", length(ids))

  ids
}

write_ids <- function(ids, file) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & ids != ""]
  fwrite(
    data.table(SNP = ids),
    file = file,
    quote = FALSE,
    sep = "\t",
    col.names = FALSE
  )
  invisible(file)
}

variant_pool <- function(bim, variant_class) {
  variant_class <- match.arg(variant_class, c("SV", "SNV_INDEL"))

  if (variant_class == "SV") {
    # Main rule: max allele length >= 50.
    # Extra rule: paragraph/deepvariant IDs containing INS/DEL/DUP/INV are SV-like.
    # Because some INS/DEL IDs already encode length, but allele length is the safer general rule.
    out <- bim %>%
      filter(max_allele_len >= 50 | id_has_sv_tag)
  } else {
    # SNV + small INDEL.
    # For Pangenie, SNV has max allele length 1 and small INDEL has max allele length < 50.
    # For paragraph/deepvariant, this captures SNVs and small indels.
    out <- bim %>%
      filter(max_allele_len < 50 & !id_has_sv_tag)
  }

  out
}

standardize_sig_ids_to_bim <- function(sig_ids, bim, label) {
  sig_ids <- unique(as.character(sig_ids))
  sig_ids <- sig_ids[!is.na(sig_ids) & sig_ids != ""]

  found <- sig_ids[sig_ids %in% bim$SNP]
  missing <- setdiff(sig_ids, found)

  message("[", label, "] input IDs: ", length(sig_ids),
          "; found in BIM: ", length(found),
          "; missing: ", length(missing))

  if (length(missing) > 0) {
    missing_file <- file.path(ids_dir, paste0(label, ".missing_from_bim.ids"))
    write_ids(missing, missing_file)
    message("[", label, "] missing IDs written to: ", missing_file)
  }

  unique(found)
}

make_null_by_chr <- function(
  bim,
  sig_ids,
  exclude_ids,
  variant_class,
  label,
  seed = 1,
  null_ratio = 10,
  pool_bim = NULL
) {
  set.seed(seed)

  if (is.null(pool_bim)) {
    pool_bim <- bim
  }

  sig_ids <- unique(as.character(sig_ids))
  exclude_ids <- unique(as.character(exclude_ids))
  null_ratio <- as.integer(null_ratio)

  if (is.na(null_ratio) || null_ratio < 1) {
    stop("[", label, "] null_ratio must be an integer >= 1.")
  }

  sig_df <- bim %>%
    filter(SNP %in% sig_ids)

  if (nrow(sig_df) == 0) {
    warning("[", label, "] no significant IDs found in BIM; return empty null.")
    return(character(0))
  }

  chr_n <- sig_df %>%
    count(CHR, name = "n_sig")

  chr_n <- chr_n %>%
    mutate(n_needed = n_sig * null_ratio)

  pool <- variant_pool(pool_bim, variant_class) %>%
    filter(!SNP %in% exclude_ids)

  if (nrow(pool) == 0) {
    stop("[", label, "] null pool is empty.")
  }

  selected <- character(0)
  shortage <- 0

  for (i in seq_len(nrow(chr_n))) {
    chr_i <- chr_n$CHR[i]
    n_i <- chr_n$n_needed[i]

    pool_i <- pool %>%
      filter(CHR == chr_i, !SNP %in% selected)

    if (nrow(pool_i) >= n_i) {
      selected_i <- sample(pool_i$SNP, n_i, replace = FALSE)
    } else {
      selected_i <- pool_i$SNP
      shortage <- shortage + (n_i - length(selected_i))
      warning(
        "[", label, "] CHR ", chr_i,
        " has fewer null candidates than needed. Needed: ", n_i,
        "; available: ", length(selected_i),
        ". Shortage will be sampled from genome-wide pool."
      )
    }

    selected <- c(selected, selected_i)
  }

  if (shortage > 0) {
    remaining_pool <- pool %>%
      filter(!SNP %in% selected)

    if (nrow(remaining_pool) < shortage) {
      stop(
        "[", label, "] not enough genome-wide null candidates. Shortage: ",
        shortage, "; remaining pool: ", nrow(remaining_pool)
      )
    }

    selected_extra <- sample(remaining_pool$SNP, shortage, replace = FALSE)
    selected <- c(selected, selected_extra)
  }

  selected <- unique(selected)
  target_n <- sum(chr_n$n_needed)

  if (length(selected) != target_n) {
    warning(
      "[", label, "] null count differs from target. sig=",
      length(sig_ids), "; null_ratio=", null_ratio,
      "; target_null=", target_n, "; actual_null=", length(selected)
    )
  }

  selected
}

add_job <- function(
  job_list,
  source_id,
  bfile,
  plot_group,
  id_file,
  out_prefix,
  maf_filter = NA_real_
) {
  if (!file.exists(id_file)) {
    stop("ID file does not exist: ", id_file)
  }

  n_ids <- nrow(fread(id_file, header = FALSE))
  if (n_ids == 0) {
    warning("[", source_id, " / ", plot_group, "] empty ID file; skip PLINK job.")
    return(job_list)
  }

  job <- data.table(
    source_id = source_id,
    bfile = bfile,
    plot_group = plot_group,
    id_file = id_file,
    out_prefix = out_prefix,
    n_index = n_ids,
    maf_filter = as.numeric(maf_filter)
  )

  append(job_list, list(job))
}

make_plink_cmd <- function(plink_bin, bfile, id_file, out_prefix, maf_filter = NA_real_) {
  parts <- c(
    shQuote(plink_bin),
    "--bfile", shQuote(bfile),
    "--threads", "1",
    "--ld-snp-list", shQuote(id_file),
    "--r2", "gz",
    "--ld-window-kb", ld_window_kb,
    "--ld-window", ld_window_variant_count,
    "--ld-window-r2", ld_window_r2
  )
  if (is.finite(maf_filter)) {
    parts <- c(parts, "--maf", format(maf_filter, scientific = FALSE, trim = TRUE))
  }
  parts <- c(parts, "--out", shQuote(out_prefix))
  paste(parts, collapse = " ")
}

# ============================================================
# Main: prepare sig/null IDs and PLINK jobs
# ============================================================

check_file(canonical_lead_file)
check_file(canonical_mapping_file)
canonical_leads <- fread(canonical_lead_file, showProgress = FALSE)
canonical_mapping <- fread(canonical_mapping_file, showProgress = FALSE)
required_canonical <- c("source_set", "lead_id", "canonical_variant_group")
if (length(setdiff(required_canonical, names(canonical_leads)))) {
  stop("Canonical lead table lacks: ", paste(setdiff(required_canonical, names(canonical_leads)), collapse = ", "))
}
if (length(setdiff(c("source_set", "lead_id"), names(canonical_mapping)))) {
  stop("Canonical mapping table lacks source_set or lead_id")
}
canonical_leads[, source_set := as.character(source_set)]
canonical_leads[, lead_id := as.character(lead_id)]
canonical_leads[, canonical_variant_group := tolower(as.character(canonical_variant_group))]
canonical_mapping[, source_set := as.character(source_set)]
canonical_mapping[, lead_id := as.character(lead_id)]
n_canonical_sv <- canonical_leads[canonical_variant_group == "sv", .N]
n_canonical_snv_indel <- canonical_leads[canonical_variant_group == "snv_indel", .N]
if (n_canonical_sv < 1L || n_canonical_snv_indel < 1L) {
  stop("Canonical lead table must contain at least one SV and one SNV/indel lead")
}
message("[CANONICAL] SV=", n_canonical_sv, "; SNV/indel=", n_canonical_snv_indel)

jobs <- list()
id_summary <- list()

for (set_name in set_names) {
  message("\n========== Processing ", set_name, " ==========")

  bfile <- file.path(pangenie_base, set_name, "NGS.QCsite.QCind")
  bim <- read_bim(bfile)
  maf_file <- pangenie_maf_lists[[set_name]]
  check_file(maf_file)
  maf_ids <- unique(as.character(fread(maf_file, header = FALSE)[[1L]]))
  maf_ids <- maf_ids[!is.na(maf_ids) & maf_ids != ""]
  before_maf_n <- nrow(bim)
  bim <- bim[SNP %in% maf_ids]
  if (nrow(bim) == 0L) {
    stop("[", set_name, "] no variants remain after MAF >= ", pangenie_maf)
  }
  message(
    "[", set_name, "] MAF >= ", pangenie_maf,
    " BIM variants: ", nrow(bim), " / ", before_maf_n
  )

  snv_sig_raw <- read_canonical_ids(canonical_leads, set_name, "SNV_INDEL")
  sv_sig_raw  <- read_canonical_ids(canonical_leads, set_name, "SV")

  snv_sig <- standardize_sig_ids_to_bim(
    snv_sig_raw, bim, paste0(set_name, ".SNV_INDEL_sig")
  )
  sv_sig <- standardize_sig_ids_to_bim(
    sv_sig_raw, bim, paste0(set_name, ".SV_sig")
  )

  all_sig_in_this_bim <- unique(c(
    snv_sig,
    sv_sig,
    canonical_mapping[source_set == set_name & lead_id %in% bim$SNP, lead_id],
    read_manual_excluded_ids(manual_exclusion_file, set_name)
  ))
  all_sig_in_this_bim <- unique(all_sig_in_this_bim[all_sig_in_this_bim %in% bim$SNP])

  snv_null <- make_null_by_chr(
    bim = bim,
    sig_ids = snv_sig,
    exclude_ids = all_sig_in_this_bim,
    variant_class = "SNV_INDEL",
    label = paste0(set_name, ".SNV_INDEL_null"),
    seed = pipeline_seed,
    null_ratio = null_ratio
  )

  sv_null <- make_null_by_chr(
    bim = bim,
    sig_ids = sv_sig,
    exclude_ids = all_sig_in_this_bim,
    variant_class = "SV",
    label = paste0(set_name, ".SV_null"),
    seed = pipeline_seed,
    null_ratio = null_ratio
  )

  id_files <- list(
    SNV_INDEL_sig  = file.path(ids_dir, paste0(set_name, ".SNV_INDEL_sig.ids")),
    SV_sig         = file.path(ids_dir, paste0(set_name, ".SV_sig.ids")),
    SNV_INDEL_null = file.path(ids_dir, paste0(set_name, ".SNV_INDEL_null.ids")),
    SV_null        = file.path(ids_dir, paste0(set_name, ".SV_null.ids"))
  )

  write_ids(snv_sig,  id_files$SNV_INDEL_sig)
  write_ids(sv_sig,   id_files$SV_sig)
  write_ids(snv_null, id_files$SNV_INDEL_null)
  write_ids(sv_null,  id_files$SV_null)

  for (pg in names(id_files)) {
    out_prefix <- file.path(ld_dir, paste0(set_name, ".", pg))
    jobs <- add_job(
      job_list = jobs,
      source_id = set_name,
      bfile = bfile,
      plot_group = pg,
      id_file = id_files[[pg]],
      out_prefix = out_prefix,
      maf_filter = pangenie_maf
    )
  }

  id_summary[[set_name]] <- data.table(
    source_id = set_name,
    group = c("SNV_INDEL_sig", "SV_sig", "SNV_INDEL_null", "SV_null"),
    n = c(length(snv_sig), length(sv_sig), length(snv_null), length(sv_null))
  )
}

# ============================================================
# Deepvariant / paragraph source
# ============================================================

message("\n========== Processing deepvariant_paragraph ==========")

paragraph_bim <- read_bim(paragraph_bfile)
check_file(deepvariant_paragraph_null_maf_list)
paragraph_null_maf_ids <- unique(as.character(fread(deepvariant_paragraph_null_maf_list, header = FALSE)[[1L]]))
paragraph_null_maf_ids <- paragraph_null_maf_ids[!is.na(paragraph_null_maf_ids) & paragraph_null_maf_ids != ""]
paragraph_null_bim <- paragraph_bim[SNP %in% paragraph_null_maf_ids]
if (nrow(paragraph_null_bim) == 0L) {
  stop("[deepvariant_paragraph] no variants remain in null pool after MAF >= ", deepvariant_paragraph_null_maf)
}
message(
  "[deepvariant_paragraph] null pool MAF >= ", deepvariant_paragraph_null_maf,
  " BIM variants: ", nrow(paragraph_null_bim), " / ", nrow(paragraph_bim),
  "; significant IDs are not filtered by this threshold"
)

paragraph_sv_sig_raw <- read_canonical_ids(canonical_leads, "paragraph", "SV")
deepvariant_snv_sig_raw <- read_canonical_ids(canonical_leads, "deepvariant", "SNV_INDEL")

paragraph_sv_sig <- standardize_sig_ids_to_bim(
  paragraph_sv_sig_raw,
  paragraph_bim,
  "deepvariant_paragraph.SV_sig"
)

deepvariant_snv_sig <- standardize_sig_ids_to_bim(
  deepvariant_snv_sig_raw,
  paragraph_bim,
  "deepvariant_paragraph.SNV_INDEL_sig"
)

all_sig_paragraph <- unique(c(
  paragraph_sv_sig,
  deepvariant_snv_sig,
  canonical_mapping[source_set %in% c("paragraph", "deepvariant") & lead_id %in% paragraph_bim$SNP, lead_id],
  read_manual_excluded_ids(manual_exclusion_file, c("paragraph", "deepvariant"))
))
all_sig_paragraph <- unique(all_sig_paragraph[all_sig_paragraph %in% paragraph_bim$SNP])

paragraph_sv_null <- make_null_by_chr(
  bim = paragraph_bim,
  sig_ids = paragraph_sv_sig,
  exclude_ids = all_sig_paragraph,
  variant_class = "SV",
  label = "deepvariant_paragraph.SV_null",
  seed = pipeline_seed,
  null_ratio = null_ratio,
  pool_bim = paragraph_null_bim
)

deepvariant_snv_null <- make_null_by_chr(
  bim = paragraph_bim,
  sig_ids = deepvariant_snv_sig,
  exclude_ids = all_sig_paragraph,
  variant_class = "SNV_INDEL",
  label = "deepvariant_paragraph.SNV_INDEL_null",
  seed = pipeline_seed,
  null_ratio = null_ratio,
  pool_bim = paragraph_null_bim
)

paragraph_id_files <- list(
  SNV_INDEL_sig  = file.path(ids_dir, "deepvariant_paragraph.SNV_INDEL_sig.ids"),
  SV_sig         = file.path(ids_dir, "deepvariant_paragraph.SV_sig.ids"),
  SNV_INDEL_null = file.path(ids_dir, "deepvariant_paragraph.SNV_INDEL_null.ids"),
  SV_null        = file.path(ids_dir, "deepvariant_paragraph.SV_null.ids")
)

write_ids(deepvariant_snv_sig,  paragraph_id_files$SNV_INDEL_sig)
write_ids(paragraph_sv_sig,     paragraph_id_files$SV_sig)
write_ids(deepvariant_snv_null, paragraph_id_files$SNV_INDEL_null)
write_ids(paragraph_sv_null,    paragraph_id_files$SV_null)

for (pg in names(paragraph_id_files)) {
  out_prefix <- file.path(ld_dir, paste0("deepvariant_paragraph.", pg))
  jobs <- add_job(
    job_list = jobs,
    source_id = "deepvariant_paragraph",
    bfile = paragraph_bfile,
    plot_group = pg,
    id_file = paragraph_id_files[[pg]],
    out_prefix = out_prefix
  )
}

id_summary[["deepvariant_paragraph"]] <- data.table(
  source_id = "deepvariant_paragraph",
  group = c("SNV_INDEL_sig", "SV_sig", "SNV_INDEL_null", "SV_null"),
  n = c(
    length(deepvariant_snv_sig),
    length(paragraph_sv_sig),
    length(deepvariant_snv_null),
    length(paragraph_sv_null)
  )
)

# ============================================================
# Write metadata and PLINK command script
# ============================================================

jobs_dt <- rbindlist(jobs, fill = TRUE)
id_summary_dt <- rbindlist(id_summary, fill = TRUE)

metadata_file <- file.path(outdir, "ld_jobs.metadata.tsv")
id_summary_file <- file.path(outdir, "id_summary.tsv")
null_config_file <- file.path(outdir, "null_matching_config.tsv")

fwrite(jobs_dt, metadata_file, sep = "\t", quote = FALSE)
fwrite(id_summary_dt, id_summary_file, sep = "\t", quote = FALSE)

null_config <- data.table(
  parameter = c(
    "null_ratio",
    "random_seed",
    "matching_level",
    "variant_class_rule",
    "exclude_ids",
    "sampling",
    "pangenie_maf_filter",
    "pangenie_maf_lists",
    "deepvariant_paragraph_null_maf_filter",
    "deepvariant_paragraph_null_maf_list",
    "deepvariant_paragraph_sig_maf_filter",
    "ld_execution"
  ),
  value = c(
    as.character(null_ratio),
    as.character(pipeline_seed),
    "chromosome + variant_class",
    "SV: max_allele_len >= 50 or INS/DEL/DUP/INV ID tag; SNV_INDEL: max_allele_len < 50 and no SV ID tag",
    "all source-level significant IDs in the same BIM/source are excluded from the null pool, including signals merged out of the canonical index set",
    "without replacement within chromosome; chromosome shortage supplemented from genome-wide remaining pool",
    as.character(pangenie_maf),
    paste(pangenie_maf_lists, collapse = ";"),
    as.character(deepvariant_paragraph_null_maf),
    deepvariant_paragraph_null_maf_list,
    "none; significant index IDs are retained if present in the deepvariant_paragraph BIM",
    paste0("maximum ", max_plink_jobs, " concurrent PLINK LD jobs; --threads 1 per job; completed gzip outputs are resumable")
  )
)

fwrite(null_config, null_config_file, sep = "\t", quote = FALSE)

cmd_file <- file.path(cmd_dir, "run_plink_ld.sh")

cmds <- vapply(
  seq_len(nrow(jobs_dt)),
  function(i) {
    make_plink_cmd(
      plink_bin = plink_bin,
      bfile = jobs_dt$bfile[i],
      id_file = jobs_dt$id_file[i],
      out_prefix = jobs_dt$out_prefix[i],
      maf_filter = jobs_dt$maf_filter[i]
    )
  },
  character(1)
)

labels <- paste(jobs_dt$source_id, jobs_dt$plot_group, sep = ".")
ld_files <- paste0(jobs_dt$out_prefix, ".ld.gz")
task_lines <- paste0(
  "  ", shQuote(labels), " ", shQuote(ld_files), " ", shQuote(jobs_dt$id_file), " ", shQuote(cmds), " \\"
)
verify_lines <- paste0("  ", shQuote(ld_files), " \\")

cmd_lines <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  paste0("mkdir -p ", shQuote(ld_dir)),
  "",
  sprintf("MAX_PLINK_JOBS=%d", max_plink_jobs),
  "ld_complete() {",
  "  local ld_file=\"$1\"",
  "  local id_file=\"$2\"",
  "  [[ -s \"$ld_file\" ]] && [[ \"$ld_file\" -nt \"$id_file\" ]] && gzip -t \"$ld_file\" 2>/dev/null",
  "}",
  "run_one() {",
  "  local label=\"$1\"",
  "  local ld_file=\"$2\"",
  "  local id_file=\"$3\"",
  "  local cmd=\"$4\"",
  "  local status",
  "  if ld_complete \"$ld_file\" \"$id_file\"; then",
  "    printf '[SKIP LD] %s | %s\\n' \"$label\" \"$ld_file\"",
  "    return 0",
  "  fi",
  "  printf '[LAUNCH LD] %s\\n' \"$label\"",
  "  if bash -c \"$cmd\"; then",
  "    if ld_complete \"$ld_file\" \"$id_file\"; then",
  "      printf '[DONE LD] %s\\n' \"$label\"",
  "      return 0",
  "    fi",
  "    printf '[ERROR LD] PLINK returned zero but gzip output is missing/invalid: %s\\n' \"$label\" >&2",
  "    return 1",
  "  else",
  "    status=$?",
  "    printf '[ERROR LD] PLINK exit %s: %s\\n' \"$status\" \"$label\" >&2",
  "    return 1",
  "  fi",
  "}",
  "export -f ld_complete run_one",
  "",
  "xargs_status=0",
  "printf '%s\\0%s\\0%s\\0%s\\0' \\",
  task_lines,
  "  | xargs -0 -n 4 -P \"$MAX_PLINK_JOBS\" bash -c 'run_one \"$1\" \"$2\" \"$3\" \"$4\"' _ || xargs_status=$?",
  "",
  "missing=0",
  "for ld_file in \\",
  verify_lines,
  "; do",
  "  if [[ ! -s \"$ld_file\" ]] || ! gzip -t \"$ld_file\" 2>/dev/null; then",
  "    printf '[MISSING LD] %s\\n' \"$ld_file\" >&2",
  "    ((missing += 1))",
  "  fi",
  "done",
  "if (( missing > 0 )); then",
  "  printf '[ERROR] %s LD job(s) remain missing/invalid; xargs status=%s\\n' \"$missing\" \"$xargs_status\" >&2",
  "  exit 1",
  "fi",
  "echo '[DONE] All PLINK LD-decay jobs finished.'"
)

writeLines(cmd_lines, cmd_file)
Sys.chmod(cmd_file, mode = "0755")

message("\nDone.")
message("Metadata: ", metadata_file)
message("ID summary: ", id_summary_file)
message("Null matching config: ", null_config_file)
message("PLINK command script: ", cmd_file)
message("\nNext step:")
message("bash ", cmd_file)
