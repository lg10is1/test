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
rdata_dir <- file.path(outdir, "rdata")
table_dir <- file.path(outdir, "tables")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(rdata_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
manual_exclusion_file <- file.path(script_dir, "manual_excluded_leads.tsv")

set_names <- c("set00")
sig_p <- 5e-6

pangenie_base <- "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas"

clump_dir <- "/path/to/EOSCZ_PROJECT/Figure3/01.GWAS_figure.version20/clumping_by_set_subtype"

read_manual_exclusions <- function(path) {
  if (!file.exists(path)) {
    message("[MANUAL EXCLUSION] No manual exclusion file found: ", path)
    return(data.table())
  }
  dt <- fread(path, showProgress = FALSE)
  required <- c("source_set", "variant_type", "lead_id")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop("Manual exclusion file missing columns: ", paste(missing, collapse = ", "))
  }
  dt <- dt[, .(
    source_set = tolower(trimws(as.character(source_set))),
    variant_type = tolower(trimws(as.character(variant_type))),
    lead_id = trimws(as.character(lead_id)),
    reason = if ("reason" %in% names(dt)) as.character(reason) else NA_character_
  )]
  dt <- dt[source_set != "" & variant_type != "" & lead_id != ""]
  unique(dt, by = c("source_set", "variant_type", "lead_id"))
}

apply_manual_exclusion_to_leads <- function(dt, exclusions, audit_file) {
  if (!nrow(exclusions) || !nrow(dt)) {
    fwrite(data.table(status = "NOOP", removed = 0L), audit_file, sep = "\t", quote = FALSE, na = "NA")
    return(dt)
  }
  required <- c("source_set", "variant_type", "lead_id")
  missing <- setdiff(required, names(dt))
  if (length(missing)) stop("Lead table missing columns for manual exclusion: ", paste(missing, collapse = ", "))

  excl_keys <- paste(exclusions$source_set, exclusions$variant_type, exclusions$lead_id, sep = "\r")
  lead_keys <- paste(
    tolower(trimws(as.character(dt$source_set))),
    tolower(trimws(as.character(dt$variant_type))),
    trimws(as.character(dt$lead_id)),
    sep = "\r"
  )
  remove <- lead_keys %in% excl_keys

  removed <- copy(dt[remove])
  if (nrow(removed)) {
    removed[, manual_exclusion_key := lead_keys[remove]]
    removed <- merge(
      removed,
      exclusions[, .(
        manual_exclusion_key = paste(source_set, variant_type, lead_id, sep = "\r"),
        manual_exclusion_reason = reason
      )],
      by = "manual_exclusion_key",
      all.x = TRUE
    )
  }

  audit <- data.table(
    status = if (nrow(removed)) "FILTERED" else "UNCHANGED",
    rows_before = nrow(dt),
    rows_after = nrow(dt) - nrow(removed),
    removed = nrow(removed),
    unique_removed_leads = if (nrow(removed)) uniqueN(removed$lead_id) else 0L
  )
  fwrite(audit, audit_file, sep = "\t", quote = FALSE, na = "NA")

  removed_file <- sub("\\.tsv$", ".removed.tsv", audit_file)
  fwrite(removed, removed_file, sep = "\t", quote = FALSE, na = "NA")

  message(
    "[MANUAL EXCLUSION] removed ", nrow(removed),
    " rows / ", audit$unique_removed_leads,
    " unique leads after extracting all sources. Audit: ", audit_file
  )
  dt[!remove]
}

# ============================================================
# Pangenie GWAS summary files
# ============================================================

manual_gwas_files <- list(
  ## Public example: set00 only. Add set01/set02 entries if available.
  "set00__sv" = "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/gwas/SCZ.mlm.ngspc.fastGWA",
  "set00__snv_indel" = "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/gwas/SCZ.mlm.ngspc.fastGWA"
)

# ============================================================
# Deepvariant / paragraph input
# ============================================================

deepvariant_snv_clump_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/deepvar/04_clumping/SCZ.deepvar.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped"
)

deepvariant_snv_gwas_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/deepvar/03_gwas/SCZ.deepvar.mlm.geno0.1.maf0.01.fastGWA"
)

paragraph_sv_clump_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/paragraph_test/04_clumping/SCZ.paragraph_test.pcsrc_deepvar_pc20_grm_deepvar_with_batch.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped"
)

paragraph_sv_gwas_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/paragraph_test/03_gwas/SCZ.paragraph_test.pcsrc_deepvar_pc20_grm_deepvar_with_batch.mlm.geno0.1.maf0.01.fastGWA"
)

# ============================================================
# Helper functions
# ============================================================

check_file <- function(x) {
  if (is.na(x) || x == "") {
    stop("File path is NA or empty.")
  }

  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }

  if (dir.exists(x)) {
    stop("Path is a directory, not a file: ", x)
  }

  invisible(TRUE)
}

first_existing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  hit <- x[file.exists(x) & !dir.exists(x)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_clumped_ids <- function(file, id_col = "SNP") {
  check_file(file)
  x <- fread(file)

  if (!id_col %in% colnames(x)) {
    stop(
      "Column ", id_col, " not found in clumped file: ", file,
      "\nColumns are: ", paste(colnames(x), collapse = ", ")
    )
  }

  chr_cols <- intersect(c("CHR", "#CHROM", "CHROM"), colnames(x))
  if (!length(chr_cols)) {
    stop("Clumped file lacks a chromosome column required for chr1-22 filtering: ", file)
  }
  chr_value <- sub("^chr", "", as.character(x[[chr_cols[1L]]]), ignore.case = TRUE)
  x <- x[chr_value %in% as.character(1:22)]

  ids <- unique(as.character(x[[id_col]]))
  ids <- ids[!is.na(ids) & ids != ""]

  message("[CLUMP] ", basename(file), " | n_ids=", length(ids))

  ids
}

read_tsv_ids <- function(file, id_col = "ID") {
  check_file(file)
  x <- fread(file)

  if (!id_col %in% colnames(x)) {
    stop(
      "Column ", id_col, " not found in file: ", file,
      "\nColumns are: ", paste(colnames(x), collapse = ", ")
    )
  }

  ids <- unique(as.character(x[[id_col]]))
  ids <- ids[!is.na(ids) & ids != ""]

  message("[TSV] ", basename(file), " | n_ids=", length(ids))

  ids
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

resolve_pangenie_gwas_file <- function(set_name, variant_type) {
  key <- paste(set_name, variant_type, sep = "__")

  manual <- manual_gwas_files[[key]]

  if (!is.null(manual) && !is.na(manual) && manual != "") {
    check_file(manual)
    return(manual)
  }

  stop(
    "Cannot resolve Pangenie GWAS summary file for source_set=", set_name,
    ", variant_type=", variant_type,
    "\nPlease fill manual_gwas_files[['", key, "']] in the config block."
  )
}

standardize_common_colnames <- function(x) {
  #                      ?
  if ("#CHROM" %in% colnames(x) && !"CHR" %in% colnames(x)) {
    setnames(x, "#CHROM", "CHR")
  }

  if ("CHROM" %in% colnames(x) && !"CHR" %in% colnames(x)) {
    setnames(x, "CHROM", "CHR")
  }

  if ("BP" %in% colnames(x) && !"POS" %in% colnames(x)) {
    x[, POS := BP]
  }

  x
}

find_id_col_in_gwas <- function(gwas, lead_ids) {
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

  candidate_cols <- intersect(candidate_cols, colnames(gwas))

  if (length(candidate_cols) == 0) {
    stop(
      "No candidate ID column found in GWAS summary.\n",
      "Columns are: ", paste(colnames(gwas), collapse = ", ")
    )
  }

  overlap_n <- sapply(candidate_cols, function(cc) {
    sum(as.character(gwas[[cc]]) %in% lead_ids)
  })

  best_col <- candidate_cols[which.max(overlap_n)]
  best_n <- max(overlap_n)

  message("[GWAS] candidate ID column overlap:")
  for (cc in candidate_cols) {
    message("  - ", cc, ": ", overlap_n[[cc]])
  }

  if (best_n == 0) {
    stop(
      "No lead IDs matched any candidate ID column in GWAS summary.\n",
      "Candidate columns checked: ", paste(candidate_cols, collapse = ", "),
      "\nExample lead IDs: ", paste(head(lead_ids, 10), collapse = ", ")
    )
  }

  best_col
}

extract_leads_from_gwas <- function(
  lead_ids,
  source_set,
  variant_type,
  lead_source,
  lead_source_file,
  gwas_file
) {
  lead_ids <- unique(as.character(lead_ids))
  lead_ids <- lead_ids[!is.na(lead_ids) & lead_ids != ""]

  check_file(lead_source_file)
  check_file(gwas_file)

  message("\n============================================================")
  message("[GWAS] source_set=", source_set,
          " | variant_type=", variant_type,
          " | n_leads=", length(lead_ids))
  message("[GWAS] lead source: ", lead_source_file)
  message("[GWAS] summary file: ", gwas_file)

  gwas <- fread(gwas_file)
  gwas <- standardize_common_colnames(gwas)
  if (!"CHR" %in% colnames(gwas)) {
    stop("GWAS file lacks CHR after column standardization: ", gwas_file)
  }
  gwas[, CHR := sub("^chr", "", as.character(CHR), ignore.case = TRUE)]
  gwas <- gwas[CHR %in% as.character(1:22)]

  id_col <- find_id_col_in_gwas(gwas, lead_ids)

  message("[GWAS] using ID column: ", id_col)

  gwas[, gwas_id_for_match := as.character(get(id_col))]
  gwas_sub <- gwas[gwas_id_for_match %in% lead_ids]

  if (nrow(gwas_sub) > 0) {
    gwas_sub[, duplicate_index_in_gwas := seq_len(.N), by = gwas_id_for_match]
  }

  lead_dt <- data.table(
    lead_id = lead_ids,
    source_set = source_set,
    variant_type = variant_type,
    lead_source = lead_source,
    lead_source_file = lead_source_file,
    gwas_file = gwas_file
  )

  if (nrow(gwas_sub) == 0) {
    out <- copy(lead_dt)
    out[, found_in_gwas := FALSE]
    return(out)
  }

  setnames(gwas_sub, "gwas_id_for_match", "lead_id")
  gwas_sub[, found_in_gwas := TRUE]

  out <- merge(
    lead_dt,
    gwas_sub,
    by = "lead_id",
    all.x = TRUE,
    allow.cartesian = TRUE
  )

  out[is.na(found_in_gwas), found_in_gwas := FALSE]

  front_cols <- c(
    "source_set",
    "variant_type",
    "lead_id",
    "lead_source",
    "lead_source_file",
    "gwas_file",
    "found_in_gwas"
  )

  other_cols <- setdiff(colnames(out), front_cols)
  setcolorder(out, c(front_cols, other_cols))

  n_found <- uniqueN(out[found_in_gwas == TRUE]$lead_id)
  n_missing <- length(setdiff(lead_ids, out[found_in_gwas == TRUE]$lead_id))

  message("[GWAS] found lead IDs: ", n_found, " / ", length(lead_ids),
          "; missing: ", n_missing)

  out
}

# ============================================================
# Build lead specs
# ============================================================

lead_specs <- list()

# -----------------------------
# Pangenie set00 / set01 / set02
# -----------------------------

for (set_name in set_names) {
  snv_clump_file <- file.path(
    clump_dir,
    paste0(set_name, ".SNV_INDEL.clump_p1_5e-06.r2_0.01.kb_1000.clumped")
  )

  sv_clump_file <- file.path(
    clump_dir,
    paste0(set_name, ".SV.clump_p1_5e-06.r2_0.01.kb_1000.clumped")
  )

  snv_ids <- read_clumped_ids(snv_clump_file, id_col = "SNP")
  sv_ids <- read_clumped_ids(sv_clump_file, id_col = "SNP")

  lead_specs[[length(lead_specs) + 1]] <- list(
    source_set = set_name,
    variant_type = "snv_indel",
    lead_source = "clumped",
    lead_source_file = snv_clump_file,
    lead_ids = snv_ids,
    gwas_file = resolve_pangenie_gwas_file(set_name, "snv_indel")
  )

  lead_specs[[length(lead_specs) + 1]] <- list(
    source_set = set_name,
    variant_type = "sv",
    lead_source = "clumped",
    lead_source_file = sv_clump_file,
    lead_ids = sv_ids,
    gwas_file = resolve_pangenie_gwas_file(set_name, "sv")
  )
}

# -----------------------------
# Paragraph SV
# lead IDs from the new source-specific Paragraph clumping
# GWAS summary from the matching source-specific Paragraph MLM
# -----------------------------

paragraph_sv_ids <- read_clumped_ids(paragraph_sv_clump_file, id_col = "SNP")

lead_specs[[length(lead_specs) + 1]] <- list(
  source_set = "paragraph",
  variant_type = "sv",
  lead_source = "clumped",
  lead_source_file = paragraph_sv_clump_file,
  lead_ids = paragraph_sv_ids,
  gwas_file = paragraph_sv_gwas_file
)

# -----------------------------
# Deepvariant SNV/INDEL
# lead IDs from deepvariant clumped
# GWAS summary from deepvariant fastGWA
# -----------------------------

deepvariant_snv_ids <- read_clumped_ids(deepvariant_snv_clump_file, id_col = "SNP")

lead_specs[[length(lead_specs) + 1]] <- list(
  source_set = "deepvariant",
  variant_type = "snv_indel",
  lead_source = "clumped",
  lead_source_file = deepvariant_snv_clump_file,
  lead_ids = deepvariant_snv_ids,
  gwas_file = deepvariant_snv_gwas_file
)

# ============================================================
# Extract all leads from GWAS summary
# ============================================================

all_leads_list <- lapply(lead_specs, function(sp) {
  extract_leads_from_gwas(
    lead_ids = sp$lead_ids,
    source_set = sp$source_set,
    variant_type = sp$variant_type,
    lead_source = sp$lead_source,
    lead_source_file = sp$lead_source_file,
    gwas_file = sp$gwas_file
  )
})

all_lead_sig_gwas <- rbindlist(all_leads_list, fill = TRUE)
manual_exclusions <- read_manual_exclusions(manual_exclusion_file)
all_lead_sig_gwas <- apply_manual_exclusion_to_leads(
  all_lead_sig_gwas,
  manual_exclusions,
  file.path(table_dir, "lead_sig_from_gwas.manual_exclusion.audit.tsv")
)

# ============================================================
# Add unified useful columns
# ============================================================

#           ?CHR/POS/SNP/P              ?
if (!"CHR" %in% colnames(all_lead_sig_gwas)) {
  if ("#CHROM" %in% colnames(all_lead_sig_gwas)) {
    all_lead_sig_gwas[, CHR := `#CHROM`]
  } else if ("CHROM" %in% colnames(all_lead_sig_gwas)) {
    all_lead_sig_gwas[, CHR := CHROM]
  }
}

if (!"POS" %in% colnames(all_lead_sig_gwas)) {
  if ("BP" %in% colnames(all_lead_sig_gwas)) {
    all_lead_sig_gwas[, POS := BP]
  }
}

if (!"SNP" %in% colnames(all_lead_sig_gwas)) {
  if ("ID" %in% colnames(all_lead_sig_gwas)) {
    all_lead_sig_gwas[, SNP := ID]
  }
}

#                       GWAS    ?
all_lead_sig_gwas[, lead_chr := NA_character_]
all_lead_sig_gwas[, lead_pos := NA_real_]
all_lead_sig_gwas[, lead_p := NA_real_]

if ("CHR" %in% colnames(all_lead_sig_gwas)) {
  all_lead_sig_gwas[, lead_chr := sub("^chr", "", as.character(CHR), ignore.case = TRUE)]
}

if ("POS" %in% colnames(all_lead_sig_gwas)) {
  all_lead_sig_gwas[, lead_pos := suppressWarnings(as.numeric(POS))]
}

if ("P" %in% colnames(all_lead_sig_gwas)) {
  all_lead_sig_gwas[, lead_p := suppressWarnings(as.numeric(P))]
}

all_lead_sig_gwas <- all_lead_sig_gwas[lead_chr %in% as.character(1:22)]

#    
setorder(all_lead_sig_gwas, source_set, variant_type, lead_p, na.last = TRUE)

missing_leads <- all_lead_sig_gwas[found_in_gwas == FALSE]

summary_dt <- all_lead_sig_gwas[, .(
  n_rows = .N,
  n_unique_leads = uniqueN(lead_id),
  n_found_unique_leads = uniqueN(lead_id[found_in_gwas == TRUE]),
  n_missing_unique_leads = uniqueN(lead_id[found_in_gwas == FALSE])
), by = .(source_set, variant_type)]

setorder(summary_dt, source_set, variant_type)

# ============================================================
# Reorder columns for readability
# ============================================================

front_cols <- c(
  "source_set",
  "variant_type",
  "lead_id",
  "lead_chr",
  "lead_pos",
  "lead_p",
  "found_in_gwas",
  "lead_source",
  "lead_source_file",
  "gwas_file"
)

front_cols <- intersect(front_cols, colnames(all_lead_sig_gwas))
other_cols <- setdiff(colnames(all_lead_sig_gwas), front_cols)
setcolorder(all_lead_sig_gwas, c(front_cols, other_cols))

# ============================================================
# Output
# ============================================================

out_all_tsv <- file.path(table_dir, "lead_sig_from_gwas.all.tsv")
out_all_csv <- file.path(table_dir, "lead_sig_from_gwas.all.csv")

out_missing_tsv <- file.path(table_dir, "lead_sig_from_gwas.missing.tsv")
out_missing_csv <- file.path(table_dir, "lead_sig_from_gwas.missing.csv")

out_summary_tsv <- file.path(table_dir, "lead_sig_from_gwas.summary.tsv")
out_summary_csv <- file.path(table_dir, "lead_sig_from_gwas.summary.csv")

out_rdata <- file.path(rdata_dir, "lead_sig_from_gwas.Rdata")

fwrite(all_lead_sig_gwas, out_all_tsv, sep = "\t", quote = FALSE, na = "NA")
fwrite(all_lead_sig_gwas, out_all_csv, sep = ",", quote = TRUE, na = "NA")

fwrite(missing_leads, out_missing_tsv, sep = "\t", quote = FALSE, na = "NA")
fwrite(missing_leads, out_missing_csv, sep = ",", quote = TRUE, na = "NA")

fwrite(summary_dt, out_summary_tsv, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary_dt, out_summary_csv, sep = ",", quote = TRUE, na = "NA")

save(
  all_lead_sig_gwas,
  missing_leads,
  summary_dt,
  lead_specs,
  file = out_rdata
)

message("\n============================================================")
message("Done.")
message("All lead GWAS summary TSV: ", out_all_tsv)
message("All lead GWAS summary CSV: ", out_all_csv)
message("Missing lead IDs TSV: ", out_missing_tsv)
message("Missing lead IDs CSV: ", out_missing_csv)
message("Summary TSV: ", out_summary_tsv)
message("Summary CSV: ", out_summary_csv)
message("RData: ", out_rdata)

message("\nQuick check:")
print(summary_dt)

