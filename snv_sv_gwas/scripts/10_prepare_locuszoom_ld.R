#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

base <- "/path/to/EOSCZ_PROJECT"
figure3 <- file.path(base, "figure_analysis")
ld_root <- file.path(figure3, "SV_SNV_LD", "LD_decay_public")
lead_file <- file.path(ld_root, "tables", "lead_sig_from_gwas.canonical_1000kb.tsv")
out_root <- file.path(ld_root, "locuszoom")
window_kb <- 500L
plink_bin <- Sys.getenv("PLINK_BIN", "plink")
max_plink_jobs <- 2L

pangenie <- file.path(base, "TGS_callset", "Pangenie_v3", "06.gwas")
dv_para <- file.path(base, "GWAS", "Deepvariant_paragraph")
merged_bfile <- file.path(dv_para, "chr_all2.strict_step2_genimi.common_samples.merged")
deepvar_gwas <- file.path(
  dv_para, "deepvar_gwas", "deepvar", "03_gwas",
  "SCZ.deepvar.mlm.geno0.1.maf0.01.fastGWA"
)
paragraph_gwas <- file.path(
  dv_para, "deepvar_gwas", "paragraph_test", "03_gwas",
  paste0(
    "SCZ.paragraph_test.pcsrc_deepvar_pc20_grm_deepvar_with_batch.mlm.",
    "geno0.1.maf0.01.fastGWA"
  )
)

source_config <- rbindlist(list(
  rbindlist(lapply(c("set00"), function(s) {
    data.table(
      source_set = s,
      bfile = file.path(pangenie, s, "NGS.QCsite.QCind"),
      regional_gwas_file = file.path(pangenie, s, "gwas", "SCZ.mlm.ngspc.fastGWA"),
      overlay_gwas_file = NA_character_
    )
  })),
  data.table(
    source_set = c("deepvariant", "paragraph"),
    bfile = merged_bfile,
    regional_gwas_file = deepvar_gwas,
    overlay_gwas_file = paragraph_gwas
  )
), use.names = TRUE)

check_file <- function(path, label = "file") {
  if (!file.exists(path) || dir.exists(path) || file.info(path)$size == 0) {
    stop("Missing or empty ", label, ": ", path)
  }
}

check_bfile <- function(prefix) {
  missing <- paste0(prefix, c(".bed", ".bim", ".fam"))
  missing <- missing[!file.exists(missing)]
  if (length(missing)) stop("Missing bfile component(s):\n", paste(missing, collapse = "\n"))
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

quote_cmd <- function(x) shQuote(as.character(x), type = "sh")

check_file(lead_file, "lead table")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(source_config))) {
  check_bfile(source_config$bfile[i])
  check_file(source_config$regional_gwas_file[i], "GWAS")
  if (!is.na(source_config$overlay_gwas_file[i])) {
    check_file(source_config$overlay_gwas_file[i], "overlay GWAS")
  }
}

leads <- fread(lead_file, showProgress = FALSE)
required <- c(
  "source_set", "variant_type", "lead_id", "lead_chr", "lead_pos",
  "lead_p", "found_in_gwas"
)
missing_cols <- setdiff(required, names(leads))
if (length(missing_cols)) stop("Lead table missing: ", paste(missing_cols, collapse = ", "))

leads[, source_set := as.character(source_set)]
leads[, variant_type := tolower(as.character(variant_type))]
leads[, lead_id := as.character(lead_id)]
leads[, lead_chr := sub("^chr", "", as.character(lead_chr), ignore.case = TRUE)]
leads[, lead_pos := suppressWarnings(as.integer(lead_pos))]
leads[, lead_p := suppressWarnings(as.numeric(lead_p))]
leads[, found_ok := tolower(as.character(found_in_gwas)) %in% c("true", "1", "t")]
leads <- leads[
  found_ok & !is.na(lead_id) & lead_id != "" &
    lead_chr %in% as.character(1:22) & !is.na(lead_pos) &
    is.finite(lead_p) & lead_p > 0 & lead_p <= 1
]
leads[, plot_class := fifelse(variant_type %in% c("snv_indel", "sgv"), "SGV", "SV")]
setorder(leads, source_set, variant_type, lead_id, lead_p)
leads <- leads[, .SD[1L], by = .(source_set, variant_type, lead_id)]

unknown_sources <- setdiff(unique(leads$source_set), source_config$source_set)
if (length(unknown_sources)) {
  stop("No locuszoom source configuration for: ", paste(unknown_sources, collapse = ", "))
}

jobs <- merge(leads, source_config, by = "source_set", all.x = TRUE, sort = FALSE)
bim_cache <- new.env(parent = emptyenv())
get_bim_ids <- function(prefix) {
  key <- prefix
  if (!exists(key, envir = bim_cache, inherits = FALSE)) {
    bim <- fread(
      paste0(prefix, ".bim"), header = FALSE,
      col.names = c("CHR", "SNP", "BP", "A1", "A2"),
      select = c(1, 2, 4, 5, 6), showProgress = FALSE
    )
    bim[, SNP := as.character(SNP)]
    assign(key, bim, envir = bim_cache)
  }
  get(key, envir = bim_cache, inherits = FALSE)
}

jobs[, lead_in_bim := FALSE]
for (i in seq_len(nrow(jobs))) {
  jobs$lead_in_bim[i] <- jobs$lead_id[i] %in% get_bim_ids(jobs$bfile[i])$SNP
}

rejected <- jobs[lead_in_bim == FALSE]
fwrite(rejected, file.path(out_root, "locuszoom_leads.rejected.tsv"), sep = "\t", na = "NA")
if (nrow(rejected)) {
  stop(
    nrow(rejected), " lead(s) are absent from their source BIM. See ",
    file.path(out_root, "locuszoom_leads.rejected.tsv")
  )
}

jobs[, stem := safe_name(paste(source_set, paste0("chr", lead_chr), lead_pos, lead_id, sep = "_"))]
jobs[, ld_prefix := file.path(out_root, plot_class, "ld", source_set, stem)]
jobs[, plot_prefix := file.path(out_root, plot_class, "plots", source_set, stem)]
jobs[, regional_table := file.path(
  out_root, plot_class, "tables", source_set, paste0(stem, ".regional_ld.tsv.gz")
)]

for (d in unique(c(
  dirname(jobs$ld_prefix), dirname(jobs$plot_prefix), dirname(jobs$regional_table)
))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

jobs[, ld_file := paste0(ld_prefix, ".ld.gz")]
jobs[, command := paste(
  quote_cmd(plink_bin),
  "--bfile", quote_cmd(bfile),
  "--threads 1",
  "--ld-snp", quote_cmd(lead_id),
  "--r2 gz",
  "--ld-window-kb", window_kb,
  "--ld-window", 999999,
  "--ld-window-r2 0",
  "--out", quote_cmd(ld_prefix)
)]

metadata_file <- file.path(out_root, "locuszoom_jobs.metadata.tsv")
command_file <- file.path(out_root, "run_locuszoom_ld.sh")
fwrite(jobs, metadata_file, sep = "\t", quote = FALSE, na = "NA")
task_lines <- paste0(
  "  ", quote_cmd(jobs$ld_file), " ", quote_cmd(jobs$command), " \\"
)
verify_lines <- paste0("  ", quote_cmd(jobs$ld_file), " \\")
writeLines(c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  sprintf("MAX_PLINK_JOBS=%d", max_plink_jobs),
  sprintf(
    "printf '[INFO] Running %d lead-specific locuszoom LD jobs (up to %%s concurrent; 1 thread each)\\n' \"$MAX_PLINK_JOBS\"",
    nrow(jobs)
  ),
  "",
  "ld_complete() {",
  "  [[ -s \"$1\" ]] && gzip -t \"$1\" 2>/dev/null",
  "}",
  "run_one() {",
  "  local ld_file=\"$1\"",
  "  local cmd=\"$2\"",
  "  local status",
  "  if ld_complete \"$ld_file\"; then",
  "    printf '[SKIP] Existing LD: %s\\n' \"$ld_file\"",
  "    return 0",
  "  fi",
  "  printf '[RUN] Missing LD: %s\\n' \"$ld_file\"",
  "  if bash -c \"$cmd\"; then",
  "    if ld_complete \"$ld_file\"; then",
  "      printf '[DONE] LD: %s\\n' \"$ld_file\"",
  "      return 0",
  "    fi",
  "    printf '[ERROR] PLINK returned zero but LD is missing/empty: %s\\n' \"$ld_file\" >&2",
  "    return 1",
  "  else",
  "    status=$?",
  "    printf '[ERROR] PLINK exit %s: %s\\n' \"$status\" \"$ld_file\" >&2",
  "    return 1",
  "  fi",
  "}",
  "export -f ld_complete run_one",
  "",
  "xargs_status=0",
  "printf '%s\\0%s\\0' \\",
  task_lines,
  "  | xargs -0 -n 2 -P \"$MAX_PLINK_JOBS\" bash -c 'run_one \"$1\" \"$2\"' _ || xargs_status=$?",
  "",
  "missing=0",
  "for ld_file in \\",
  verify_lines,
  "; do",
  "  if ! ld_complete \"$ld_file\"; then",
  "    printf '[MISSING] %s\\n' \"$ld_file\" >&2",
  "    ((missing += 1))",
  "  fi",
  "done",
  "if (( missing > 0 )); then",
  "  printf '[ERROR] %s LD job(s) remain missing; xargs status=%s\\n' \"$missing\" \"$xargs_status\" >&2",
  "  exit 1",
  "fi",
  "echo '[DONE] All locuszoom LD jobs finished.'"
), command_file, useBytes = TRUE)
Sys.chmod(command_file, mode = "0755")

class_summary <- jobs[, .(n_leads = .N), by = .(plot_class, source_set)]
fwrite(class_summary, file.path(out_root, "locuszoom_jobs.summary.tsv"), sep = "\t")
fwrite(source_config, file.path(out_root, "locuszoom_source_config.tsv"), sep = "\t", na = "NA")

message("[DONE] Locuszoom metadata: ", metadata_file)
message("[DONE] PLINK LD commands: ", command_file)
print(class_summary)
