#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || identical(x, "")) y else x

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Invalid argument near: ", x[[i]])
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
sources <- strsplit(args[["sources"]] %||% "set00,set01,set02", ",", fixed = TRUE)[[1]]
gwas_base <- args[["gwas-base"]] %||% Sys.getenv("GWAS_BASE", "")
if (!nzchar(gwas_base)) stop("Set --gwas-base or GWAS_BASE")
out_dir <- args[["out-dir"]] %||% stop("Missing --out-dir")
small_mode <- args[["small-variant-mode"]] %||% "SNV_INDEL_LT50"
if (!small_mode %in% c("SNV_INDEL_LT50", "STRICT_SNV")) stop("Invalid small-variant mode: ", small_mode)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bfile_for <- function(source) file.path(gwas_base, source, "NGS.QCsite.QCind")
suffix_for <- function(source) {
  if (grepl("^set[0-9][0-9]$", source)) return(sub("^set", "_", source))
  paste0("_", source)
}

read_bim <- function(source) {
  prefix <- bfile_for(source)
  files <- paste0(prefix, c(".bed", ".bim", ".fam"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing PLINK file(s): ", paste(missing, collapse = ", "))
  x <- fread(paste0(prefix, ".bim"), header = FALSE)
  if (ncol(x) < 6L) stop("BIM has fewer than six columns: ", prefix, ".bim")
  setnames(x, 1:6, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  x[, `:=`(
    source = source,
    bfile = prefix,
    SNP = as.character(SNP),
    A1 = as.character(A1),
    A2 = as.character(A2)
  )]
  x[, max_allele_len := pmax(nchar(A1), nchar(A2))]
  x[, sv_tag := grepl("(^|[-_])(INS|DEL|DUP|INV)([-_]|$)", SNP, ignore.case = TRUE)]
  x[, strict_snv := nchar(A1) == 1L & nchar(A2) == 1L &
      grepl("^[ACGT]$", A1, ignore.case = TRUE) & grepl("^[ACGT]$", A2, ignore.case = TRUE)]
  x[, allele_lo := fifelse(A1 <= A2, A1, A2)]
  x[, allele_hi := fifelse(A1 <= A2, A2, A1)]
  x[, variant_key := paste(CHR, BP, allele_lo, allele_hi, sep = ":")]
  x
}

candidate_by_component <- function(bim, source, component) {
  if (component == "SV") {
    out <- bim[max_allele_len >= 50L | sv_tag]
    definition <- "BIM fallback: maximum allele length >=50 bp or INS/DEL/DUP/INV-tagged ID"
  } else {
    if (small_mode == "STRICT_SNV") {
      out <- bim[strict_snv & !sv_tag]
      definition <- "BIM fallback: single-base A/C/G/T alleles excluding SV-tagged IDs"
    } else {
      out <- bim[max_allele_len < 50L & !sv_tag]
      definition <- "BIM fallback: maximum allele length <50 bp excluding SV-tagged IDs"
    }
  }
  out[, selection_source := definition]
  out
}

selected_list <- list()
summary_list <- list()
source_bfiles <- list()
seen_keys <- list(SNV_INDEL = character(), SV = character())

for (source in sources) {
  message("[INFO] Reading BIM and selecting variants for ", source)
  bim <- read_bim(source)
  source_bfiles[[source]] <- unique(bim[, .(source, bfile)])
  for (component in c("SNV_INDEL", "SV")) {
    candidates <- candidate_by_component(bim, source, component)
    n_input <- nrow(candidates)
    if (n_input) {
      candidates <- candidates[!duplicated(variant_key)]
    }
    n_unique <- nrow(candidates)
    if (n_unique) {
      keep <- !candidates$variant_key %chin% seen_keys[[component]]
      chosen <- candidates[keep]
      if (nrow(chosen)) {
        chosen[, `:=`(
          component = component,
          old_id = SNP,
          new_id_base = paste0(SNP, suffix_for(source))
        )]
        selected_list[[paste(source, component, sep = ".")]] <- chosen
        seen_keys[[component]] <- c(seen_keys[[component]], chosen$variant_key)
      }
    } else {
      chosen <- candidates[0]
    }
    summary_list[[paste(source, component, sep = ".")]] <- data.table(
      component = component,
      source = source,
      bfile = bfile_for(source),
      n_input_variants = n_input,
      n_unique_variant_keys_within_source = n_unique,
      n_selected_after_cross_source_dedup = nrow(chosen),
      n_dropped_as_duplicate_of_prior_source = n_unique - nrow(chosen),
      selection_source = if (n_input) candidates$selection_source[[1]] else "BIM classification"
    )
  }
  rm(bim)
  gc(verbose = FALSE)
}

selected <- rbindlist(selected_list, use.names = TRUE, fill = TRUE)
if (!nrow(selected)) stop("No variants selected after merging/deduplication")
selected[, new_id := make.unique(new_id_base, sep = "_dup")]
summary <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
setorder(summary, component, source)

fwrite(selected[, .(component, source, old_id, new_id, CHR, BP, A1, A2, variant_key, bfile,
  selection_source)],
  file.path(out_dir, "pangenie.variant_map.tsv"), sep = "\t", quote = FALSE)
fwrite(summary, file.path(out_dir, "pangenie.variant_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(summary, file.path(out_dir, "pangenie.variant_summary.csv"))
fwrite(rbindlist(source_bfiles, use.names = TRUE, fill = TRUE),
  file.path(out_dir, "source_bfiles.tsv"), sep = "\t", quote = FALSE)

for (comp in c("SNV_INDEL", "SV")) {
  component_dir <- file.path(out_dir, comp)
  dir.create(component_dir, recursive = TRUE, showWarnings = FALSE)
  for (source in sources) {
    src <- source
    x <- selected[component == comp & source == src]
    extract_file <- file.path(component_dir, paste0(src, ".", comp, ".old_ids.extract"))
    rename_file <- file.path(component_dir, paste0(src, ".", comp, ".rename.tsv"))
    fwrite(x[, .(old_id)], extract_file, sep = "\t", col.names = FALSE, quote = FALSE)
    fwrite(x[, .(old_id, new_id)], rename_file, sep = "\t", col.names = FALSE, quote = FALSE)
  }
}

message("[DONE] Pangenie variant maps: ", out_dir)
