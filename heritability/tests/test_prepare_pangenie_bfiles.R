#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
script_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile("pangenie_prepare_test_")
gwas_base <- file.path(tmp, "gwas")
out_dir <- file.path(tmp, "out")
dir.create(gwas_base, recursive = TRUE)

write_source <- function(source, bim) {
  prefix <- file.path(gwas_base, source, "NGS.QCsite.QCind")
  dir.create(dirname(prefix), recursive = TRUE)
  writeBin(as.raw(1L), paste0(prefix, ".bed"))
  fwrite(data.table("F1", "I1", 0, 0, 1, -9), paste0(prefix, ".fam"),
    sep = "\t", col.names = FALSE)
  fwrite(bim, paste0(prefix, ".bim"), sep = "\t", col.names = FALSE)
}

long_alt <- paste(rep("A", 50), collapse = "")
write_source("set00", data.table(
  CHR = 1, SNP = c("snv_a", "sv_a"), CM = 0, BP = c(100, 200),
  A1 = "A", A2 = c("G", long_alt)
))
write_source("set01", data.table(
  CHR = 1, SNP = c("snv_duplicate", "indel_b", "sv_b"), CM = 0,
  BP = c(100, 110, 210), A1 = "A", A2 = c("G", "AT", long_alt)
))
write_source("set02", data.table(
  CHR = 1, SNP = c("snv_c", "sv_c"), CM = 0, BP = c(120, 220),
  A1 = "C", A2 = c("T", long_alt)
))

cmd <- c(
  file.path(script_dir, "prepare_pangenie_bfiles.R"),
  "--gwas-base", gwas_base,
  "--out-dir", out_dir
)
if (system2(rscript, shQuote(cmd)) != 0L) stop("prepare_pangenie_bfiles.R failed")

variant_map <- fread(file.path(out_dir, "pangenie.variant_map.tsv"))
summary <- fread(file.path(out_dir, "pangenie.variant_summary.tsv"))
stopifnot(
  nrow(summary) == 6L,
  !anyDuplicated(variant_map[, .(component, variant_key)]),
  !"snv_duplicate" %in% variant_map$old_id,
  all(c("SNV_INDEL", "SV") %in% variant_map$component)
)
for (component in c("SNV_INDEL", "SV")) {
  for (source in c("set00", "set01", "set02")) {
    stopifnot(
      file.exists(file.path(out_dir, component, paste0(source, ".", component, ".old_ids.extract"))),
      file.exists(file.path(out_dir, component, paste0(source, ".", component, ".rename.tsv")))
    )
  }
}

message("[PASS] Pangenie BIM classification, deduplication, and map generation | tmp=", tmp)
