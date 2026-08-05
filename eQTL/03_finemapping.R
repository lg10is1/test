#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 5) {
  stop(paste(
    "Usage: Rscript 03_finemapping.R QTL_RESULTS REFERENCE_BFILE OUTPUT_DIR",
    "[N_GWAS=139] [THREADS=16]"
  ))
}

suppressPackageStartupMessages(library(data.table))

qtl_file <- args[1]
reference_bfile <- args[2]
output_dir <- args[3]
n_gwas <- if (length(args) >= 4) as.integer(args[4]) else 139L
threads <- if (length(args) >= 5) as.integer(args[5]) else 16L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qtl <- fread(qtl_file, header = FALSE)
if (ncol(qtl) >= 15) {
  beta_column <- 14
  best_column <- 15
} else {
  beta_column <- 13
  best_column <- 14
}

associations <- data.table(
  phenotype = qtl[[1]],
  variant = qtl[[8]],
  chromosome = qtl[[9]],
  position = qtl[[10]],
  p = qtl[[12]],
  beta = qtl[[beta_column]],
  best = qtl[[best_column]]
)

targets <- unique(
  associations[best == 1 & p < 5e-3 & chromosome != "chrX", phenotype]
)

for (phenotype_id in targets) {
  locus <- associations[phenotype == phenotype_id & !is.na(p) & !is.na(beta)]
  safe_id <- gsub("[^A-Za-z0-9_.-]", "_", phenotype_id)
  locus_dir <- file.path(output_dir, safe_id)
  dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)

  extract_file <- file.path(locus_dir, "variants.txt")
  fwrite(locus[, .(variant)], extract_file, col.names = FALSE)

  prefix <- file.path(locus_dir, safe_id)
  plink_status <- system2(
    "plink",
    c(
      "--bfile", reference_bfile,
      "--extract", extract_file,
      "--keep-allele-order",
      "--make-bed",
      "--r", "square",
      "--out", prefix
    )
  )
  if (plink_status != 0) next

  bim <- fread(paste0(prefix, ".bim"), header = FALSE)
  setnames(bim, c("chr", "variant", "cm", "bp", "A1", "A2"))
  summary <- merge(bim, locus, by = "variant")
  if (nrow(summary) == 0) next

  summary[, z := -qnorm(pmax(p, .Machine$double.xmin) / 2) * sign(beta)]
  summary[, se := beta / z]
  summary[, bp2 := seq_len(.N)]

  sst <- summary[, .(
    chr = sub("^chr", "", chr),
    snp = variant,
    bp,
    A1,
    A2,
    beta,
    se,
    stat = z,
    p,
    bp2
  )]

  bim[, bp := seq_len(.N)]
  fwrite(bim, paste0(prefix, ".bim"), sep = "\t", col.names = FALSE)
  sst_file <- paste0(prefix, ".sst")
  fwrite(sst, sst_file, sep = "\t")

  system2(
    "SuSiEx",
    c(
      paste0("--sst_file=", sst_file),
      paste0("--n_gwas=", n_gwas),
      paste0("--ref_file=", prefix),
      paste0("--ld_file=", prefix),
      paste0("--out_dir=", locus_dir),
      paste0("--out_name=", safe_id),
      paste0("--chr=", sst$chr[1]),
      "--bp=1,1000000000",
      "--snp_col=2", "--chr_col=1", "--bp_col=10",
      "--a1_col=4", "--a2_col=5", "--eff_col=6",
      "--se_col=7", "--pval_col=9", "--plink=plink",
      "--mult-step=True", "--keep-ambig=True",
      paste0("--threads=", threads),
      "--pval_thresh=0.99", "--maf=0.00001",
      "--n_sig=10", "--min_purity=0.0"
    )
  )
}

