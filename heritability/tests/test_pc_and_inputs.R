#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
script_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile("figure7_pc_inputs_test_"); dir.create(tmp, recursive = TRUE)
set.seed(7)
n <- 120L
ids <- data.table(FID = paste0("F", seq_len(n)), IID = paste0("I", seq_len(n)))
pheno <- copy(ids)[, PHENO := rep(c(1, 2), length.out = n)]
batch <- copy(ids)[, BATCH := rep(letters[1:4], length.out = n)]
fwrite(pheno, file.path(tmp, "pheno.txt"), sep = "\t", col.names = FALSE)
fwrite(batch, file.path(tmp, "batch.txt"), sep = "\t", col.names = FALSE)
for (nm in c("snv", "sv", "tr")) fwrite(ids, file.path(tmp, paste0(nm, ".grm.id")), sep = "\t", col.names = FALSE)
fwrite(ids[1:2], file.path(tmp, "remove_sample.tsv"), sep = "\t", col.names = FALSE)

make_pc <- function(component) {
  x <- copy(ids)
  for (j in 1:20) x[, (paste0("PC", j)) := rnorm(n)]
  fwrite(x[, c("FID", "IID", paste0("PC", 1:10)), with = FALSE], file.path(tmp, paste0(component, ".pc10.tsv")), sep = "\t")
  fwrite(x, file.path(tmp, paste0(component, ".pc20.tsv")), sep = "\t")
}
for (component in c("snv", "sv", "tr")) make_pc(component)

out_dir <- file.path(tmp, "inputs")
cmd <- c(file.path(script_dir, "build_analysis_inputs.R"), "--source", "pangenie",
  "--pheno", file.path(tmp, "pheno.txt"), "--batch", file.path(tmp, "batch.txt"),
  "--snv-pc10", file.path(tmp, "snv.pc10.tsv"), "--snv-pc20", file.path(tmp, "snv.pc20.tsv"),
  "--sv-pc10", file.path(tmp, "sv.pc10.tsv"), "--sv-pc20", file.path(tmp, "sv.pc20.tsv"),
  "--tr-pc10", file.path(tmp, "tr.pc10.tsv"), "--tr-pc20", file.path(tmp, "tr.pc20.tsv"),
  "--snv-grm-id", file.path(tmp, "snv.grm.id"), "--sv-grm-id", file.path(tmp, "sv.grm.id"),
  "--tr-grm-id", file.path(tmp, "tr.grm.id"), "--remove-file", file.path(tmp, "remove_sample.tsv"), "--out-dir", out_dir)
if (system2(rscript, shQuote(cmd)) != 0L) stop("build_analysis_inputs.R failed")
manifest <- fread(file.path(out_dir, "pangenie.qcovar_manifest.tsv"))
stopifnot(nrow(manifest) == 14L, all(manifest$valid), all(manifest$n_samples == n - 2L))
qcovar_first_line <- readLines(file.path(out_dir, "pangenie.pc10.snv_pc.qcovar"), n = 1L)
stopifnot(!grepl('"', qcovar_first_line, fixed = TRUE))

sscore <- copy(ids)
sscore[, `:=`(ALLELE_CT = 100, DENOM = 100)]
for (j in 1:20) sscore[, (paste0("PC", j, "_AVG")) := rnorm(n)]
setnames(sscore, "FID", "#FID")
fwrite(sscore, file.path(tmp, "projected.sscore"), sep = "\t")
fwrite(data.table(value = rev(seq_len(20))), file.path(tmp, "pca.eigenval"), col.names = FALSE)
fwrite(ids[1:90], file.path(tmp, "unrelated.keep"), sep = "\t", col.names = FALSE)
pc_prefix <- file.path(tmp, "pangenie.SNV_INDEL")
cmd <- c(file.path(script_dir, "prepare_pc_outputs.R"), "--source", "pangenie", "--component", "SNV_INDEL",
  "--method", "unrelated_projected", "--scores", file.path(tmp, "projected.sscore"),
  "--eigenval", file.path(tmp, "pca.eigenval"), "--unrelated", file.path(tmp, "unrelated.keep"),
  "--pheno", file.path(tmp, "pheno.txt"), "--batch", file.path(tmp, "batch.txt"), "--out-prefix", pc_prefix)
if (system2(rscript, shQuote(cmd)) != 0L) stop("prepare_pc_outputs.R failed")
stopifnot(nrow(fread(paste0(pc_prefix, ".pc10.qcovar.tsv"))) == n)

plot_prefix <- file.path(tmp, "plots", "pangenie.SNV_INDEL")
cmd <- c(file.path(script_dir, "plot_pca.R"), "--plot-data", paste0(pc_prefix, ".plot.tsv"),
  "--eigenval", paste0(pc_prefix, ".eigenval.tsv"), "--out-prefix", plot_prefix)
if (system2(rscript, shQuote(cmd)) != 0L) stop("plot_pca.R failed")
stopifnot(
  file.exists(paste0(plot_prefix, ".pc1_pc2_case.png")),
  file.exists(paste0(plot_prefix, ".pc9_pc10_batch.png")),
  file.exists(paste0(plot_prefix, ".case_control_pc1_10.pdf")),
  file.exists(paste0(plot_prefix, ".batch_pc1_10.pdf")),
  file.exists(paste0(plot_prefix, ".all_pca_plots.pdf"))
)
message("[PASS] common sample keep, 14 qcovars, PC conversion, and PCA rendering | tmp=", tmp)
