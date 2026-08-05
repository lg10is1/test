#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
script_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile("figure7_heritability_test_")
dir.create(tmp, recursive = TRUE)
sources <- "pangenie"

write_hsq <- function(path, components) {
  n <- length(components)
  if (n == 1L) {
    x <- data.table(Source = c("V(G)", "V(e)", "Vp", "V(G)/Vp", "V(G)/Vp_L", "Pval", "n"),
      Variance = c(.25, .75, 1, .25, .14, .01, 1000), SE = c(.03, .03, .01, .03, .02, NA, NA))
  } else {
    src <- c(paste0("V(G", seq_len(n), ")"), "V(e)", "Vp",
      paste0("V(G", seq_len(n), ")/Vp"), "Sum of V(G)/Vp",
      paste0("V(G", seq_len(n), ")/Vp_L"), "Sum of V(G)_L/Vp", "Pval", "n")
    vals <- c(rep(.1, n), .7, 1, rep(.1, n), .1 * n, rep(.056, n), .056 * n, .02, 1000)
    ses <- c(rep(.02, n), .03, .01, rep(.02, n), .03, rep(.012, n), .018, NA, NA)
    x <- data.table(Source = src, Variance = vals, SE = ses)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, path, sep = "\t", na = "NA")
}

for (source in sources) {
  input_dir <- file.path(tmp, source, "inputs")
  reml_dir <- file.path(tmp, source, "reml")
  log_dir <- file.path(tmp, source, "logs", "reml")
  dir.create(input_dir, recursive = TRUE)
  cmd <- c(file.path(script_dir, "make_reml_tasks.R"),
    "--source", source, "--input-dir", input_dir, "--reml-dir", reml_dir,
    "--log-dir", log_dir, "--snv-grm", "/fake/snv", "--sv-grm", "/fake/sv",
    "--tr-grm", "/fake/tr", "--output", file.path(input_dir, paste0(source, ".reml_tasks.tsv")))
  rc <- system2(rscript, shQuote(cmd))
  if (rc != 0L) stop("make_reml_tasks.R failed for ", source)
  tasks <- fread(file.path(input_dir, paste0(source, ".reml_tasks.tsv")))
  stopifnot(nrow(tasks) == 34L, !anyDuplicated(tasks$task_id))
  statuses <- tasks[, .(source, task_id, model, pc_n, adjustment)]
  statuses[, `:=`(status = "success", exit_code = 0L, error_reason = "", log_file = tasks$log_file,
    hsq_file = paste0(tasks$out_prefix, ".hsq"), finished_at = "test")]
  statuses[seq(2L, nrow(statuses), by = 7L), `:=`(status = "failed", exit_code = 1L, error_reason = "matrix not invertible")]
  statuses[seq(3L, nrow(statuses), by = 7L), `:=`(status = "warning", exit_code = 0L, error_reason = "REML did not converge")]
  statuses[seq(4L, nrow(statuses), by = 7L), `:=`(status = "failed", exit_code = 0L, error_reason = "missing .hsq")]
  for (i in seq_len(nrow(tasks))) {
    if (statuses$status[[i]] %in% c("success", "warning")) write_hsq(paste0(tasks$out_prefix[[i]], ".hsq"), strsplit(tasks$components[[i]], ",", fixed = TRUE)[[1]])
  }
  status_dir <- file.path(tmp, source, "status"); dir.create(status_dir, recursive = TRUE)
  fwrite(statuses, file.path(status_dir, paste0(source, ".reml_status.tsv")), sep = "\t")
}

summary_dir <- file.path(tmp, "summary"); plot_dir <- file.path(tmp, "plots")
cmd <- c(file.path(script_dir, "collect_results.R"), "--sources", paste(sources, collapse = ","),
  "--out-base", tmp, "--summary-dir", summary_dir, "--plot-dir", plot_dir)
rc <- system2(rscript, shQuote(cmd))
if (rc != 0L) stop("collect_results.R failed")

totals <- fread(file.path(summary_dir, "heritability_model_totals.tsv"))
focus <- fread(file.path(summary_dir, "heritability_focus_snv_models.long.tsv"))
errors <- fread(file.path(summary_dir, "heritability_errors.tsv"))
stopifnot(nrow(totals) == 34L, nrow(focus) == 16L, nrow(errors) > 0L)
stopifnot(any(grepl("not invertible", errors$error_reason, fixed = TRUE)))
stopifnot(file.exists(file.path(summary_dir, "heritability_model_totals.csv")))
stopifnot(file.exists(file.path(summary_dir, "heritability_focus_snv_pc_only.wide.csv")))
stopifnot(file.exists(file.path(summary_dir, "heritability_liability_components_for_plot.tsv")))
stopifnot(file.exists(file.path(summary_dir, "heritability_liability_selected5_for_plot.tsv")))
stopifnot(file.exists(file.path(summary_dir, "heritability_liability_overview_unfaceted_for_plot.tsv")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_components_pc10.png")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_components_pc20.pdf")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_components_all.pdf")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_selected5_components.png")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_selected5_components.pdf")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_overview_unfaceted_pc10.png")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_overview_unfaceted_pc20.pdf")))
stopifnot(file.exists(file.path(plot_dir, "heritability_liability_overview_unfaceted_all.pdf")))
message("[PASS] task matrix, failure retention, CSV/TSV outputs, and plots | tmp=", tmp)
