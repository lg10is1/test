#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || identical(x, "")) y else x
parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Invalid argument near: ", x[[i]])
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
source_name <- args[["source"]] %||% stop("Missing --source")
input_dir <- args[["input-dir"]] %||% stop("Missing --input-dir")
reml_dir <- args[["reml-dir"]] %||% stop("Missing --reml-dir")
log_dir <- args[["log-dir"]] %||% stop("Missing --log-dir")
snv_grm <- args[["snv-grm"]] %||% stop("Missing --snv-grm")
sv_grm <- args[["sv-grm"]] %||% stop("Missing --sv-grm")
tr_grm <- args[["tr-grm"]] %||% stop("Missing --tr-grm")
output <- args[["output"]] %||% stop("Missing --output")

spec <- rbindlist(list(
  data.table(model = "SNV_INDEL.marginal", input_type = "grm", input_path = snv_grm,
    components = "SNV_INDEL", adjustment = "snv_pc"),
  data.table(model = "SV.marginal", input_type = "grm", input_path = sv_grm,
    components = "SV", adjustment = c("snv_pc", "sv_pc", "snv_sv_pc")),
  data.table(model = "TR.marginal", input_type = "grm", input_path = tr_grm,
    components = "TR", adjustment = c("snv_pc", "tr_pc", "snv_tr_pc")),
  data.table(model = "SNV_INDEL_SV.joint", input_type = "mgrm",
    input_path = file.path(input_dir, paste0(source_name, ".SNV_INDEL_SV.mgrm")),
    components = "SNV_INDEL,SV", adjustment = c("snv_pc", "snv_sv_pc")),
  data.table(model = "SNV_INDEL_TR.joint", input_type = "mgrm",
    input_path = file.path(input_dir, paste0(source_name, ".SNV_INDEL_TR.mgrm")),
    components = "SNV_INDEL,TR", adjustment = c("snv_pc", "snv_tr_pc")),
  data.table(model = "SV_TR.joint", input_type = "mgrm",
    input_path = file.path(input_dir, paste0(source_name, ".SV_TR.mgrm")),
    components = "SV,TR", adjustment = c("snv_pc", "sv_tr_pc", "snv_sv_tr_pc")),
  data.table(model = "SNV_INDEL_SV_TR.joint", input_type = "mgrm",
    input_path = file.path(input_dir, paste0(source_name, ".SNV_INDEL_SV_TR.mgrm")),
    components = "SNV_INDEL,SV,TR", adjustment = c("snv_pc", "snv_tr_pc", "snv_sv_tr_pc"))
))

tasks <- rbindlist(lapply(c(10L, 20L), function(n_pc) copy(spec)[, pc_n := n_pc]))
tasks[, `:=`(
  source = source_name,
  task_id = sprintf("%s.pc%d.%s", model, pc_n, adjustment),
  qcovar = file.path(input_dir, sprintf("%s.pc%d.%s.qcovar", source_name, pc_n, adjustment))
)]
tasks[, out_prefix := file.path(reml_dir, task_id)]
tasks[, log_file := file.path(log_dir, paste0(task_id, ".log"))]
setcolorder(tasks, c("source", "task_id", "model", "pc_n", "adjustment", "input_type",
  "input_path", "qcovar", "out_prefix", "log_file", "components"))
if (nrow(tasks) != 34L) stop("Internal error: expected 34 tasks for the analysis, got ", nrow(tasks))
if (anyDuplicated(tasks$task_id)) stop("Internal error: duplicated task IDs")

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
dir.create(reml_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(tasks, output, sep = "\t")
fwrite(tasks, sub("\\.tsv$", ".csv", output))
message("[DONE] REML task manifest: ", output, " | tasks=", nrow(tasks))
