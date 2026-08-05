#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || identical(x, "")) y else x

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[sub("^--", "", key)]] <- TRUE
      i <- i + 1L
    } else {
      out[[sub("^--", "", key)]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
default_base <- Sys.getenv(
  "PANGENIE_HERITABILITY_OUT",
  file.path(dirname(script_dir), "results", "heritability", "pangenie")
)
base_dir <- args[["base-dir"]] %||% default_base
input_file <- args[["input"]] %||%
  file.path(base_dir, "summary", "heritability_liability_components_for_plot.tsv")
out_prefix <- args[["out-prefix"]] %||%
  file.path(base_dir, "plots", "heritability_liability_selected5_snvpc_pc20")
pc_n <- as.integer(args[["pc-n"]] %||% "20")
source_filter <- args[["source"]] %||% ""

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file,
       "\nRun summarise first, e.g. bash run_pangenie_heritability.sh --mode summarise")
}

d <- fread(input_file)
required <- c("source", "pc_n", "model_label", "adjustment_label",
  "component_label", "h2_liability_component")
missing <- setdiff(required, names(d))
if (length(missing)) stop("Input lacks required columns: ", paste(missing, collapse = ", "))

model_levels <- c("SNV", "SV", "TR", "SNV+SV", "SNV+SV+TR")
component_levels <- c("SNV_INDEL", "SV", "TR", "Unknown")
component_colours <- c(SNV_INDEL = "#4C78A8", SV = "#F58518", TR = "#54A24B", Unknown = "grey65")

pc_value <- pc_n
plot_dt <- d[
  pc_n == pc_value &
    adjustment_label == "SNV PC" &
    model_label %chin% model_levels
]
if (nzchar(source_filter)) plot_dt <- plot_dt[source == source_filter]
if (!nrow(plot_dt)) {
  stop("No rows after filtering: pc_n=", pc_n,
       ", adjustment_label='SNV PC', models=", paste(model_levels, collapse = ", "),
       if (nzchar(source_filter)) paste0(", source=", source_filter) else "")
}

plot_dt[, `:=`(
  model_label = factor(as.character(model_label), levels = model_levels),
  component_label = factor(as.character(component_label), levels = component_levels),
  h2_liability_component = as.numeric(h2_liability_component)
)]
setorder(plot_dt, source, model_label, component_label)

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(plot_dt, paste0(out_prefix, ".plot_data.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(plot_dt, paste0(out_prefix, ".plot_data.csv"), quote = TRUE, na = "NA")

p <- ggplot(plot_dt, aes(x = model_label, y = h2_liability_component, fill = component_label)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.15, na.rm = TRUE) +
  scale_fill_manual(values = component_colours, drop = FALSE) +
  labs(
    title = sprintf("Liability-scale heritability, SNV PC adjustment (PC%d)", pc_n),
    subtitle = "Selected models: SNV, SV, TR, SNV+SV, SNV+SV+TR; joint models are stacked by GRM component",
    x = NULL,
    y = "Liability-scale heritability",
    fill = "GRM component"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

if (uniqueN(plot_dt$source) > 1L) {
  p <- p + facet_grid(source ~ ., scales = "free_y")
}

height <- if (uniqueN(plot_dt$source) > 1L) 10 else 5.5
width <- 8.5
ggsave(paste0(out_prefix, ".png"), p, width = width, height = height, dpi = 220, limitsize = FALSE)
ggsave(paste0(out_prefix, ".pdf"), p, width = width, height = height, limitsize = FALSE)

message("[DONE] Plot data: ", out_prefix, ".plot_data.tsv")
message("[DONE] Plot PNG:  ", out_prefix, ".png")
message("[DONE] Plot PDF:  ", out_prefix, ".pdf")
