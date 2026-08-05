#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

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
plot_file <- args[["plot-data"]] %||% stop("Missing --plot-data")
eigen_file <- args[["eigenval"]] %||% stop("Missing --eigenval")
out_prefix <- args[["out-prefix"]] %||% stop("Missing --out-prefix")
for (f in c(plot_file, eigen_file)) if (!file.exists(f)) stop("Missing file: ", f)

d <- fread(plot_file)
e <- fread(eigen_file)
needed <- c(paste0("PC", 1:10), "pca_reference", "phenotype", "batch")
if (!all(needed %in% names(d))) stop("PCA plot table lacks: ", paste(setdiff(needed, names(d)), collapse = ", "))
title_base <- paste(unique(d$source), unique(d$component), unique(d$pca_method), sep = " | ")

d[, phenotype_label := fifelse(as.character(phenotype) %in% c("2", "case", "CASE"), "Case",
  fifelse(as.character(phenotype) %in% c("1", "control", "CONTROL"), "Control", as.character(phenotype)))]
d[, batch_label := as.factor(batch)]
d[, reference_label := as.factor(pca_reference)]

pc_pct <- setNames(e$percent, e$PC)
axis_label <- function(pc) {
  pct <- pc_pct[[pc]]
  if (!is.null(pct) && is.finite(pct)) sprintf("%s (%.2f%%)", pc, pct) else pc
}

make_scatter <- function(x, y, colour_col, colour_label) {
  ggplot(d, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour_col]], shape = reference_label)) +
    geom_point(size = 1.8, alpha = 0.78, na.rm = TRUE) +
    labs(title = title_base, subtitle = paste(colour_label, "colour; reference/projection shape"),
      x = axis_label(x), y = axis_label(y), colour = colour_label, shape = "PCA sample") +
    theme_bw(base_size = 11) + theme(legend.position = "right")
}

pc_pairs <- data.table(
  x = paste0("PC", seq(1L, 9L, by = 2L)),
  y = paste0("PC", seq(2L, 10L, by = 2L))
)
case_plots <- list()
batch_plots <- list()
for (i in seq_len(nrow(pc_pairs))) {
  x <- pc_pairs$x[[i]]
  y <- pc_pairs$y[[i]]
  pair_name <- tolower(paste(x, y, sep = "_"))
  case_plots[[paste0(pair_name, "_case")]] <- make_scatter(x, y, "phenotype_label", "Phenotype")
  batch_plots[[paste0(pair_name, "_batch")]] <- make_scatter(x, y, "batch_label", "Batch")
}

scree <- e[seq_len(min(20L, .N))]
scree_plot <- ggplot(scree, aes(x = factor(PC, levels = PC), y = percent)) +
  geom_col(fill = "#356AA0") +
  labs(title = title_base, subtitle = "PCA scree plot", x = "Principal component", y = "Variance explained (%)") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plots <- c(case_plots, batch_plots, list(scree = scree_plot))

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
for (nm in names(plots)) {
  ggsave(paste0(out_prefix, ".", nm, ".png"), plots[[nm]], width = 8.2, height = 6.2, dpi = 180)
  ggsave(paste0(out_prefix, ".", nm, ".pdf"), plots[[nm]], width = 8.2, height = 6.2)
}

pdf(paste0(out_prefix, ".all_pca_plots.pdf"), width = 8.2, height = 6.2, onefile = TRUE)
for (p in plots) print(p)
dev.off()
pdf(paste0(out_prefix, ".case_control_pc1_10.pdf"), width = 8.2, height = 6.2, onefile = TRUE)
for (p in case_plots) print(p)
dev.off()
pdf(paste0(out_prefix, ".batch_pc1_10.pdf"), width = 8.2, height = 6.2, onefile = TRUE)
for (p in batch_plots) print(p)
dev.off()
pdf(paste0(out_prefix, ".scree.pdf"), width = 8.2, height = 6.2, onefile = TRUE)
print(scree_plot)
dev.off()
message("[DONE] PCA plots: ", out_prefix)
