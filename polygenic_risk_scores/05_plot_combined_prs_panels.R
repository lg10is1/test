#!/usr/bin/env Rscript
############################################################
## Combined PRS summary figure.
##
## Layout:
##   rows 1-2: LRS (Nagelkerke + five PRS violin panels)
##   rows 3-4: NGS (Nagelkerke + five PRS violin panels)
##   3 columns x 4 rows = 12 panels
############################################################

cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(script_dir, "prs_common.R"))

base_dir <- Sys.getenv("RESULT_ROOT", file.path(dirname(script_dir), "results", "prs"))
n_pcs <- as.integer(Sys.getenv("N_PCS", "20"))

lrs_mark <- paste0("lrs_pc", n_pcs)
ngs_mark <- paste0("ngs_deepvar_pc", n_pcs, "_batch")
lrs_table_dir <- file.path(base_dir, "figure", lrs_mark, "tables")
ngs_table_dir <- file.path(base_dir, "figure", ngs_mark, "tables")
outdir <- file.path(base_dir, "figure", "combined_lrs_ngs_pc20")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

trait_order <- c("SCZ", "BIP", "ADHD", "MDD", "ASD")

read_dataset_tables <- function(table_dir, mark, dataset_label) {
  residual_file <- file.path(table_dir, paste0("merged_prs_residual_", mark, ".tsv"))
  result_file <- file.path(table_dir, paste0("continuous_or_auc_pseudor2_", mark, ".tsv"))
  if (!file.exists(residual_file)) stop("Missing residual table: ", residual_file)
  if (!file.exists(result_file)) stop("Missing PRS result table: ", result_file)

  dat <- fread(residual_file)
  prs_result <- fread(result_file)
  if (!all(c("sample", "type", "y", trait_order) %in% names(dat))) {
    stop("Residual table lacks required columns: ", residual_file)
  }
  if (!all(c("Trait", "Nagelkerke", "P", "OR_per_1SD_adjusted_residual", "CI_low", "CI_high") %in% names(prs_result))) {
    stop("PRS result table lacks required columns: ", result_file)
  }
  dat[, dataset := dataset_label]
  dat[, type := factor(as.character(type), levels = c("HC", "Case"))]
  prs_result[, dataset := dataset_label]
  prs_result[, Trait := factor(as.character(Trait), levels = trait_order)]
  list(dat = dat, result = prs_result, residual_file = residual_file, result_file = result_file)
}

format_p_label <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(p < 1e-4, paste0("P = ", formatC(p, format = "e", digits = 2)), paste0("P = ", signif(p, 3)))
  )
}

make_nagelkerke_panel <- function(prs_result, dataset_label) {
  plot_dt <- copy(prs_result)
  plot_dt[, Trait := factor(as.character(Trait), levels = trait_order)]
  ymax <- max(plot_dt$Nagelkerke, na.rm = TRUE) * 1.20
  if (!is.finite(ymax) || ymax <= 0) ymax <- 0.01
  ggplot(plot_dt, aes(x = Trait, y = Nagelkerke, fill = Trait)) +
    geom_col(width = 0.72, colour = "grey30", linewidth = 0.15) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * Nagelkerke)), vjust = -0.28, size = 2.6) +
    scale_fill_manual(values = trait_cols) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    labs(title = paste0(dataset_label, ": Nagelkerke pseudo-R2"), x = NULL, y = "Nagelkerke pseudo-R2") +
    theme_prs(9) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", size = 9)
    )
}

make_violin_panel <- function(dat, prs_result, trait, dataset_label) {
  plot_dt <- dat[, .(sample, type, y, PRS_residual = as.numeric(get(trait)))]
  plot_dt <- plot_dt[!is.na(PRS_residual) & !is.na(type)]
  stat_dt <- prs_result[as.character(Trait) == trait]
  if (!nrow(stat_dt)) stop("Missing PRS result row for ", dataset_label, " ", trait)
  yr <- range(plot_dt$PRS_residual, na.rm = TRUE)
  ypad <- diff(yr) * 0.10
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.25
  p_anno <- data.frame(
    group1 = "HC",
    group2 = "Case",
    y.position = yr[2] + ypad,
    p_label = format_p_label(stat_dt$P[1])
  )
  ggplot(plot_dt, aes(x = type, y = PRS_residual)) +
    geom_half_violin(aes(fill = type), alpha = 0.72, side = "r", trim = FALSE) +
    geom_half_boxplot(aes(fill = type), side = "r", errorbar.draw = FALSE, width = 0.20, outlier.shape = NA) +
    geom_half_point(aes(color = type), side = "l", alpha = 0.55, size = 0.42, show.legend = FALSE) +
    stat_pvalue_manual(
      p_anno,
      label = "p_label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.25,
      size = 2.6,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = case_control_cols) +
    scale_color_manual(values = case_control_cols) +
    labs(title = paste0(dataset_label, ": ", trait, " PRS"), x = NULL, y = "Adjusted PRS residual") +
    theme_prs(9) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 9)
    )
}

lrs <- read_dataset_tables(lrs_table_dir, lrs_mark, "LRS")
ngs <- read_dataset_tables(ngs_table_dir, ngs_mark, "NGS")

panel_list <- list(
  make_nagelkerke_panel(lrs$result, "LRS"),
  make_violin_panel(lrs$dat, lrs$result, "SCZ", "LRS"),
  make_violin_panel(lrs$dat, lrs$result, "BIP", "LRS"),
  make_violin_panel(lrs$dat, lrs$result, "ADHD", "LRS"),
  make_violin_panel(lrs$dat, lrs$result, "MDD", "LRS"),
  make_violin_panel(lrs$dat, lrs$result, "ASD", "LRS"),
  make_nagelkerke_panel(ngs$result, "NGS"),
  make_violin_panel(ngs$dat, ngs$result, "SCZ", "NGS"),
  make_violin_panel(ngs$dat, ngs$result, "BIP", "NGS"),
  make_violin_panel(ngs$dat, ngs$result, "ADHD", "NGS"),
  make_violin_panel(ngs$dat, ngs$result, "MDD", "NGS"),
  make_violin_panel(ngs$dat, ngs$result, "ASD", "NGS")
)

combined <- plot_grid(
  plotlist = panel_list,
  ncol = 3,
  labels = LETTERS[seq_along(panel_list)],
  label_size = 10,
  align = "hv",
  axis = "tblr"
)

out_prefix <- file.path(outdir, paste0("combined_nagelkerke_violin_lrs_ngs_pc", n_pcs))
ggsave(paste0(out_prefix, ".pdf"), combined, width = 24, height = 28, units = "cm", bg = "white")
ggsave(paste0(out_prefix, ".tiff"), combined, width = 24, height = 28, units = "cm",
       bg = "white", dpi = 300, compression = "lzw")

manifest <- data.table(
  dataset = c("LRS", "LRS", "NGS", "NGS"),
  file_type = c("residual", "continuous_prs_result", "residual", "continuous_prs_result"),
  file = c(lrs$residual_file, lrs$result_file, ngs$residual_file, ngs$result_file)
)
fwrite(manifest, paste0(out_prefix, ".input_manifest.tsv"), sep = "\t", quote = FALSE)

message("[DONE] Combined PRS panel PDF:  ", paste0(out_prefix, ".pdf"))
message("[DONE] Combined PRS panel TIFF: ", paste0(out_prefix, ".tiff"))
