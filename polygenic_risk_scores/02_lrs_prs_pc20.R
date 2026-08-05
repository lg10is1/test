#!/usr/bin/env Rscript
############################################################
## LRS previous PRS case-control analysis.
## Adjustment: 20 LRS projected PCs only.
## Samples listed in exclude_samples.tsv are removed before PC adjustment.
############################################################

cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(script_dir, "prs_common.R"))

base_dir <- Sys.getenv("RESULT_ROOT", file.path(dirname(script_dir), "results", "prs"))
pc_file <- env_required("LRS_PC_FILE")
exclude_file <- env_required("LRS_EXCLUDE_FILE")
n_pcs <- as.integer(Sys.getenv("N_PCS", "20"))
n_quantile <- as.integer(Sys.getenv("N_QUANTILE", "5"))

mark <- paste0("lrs_pc", n_pcs)
outdir <- file.path(base_dir, "figure", mark)
table_dir <- file.path(outdir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("SCZ", "BIP", "ADHD", "ASD", "MDD")
trait_order <- c("SCZ", "BIP", "ADHD", "MDD", "ASD")

score_files <- list(
  SCZ = file.path(base_dir, "SCZ/prscsx_eur.sscore"),
  BIP = file.path(base_dir, "BIP/prscsx_eur.sscore"),
  ADHD = file.path(base_dir, "ADHD/prscsx_eur.sscore"),
  ASD = file.path(base_dir, "ASD/prscsx_eur.sscore"),
  MDD = file.path(base_dir, "MDD/prscsx_eur.sscore")
)
missing_scores <- score_files[!file.exists(unlist(score_files))]
if (length(missing_scores)) stop("Missing PRS score file(s): ", paste(unlist(missing_scores), collapse = ", "))

read_exclude_samples <- function(file) {
  if (!file.exists(file)) stop("Missing LRS exclude sample file: ", file)
  x <- fread(file, header = "auto", fill = TRUE)
  sample_col <- if ("sample" %in% names(x)) "sample" else names(x)[1]
  out <- unique(as.character(x[[sample_col]]))
  out <- out[!is.na(out) & nzchar(out) & out != sample_col]
  out
}

prs <- Reduce(
  function(x, y) merge(x, y, by = "sample", all = FALSE),
  lapply(names(score_files), function(trait) read_prs(score_files[[trait]], trait))
)

pc <- read_pc(pc_file, n_pcs = n_pcs, has_header = TRUE)
pc_cols <- paste0("PC", seq_len(n_pcs))

dat <- merge(prs, pc, by = "sample", all = FALSE) %>% as.data.frame()
exclude_samples <- read_exclude_samples(exclude_file)
matched_exclude <- intersect(dat$sample, exclude_samples)
sample_filter_summary <- data.table(
  analysis = mark,
  exclude_file = exclude_file,
  n_exclude_file = length(exclude_samples),
  n_matched_removed = length(matched_exclude),
  n_before_remove = nrow(dat),
  n_after_remove = nrow(dat) - length(matched_exclude)
)
fwrite(sample_filter_summary, file.path(table_dir, paste0("excluded_samples_summary_", mark, ".tsv")), sep = "\t")
fwrite(
  data.table(sample = matched_exclude),
  file.path(table_dir, paste0("excluded_samples_applied_", mark, ".tsv")),
  sep = "\t",
  quote = FALSE
)
dat <- dat %>% filter(!sample %in% matched_exclude)

dat <- adjust_prs_by_covariates(dat, traits, covar_cols = pc_cols, mark = paste0("pc", n_pcs))
used_covariates <- attr(dat, "adjustment_covariates")
dat <- dat %>%
  mutate(
    type = ifelse(grepl("^C", sample), "HC", "Case"),
    type = factor(type, levels = c("HC", "Case")),
    y = ifelse(type == "Case", 1L, 0L)
  )

fwrite(
  dat,
  file.path(table_dir, paste0("merged_prs_residual_", mark, ".tsv")),
  sep = "\t",
  quote = FALSE
)

covar_manifest <- data.table(
  analysis = mark,
  adjustment = "LRS projected PCs only",
  n_pcs = n_pcs,
  pc_file = pc_file,
  exclude_file = exclude_file,
  covariate = used_covariates
)
fwrite(covar_manifest, file.path(table_dir, paste0("adjustment_covariates_", mark, ".tsv")), sep = "\t")

prs_long <- dat %>%
  pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "PRS_residual") %>%
  mutate(Trait = factor(Trait, levels = trait_order))

prs_result <- do.call(rbind, lapply(traits, function(x) calc_prs_result(dat, x)))
prs_result$Trait <- factor(prs_result$Trait, levels = trait_order)
fwrite(prs_result, file.path(table_dir, paste0("continuous_or_auc_pseudor2_", mark, ".tsv")), sep = "\t", quote = FALSE)

format_p_label <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(p < 1e-4, paste0("P = ", formatC(p, format = "e", digits = 2)), paste0("P = ", signif(p, 3)))
  )
}

p_annotation_y <- prs_long %>%
  group_by(Trait) %>%
  summarise(
    ymin = min(PRS_residual, na.rm = TRUE),
    ymax = max(PRS_residual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y.position = ymax + pmax(0.25, 0.10 * (ymax - ymin))) %>%
  select(Trait, y.position)

p_annotation <- prs_result %>%
  transmute(
    Trait = factor(as.character(Trait), levels = trait_order),
    group1 = "HC",
    group2 = "Case",
    p_label = format_p_label(P),
    logistic_P = P,
    OR_per_1SD_adjusted_residual = OR_per_1SD_adjusted_residual,
    CI_low = CI_low,
    CI_high = CI_high
  ) %>%
  left_join(p_annotation_y, by = "Trait")
fwrite(p_annotation, file.path(table_dir, paste0("case_control_logistic_p_annotation_", mark, ".tsv")), sep = "\t", quote = FALSE)

p_violin <- ggplot(prs_long, aes(x = type, y = PRS_residual)) +
  geom_half_violin(aes(fill = type), alpha = 0.7, side = "r", trim = FALSE) +
  geom_half_boxplot(aes(fill = type), side = "r", errorbar.draw = FALSE, width = 0.2, outlier.shape = NA) +
  geom_half_point(aes(color = type), side = "l", alpha = 0.55, size = 0.55, show.legend = FALSE) +
  stat_pvalue_manual(
    p_annotation,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = case_control_cols) +
  scale_color_manual(values = case_control_cols) +
  facet_wrap(~Trait, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = paste0("Standardized PRS residual after LRS ", n_pcs, " PCs")) +
  theme_prs(10) +
  theme(legend.position = "none")
save_plot_pair(p_violin, file.path(outdir, paste0("prs_case_control_", mark)), 31.5, 6.3)
p_nagelkerke <- ggplot(prs_result, aes(x = Trait, y = Nagelkerke)) +
  geom_col(aes(fill = Trait), width = 0.7) +
  scale_fill_manual(values = trait_cols) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Nagelkerke pseudo-R2") +
  theme_prs(10) +
  theme(legend.position = "none")
save_plot_pair(p_nagelkerke, file.path(outdir, paste0("nagelkerke_pseudor2_", mark)), 31.5, 6.3)

r2_long <- prs_result %>%
  select(Trait, McFadden, CoxSnell, Nagelkerke, Tjur) %>%
  pivot_longer(cols = c(McFadden, CoxSnell, Nagelkerke, Tjur), names_to = "Metric", values_to = "R2") %>%
  mutate(Metric = factor(Metric, levels = c("McFadden", "CoxSnell", "Nagelkerke", "Tjur")))

p_r2_all <- ggplot(r2_long, aes(x = Trait, y = R2)) +
  geom_col(aes(fill = Trait), width = 0.7) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = trait_cols) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Pseudo-R2") +
  theme_prs(10) +
  theme(legend.position = "none")
save_plot_pair(p_r2_all, file.path(outdir, paste0("all_pseudor2_metrics_", mark)), 31.5, 6.3)

quantile_result <- do.call(
  rbind,
  lapply(traits, function(x) make_quantile_result(dat, x, n_quantile = n_quantile, ref_quantile = "Q1"))
)
quantile_result$Trait <- factor(quantile_result$Trait, levels = trait_order)
fwrite(
  quantile_result,
  file.path(table_dir, paste0("quantile_or_Q", n_quantile, "_", mark, ".tsv")),
  sep = "\t",
  quote = FALSE
)

p_quantile_or <- ggplot(quantile_result, aes(x = Quantile_num, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high, color = Trait), width = 0.12, linewidth = 0.4) +
  geom_point(aes(color = Trait), size = 1.8) +
  facet_wrap(~Trait, nrow = 1, scales = "free_y") +
  scale_color_manual(values = trait_cols) +
  scale_x_continuous(breaks = seq_len(n_quantile), labels = paste0("Q", seq_len(n_quantile))) +
  scale_y_log10() +
  labs(x = "Adjusted PRS residual quantile", y = "Odds ratio vs Q1") +
  theme_prs(10) +
  theme(legend.position = "none")
save_plot_pair(p_quantile_or, file.path(outdir, paste0("quantile_or_Q", n_quantile, "_", mark)), 31.5, 6.3)

p_quantile_case <- ggplot(quantile_result, aes(x = Quantile, y = case_fraction)) +
  geom_col(aes(fill = Trait), width = 0.7) +
  geom_text(aes(label = paste0(n_case, "/", n)), vjust = -0.35, size = 2.5) +
  facet_wrap(~Trait, nrow = 1) +
  scale_fill_manual(values = trait_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Adjusted PRS residual quantile", y = "Case proportion") +
  theme_prs(10) +
  theme(legend.position = "none")
save_plot_pair(p_quantile_case, file.path(outdir, paste0("quantile_case_fraction_Q", n_quantile, "_", mark)), 31.5, 6.3)

message("[DONE] LRS PRS PC20 analysis: ", outdir)
