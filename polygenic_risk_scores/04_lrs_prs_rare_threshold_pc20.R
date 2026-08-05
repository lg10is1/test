#!/usr/bin/env Rscript
############################################################
## LRS SCZ PRS vs rare burden threshold groups.
## Adjustment: 20 LRS projected PCs only.
## Samples listed in exclude_samples.tsv are removed before PC adjustment.
##
## For each threshold k in 1, 2, 3:
##   HC
##   Case rare burden < k
##   Case rare burden >= k
############################################################

cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(script_dir, "prs_common.R"))
suppressPackageStartupMessages(library(ggridges))

base_dir <- Sys.getenv("RESULT_ROOT", file.path(dirname(script_dir), "results", "prs"))
prs_file <- Sys.getenv("SCZ_PRS_FILE", file.path(base_dir, "SCZ/prscsx_eur.sscore"))
burden_file <- Sys.getenv("RARE_BURDEN_FILE", file.path(base_dir, "data/rare_variants_clean_burden.tsv"))
pc_file <- env_required("LRS_PC_FILE")
exclude_file <- env_required("LRS_EXCLUDE_FILE")
n_pcs <- as.integer(Sys.getenv("N_PCS", "20"))
thresholds <- as.integer(split_csv(Sys.getenv("RARE_THRESHOLDS", "1,2,3")))

mark <- paste0("lrs_scz_prs_pc", n_pcs, "_rare_thresholds")
outdir <- file.path(base_dir, "figure", "rare_threshold_lrs_pc20")
table_dir <- file.path(outdir, "tables")
single_plot_dir <- file.path(outdir, "threshold_panels")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(single_plot_dir, recursive = TRUE, showWarnings = FALSE)

read_exclude_samples <- function(file) {
  if (!file.exists(file)) stop("Missing LRS exclude sample file: ", file)
  x <- fread(file, header = "auto", fill = TRUE)
  sample_col <- if ("sample" %in% names(x)) "sample" else names(x)[1]
  out <- unique(as.character(x[[sample_col]]))
  out <- out[!is.na(out) & nzchar(out) & out != sample_col]
  out
}

prs <- read_prs(prs_file, "SCZ") %>%
  mutate(is_HC = grepl("^C", sample))
pc <- read_pc(pc_file, n_pcs = n_pcs, has_header = TRUE)
pc_cols <- paste0("PC", seq_len(n_pcs))

dat <- prs %>%
  inner_join(pc, by = "sample") %>%
  as.data.frame()
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
dat <- dat %>%
  filter(!sample %in% matched_exclude) %>%
  adjust_prs_by_covariates("SCZ", covar_cols = pc_cols, mark = paste0("pc", n_pcs))
used_covariates <- attr(dat, "adjustment_covariates")
dat$SCZ_PRS <- dat$SCZ

burden <- fread(burden_file)
if (!all(c("sample", "sum_rare") %in% names(burden))) {
  stop("Burden file must contain sample and sum_rare: ", burden_file)
}
for (cc in c("RareSV", "RareTR", "RareSNV")) {
  if (!cc %in% names(burden)) burden[, (cc) := 0L]
}
burden <- burden %>%
  transmute(
    sample = as.character(sample),
    RareSV = as.integer(RareSV),
    RareTR = as.integer(RareTR),
    RareSNV = as.integer(RareSNV),
    sum_rare = as.integer(sum_rare)
  )

dat <- dat %>%
  left_join(burden, by = "sample") %>%
  mutate(
    RareSV = ifelse(is.na(RareSV), 0L, RareSV),
    RareTR = ifelse(is.na(RareTR), 0L, RareTR),
    RareSNV = ifelse(is.na(RareSNV), 0L, RareSNV),
    sum_rare = ifelse(is.na(sum_rare), 0L, sum_rare)
  )

fwrite(
  dat,
  file.path(table_dir, paste0("scz_prs_pc_adjusted_with_burden_", mark, ".tsv")),
  sep = "\t",
  quote = FALSE
)

covar_manifest <- data.table(
  analysis = mark,
  adjustment = "LRS projected PCs only",
  n_pcs = n_pcs,
  pc_file = pc_file,
  burden_file = burden_file,
  exclude_file = exclude_file,
  covariate = used_covariates
)
fwrite(covar_manifest, file.path(table_dir, paste0("adjustment_covariates_", mark, ".tsv")), sep = "\t")

make_threshold_data <- function(k) {
  out <- dat %>%
    mutate(
      threshold = k,
      threshold_label = paste0("Threshold >= ", k, " rare"),
      threshold_group = case_when(
        is_HC ~ "HC",
        !is_HC & sum_rare < k ~ paste0("Case <", k),
        !is_HC & sum_rare >= k ~ paste0("Case >=", k)
      )
    )
  out$threshold_group <- factor(out$threshold_group, levels = c("HC", paste0("Case <", k), paste0("Case >=", k)))
  out
}

threshold_dat <- bind_rows(lapply(thresholds, make_threshold_data))
fwrite(
  threshold_dat,
  file.path(table_dir, paste0("threshold_grouped_data_", mark, ".tsv")),
  sep = "\t",
  quote = FALSE
)

group_summary <- threshold_dat %>%
  group_by(threshold, threshold_label, threshold_group) %>%
  summarise(
    n = n(),
    mean_SCZ_PRS = mean(SCZ_PRS, na.rm = TRUE),
    sd_SCZ_PRS = sd(SCZ_PRS, na.rm = TRUE),
    se_SCZ_PRS = sd_SCZ_PRS / sqrt(n),
    ci95_low = mean_SCZ_PRS - 1.96 * se_SCZ_PRS,
    ci95_high = mean_SCZ_PRS + 1.96 * se_SCZ_PRS,
    median_SCZ_PRS = median(SCZ_PRS, na.rm = TRUE),
    q25_SCZ_PRS = quantile(SCZ_PRS, 0.25, na.rm = TRUE),
    q75_SCZ_PRS = quantile(SCZ_PRS, 0.75, na.rm = TRUE),
    mean_sum_rare = mean(sum_rare, na.rm = TRUE),
    median_sum_rare = median(sum_rare, na.rm = TRUE),
    max_sum_rare = max(sum_rare, na.rm = TRUE),
    .groups = "drop"
  )
fwrite(group_summary, file.path(table_dir, paste0("group_summary_", mark, ".tsv")), sep = "\t", quote = FALSE)

omnibus_results <- bind_rows(lapply(thresholds, function(k) {
  dd <- threshold_dat %>% filter(threshold == k)
  aov_fit <- aov(SCZ_PRS ~ threshold_group, data = dd)
  kw <- kruskal.test(SCZ_PRS ~ threshold_group, data = dd)
  data.frame(
    threshold = k,
    threshold_label = paste0("Threshold >= ", k, " rare"),
    anova_p = summary(aov_fit)[[1]][["Pr(>F)"]][1],
    kruskal_chisq = as.numeric(kw$statistic),
    kruskal_df = as.numeric(kw$parameter),
    kruskal_p = kw$p.value
  )
}))
omnibus_results$anova_p_FDR <- p.adjust(omnibus_results$anova_p, method = "BH")
omnibus_results$kruskal_p_FDR <- p.adjust(omnibus_results$kruskal_p, method = "BH")
fwrite(omnibus_results, file.path(table_dir, paste0("omnibus_tests_", mark, ".tsv")), sep = "\t", quote = FALSE)

pairwise_results <- bind_rows(lapply(thresholds, function(k) {
  dd <- threshold_dat %>% filter(threshold == k)
  groups <- c("HC", paste0("Case <", k), paste0("Case >=", k))
  pairs <- list(c(groups[1], groups[2]), c(groups[1], groups[3]), c(groups[2], groups[3]))
  res <- bind_rows(lapply(pairs, function(z) pairwise_score_test(dd, "threshold_group", "SCZ_PRS", z[1], z[2])))
  res$threshold <- k
  res$threshold_label <- paste0("Threshold >= ", k, " rare")
  res
})) %>%
  group_by(threshold) %>%
  mutate(
    welch_p_FDR_within_threshold = p.adjust(welch_p, method = "BH"),
    wilcox_p_FDR_within_threshold = p.adjust(wilcox_p, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    welch_p_FDR_all = p.adjust(welch_p, method = "BH"),
    wilcox_p_FDR_all = p.adjust(wilcox_p, method = "BH")
  ) %>%
  select(threshold, threshold_label, everything())
fwrite(pairwise_results, file.path(table_dir, paste0("pairwise_welch_wilcox_effects_", mark, ".tsv")), sep = "\t", quote = FALSE)

format_p_label <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(p < 1e-4, "P < 1e-4", paste0("P = ", signif(p, 3)))
  )
}

p_annotation <- pairwise_results %>%
  group_by(threshold, threshold_label) %>%
  mutate(
    p_label = format_p_label(welch_p),
    y.position = max(threshold_dat$SCZ_PRS[threshold_dat$threshold == first(threshold)], na.rm = TRUE) +
      0.25 + 0.30 * (row_number() - 1L)
  ) %>%
  ungroup() %>%
  select(threshold, threshold_label, group1, group2, y.position, p_label)

case_cor_results <- bind_rows(lapply(thresholds, function(k) {
  dd <- threshold_dat %>% filter(threshold == k, !is_HC)
  sp <- cor.test(dd$SCZ_PRS, dd$sum_rare, method = "spearman", exact = FALSE)
  data.frame(
    threshold = k,
    threshold_label = paste0("Threshold >= ", k, " rare"),
    test = "Spearman correlation among cases",
    rho = as.numeric(sp$estimate),
    p_value = sp$p.value,
    n = nrow(dd)
  )
}))
case_cor_results$p_value_FDR <- p.adjust(case_cor_results$p_value, method = "BH")
fwrite(case_cor_results, file.path(table_dir, paste0("case_sum_rare_vs_scz_prs_spearman_", mark, ".tsv")), sep = "\t", quote = FALSE)

group_cols <- c(
  "HC" = "#da7271",
  "Case <1" = "#47a1a2",
  "Case >=1" = "#1f78b4",
  "Case <2" = "#47a1a2",
  "Case >=2" = "#1f78b4",
  "Case <3" = "#47a1a2",
  "Case >=3" = "#1f78b4"
)

p_violin <- ggplot(threshold_dat, aes(x = threshold_group, y = SCZ_PRS)) +
  geom_half_violin(aes(fill = threshold_group), side = "r", alpha = 0.7, trim = FALSE) +
  geom_half_boxplot(aes(fill = threshold_group), side = "r", width = 0.22, errorbar.draw = FALSE, outlier.shape = NA) +
  geom_half_point(aes(color = threshold_group), side = "l", alpha = 0.55, size = 0.65, show.legend = FALSE) +
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
  facet_wrap(~threshold_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(x = NULL, y = paste0("Standardized SCZ PRS residual after LRS ", n_pcs, " PCs")) +
  theme_prs(10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
save_plot_pair(p_violin, file.path(outdir, paste0("scz_prs_threshold_group_violin_", mark)), 30, 10)

p_mean <- ggplot(group_summary, aes(x = threshold_group, y = mean_SCZ_PRS)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high, color = threshold_group), width = 0.15, linewidth = 0.5) +
  geom_point(aes(color = threshold_group), size = 2.4) +
  facet_wrap(~threshold_label, nrow = 1, scales = "free_x") +
  scale_color_manual(values = group_cols) +
  labs(x = NULL, y = "Mean adjusted SCZ PRS residual +/- 95% CI") +
  theme_prs(10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
save_plot_pair(p_mean, file.path(outdir, paste0("scz_prs_threshold_group_mean_ci_", mark)), 30, 9)

p_ridge <- ggplot(threshold_dat, aes(x = SCZ_PRS, y = threshold_group, fill = threshold_group)) +
  geom_density_ridges(alpha = 0.8, scale = 1.1, rel_min_height = 0.01) +
  facet_wrap(~threshold_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  labs(x = paste0("Standardized SCZ PRS residual after LRS ", n_pcs, " PCs"), y = NULL) +
  theme_prs(10) +
  theme(legend.position = "none")
save_plot_pair(p_ridge, file.path(outdir, paste0("scz_prs_threshold_group_density_ridge_", mark)), 30, 9)

for (k in thresholds) {
  dd <- threshold_dat %>% filter(threshold == k)
  label <- paste0("threshold", k)
  p_anno_one <- p_annotation %>% filter(threshold == k)
  p_one <- ggplot(dd, aes(x = threshold_group, y = SCZ_PRS)) +
    geom_half_violin(aes(fill = threshold_group), side = "r", alpha = 0.7, trim = FALSE) +
    geom_half_boxplot(aes(fill = threshold_group), side = "r", width = 0.22, errorbar.draw = FALSE, outlier.shape = NA) +
    geom_half_point(aes(color = threshold_group), side = "l", alpha = 0.55, size = 0.75, show.legend = FALSE) +
    stat_pvalue_manual(
      p_anno_one,
      label = "p_label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.3,
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = group_cols) +
    scale_color_manual(values = group_cols) +
    labs(x = NULL, y = paste0("Standardized SCZ PRS residual after LRS ", n_pcs, " PCs"), title = paste0("Rare burden threshold >= ", k)) +
    theme_prs(10) +
    theme(legend.position = "none")
  save_plot_pair(p_one, file.path(single_plot_dir, paste0("scz_prs_rare_", label, "_violin_", mark)), 14, 10)
}

case_dat <- dat %>% filter(!is_HC)
p_burden <- ggplot(case_dat, aes(x = sum_rare, y = SCZ_PRS)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1.1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
  labs(
    x = "Rare variant burden count among cases",
    y = paste0("Standardized SCZ PRS residual after LRS ", n_pcs, " PCs")
  ) +
  theme_prs(10)
save_plot_pair(p_burden, file.path(outdir, paste0("case_sum_rare_vs_scz_prs_", mark)), 10, 8)

message("[DONE] LRS rare-threshold PRS analysis: ", outdir)
