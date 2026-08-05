#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(qqman)
})

# ============================================================
# Config
# ============================================================

workdir <- "/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD"
outdir <- file.path(workdir, "LD_decay_public")

ids_dir <- file.path(outdir, "ids")
ld_dir  <- file.path(outdir, "ld")
rdata_dir <- file.path(outdir, "rdata")
fig_dir <- file.path(outdir, "figures")

metadata_file <- file.path(outdir, "ld_jobs.metadata.tsv")

max_dist <- 1000000
bin_size <- 1000
n_boot <- 500
seed <- 1L
xmax_decay <- 250000
max_dist_for_maxr2 <- 1000000

parse_seed_arg <- function(args, default = 1L) {
  seed <- as.integer(default)
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--seed") {
      i <- i + 1L
      if (i > length(args)) stop("--seed requires a value")
      seed <- suppressWarnings(as.integer(args[[i]]))
    } else {
      stop("Unknown argument: ", args[[i]])
    }
    i <- i + 1L
  }
  if (is.na(seed)) stop("--seed must be an integer")
  seed
}

seed <- parse_seed_arg(commandArgs(trailingOnly = TRUE), default = seed)

dir.create(rdata_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper functions
# ============================================================

check_file <- function(x) {
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }
  invisible(TRUE)
}

read_index_ids <- function(id_file) {
  check_file(id_file)
  raw <- readLines(id_file, warn = FALSE)
  raw <- trimws(raw)
  raw <- raw[!is.na(raw) & raw != ""]
  if (!length(raw)) {
    return(character())
  }
  parts <- strsplit(raw, "[[:space:]]+")
  ids <- vapply(parts, function(x) x[[1L]], character(1L))
  ids <- unique(as.character(ids))
  ids[!is.na(ids) & ids != ""]
}

find_first_col <- function(candidates, x) {
  out <- intersect(candidates, colnames(x))
  if (length(out) == 0) {
    return(NA_character_)
  }
  out[1]
}
summarise_one_ld_file <- function(
  ld_file,
  id_file,
  source_id,
  plot_group,
  max_dist = 1000000,
  bin_size = 1000
) {
  check_file(ld_file)
  index_ids <- read_index_ids(id_file)

  message("[READ] ", basename(ld_file), " | group=", plot_group, " | source=", source_id)

  ld <- fread(ld_file)

  if (nrow(ld) == 0) {
    warning("Empty LD file: ", ld_file)
    return(NULL)
  }

  r2_col <- find_first_col(c("R2", "R^2", "UNPHASED_R2", "PHASED_R2"), ld)
  if (is.na(r2_col)) {
    stop("Cannot find R2 column in: ", ld_file,
         "\nColumns are: ", paste(colnames(ld), collapse = ", "))
  }

  required_cols <- c("SNP_A", "SNP_B", "BP_A", "BP_B")
  missing_cols <- setdiff(required_cols, colnames(ld))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in LD file: ", ld_file,
      "\nMissing: ", paste(missing_cols, collapse = ", "),
      "\nColumns are: ", paste(colnames(ld), collapse = ", ")
    )
  }

  ld[, SNP_A := as.character(SNP_A)]
  ld[, SNP_B := as.character(SNP_B)]

  ld[, is_a_index := SNP_A %in% index_ids]
  ld[, is_b_index := SNP_B %in% index_ids]
  ld <- ld[is_a_index | is_b_index]

  if (nrow(ld) == 0) {
    warning("No rows contain index variants in LD file: ", ld_file)
    return(NULL)
  }

  # Robust orientation:
  # Usually --ld-snp-list makes SNP_A the index, but this also handles SNP_B.
  ld[, index := fifelse(is_a_index, SNP_A, SNP_B)]
  ld[, target := fifelse(is_a_index, SNP_B, SNP_A)]
  ld[, bp_index := fifelse(is_a_index, BP_A, BP_B)]
  ld[, bp_target := fifelse(is_a_index, BP_B, BP_A)]

  ld[, dis := abs(as.numeric(bp_target) - as.numeric(bp_index))]
  ld[, ld_r2 := as.numeric(get(r2_col))]

  #     self-LD     0-1 kb bin          ?bin=0, LD=1
  ld <- ld[dis > 0 & dis <= max_dist]

  if (nrow(ld) == 0) {
    warning("No LD pairs after distance filtering: ", ld_file)
    return(NULL)
  }

  ld[, bin := ceiling(dis / bin_size) * bin_size]
  ld[, index_key := paste(source_id, index, sep = "::")]

  #           source_id  ?plot_group        ?
  source_id_value <- source_id
  plot_group_value <- plot_group

  ld[, source_id := source_id_value]
  ld[, plot_group := plot_group_value]

  per_index_bin <- ld[, .(
    mean_ld_index = mean(ld_r2, na.rm = TRUE),
    n_pair = .N
  ), by = .(plot_group, source_id, index, index_key, bin)]

  per_index_bin
}

read_ld_for_maxr2 <- function(ld_file, id_file, source_id, plot_group, max_dist = 1000000) {
  check_file(ld_file)
  check_file(id_file)
  
  index_ids <- read_index_ids(id_file)
  
  message("[READ maxR2] ", basename(ld_file), " | ", source_id, " | ", plot_group)
  
  ld <- fread(ld_file)
  
  if (nrow(ld) == 0) {
    warning("Empty LD file: ", ld_file)
    return(NULL)
  }
  
  r2_col <- find_first_col(c("R2", "R^2", "UNPHASED_R2", "PHASED_R2"), ld)
  if (is.na(r2_col)) {
    stop("Cannot find R2 column in: ", ld_file,
         "\nColumns: ", paste(colnames(ld), collapse = ", "))
  }
  
  required_cols <- c("SNP_A", "SNP_B", "BP_A", "BP_B")
  missing_cols <- setdiff(required_cols, colnames(ld))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", ld_file, ": ",
         paste(missing_cols, collapse = ", "))
  }
  
  ld[, SNP_A := as.character(SNP_A)]
  ld[, SNP_B := as.character(SNP_B)]
  
  ld[, is_a_index := SNP_A %in% index_ids]
  ld[, is_b_index := SNP_B %in% index_ids]
  ld <- ld[is_a_index | is_b_index]
  
  if (nrow(ld) == 0) {
    warning("No LD rows contain index variants: ", ld_file)
    return(NULL)
  }
  
  ld[, index := fifelse(is_a_index, SNP_A, SNP_B)]
  ld[, target := fifelse(is_a_index, SNP_B, SNP_A)]
  ld[, bp_index := fifelse(is_a_index, BP_A, BP_B)]
  ld[, bp_target := fifelse(is_a_index, BP_B, BP_A)]
  
  ld[, dis := abs(as.numeric(bp_target) - as.numeric(bp_index))]
  ld[, ld_r2 := as.numeric(get(r2_col))]
  
  ld <- ld[dis > 0 & dis <= max_dist]
  
  if (nrow(ld) == 0) {
    warning("No LD rows after distance filtering: ", ld_file)
    return(NULL)
  }
  
  ld[, source_id := source_id]
  ld[, plot_group := plot_group]
  ld[, index_key := paste(source_id, index, sep = "::")]
  
  out <- ld[, .(
    max_r2 = max(ld_r2, na.rm = TRUE),
    n_pair = .N
  ), by = .(plot_group, source_id, index, index_key)]
  
  out
}

make_observed_decay <- function(per_index_bin) {
  dt <- as.data.table(per_index_bin)

  obs <- dt[, .(
    ld = mean(mean_ld_index, na.rm = TRUE),
    sd = sd(mean_ld_index, na.rm = TRUE),
    n_index = uniqueN(index_key),
    n_pair_bin = sum(n_pair, na.rm = TRUE),
    se = sd(mean_ld_index, na.rm = TRUE) / sqrt(uniqueN(index_key))
  ), by = .(plot_group, bin)]

  start <- dt[, .(
    bin = 0,
    ld = 1,
    sd = 0,
    n_index = uniqueN(index_key),
    n_pair_bin = NA_real_,
    se = 0
  ), by = .(plot_group)]

  out <- rbindlist(list(start, obs), fill = TRUE)
  setorder(out, plot_group, bin)
  out
}

bootstrap_decay <- function(per_index_bin, n_boot = 500, seed = 1) {
  set.seed(seed)

  dt <- as.data.table(per_index_bin)
  groups <- sort(unique(dt$plot_group))

  boot_results <- vector("list", length(groups) * n_boot)
  k <- 1

  for (g in groups) {
    dt_g <- dt[plot_group == g]
    ids <- unique(dt_g$index_key)

    message("[BOOT] ", g, " | n_index=", length(ids), " | n_boot=", n_boot)

    if (length(ids) < 2) {
      warning("[BOOT] group ", g, " has fewer than 2 index variants; CI may be NA.")
    }

    for (b in seq_len(n_boot)) {
      sampled_ids <- sample(ids, length(ids), replace = TRUE)
      weight_dt <- data.table(index_key = sampled_ids)[, .(w = .N), by = index_key]

      x <- merge(dt_g, weight_dt, by = "index_key", allow.cartesian = TRUE)

      z <- x[, .(
        ld = sum(mean_ld_index * w, na.rm = TRUE) / sum(w[!is.na(mean_ld_index)])
      ), by = .(plot_group, bin)]

      z[, boot := b]
      boot_results[[k]] <- z
      k <- k + 1
    }
  }

  boot_all <- rbindlist(boot_results, fill = TRUE)

  ci <- boot_all[, .(
    ld_boot_mean = mean(ld, na.rm = TRUE),
    lower = quantile(ld, probs = 0.025, na.rm = TRUE),
    upper = quantile(ld, probs = 0.975, na.rm = TRUE)
  ), by = .(plot_group, bin)]

  start <- dt[, .(
    bin = 0,
    ld_boot_mean = 1,
    lower = 1,
    upper = 1
  ), by = .(plot_group)]

  ci <- rbindlist(list(start, ci), fill = TRUE)
  setorder(ci, plot_group, bin)

  list(
    boot_all = boot_all,
    ci = ci
  )
}

calc_auc_on_grid <- function(x, y, max_dist = 250000, bin_size = 1000) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- !is.na(x) & !is.na(y) & x >= 0 & x <= max_dist
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2) {
    return(NA_real_)
  }

  dt <- data.table(x = x, y = y)
  dt <- dt[, .(y = mean(y, na.rm = TRUE)), by = x]
  setorder(dt, x)

  grid <- seq(0, max_dist, by = bin_size)
  interp <- approx(
    x = dt$x,
    y = dt$y,
    xout = grid,
    rule = 2,
    ties = mean
  )

  sum(diff(interp$x) * (head(interp$y, -1) + tail(interp$y, -1)) / 2)
}

summarise_auc <- function(plot_dt, boot_all, max_dist = 250000, bin_size = 1000) {
  plot_dt <- as.data.table(plot_dt)
  boot_all <- as.data.table(boot_all)

  observed_auc <- plot_dt[, .(
    auc = calc_auc_on_grid(bin, ld, max_dist = max_dist, bin_size = bin_size)
  ), by = plot_group]

  observed_auc[, auc_mean_ld := auc / max_dist]
  observed_auc[, max_dist := max_dist]

  boot_start <- unique(boot_all[, .(plot_group, boot)])
  boot_start[, `:=`(
    bin = 0,
    ld = 1
  )]

  boot_auc_input <- rbindlist(
    list(
      boot_all[, .(plot_group, boot, bin, ld)],
      boot_start[, .(plot_group, boot, bin, ld)]
    ),
    fill = TRUE
  )

  bootstrap_auc <- boot_auc_input[, .(
    auc = calc_auc_on_grid(bin, ld, max_dist = max_dist, bin_size = bin_size)
  ), by = .(plot_group, boot)]

  bootstrap_auc[, auc_mean_ld := auc / max_dist]
  bootstrap_auc[, max_dist := max_dist]

  pair_list <- list(
    c("SV_sig", "SV_null"),
    c("SNV_INDEL_sig", "SNV_INDEL_null"),
    c("SV_sig", "SNV_INDEL_sig"),
    c("SV_null", "SNV_INDEL_null")
  )

  observed_wide <- dcast(
    observed_auc,
    max_dist ~ plot_group,
    value.var = "auc_mean_ld"
  )

  bootstrap_wide <- dcast(
    bootstrap_auc,
    boot + max_dist ~ plot_group,
    value.var = "auc_mean_ld"
  )

  pairwise_auc_diff <- rbindlist(lapply(pair_list, function(pair) {
    a <- pair[1]
    b <- pair[2]

    if (!all(c(a, b) %in% colnames(bootstrap_wide))) {
      return(NULL)
    }

    boot_diff <- bootstrap_wide[[a]] - bootstrap_wide[[b]]
    obs_diff <- if (all(c(a, b) %in% colnames(observed_wide))) {
      observed_wide[[a]][1] - observed_wide[[b]][1]
    } else {
      NA_real_
    }

    p_two_sided <- 2 * min(
      mean(boot_diff <= 0, na.rm = TRUE),
      mean(boot_diff >= 0, na.rm = TRUE)
    )
    p_two_sided <- min(p_two_sided, 1)

    data.table(
      group_a = a,
      group_b = b,
      statistic = "AUC_mean_LD_0_250kb_difference",
      interpretation = "positive means group_a has higher LD / slower decay than group_b",
      observed_diff = obs_diff,
      boot_mean_diff = mean(boot_diff, na.rm = TRUE),
      boot_median_diff = median(boot_diff, na.rm = TRUE),
      ci_lower = quantile(boot_diff, 0.025, na.rm = TRUE),
      ci_upper = quantile(boot_diff, 0.975, na.rm = TRUE),
      p_two_sided = p_two_sided,
      n_boot = sum(!is.na(boot_diff)),
      max_dist = max_dist
    )
  }), fill = TRUE)

  list(
    observed_auc = observed_auc,
    bootstrap_auc = bootstrap_auc,
    pairwise_auc_diff = pairwise_auc_diff
  )
}

# ============================================================
# Main
# ============================================================

check_file(metadata_file)
metadata <- fread(metadata_file)

required_meta_cols <- c("source_id", "bfile", "plot_group", "id_file", "out_prefix", "n_index")
missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata))
if (length(missing_meta_cols) > 0) {
  stop("Metadata missing columns: ", paste(missing_meta_cols, collapse = ", "))
}

per_file_list <- list()

for (i in seq_len(nrow(metadata))) {
  out_prefix <- metadata$out_prefix[i]

  # PLINK --r2 gz normally writes .ld.gz.
  ld_file <- paste0(out_prefix, ".ld.gz")

  if (!file.exists(ld_file)) {
    warning("LD file does not exist, skip: ", ld_file)
    next
  }

  tmp <- summarise_one_ld_file(
    ld_file = ld_file,
    id_file = metadata$id_file[i],
    source_id = metadata$source_id[i],
    plot_group = metadata$plot_group[i],
    max_dist = max_dist,
    bin_size = bin_size
  )

  if (!is.null(tmp) && nrow(tmp) > 0) {
    per_file_list[[length(per_file_list) + 1]] <- tmp
  }
}

if (length(per_file_list) == 0) {
  stop("No LD data were loaded. Check PLINK outputs.")
}

per_index_bin <- rbindlist(per_file_list, fill = TRUE)

# Keep final groups in desired order.
group_levels <- c("SV_sig", "SNV_INDEL_sig", "SV_null", "SNV_INDEL_null")
per_index_bin[, plot_group := factor(plot_group, levels = group_levels)]

observed_decay <- make_observed_decay(per_index_bin)
boot <- bootstrap_decay(per_index_bin, n_boot = n_boot, seed = seed)
decay_ci <- boot$ci

observed_decay[, plot_group := factor(plot_group, levels = group_levels)]
decay_ci[, plot_group := factor(plot_group, levels = group_levels)]

plot_dt <- merge(
  observed_decay,
  decay_ci,
  by = c("plot_group", "bin"),
  all.x = TRUE
)

# Use observed line, bootstrap CI ribbon.
# If you prefer bootstrap mean line, change y=ld to y=ld_boot_mean.
plot_dt <- plot_dt[order(plot_group, bin)]

auc_results <- summarise_auc(
  plot_dt = plot_dt,
  boot_all = boot$boot_all,
  max_dist = xmax_decay,
  bin_size = bin_size
)

observed_auc <- auc_results$observed_auc
bootstrap_auc <- auc_results$bootstrap_auc
pairwise_auc_diff <- auc_results$pairwise_auc_diff

save(
  per_index_bin,
  observed_decay,
  decay_ci,
  plot_dt,
  observed_auc,
  bootstrap_auc,
  pairwise_auc_diff,
  file = file.path(rdata_dir, "ld_decay_v2_summary.Rdata")
)

fwrite(per_index_bin, file.path(rdata_dir, "per_index_bin.tsv.gz"), sep = "\t")
fwrite(observed_decay, file.path(rdata_dir, "observed_decay.tsv"), sep = "\t")
fwrite(decay_ci, file.path(rdata_dir, "bootstrap_ci.tsv"), sep = "\t")
fwrite(plot_dt, file.path(rdata_dir, "plot_data.tsv"), sep = "\t")
fwrite(observed_auc, file.path(rdata_dir, "ld_decay_auc_observed.tsv"), sep = "\t")
fwrite(bootstrap_auc, file.path(rdata_dir, "ld_decay_auc_bootstrap.tsv.gz"), sep = "\t")
fwrite(pairwise_auc_diff, file.path(rdata_dir, "ld_decay_auc_pairwise_diff.tsv"), sep = "\t")
fwrite(
  data.table(parameter = c("bootstrap_seed", "bootstrap_replicates"), value = c(seed, n_boot)),
  file.path(rdata_dir, "ld_decay_bootstrap_config.tsv"),
  sep = "\t"
)

# ============================================================
# Part 1. Old-style LD decay plot, 0-250 kb
# ============================================================

plot_dt2 <- as.data.table(plot_dt)

#     group                  
plot_dt2[, type := as.character(plot_group)]
plot_dt2[type == "SV_sig", type := "sv_sig"]
plot_dt2[type == "SNV_INDEL_sig", type := "snv_sig"]
plot_dt2[type == "SV_null", type := "sv_null"]
plot_dt2[type == "SNV_INDEL_null", type := "snv_null"]

type_levels <- c("sv_null", "snv_null", "snv_sig", "sv_sig")
plot_dt2[, type := factor(type, levels = type_levels)]

#        ?# scale_color_manual(values=c('#47a1a2',"#f9bfcb",'#1f78b4','#da7271'))
color_values_old <- c(
  "sv_null" = "#47a1a2",
  "snv_null" = "#f9bfcb",
  "snv_sig" = "#1f78b4",
  "sv_sig" = "#da7271"
)

gg_decay <- ggplot(plot_dt2, aes(x = bin, y = ld)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, group = type),
    fill = "grey70",
    alpha = 0.18,
    color = NA
  ) +
  geom_line(
    aes(y = lower, group = type),
    color = "grey70",
    linewidth = 0.20,
    alpha = 0.55
  ) +
  geom_line(
    aes(y = upper, group = type),
    color = "grey70",
    linewidth = 0.20,
    alpha = 0.55
  ) +
  geom_line(aes(color = type), linewidth = 0.35) +
  theme_cowplot(8) +
  scale_x_continuous(
    breaks = seq(0, xmax_decay, length.out = 5),
    labels = round(seq(0, xmax_decay, length.out = 5) / 1000, 0),
    limits = c(0, xmax_decay),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00)
  ) +
  labs(
    x = "Distance (KB)",
    y = expression(LD~(R^2)),
    color = NULL
  ) +
  scale_color_manual(values = color_values_old, drop = FALSE) +
  theme_cowplot(8) +
  theme(
    legend.position = c(0.70, 0.70),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.30, "cm"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )

ggsave(
  gg_decay,
  file = file.path(fig_dir, "f3b_2_v2.svg"),
  width = 6,
  height = 5,
  units = "cm"
)

ggsave(
  gg_decay,
  file = file.path(fig_dir, "f3b_2_v2.pdf"),
  width = 6,
  height = 5,
  units = "cm"
)

ggsave(
  gg_decay,
  file = file.path(fig_dir, "f3b_2_v2.tiff"),
  width = 6,
  height = 5,
  units = "cm",
  dpi = 300,
  bg = "white"
)

message("[DONE] LD decay old-style plot:")
message("  ", file.path(fig_dir, "f3b_2_v2.svg"))
message("[DONE] LD decay AUC statistics:")
message("  ", file.path(rdata_dir, "ld_decay_auc_observed.tsv"))
message("  ", file.path(rdata_dir, "ld_decay_auc_bootstrap.tsv.gz"))
message("  ", file.path(rdata_dir, "ld_decay_auc_pairwise_diff.tsv"))
print(pairwise_auc_diff)

# ============================================================
# Part 2. SV_sig vs SV_null max R2 histogram
# ============================================================

metadata_sv <- metadata[plot_group %in% c("SV_sig", "SV_null")]

if (nrow(metadata_sv) == 0) {
  stop("No SV_sig / SV_null rows found in metadata: ", metadata_file)
}

maxr2_list <- list()

for (i in seq_len(nrow(metadata_sv))) {
  ld_file <- paste0(metadata_sv$out_prefix[i], ".ld.gz")
  
  if (!file.exists(ld_file)) {
    warning("LD file missing, skip: ", ld_file)
    next
  }
  
  tmp <- read_ld_for_maxr2(
    ld_file = ld_file,
    id_file = metadata_sv$id_file[i],
    source_id = metadata_sv$source_id[i],
    plot_group = metadata_sv$plot_group[i],
    max_dist = max_dist_for_maxr2
  )
  
  if (!is.null(tmp) && nrow(tmp) > 0) {
    maxr2_list[[length(maxr2_list) + 1]] <- tmp
  }
}

if (length(maxr2_list) == 0) {
  stop("No maxR2 data loaded.")
}

sv_maxr2 <- rbindlist(maxr2_list, fill = TRUE)

sv_maxr2[, type := fifelse(plot_group == "SV_sig", "sv_sig", "sv_null")]
sv_maxr2[, type := factor(type, levels = c("sv_null", "sv_sig"))]

fwrite(
  sv_maxr2,
  file.path(rdata_dir, "sv_sig_vs_null.maxR2.tsv"),
  sep = "\t",
  quote = FALSE
)

save(
  sv_maxr2,
  file = file.path(rdata_dir, "sv_sig_vs_null.maxR2.Rdata")
)

gg_maxr2_hist <- ggplot(sv_maxr2, aes(x = max_r2, fill = type)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 10,
    position = "identity",
    alpha = 0.65,
    boundary = 0
  ) +
  theme_cowplot(8) +
  scale_fill_manual(
    values = c(
      "sv_null" = "#708090",
      "sv_sig" = "#da7271"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 1.01),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = expression(Max~R^2),
    y = "Density",
    fill = NULL
  ) +
  theme(
    legend.position = c(0.72, 0.78),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.30, "cm"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )

ggsave(
  gg_maxr2_hist,
  file = file.path(fig_dir, "sv_maxR2_pair_observed_hist_v2.svg"),
  width = 5,
  height = 5,
  units = "cm"
)

ggsave(
  gg_maxr2_hist,
  file = file.path(fig_dir, "sv_maxR2_pair_observed_hist_v2.pdf"),
  width = 5,
  height = 5,
  units = "cm"
)

ggsave(
  gg_maxr2_hist,
  file = file.path(fig_dir, "sv_maxR2_pair_observed_hist_v2.tiff"),
  width = 5,
  height = 5,
  units = "cm",
  dpi = 300,
  bg = "white"
)

message("[DONE] SV maxR2 histogram:")
message("  ", file.path(fig_dir, "sv_maxR2_pair_observed_hist_v2.svg"))

# ============================================================
# Optional: density plot version
# ============================================================

gg_maxr2_density <- ggplot(sv_maxr2, aes(x = max_r2, fill = type, color = type)) +
  geom_density(alpha = 0.25, linewidth = 0.35) +
  theme_cowplot(8) +
  scale_fill_manual(
    values = c(
      "sv_null" = "#708090",
      "sv_sig" = "#da7271"
    ),
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(
      "sv_null" = "#708090",
      "sv_sig" = "#da7271"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 1.01),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = expression(Max~R^2),
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme(
    legend.position = c(0.72, 0.78),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.30, "cm"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )

ggsave(
  gg_maxr2_density,
  file = file.path(fig_dir, "sv_maxR2_pair_observed_density_v2.svg"),
  width = 5,
  height = 5,
  units = "cm"
)

ggsave(
  gg_maxr2_density,
  file = file.path(fig_dir, "sv_maxR2_pair_observed_density_v2.pdf"),
  width = 5,
  height = 5,
  units = "cm"
)

message("[DONE] SV maxR2 density:")
message("  ", file.path(fig_dir, "sv_maxR2_pair_observed_density_v2.svg"))

# ============================================================
# Summary
# ============================================================

message("\nSummary of max R2:")
print(
  sv_maxr2[, .(
    n_index = .N,
    mean_max_r2 = mean(max_r2, na.rm = TRUE),
    median_max_r2 = median(max_r2, na.rm = TRUE),
    q25 = quantile(max_r2, 0.25, na.rm = TRUE),
    q75 = quantile(max_r2, 0.75, na.rm = TRUE)
  ), by = type]
)

message("\nAll done.")

