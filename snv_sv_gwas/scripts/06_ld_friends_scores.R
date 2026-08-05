suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# LD score / LD friends comparison for the four existing LD-decay groups.
workdir <- "/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD"
ld_decay_dir <- file.path(workdir, "LD_decay_public")
metadata_file <- file.path(ld_decay_dir, "ld_jobs.metadata.tsv")
analysis_dir <- file.path(ld_decay_dir, "ld_friends_scores")
figure_dir <- file.path(ld_decay_dir, "figures")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse the existing 1 Mb PLINK LD files and report both the original 1 Mb
# analysis and a 250 kb sensitivity analysis matching the LD-decay figure.
analysis_windows_bp <- c(1000000L, 250000L)
max_input_dist_bp <- max(analysis_windows_bp)
ld_friend_p_threshold <- 0.05

group_levels <- c(
  "SNV_INDEL_null",
  "SNV_INDEL_sig",
  "SV_null",
  "SV_sig"
)
partner_scope_levels <- c("ALL", "SNV_INDEL", "SV")
group_colors <- c(
  "SNV_INDEL_null" = "#f9bfcb",
  "SNV_INDEL_sig" = "#1f78b4",
  "SV_null" = "#47a1a2",
  "SV_sig" = "#da7271"
)

check_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Required file is missing or empty: ", path)
  }
  invisible(TRUE)
}

read_ids <- function(path) {
  check_file(path)
  unique(as.character(fread(path, header = FALSE)[[1]]))
}

read_bim_classes <- function(bfile) {
  bim_file <- paste0(bfile, ".bim")
  check_file(bim_file)
  bim <- fread(
    bim_file,
    header = FALSE,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2")
  )
  bim[, `:=`(
    SNP = as.character(SNP),
    A1 = as.character(A1),
    A2 = as.character(A2)
  )]
  bim[, max_allele_len := pmax(nchar(A1), nchar(A2))]
  bim[, id_has_sv_tag := grepl("(^|[-_])(INS|DEL|DUP|INV)([-_]|$)", SNP)]
  bim[, partner_class := fifelse(
    max_allele_len >= 50 | id_has_sv_tag,
    "SV",
    "SNV_INDEL"
  )]
  unique(bim[, .(SNP, partner_class)])
}

read_sample_size <- function(bfile) {
  fam_file <- paste0(bfile, ".fam")
  check_file(fam_file)
  n <- nrow(fread(fam_file, header = FALSE, select = 1))
  if (is.na(n) || n <= 2) {
    stop("Invalid sample size from FAM file: ", fam_file, " (N=", n, ")")
  }
  n
}

find_r2_column <- function(x) {
  found <- intersect(
    c("R2", "R^2", "UNPHASED_R2", "PHASED_R2"),
    names(x)
  )
  if (length(found) == 0) {
    stop("No R2 column found. Columns: ", paste(names(x), collapse = ", "))
  }
  found[1]
}

pearson_p_from_r2 <- function(r2, n) {
  r2 <- pmin(pmax(as.numeric(r2), 0), 1)
  out <- numeric(length(r2))
  perfect <- r2 >= 1
  out[perfect] <- 0
  if (any(!perfect)) {
    t_stat <- sqrt(r2[!perfect] * (n - 2) / (1 - r2[!perfect]))
    out[!perfect] <- 2 * pt(-abs(t_stat), df = n - 2)
  }
  out
}

read_target_pairs <- function(meta_row) {
  source_id_value <- as.character(meta_row$source_id)
  plot_group_value <- as.character(meta_row$plot_group)
  bfile <- as.character(meta_row$bfile)
  id_file <- as.character(meta_row$id_file)
  ld_file <- paste0(as.character(meta_row$out_prefix), ".ld.gz")

  check_file(ld_file)
  target_ids <- read_ids(id_file)
  variant_classes <- read_bim_classes(bfile)
  sample_n_value <- read_sample_size(bfile)

  message(
    "[READ] ", basename(ld_file),
    " | source=", source_id_value,
    " | group=", plot_group_value,
    " | targets=", length(target_ids),
    " | N=", sample_n_value
  )

  targets <- data.table(
    source_id = source_id_value,
    plot_group = plot_group_value,
    index = target_ids,
    sample_n = sample_n_value
  )

  ld <- fread(ld_file)
  if (nrow(ld) == 0) {
    warning("Empty LD file: ", ld_file)
    return(list(targets = targets, pairs = NULL))
  }

  required <- c("SNP_A", "SNP_B", "BP_A", "BP_B")
  missing <- setdiff(required, names(ld))
  if (length(missing) > 0) {
    stop("Missing columns in ", ld_file, ": ", paste(missing, collapse = ", "))
  }

  r2_col <- find_r2_column(ld)
  ld[, `:=`(
    SNP_A = as.character(SNP_A),
    SNP_B = as.character(SNP_B),
    BP_A = as.numeric(BP_A),
    BP_B = as.numeric(BP_B),
    ld_r2 = as.numeric(get(r2_col))
  )]

  # Expand both orientations. This preserves both targets when a pair contains
  # two variants from the target list.
  from_a <- ld[SNP_A %in% target_ids, .(
    index = SNP_A,
    partner = SNP_B,
    bp_index = BP_A,
    bp_partner = BP_B,
    ld_r2
  )]
  from_b <- ld[SNP_B %in% target_ids, .(
    index = SNP_B,
    partner = SNP_A,
    bp_index = BP_B,
    bp_partner = BP_A,
    ld_r2
  )]

  pairs <- rbindlist(list(from_a, from_b), use.names = TRUE)
  pairs[, distance_bp := abs(bp_partner - bp_index)]
  pairs <- pairs[
    index != partner &
      distance_bp > 0 &
      distance_bp <= max_input_dist_bp &
      !is.na(ld_r2)
  ]
  setorder(pairs, index, partner, -ld_r2)
  pairs <- unique(pairs, by = c("index", "partner"))

  pairs <- merge(
    pairs,
    variant_classes,
    by.x = "partner",
    by.y = "SNP",
    all.x = TRUE
  )

  if (anyNA(pairs$partner_class)) {
    missing_partner <- unique(pairs[is.na(partner_class), partner])
    stop(
      "Partner variants not found in BIM for source ", source_id_value,
      ". Examples: ", paste(head(missing_partner, 10), collapse = ", ")
    )
  }

  pairs[, `:=`(
    source_id = source_id_value,
    plot_group = plot_group_value,
    sample_n = sample_n_value,
    ld_p = pearson_p_from_r2(ld_r2, sample_n_value)
  )]

  list(targets = targets, pairs = pairs)
}

summarise_scope <- function(targets, pairs, scope_value, window_bp_value) {
  scope_pairs <- if (scope_value == "ALL") {
    if (is.null(pairs)) NULL else pairs[distance_bp <= window_bp_value]
  } else if (is.null(pairs)) {
    NULL
  } else {
    pairs[
      partner_class == scope_value &
        distance_bp <= window_bp_value
    ]
  }

  if (is.null(scope_pairs) || nrow(scope_pairs) == 0) {
    return(targets[, .(
      source_id,
      plot_group,
      index,
      sample_n,
      window_bp = window_bp_value,
      partner_scope = scope_value,
      n_partners = 0L,
      sum_r2 = 0,
      ld_score = 1,
      n_ld_friends = 0L,
      friends_rate = NA_real_,
      mean_r2 = NA_real_
    )])
  }

  observed <- scope_pairs[, .(
    n_partners = .N,
    sum_r2 = sum(ld_r2, na.rm = TRUE),
    n_ld_friends = sum(ld_p < ld_friend_p_threshold, na.rm = TRUE)
  ), by = .(source_id, plot_group, index, sample_n)]
  observed[, `:=`(
    window_bp = window_bp_value,
    partner_scope = scope_value,
    ld_score = 1 + sum_r2,
    friends_rate = n_ld_friends / n_partners,
    mean_r2 = sum_r2 / n_partners
  )]

  out <- merge(
    targets,
    observed,
    by = c("source_id", "plot_group", "index", "sample_n"),
    all.x = TRUE
  )
  out[is.na(n_partners), `:=`(
    window_bp = window_bp_value,
    partner_scope = scope_value,
    n_partners = 0L,
    sum_r2 = 0,
    ld_score = 1,
    n_ld_friends = 0L,
    friends_rate = NA_real_,
    mean_r2 = NA_real_
  )]
  out
}

rank_biserial <- function(w_stat, n_sig, n_null) {
  if (n_sig == 0 || n_null == 0) {
    return(NA_real_)
  }
  2 * as.numeric(w_stat) / (n_sig * n_null) - 1
}

run_one_wilcoxon <- function(
  dt,
  metric_value,
  target_class_value,
  partner_scope_value,
  source_id_value
) {
  sig_group <- paste0(target_class_value, "_sig")
  null_group <- paste0(target_class_value, "_null")

  x <- dt[
    plot_group == sig_group &
      partner_scope == partner_scope_value &
      source_id == source_id_value,
    get(metric_value)
  ]
  y <- dt[
    plot_group == null_group &
      partner_scope == partner_scope_value &
      source_id == source_id_value,
    get(metric_value)
  ]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) == 0 || length(y) == 0) {
    return(NULL)
  }

  wt <- suppressWarnings(
    wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
  )

  data.table(
    analysis_level = "source",
    source_id = source_id_value,
    target_class = target_class_value,
    partner_scope = partner_scope_value,
    metric = metric_value,
    n_sig = length(x),
    n_null = length(y),
    median_sig = as.numeric(median(x)),
    median_null = as.numeric(median(y)),
    median_difference = as.numeric(median(x) - median(y)),
    wilcoxon_w = as.numeric(wt$statistic),
    rank_biserial = rank_biserial(wt$statistic, length(x), length(y)),
    p_value = wt$p.value
  )
}

run_tests <- function(metrics) {
  test_metrics <- c("ld_score", "n_ld_friends", "friends_rate", "mean_r2")
  target_classes <- c("SNV_INDEL", "SV")
  out <- list()

  for (source_value in unique(metrics$source_id)) {
    for (target_value in target_classes) {
      for (scope_value in partner_scope_levels) {
        for (metric_value in test_metrics) {
          z <- run_one_wilcoxon(
            metrics,
            metric_value,
            target_value,
            scope_value,
            source_value
          )
          if (!is.null(z)) {
            out[[length(out) + 1]] <- z
          }
        }
      }
    }
  }

  pooled <- copy(metrics)
  pooled[, source_id := "POOLED"]
  for (target_value in target_classes) {
    for (scope_value in partner_scope_levels) {
      for (metric_value in test_metrics) {
        z <- run_one_wilcoxon(
          pooled,
          metric_value,
          target_value,
          scope_value,
          "POOLED"
        )
        if (!is.null(z)) {
          z[, analysis_level := "pooled"]
          out[[length(out) + 1]] <- z
        }
      }
    }
  }

  tests <- rbindlist(out, fill = TRUE)
  tests[, p_adj_bh := p.adjust(p_value, method = "BH"), by = analysis_level]
  setorder(tests, analysis_level, source_id, metric, partner_scope, target_class)
  tests
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 0.001,
      format(p, scientific = TRUE, digits = 2),
      sprintf("%.3f", p)
    )
  )
}

make_plot_data <- function(metrics, scopes) {
  out <- copy(metrics[partner_scope %in% scopes])
  out[, plot_group := factor(plot_group, levels = group_levels)]
  melt(
    out,
    id.vars = c(
      "source_id", "plot_group", "index", "window_bp", "partner_scope"
    ),
    measure.vars = c("ld_score", "n_ld_friends"),
    variable.name = "metric",
    value.name = "value"
  )
}

make_annotations <- function(plot_long, tests, scopes) {
  ann <- copy(tests[
    analysis_level == "pooled" &
      partner_scope %in% scopes &
      metric %in% c("ld_score", "n_ld_friends")
  ])
  if (nrow(ann) == 0) {
    return(NULL)
  }

  ann[, label := paste0("BH-FDR=", format_p(p_adj_bh))]
  ymax <- plot_long[, .(
    y = max(value, na.rm = TRUE)
  ), by = .(partner_scope, metric)]
  ann <- merge(ann, ymax, by = c("partner_scope", "metric"), all.x = TRUE)
  ann[, `:=`(
    x = fifelse(target_class == "SNV_INDEL", 1.5, 3.5),
    y = fifelse(y <= 0, 1, y * 1.15)
  )]
  ann
}

draw_metric_plot <- function(plot_long, annotations, split_partner, output_stem) {
  plot_long[, metric_label := factor(
    metric,
    levels = c("ld_score", "n_ld_friends"),
    labels = c("LD score", "Number of LD friends")
  )]

  p <- ggplot(plot_long, aes(x = plot_group, y = value, fill = plot_group)) +
    geom_boxplot(
      width = 0.68,
      outlier.shape = NA,
      linewidth = 0.25
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_x_discrete(
      labels = c(
        "SNV_INDEL_null" = "SNV/INDEL\nnull",
        "SNV_INDEL_sig" = "SNV/INDEL\nsig",
        "SV_null" = "SV\nnull",
        "SV_sig" = "SV\nsig"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1)
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 8),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.margin = margin(4, 4, 4, 4, unit = "pt")
    )

  if (split_partner) {
    p <- p + facet_grid(
      metric_label ~ partner_scope,
      scales = "free_y",
      labeller = labeller(
        partner_scope = c(
          "SNV_INDEL" = "SNV/INDEL partners",
          "SV" = "SV partners"
        )
      )
    )
  } else {
    p <- p + facet_wrap(~metric_label, scales = "free_y", nrow = 1)
  }

  if (!is.null(annotations) && nrow(annotations) > 0) {
    annotations[, metric_label := factor(
      metric,
      levels = c("ld_score", "n_ld_friends"),
      labels = c("LD score", "Number of LD friends")
    )]
    p <- p + geom_text(
      data = annotations,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.2,
      vjust = 0
    )
  }

  width_cm <- if (split_partner) 12 else 10
  height_cm <- if (split_partner) 10 else 5.5
  for (extension in c("pdf", "svg")) {
    ggsave(
      file.path(figure_dir, paste0(output_stem, ".", extension)),
      p,
      width = width_cm,
      height = height_cm,
      units = "cm"
    )
  }
  ggsave(
    file.path(figure_dir, paste0(output_stem, ".tiff")),
    p,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
}

# Main
check_file(metadata_file)
metadata <- fread(metadata_file)
required_metadata <- c(
  "source_id", "bfile", "plot_group", "id_file", "out_prefix"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop("Metadata is missing columns: ", paste(missing_metadata, collapse = ", "))
}

metadata <- metadata[plot_group %in% group_levels]
if (nrow(metadata) == 0) {
  stop("No rows for the four LD-decay groups in: ", metadata_file)
}

metric_list <- list()
for (i in seq_len(nrow(metadata))) {
  loaded <- read_target_pairs(metadata[i])
  for (window_bp_value in analysis_windows_bp) {
    for (scope_value in partner_scope_levels) {
      metric_list[[length(metric_list) + 1]] <- summarise_scope(
        loaded$targets,
        loaded$pairs,
        scope_value,
        window_bp_value
      )
    }
  }
  rm(loaded)
  invisible(gc())
}

target_metrics <- rbindlist(metric_list, fill = TRUE)
target_metrics[, plot_group := factor(plot_group, levels = group_levels)]
target_metrics[, partner_scope := factor(
  partner_scope,
  levels = partner_scope_levels
)]

make_group_summary <- function(metrics) {
  metrics[, .(
    n_target = as.integer(.N),
    median_ld_score = as.numeric(median(ld_score, na.rm = TRUE)),
    q25_ld_score = as.numeric(quantile(ld_score, 0.25, na.rm = TRUE)),
    q75_ld_score = as.numeric(quantile(ld_score, 0.75, na.rm = TRUE)),
    median_ld_friends = as.numeric(median(n_ld_friends, na.rm = TRUE)),
    q25_ld_friends = as.numeric(
      quantile(n_ld_friends, 0.25, na.rm = TRUE)
    ),
    q75_ld_friends = as.numeric(
      quantile(n_ld_friends, 0.75, na.rm = TRUE)
    ),
    median_friends_rate = as.numeric(median(friends_rate, na.rm = TRUE)),
    median_mean_r2 = as.numeric(median(mean_r2, na.rm = TRUE))
  ), by = .(source_id, plot_group, partner_scope)]
}

metrics_1mb <- target_metrics[window_bp == 1000000L]
metrics_250kb <- target_metrics[window_bp == 250000L]

group_summary_1mb <- make_group_summary(metrics_1mb)
group_summary_250kb <- make_group_summary(metrics_250kb)
wilcoxon_tests_1mb <- run_tests(metrics_1mb)
wilcoxon_tests_250kb <- run_tests(metrics_250kb)

config <- data.table(
  parameter = c(
    "target_groups",
    "null_definition",
    "analysis_windows_bp",
    "primary_window_bp",
    "sensitivity_window_bp",
    "ld_score_definition",
    "ld_friend_definition",
    "ld_friend_sample_size",
    "primary_test",
    "multiple_testing",
    "partner_scopes"
  ),
  value = c(
    paste(group_levels, collapse = ","),
    "Existing LD-decay nulls: chromosome + variant_class, null_ratio=10",
    paste(analysis_windows_bp, collapse = ","),
    "1000000",
    "250000",
    "1 + sum(r2), excluding self-pairs from the sum",
    paste0("Two-sided Pearson correlation P < ", ld_friend_p_threshold),
    "Total N from source .fam; pair-specific missingness unavailable",
    "Two-sided Wilcoxon rank-sum test: sig versus null within target class",
    "Benjamini-Hochberg FDR",
    paste(partner_scope_levels, collapse = ",")
  )
)

fwrite(
  metrics_1mb,
  file.path(analysis_dir, "ld_friends_scores.per_target.tsv.gz"),
  sep = "\t"
)
fwrite(
  group_summary_1mb,
  file.path(analysis_dir, "ld_friends_scores.group_summary.tsv"),
  sep = "\t"
)
fwrite(
  wilcoxon_tests_1mb,
  file.path(analysis_dir, "ld_friends_scores.wilcoxon.tsv"),
  sep = "\t"
)
fwrite(
  config,
  file.path(analysis_dir, "ld_friends_scores.config.tsv"),
  sep = "\t"
)

fwrite(
  metrics_250kb,
  file.path(analysis_dir, "ld_friends_scores.per_target.250kb.tsv.gz"),
  sep = "\t"
)
fwrite(
  group_summary_250kb,
  file.path(analysis_dir, "ld_friends_scores.group_summary.250kb.tsv"),
  sep = "\t"
)
fwrite(
  wilcoxon_tests_250kb,
  file.path(analysis_dir, "ld_friends_scores.wilcoxon.250kb.tsv"),
  sep = "\t"
)

all_plot_data <- make_plot_data(metrics_1mb, "ALL")
draw_metric_plot(
  all_plot_data,
  make_annotations(all_plot_data, wilcoxon_tests_1mb, "ALL"),
  split_partner = FALSE,
  output_stem = "ld_friends_scores_sig_vs_null.all_partners"
)

split_plot_data <- make_plot_data(metrics_1mb, c("SNV_INDEL", "SV"))
draw_metric_plot(
  split_plot_data,
  make_annotations(
    split_plot_data,
    wilcoxon_tests_1mb,
    c("SNV_INDEL", "SV")
  ),
  split_partner = TRUE,
  output_stem = "ld_friends_scores_sig_vs_null.by_partner_class"
)

all_plot_data_250kb <- make_plot_data(metrics_250kb, "ALL")
draw_metric_plot(
  all_plot_data_250kb,
  make_annotations(
    all_plot_data_250kb,
    wilcoxon_tests_250kb,
    "ALL"
  ),
  split_partner = FALSE,
  output_stem = "ld_friends_scores_sig_vs_null.all_partners.250kb"
)

split_plot_data_250kb <- make_plot_data(
  metrics_250kb,
  c("SNV_INDEL", "SV")
)
draw_metric_plot(
  split_plot_data_250kb,
  make_annotations(
    split_plot_data_250kb,
    wilcoxon_tests_250kb,
    c("SNV_INDEL", "SV")
  ),
  split_partner = TRUE,
  output_stem = "ld_friends_scores_sig_vs_null.by_partner_class.250kb"
)

message("\n[DONE] LD friends / LD score analysis (1 Mb + 250 kb)")
message("Tables: ", analysis_dir)
message(
  "Main figure: ",
  file.path(
    figure_dir,
    "ld_friends_scores_sig_vs_null.all_partners.pdf"
  )
)
message(
  "Partner-class figure: ",
  file.path(
    figure_dir,
    "ld_friends_scores_sig_vs_null.by_partner_class.pdf"
  )
)
message(
  "250 kb main figure: ",
  file.path(
    figure_dir,
    "ld_friends_scores_sig_vs_null.all_partners.250kb.pdf"
  )
)
message(
  "250 kb partner-class figure: ",
  file.path(
    figure_dir,
    "ld_friends_scores_sig_vs_null.by_partner_class.250kb.pdf"
  )
)
message("\nPooled primary tests (1 Mb):")
print(
  wilcoxon_tests_1mb[
    analysis_level == "pooled" &
      partner_scope == "ALL" &
      metric %in% c("ld_score", "n_ld_friends")
  ]
)
message("\nPooled sensitivity tests (250 kb):")
print(
  wilcoxon_tests_250kb[
    analysis_level == "pooled" &
      partner_scope == "ALL" &
      metric %in% c("ld_score", "n_ld_friends")
  ]
)
