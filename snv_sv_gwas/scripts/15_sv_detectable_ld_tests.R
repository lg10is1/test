#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

# Standalone post-processing analysis for the LD_decay_public PLINK LD files.
# It does not run PLINK and does not use LD-friends or LD-score results.

default_outdir <- paste0(
  "/path/to/EOSCZ_PROJECT/",
  "figure_analysis/SV_SNV_LD/LD_decay_public"
)

parse_args <- function(args) {
  out <- list(
    outdir = default_outdir,
    threshold = 0.1,
    max_dist_bp = 1000000L,
    seed = 1L
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--outdir") {
      i <- i + 1L
      if (i > length(args)) stop("--outdir requires a value")
      out$outdir <- args[[i]]
    } else if (key == "--threshold") {
      i <- i + 1L
      if (i > length(args)) stop("--threshold requires a value")
      out$threshold <- as.numeric(args[[i]])
    } else if (key == "--max-dist-bp") {
      i <- i + 1L
      if (i > length(args)) stop("--max-dist-bp requires a value")
      out$max_dist_bp <- as.integer(args[[i]])
    } else if (key == "--seed") {
      i <- i + 1L
      if (i > length(args)) stop("--seed requires a value")
      out$seed <- suppressWarnings(as.integer(args[[i]]))
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 1L
  }

  if (!is.finite(out$threshold) || out$threshold < 0 || out$threshold > 1) {
    stop("--threshold must be between 0 and 1")
  }
  if (is.na(out$max_dist_bp) || out$max_dist_bp < 1L) {
    stop("--max-dist-bp must be a positive integer")
  }
  if (is.na(out$seed)) {
    stop("--seed must be an integer")
  }
  out
}

check_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) stop("Empty file: ", path)
  invisible(TRUE)
}

find_r2_column <- function(dt) {
  hit <- intersect(c("R2", "R^2", "UNPHASED_R2", "PHASED_R2"), names(dt))
  if (length(hit) == 0L) {
    stop("No r2 column found. Columns: ", paste(names(dt), collapse = ", "))
  }
  hit[[1L]]
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
  bim[, id_has_sv_tag := grepl(
    "(^|[-_])(INS|DEL|DUP|INV)([-_]|$)",
    SNP,
    ignore.case = TRUE
  )]
  bim[, variant_class := fifelse(
    max_allele_len >= 50L | id_has_sv_tag,
    "SV",
    "SNV_INDEL"
  )]
  unique(bim[, .(SNP, CHR, BP, variant_class)])
}

read_index_table <- function(meta_row) {
  id_file <- as.character(meta_row$id_file)
  check_file(id_file)
  ids <- unique(as.character(fread(id_file, header = FALSE)[[1L]]))
  ids <- ids[!is.na(ids) & ids != ""]
  data.table(
    source_id = as.character(meta_row$source_id),
    plot_group = as.character(meta_row$plot_group),
    index = ids,
    index_key = paste(as.character(meta_row$source_id), ids, sep = "::")
  )
}

read_sv_to_snv_indel_max_r2 <- function(meta_row, max_dist_bp) {
  source_id <- as.character(meta_row$source_id)
  plot_group <- as.character(meta_row$plot_group)
  ld_file <- paste0(as.character(meta_row$out_prefix), ".ld.gz")
  check_file(ld_file)

  index_dt <- read_index_table(meta_row)
  if (nrow(index_dt) == 0L) stop("No index IDs for ", source_id, " / ", plot_group)

  bim <- read_bim_classes(as.character(meta_row$bfile))
  missing_index <- setdiff(index_dt$index, bim$SNP)
  if (length(missing_index) > 0L) {
    stop(
      "Index IDs missing from BIM for ", source_id, " / ", plot_group,
      ": ", paste(head(missing_index, 10L), collapse = ", ")
    )
  }
  non_sv_index <- bim[SNP %in% index_dt$index & variant_class != "SV", SNP]
  if (length(non_sv_index) > 0L) {
    stop(
      "SV analysis contains index IDs not classified as SV for ", source_id,
      " / ", plot_group, ": ",
      paste(head(non_sv_index, 10L), collapse = ", ")
    )
  }

  index_position <- bim[SNP %in% index_dt$index, .(
    index = SNP,
    index_chr = as.character(CHR),
    index_bp = as.numeric(BP)
  )]
  index_dt <- merge(index_dt, index_position, by = "index", all.x = TRUE, sort = FALSE)
  if (anyNA(index_dt$index_chr) || anyNA(index_dt$index_bp)) {
    stop("Missing index chromosome or position for ", source_id, " / ", plot_group)
  }

  target_ids <- bim[variant_class == "SNV_INDEL", SNP]
  if (length(target_ids) == 0L) {
    stop("No SNV/INDEL targets in BIM for ", source_id)
  }

  message(
    "[READ] ", source_id, " / ", plot_group,
    " | index SVs=", nrow(index_dt),
    " | SNV/INDEL targets=", length(target_ids),
    " | ", ld_file
  )

  ld <- fread(ld_file)
  required <- c("SNP_A", "SNP_B", "BP_A", "BP_B")
  missing_cols <- setdiff(required, names(ld))
  if (length(missing_cols) > 0L) {
    stop("LD file missing columns: ", paste(missing_cols, collapse = ", "))
  }
  r2_col <- find_r2_column(ld)

  ld[, `:=`(
    SNP_A = as.character(SNP_A),
    SNP_B = as.character(SNP_B),
    BP_A = as.numeric(BP_A),
    BP_B = as.numeric(BP_B),
    ld_r2 = as.numeric(get(r2_col))
  )]

  index_ids <- index_dt$index
  from_a <- ld[SNP_A %in% index_ids, .(
    index = SNP_A,
    target = SNP_B,
    bp_index = BP_A,
    bp_target = BP_B,
    ld_r2
  )]
  from_b <- ld[SNP_B %in% index_ids, .(
    index = SNP_B,
    target = SNP_A,
    bp_index = BP_B,
    bp_target = BP_A,
    ld_r2
  )]
  pairs <- rbindlist(list(from_a, from_b), use.names = TRUE)

  pairs[, distance_bp := abs(bp_target - bp_index)]
  pairs <- pairs[
    target %in% target_ids &
      is.finite(distance_bp) &
      distance_bp > 0 &
      distance_bp <= max_dist_bp &
      is.finite(ld_r2) &
      ld_r2 >= 0 &
      ld_r2 <= 1
  ]

  if (nrow(pairs) > 0L) {
    setorder(pairs, index, target, -ld_r2)
    pairs <- unique(pairs, by = c("index", "target"))
    observed <- pairs[, .(
      max_r2 = max(ld_r2),
      n_proximal_snv_indel_pairs = .N,
      nearest_partner_distance_bp = min(distance_bp),
      max_r2_partner = target[which.max(ld_r2)][1L],
      max_r2_partner_distance_bp = distance_bp[which.max(ld_r2)][1L]
    ), by = index]
  } else {
    observed <- data.table(
      index = character(),
      max_r2 = numeric(),
      n_proximal_snv_indel_pairs = integer(),
      nearest_partner_distance_bp = numeric(),
      max_r2_partner = character(),
      max_r2_partner_distance_bp = numeric()
    )
  }

  out <- merge(index_dt, observed, by = "index", all.x = TRUE, sort = FALSE)
  out[is.na(n_proximal_snv_indel_pairs), `:=`(
    max_r2 = 0,
    n_proximal_snv_indel_pairs = 0L
  )]
  out[, `:=`(
    max_dist_bp = as.integer(max_dist_bp),
    partner_class = "SNV_INDEL",
    no_eligible_partner = n_proximal_snv_indel_pairs == 0L
  )]
  out
}

rank_biserial_from_w <- function(w, n_sig, n_null) {
  if (n_sig == 0L || n_null == 0L) return(NA_real_)
  2 * as.numeric(w) / (n_sig * n_null) - 1
}

two_sample_ecdf_metrics <- function(x, y) {
  x <- as.numeric(x[is.finite(x)])
  y <- as.numeric(y[is.finite(y)])
  if (length(x) == 0L || length(y) == 0L) {
    return(c(ks = NA_real_, cvm = NA_real_, wasserstein = NA_real_))
  }
  grid <- sort(unique(c(x, y)))
  fx <- findInterval(grid, sort(x)) / length(x)
  fy <- findInterval(grid, sort(y)) / length(y)
  ks <- max(abs(fx - fy))
  pooled_counts <- tabulate(match(c(x, y), grid), nbins = length(grid))
  pooled_weights <- pooled_counts / (length(x) + length(y))
  cvm <- sum((fx - fy)^2 * pooled_weights)
  wasserstein <- if (length(grid) < 2L) {
    0
  } else {
    sum(abs(head(fx - fy, -1L)) * diff(grid))
  }
  c(ks = ks, cvm = cvm, wasserstein = wasserstein)
}

calculate_distribution_metrics <- function(values, labels) {
  x <- values[labels == "SV_sig"]
  y <- values[labels == "SV_null"]
  ecdf_metrics <- two_sample_ecdf_metrics(x, y)
  c(
    mean_diff_sig_minus_null = mean(x) - mean(y),
    median_diff_sig_minus_null = median(x) - median(y),
    ks_distance = unname(ecdf_metrics[["ks"]]),
    cvm_ecdf_distance = unname(ecdf_metrics[["cvm"]]),
    wasserstein_distance = unname(ecdf_metrics[["wasserstein"]])
  )
}

permute_within_strata <- function(labels, strata) {
  out <- labels
  for (idx in split(seq_along(labels), strata)) {
    out[idx] <- sample(labels[idx], length(idx), replace = FALSE)
  }
  out
}

run_stratified_permutation_tests <- function(
  dt,
  thresholds,
  n_perm = 10000L,
  seed = 1L
) {
  z <- copy(as.data.table(dt))
  z[, label := as.character(plot_group)]
  z[, stratum := paste(source_id, index_chr, sep = "::")]
  stratum_counts <- z[, .(
    n_sig = sum(label == "SV_sig"),
    n_null = sum(label == "SV_null")
  ), by = stratum]
  informative <- stratum_counts[n_sig > 0L & n_null > 0L, stratum]
  excluded <- z[!stratum %in% informative, .N]
  z <- z[stratum %in% informative]
  if (z[label == "SV_sig", .N] < 2L || z[label == "SV_null", .N] < 2L) {
    stop("Too few SVs in informative source-by-chromosome permutation strata")
  }

  values <- z$max_r2
  labels <- z$label
  strata <- z$stratum
  observed_metrics <- calculate_distribution_metrics(values, labels)
  observed_threshold_diff <- vapply(
    thresholds,
    function(threshold) {
      mean(values[labels == "SV_sig"] > threshold) -
        mean(values[labels == "SV_null"] > threshold)
    },
    numeric(1L)
  )

  set.seed(seed)
  metric_names <- names(observed_metrics)
  perm_metrics <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(metric_names),
    dimnames = list(NULL, metric_names)
  )
  perm_threshold_diff <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(thresholds)
  )

  for (b in seq_len(n_perm)) {
    permuted <- permute_within_strata(labels, strata)
    perm_metrics[b, ] <- calculate_distribution_metrics(values, permuted)
    perm_threshold_diff[b, ] <- vapply(
      thresholds,
      function(threshold) {
        mean(values[permuted == "SV_sig"] > threshold) -
          mean(values[permuted == "SV_null"] > threshold)
      },
      numeric(1L)
    )
  }

  signed_metrics <- c("mean_diff_sig_minus_null", "median_diff_sig_minus_null")
  distribution_tests <- rbindlist(lapply(metric_names, function(metric) {
    observed <- observed_metrics[[metric]]
    null_values <- perm_metrics[, metric]
    p_value <- if (metric %in% signed_metrics) {
      (1 + sum(abs(null_values) >= abs(observed), na.rm = TRUE)) /
        (1 + sum(is.finite(null_values)))
    } else {
      (1 + sum(null_values >= observed, na.rm = TRUE)) /
        (1 + sum(is.finite(null_values)))
    }
    data.table(
      test = paste0("stratified_permutation_", metric),
      statistic = observed,
      p_value = p_value,
      n_perm = n_perm,
      n_sig = sum(labels == "SV_sig"),
      n_null = sum(labels == "SV_null"),
      n_excluded_noninformative_strata = excluded,
      strata = "source_id_by_chromosome",
      alternative = "two.sided"
    )
  }))
  distribution_tests[, p_adj_bh := p.adjust(p_value, method = "BH")]

  threshold_permutation <- rbindlist(lapply(seq_along(thresholds), function(j) {
    observed <- observed_threshold_diff[[j]]
    null_values <- perm_threshold_diff[, j]
    data.table(
      threshold = thresholds[[j]],
      stratified_proportion_diff_sig_minus_null = observed,
      permutation_p_two_sided =
        (1 + sum(abs(null_values) >= abs(observed), na.rm = TRUE)) /
        (1 + sum(is.finite(null_values))),
      n_perm = n_perm,
      n_sig = sum(labels == "SV_sig"),
      n_null = sum(labels == "SV_null"),
      n_excluded_noninformative_strata = excluded,
      strata = "source_id_by_chromosome"
    )
  }))
  threshold_permutation[, permutation_p_adj_bh := p.adjust(
    permutation_p_two_sided,
    method = "BH"
  )]

  list(
    distribution_tests = distribution_tests,
    threshold_permutation = threshold_permutation,
    stratum_counts = stratum_counts,
    n_informative = nrow(z),
    n_excluded = excluded
  )
}

run_threshold_fisher_grid <- function(dt, thresholds) {
  out <- rbindlist(lapply(thresholds, function(threshold) {
    detectable <- dt$max_r2 > threshold
    tab <- table(
      factor(dt$plot_group, levels = c("SV_null", "SV_sig")),
      factor(detectable, levels = c(FALSE, TRUE))
    )
    ft <- fisher.test(tab, alternative = "two.sided")
    data.table(
      threshold = threshold,
      sig_detectable = unname(tab["SV_sig", "TRUE"]),
      sig_not_detectable = unname(tab["SV_sig", "FALSE"]),
      null_detectable = unname(tab["SV_null", "TRUE"]),
      null_not_detectable = unname(tab["SV_null", "FALSE"]),
      proportion_sig = mean(detectable[dt$plot_group == "SV_sig"]),
      proportion_null = mean(detectable[dt$plot_group == "SV_null"]),
      pooled_proportion_diff_sig_minus_null =
        mean(detectable[dt$plot_group == "SV_sig"]) -
        mean(detectable[dt$plot_group == "SV_null"]),
      odds_ratio_detectable_sig_vs_null = unname(ft$estimate),
      conf_low = unname(ft$conf.int[1L]),
      conf_high = unname(ft$conf.int[2L]),
      fisher_p_two_sided = ft$p.value
    )
  }))
  out[, fisher_p_adj_bh := p.adjust(fisher_p_two_sided, method = "BH")]
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
outdir <- normalizePath(args$outdir, winslash = "/", mustWork = FALSE)
rdata_dir <- file.path(outdir, "rdata")
fig_dir <- file.path(outdir, "figures")
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(outdir, "ld_jobs.metadata.tsv")
check_file(metadata_file)
metadata <- fread(metadata_file)
required_meta <- c("source_id", "bfile", "plot_group", "id_file", "out_prefix")
missing_meta <- setdiff(required_meta, names(metadata))
if (length(missing_meta) > 0L) {
  stop("Metadata missing columns: ", paste(missing_meta, collapse = ", "))
}

metadata_sv <- metadata[plot_group %in% c("SV_sig", "SV_null")]
if (nrow(metadata_sv) == 0L) {
  stop("No SV_sig or SV_null jobs in metadata: ", metadata_file)
}
if (!all(c("SV_sig", "SV_null") %in% metadata_sv$plot_group)) {
  stop("Both SV_sig and SV_null groups are required")
}

result_list <- lapply(seq_len(nrow(metadata_sv)), function(i) {
  read_sv_to_snv_indel_max_r2(metadata_sv[i], args$max_dist_bp)
})
sv_test <- rbindlist(result_list, fill = TRUE)
sv_test[, plot_group := factor(plot_group, levels = c("SV_null", "SV_sig"))]
sv_test[, detectable_ld := max_r2 > args$threshold]
setorder(sv_test, plot_group, source_id, index)

group_summary <- sv_test[, .(
  n_sv = .N,
  n_detectable = sum(detectable_ld),
  proportion_detectable = mean(detectable_ld),
  n_without_eligible_partner = sum(no_eligible_partner),
  median_max_r2 = median(max_r2),
  q25_max_r2 = as.numeric(quantile(max_r2, 0.25)),
  q75_max_r2 = as.numeric(quantile(max_r2, 0.75)),
  mean_max_r2 = mean(max_r2)
), by = plot_group]

contingency <- sv_test[, .N, by = .(plot_group, detectable_ld)]
contingency <- dcast(
  contingency,
  plot_group ~ detectable_ld,
  value.var = "N",
  fill = 0
)
if (!"FALSE" %in% names(contingency)) contingency[, `FALSE` := 0L]
if (!"TRUE" %in% names(contingency)) contingency[, `TRUE` := 0L]
setcolorder(contingency, c("plot_group", "FALSE", "TRUE"))

tab <- table(
  factor(sv_test$plot_group, levels = c("SV_null", "SV_sig")),
  factor(sv_test$detectable_ld, levels = c(FALSE, TRUE))
)
fisher_res <- fisher.test(tab, alternative = "two.sided")
wilcox_res <- wilcox.test(
  max_r2 ~ plot_group,
  data = sv_test,
  alternative = "two.sided",
  exact = FALSE,
  conf.int = FALSE
)

n_sig <- sv_test[plot_group == "SV_sig", .N]
n_null <- sv_test[plot_group == "SV_null", .N]
med_sig <- sv_test[plot_group == "SV_sig", median(max_r2)]
med_null <- sv_test[plot_group == "SV_null", median(max_r2)]
sig_not_detectable <- unname(tab["SV_sig", "FALSE"])
sig_detectable <- unname(tab["SV_sig", "TRUE"])
null_not_detectable <- unname(tab["SV_null", "FALSE"])
null_detectable <- unname(tab["SV_null", "TRUE"])

formal_tests <- rbindlist(list(
  data.table(
    test = "Fisher_exact_detectable_LD",
    comparison = "SV_sig_vs_SV_null",
    outcome = paste0("max_r2_gt_", format(args$threshold, trim = TRUE)),
    alternative = "two.sided",
    statistic_name = "odds_ratio_detectable_SV_sig_vs_SV_null",
    statistic = unname(fisher_res$estimate),
    conf_low = unname(fisher_res$conf.int[1L]),
    conf_high = unname(fisher_res$conf.int[2L]),
    p_value = fisher_res$p.value,
    n_sig = n_sig,
    n_null = n_null,
    sig_detectable = sig_detectable,
    sig_not_detectable = sig_not_detectable,
    null_detectable = null_detectable,
    null_not_detectable = null_not_detectable,
    threshold = args$threshold,
    max_dist_bp = args$max_dist_bp
  ),
  data.table(
    test = "Wilcoxon_rank_sum_max_r2",
    comparison = "SV_sig_vs_SV_null",
    outcome = "per_SV_max_r2_with_proximal_SNV_INDEL",
    alternative = "two.sided",
    statistic_name = "W",
    statistic = unname(wilcox_res$statistic),
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = wilcox_res$p.value,
    n_sig = n_sig,
    n_null = n_null,
    sig_detectable = sig_detectable,
    sig_not_detectable = sig_not_detectable,
    null_detectable = null_detectable,
    null_not_detectable = null_not_detectable,
    threshold = NA_real_,
    max_dist_bp = args$max_dist_bp
  )
), fill = TRUE)
formal_tests[test == "Wilcoxon_rank_sum_max_r2", `:=`(
  median_sig = med_sig,
  median_null = med_null,
  rank_biserial_null_vs_sig = rank_biserial_from_w(statistic, n_sig, n_null)
)]

threshold_grid <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
threshold_sensitivity <- run_threshold_fisher_grid(sv_test, threshold_grid)

ks_classical <- suppressWarnings(ks.test(
  sv_test[plot_group == "SV_sig", max_r2],
  sv_test[plot_group == "SV_null", max_r2],
  alternative = "two.sided",
  exact = FALSE
))
classical_distribution_tests <- data.table(
  test = c("Wilcoxon_rank_sum", "Kolmogorov_Smirnov_asymptotic"),
  statistic = c(unname(wilcox_res$statistic), unname(ks_classical$statistic)),
  p_value = c(wilcox_res$p.value, ks_classical$p.value),
  note = c(
    "per-SV max r2; exact=FALSE",
    "per-SV max r2; asymptotic P because ties are present"
  )
)
classical_distribution_tests[, p_adj_bh := p.adjust(p_value, method = "BH")]

permutation_results <- run_stratified_permutation_tests(
  sv_test,
  thresholds = threshold_grid,
  n_perm = 10000L,
  seed = args$seed
)
distribution_tests <- rbindlist(list(
  classical_distribution_tests[, .(
    test,
    statistic,
    p_value,
    p_adj_bh,
    n_perm = NA_integer_,
    n_sig = n_sig,
    n_null = n_null,
    n_excluded_noninformative_strata = NA_integer_,
    strata = "none",
    alternative = "two.sided",
    note
  )],
  permutation_results$distribution_tests[, .(
    test,
    statistic,
    p_value,
    p_adj_bh,
    n_perm,
    n_sig,
    n_null,
    n_excluded_noninformative_strata,
    strata,
    alternative,
    note = "labels permuted within informative source-by-chromosome strata"
  )]
), fill = TRUE)

threshold_sensitivity <- merge(
  threshold_sensitivity,
  permutation_results$threshold_permutation,
  by = "threshold",
  all.x = TRUE,
  sort = TRUE
)

prefix <- file.path(rdata_dir, "sv_sig_vs_null.proximal_snv_indel")
fwrite(sv_test, paste0(prefix, ".maxR2.tsv"), sep = "\t")
fwrite(group_summary, paste0(prefix, ".summary.tsv"), sep = "\t")
fwrite(contingency, paste0(prefix, ".contingency.tsv"), sep = "\t")
fwrite(formal_tests, paste0(prefix, ".formal_tests.tsv"), sep = "\t")
fwrite(threshold_sensitivity, paste0(prefix, ".threshold_sensitivity.tsv"), sep = "\t")
fwrite(distribution_tests, paste0(prefix, ".distribution_tests.tsv"), sep = "\t")
fwrite(
  permutation_results$stratum_counts,
  paste0(prefix, ".permutation_strata.tsv"),
  sep = "\t"
)

plot_dt <- copy(sv_test)
plot_dt[, type := fifelse(plot_group == "SV_sig", "sv_sig", "sv_null")]
plot_dt[, type := factor(type, levels = c("sv_null", "sv_sig"))]

gg_maxr2_hist <- ggplot(plot_dt, aes(x = max_r2, fill = type)) +
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
    drop = FALSE,
    labels = c(
      "sv_null" = "SV null",
      "sv_sig" = "SV sig"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = expression(Maximum~LD~italic(r)^2), y = "Density", fill = NULL) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.82),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  gg_maxr2_hist,
  file = file.path(fig_dir, "sv_maxR2_hist_v2.svg"),
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  gg_maxr2_hist,
  file = file.path(fig_dir, "sv_maxR2_hist_v2.pdf"),
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  gg_maxr2_hist,
  file = file.path(fig_dir, "sv_maxR2_hist_v2.tiff"),
  width = 5,
  height = 5,
  units = "cm",
  dpi = 300,
  compression = "lzw"
)

run_config <- data.table(
  parameter = c(
    "analysis_version",
    "input_metadata",
    "partner_class",
    "distance_definition",
    "max_dist_bp",
    "detectable_ld_definition",
    "missing_eligible_pair_handling",
    "fisher_alternative",
    "wilcoxon_alternative"
  ),
  value = c(
    "ld_detectable_version1",
    metadata_file,
    "SNV_INDEL only",
    "absolute difference between PLINK BIM BP coordinates; SV BP is the index reference",
    as.character(args$max_dist_bp),
    paste0("per-SV max r2 > ", args$threshold),
    "max_r2=0 and n_pairs=0",
    "two.sided",
    "two.sided"
  )
)
run_config <- rbind(
  run_config,
  data.table(
    parameter = c(
      "threshold_sensitivity_grid",
      "threshold_multiple_testing",
      "distribution_sensitivity_tests",
      "permutation_replicates",
      "permutation_seed",
      "permutation_strata",
      "permutation_noninformative_strata"
    ),
    value = c(
      paste(threshold_grid, collapse = ","),
      "BH across the fixed threshold grid",
      "Wilcoxon, asymptotic KS, stratified-permutation mean, median, KS, CvM-type ECDF, and Wasserstein tests",
      "10000",
      as.character(args$seed),
      "source_id by chromosome; lead/null counts preserved within stratum",
      "excluded from stratified permutation tests and reported"
    )
  )
)
fwrite(run_config, paste0(prefix, ".run_config.tsv"), sep = "\t")

message("[DONE] SV detectable-LD formal tests")
message("Per-SV table: ", paste0(prefix, ".maxR2.tsv"))
message("Group summary: ", paste0(prefix, ".summary.tsv"))
message("Contingency table: ", paste0(prefix, ".contingency.tsv"))
message("Formal tests: ", paste0(prefix, ".formal_tests.tsv"))
message("Threshold sensitivity: ", paste0(prefix, ".threshold_sensitivity.tsv"))
message("Distribution tests: ", paste0(prefix, ".distribution_tests.tsv"))
message("Formal maxR2 histogram: ", file.path(fig_dir, "sv_maxR2_hist_v2.pdf"))
print(group_summary)
print(formal_tests)
