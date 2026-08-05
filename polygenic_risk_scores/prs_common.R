############################################################
## Shared helpers for previous PRS downstream analyses.
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(gghalves)
  library(cowplot)
})

split_csv <- function(x) {
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

env_required <- function(name) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) stop("Required environment variable is not set: ", name)
  value
}

safe_name <- function(x) make.names(x, unique = TRUE)

theme_prs <- function(base_size = 10) {
  theme_cowplot(base_size) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = base_size),
      legend.title = element_blank()
    )
}

save_plot_pair <- function(plot, prefix, width_cm, height_cm) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    plot,
    file = paste0(prefix, ".pdf"),
    width = width_cm,
    height = height_cm,
    units = "cm",
    bg = "white"
  )
  ggsave(
    plot,
    file = paste0(prefix, ".tiff"),
    width = width_cm,
    height = height_cm,
    units = "cm",
    bg = "white",
    dpi = 300,
    compression = "lzw"
  )
}

read_prs <- function(file, trait) {
  x <- fread(file)
  if ("#FID" %in% names(x)) setnames(x, "#FID", "FID")
  if (!all(c("IID", "SCORE1_AVG") %in% names(x))) {
    stop("PRS file must contain IID and SCORE1_AVG: ", file)
  }
  x <- x[, .(sample = as.character(IID), score = as.numeric(SCORE1_AVG))]
  setnames(x, "score", paste0(trait, "_raw"))
  x
}

read_pc <- function(file, n_pcs, pc_cols = "", has_header = TRUE) {
  if (isTRUE(has_header)) {
    x <- fread(file, fill = TRUE)
    if ("#FID" %in% names(x)) setnames(x, "#FID", "FID")
    if (!all(c("FID", "IID") %in% names(x))) {
      if ("IID" %in% names(x)) {
        x[, FID := IID]
      } else {
        stop("Header PC file must contain IID: ", file)
      }
    }
  } else {
    x <- fread(file, header = FALSE, fill = TRUE)
    if (ncol(x) < n_pcs + 2L) stop("PC file has too few columns: ", file)
    setnames(x, 1:2, c("FID", "IID"))
    setnames(x, 3:ncol(x), paste0("PC", seq_len(ncol(x) - 2L)))
  }

  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  requested <- split_csv(pc_cols)
  if (!length(requested)) {
    requested <- paste0("PC", seq_len(n_pcs))
    if (!all(requested %in% names(x))) {
      requested_avg <- paste0("PC", seq_len(n_pcs), "_AVG")
      if (all(requested_avg %in% names(x))) requested <- requested_avg
    }
  }
  missing_cols <- setdiff(requested, names(x))
  if (length(missing_cols)) {
    stop("PC file missing columns: ", paste(missing_cols, collapse = ", "), " in ", file)
  }

  out <- x[, c("FID", "IID", requested), with = FALSE]
  setnames(out, requested, paste0("PC", seq_len(length(requested))))
  for (nm in paste0("PC", seq_len(length(requested)))) {
    out[, (nm) := as.numeric(get(nm))]
  }
  out[, sample := as.character(IID)]
  out[, c("sample", paste0("PC", seq_len(length(requested)))), with = FALSE]
}

read_batch <- function(file) {
  x <- fread(file, header = "auto", fill = TRUE)
  header_ok <- ncol(x) >= 3L && all(toupper(names(x)[1:2]) %in% c("FID", "IID", "#FID", "#IID"))
  if (!header_ok) x <- fread(file, header = FALSE, fill = TRUE)
  if (ncol(x) < 3L) stop("Batch file must have FID, IID and at least one batch column: ", file)
  setnames(x, 1:2, c("FID", "IID"))
  batch_cols <- names(x)[-(1:2)]
  if (any(grepl("^V[0-9]+$", batch_cols))) {
    setnames(x, batch_cols, paste0("batch", seq_along(batch_cols)))
    batch_cols <- names(x)[-(1:2)]
  }
  batch_cols <- safe_name(paste0("batch_", batch_cols))
  setnames(x, names(x)[-(1:2)], batch_cols)
  x[, sample := as.character(IID)]
  for (nm in batch_cols) x[, (nm) := as.factor(get(nm))]
  x[, c("sample", batch_cols), with = FALSE]
}

drop_bad_covariates <- function(dat, covar_cols) {
  keep <- character()
  for (cc in covar_cols) {
    z <- dat[[cc]]
    non_missing <- z[!is.na(z)]
    if (!length(non_missing)) next
    if (is.numeric(z)) {
      if (isTRUE(sd(non_missing) > 0)) keep <- c(keep, cc)
    } else {
      if (length(unique(non_missing)) > 1L) keep <- c(keep, cc)
    }
  }
  keep
}

adjust_prs_by_covariates <- function(dat, traits, covar_cols, mark) {
  covar_cols <- drop_bad_covariates(dat, covar_cols)
  if (!length(covar_cols)) stop("No usable PC/covariate columns after filtering.")
  for (trait in traits) {
    raw_col <- paste0(trait, "_raw")
    resid_col <- paste0(trait, "_resid_", mark)
    fml <- reformulate(covar_cols, response = raw_col)
    ok <- complete.cases(dat[, c(raw_col, covar_cols), drop = FALSE])
    fit <- lm(fml, data = dat[ok, , drop = FALSE])
    dat[[resid_col]] <- NA_real_
    dat[[resid_col]][which(ok)] <- resid(fit)
    dat[[trait]] <- as.numeric(scale(dat[[resid_col]]))
  }
  attr(dat, "adjustment_covariates") <- covar_cols
  dat
}

calc_auc <- function(score, label) {
  ok <- !is.na(score) & !is.na(label)
  score <- score[ok]
  label <- label[ok]
  n1 <- sum(label == 1)
  n0 <- sum(label == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_pseudo_r2_manual <- function(fit, data, y_col) {
  y <- data[[y_col]]
  pred <- predict(fit, type = "response")
  fit0 <- glm(reformulate("1", response = y_col), data = data, family = binomial)
  ll_full <- as.numeric(logLik(fit))
  ll_null <- as.numeric(logLik(fit0))
  n <- nobs(fit)
  coxsnell <- 1 - exp((2 / n) * (ll_null - ll_full))
  data.frame(
    McFadden = 1 - ll_full / ll_null,
    CoxSnell = coxsnell,
    Nagelkerke = coxsnell / (1 - exp((2 / n) * ll_null)),
    Tjur = mean(pred[y == 1], na.rm = TRUE) - mean(pred[y == 0], na.rm = TRUE)
  )
}

calc_prs_result <- function(dat, trait) {
  dd <- dat %>% filter(!is.na(.data[[trait]]), !is.na(y))
  fit <- glm(reformulate(trait, response = "y"), data = dd, family = binomial)
  coef_tab <- summary(fit)$coefficients
  beta <- coef_tab[trait, "Estimate"]
  se <- coef_tab[trait, "Std. Error"]
  r2 <- calc_pseudo_r2_manual(fit, dd, "y")
  data.frame(
    Trait = trait,
    n_HC = sum(dd$y == 0),
    n_Case = sum(dd$y == 1),
    Beta_per_1SD_adjusted_residual = beta,
    SE = se,
    OR_per_1SD_adjusted_residual = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    P = coef_tab[trait, "Pr(>|z|)"],
    AUC = calc_auc(dd[[trait]], dd$y),
    r2
  )
}

make_quantile_result <- function(dat, trait, n_quantile = 5, ref_quantile = "Q1") {
  q_levels <- paste0("Q", seq_len(n_quantile))
  dat_q <- dat %>%
    select(sample, type, y, all_of(trait)) %>%
    filter(!is.na(.data[[trait]]), !is.na(y)) %>%
    mutate(
      PRS_value = .data[[trait]],
      Quantile_num = dplyr::ntile(PRS_value, n_quantile),
      Quantile = factor(paste0("Q", Quantile_num), levels = q_levels)
    )
  q_summary <- dat_q %>%
    group_by(Quantile, Quantile_num) %>%
    summarise(
      n = n(),
      n_case = sum(y == 1),
      n_HC = sum(y == 0),
      case_fraction = mean(y == 1),
      PRS_mean = mean(PRS_value, na.rm = TRUE),
      PRS_median = median(PRS_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Trait = trait)
  dat_q$Quantile <- relevel(dat_q$Quantile, ref = ref_quantile)
  fit <- glm(y ~ Quantile, data = dat_q, family = binomial)
  coef_tab <- summary(fit)$coefficients
  or_tab <- data.frame(
    Trait = trait,
    Quantile = factor(q_levels, levels = q_levels),
    Quantile_num = seq_len(n_quantile),
    OR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P = NA_real_
  )
  or_tab[or_tab$Quantile == ref_quantile, c("OR", "CI_low", "CI_high")] <- 1
  for (q in setdiff(q_levels, ref_quantile)) {
    rn <- paste0("Quantile", q)
    if (rn %in% rownames(coef_tab)) {
      beta <- coef_tab[rn, "Estimate"]
      se <- coef_tab[rn, "Std. Error"]
      or_tab$OR[or_tab$Quantile == q] <- exp(beta)
      or_tab$CI_low[or_tab$Quantile == q] <- exp(beta - 1.96 * se)
      or_tab$CI_high[or_tab$Quantile == q] <- exp(beta + 1.96 * se)
      or_tab$P[or_tab$Quantile == q] <- coef_tab[rn, "Pr(>|z|)"]
    }
  }
  left_join(or_tab, q_summary, by = c("Trait", "Quantile", "Quantile_num"))
}

calc_cohens_d <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx < 2L || ny < 2L) return(NA_real_)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(y) - mean(x)) / sp
}

calc_cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  nx <- length(x)
  ny <- length(y)
  if (!nx || !ny) return(NA_real_)
  r <- rank(c(x, y), ties.method = "average")
  u_y <- sum(r[(nx + 1):(nx + ny)]) - ny * (ny + 1) / 2
  (2 * u_y / (nx * ny)) - 1
}

pairwise_score_test <- function(dat, group_col, score_col, g1, g2) {
  x <- dat %>% filter(.data[[group_col]] == g1) %>% pull(all_of(score_col))
  y <- dat %>% filter(.data[[group_col]] == g2) %>% pull(all_of(score_col))
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  tt <- t.test(x, y)
  wt <- wilcox.test(x, y, exact = FALSE, alternative = "two.sided")
  data.frame(
    group1 = g1,
    group2 = g2,
    n1 = length(x),
    n2 = length(y),
    mean1 = mean(x),
    mean2 = mean(y),
    median1 = median(x),
    median2 = median(y),
    mean_diff_group2_minus_group1 = mean(y) - mean(x),
    welch_t = as.numeric(tt$statistic),
    welch_p = tt$p.value,
    wilcox_w = as.numeric(wt$statistic),
    wilcox_p = wt$p.value,
    cohen_d_group2_minus_group1 = calc_cohens_d(x, y),
    cliffs_delta_group2_vs_group1 = calc_cliffs_delta(x, y)
  )
}

trait_cols <- c(
  SCZ = "#47a1a2",
  BIP = "#708090",
  ADHD = "#c65f7f",
  MDD = "#1f78b4",
  ASD = "#da7271"
)

case_control_cols <- c(HC = "#da7271", Case = "#1f78b4")
