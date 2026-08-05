options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(jsonlite)
})

analysis_date <- "2026-07-17"
args_cli <- commandArgs(trailingOnly = TRUE)
get_cli_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args_cli[startsWith(args_cli, prefix)]
  if (length(hit) > 0) return(sub(prefix, "", hit[1], fixed = TRUE))
  default
}
project_root <- normalizePath(
  get_cli_arg("project-root", Sys.getenv("EOSCZ_PROJECT_ROOT", getwd())),
  winslash = "/",
  mustWork = FALSE
)
cnv_root <- get_cli_arg("cnv-root", file.path(project_root, "cnv_analysis/protein_coding_genes"))
burden_root <- get_cli_arg("burden-root", file.path(project_root, "cnv_analysis/pathway_burden"))
out_dir <- file.path(burden_root, "controlRare_noSCZrecurrent_comparison_and_public_reference_logistic_case_cohort_26-7-17")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PUBLIC_BINARY_SAMPLE_CN <- "sample_CN_binary.xlsx"
LEGACY_BINARY_SAMPLE_CN <- "sample_copy_number_presence.xlsx"

sample_binary_path <- function(cohort_dir) {
  public <- file.path(cohort_dir, PUBLIC_BINARY_SAMPLE_CN)
  if (file.exists(public)) return(public)
  legacy <- file.path(cohort_dir, LEGACY_BINARY_SAMPLE_CN)
  if (file.exists(legacy)) return(legacy)
  public
}

gene_set_json <- file.path(burden_root, "original_36_gene_sets_symbols_for_local_burden_unique_26-7-14.json")
scz_comparison_cohort_audit <- file.path(project_root, "eoscz_ngs_hifi_rm_parents_sex_age_audit_26-7-16/SCZ_comparison_cohort_CNV_sample_haplotype_sex_age_match_all_records_26-7-16.tsv")
public_reference_sample_meta_file <- file.path(cnv_root, "slurm_scripts_public_reference/public_reference_metadata_26-7-15/public_reference_sample_sex_age_metadata_final_26-7-15.tsv")
public_reference_haplotype_meta_file <- file.path(cnv_root, "slurm_scripts_public_reference/public_reference_metadata_26-7-15/public_reference_haplotype_sex_age_metadata_final_26-7-15.tsv")

inputs <- list(
  sample = list(
    SCZ = sample_binary_path(file.path(cnv_root, "slurm_scripts_case_cohort")),
    comparison_cohort = sample_binary_path(file.path(cnv_root, "slurm_scripts_comparison_site")),
    public_reference = sample_binary_path(file.path(cnv_root, "slurm_scripts_public_reference"))
  ),
  haplotype = list(
    SCZ = file.path(cnv_root, "slurm_scripts_case_cohort/haplotype_CN.xlsx"),
    comparison_cohort = file.path(cnv_root, "slurm_scripts_comparison_site/haplotype_CN.xlsx"),
    public_reference = file.path(cnv_root, "slurm_scripts_public_reference/haplotype_CN.xlsx")
  )
)

comparisons <- list(
  filtered_cases_vs_comparison_cohort_public_reference = c("comparison_cohort", "public_reference")
)

empty_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == "" | toupper(trimws(as.character(x))) %in% c(
    "NA", "N/A", "UNKNOWN", "NO DATA", "MISSING", "RESTRICTED ACCESS",
    "NOT PROVIDED", "NOT COLLECTED", "NOT APPLICABLE", "NONE", "AGE", "unknown", "unspecified"
  )
}

clean_gene <- function(x) toupper(trimws(as.character(x)))

normalize_scz_comparison_cohort_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.scaffold$", "", x)
  x <- gsub("_2_(Pat|Mat)\\.v[0-9.]+$", "", x, ignore.case = TRUE)
  x <- gsub("_(Pat|Mat)\\.v[0-9.]+$", "", x, ignore.case = TRUE)
  x <- gsub("_(Pat|Mat)$", "", x, ignore.case = TRUE)
  x <- gsub("\\.[12][._](mat|pat)[._]R1_t2t(\\.[12])?$", "", x, ignore.case = TRUE)
  x <- gsub("\\.[12]$", "", x)
  x
}

normalise_sex <- function(x) {
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[y %in% c("male", "m", "1")] <- "male"
  out[y %in% c("female", "f", "2")] <- "female"
  out
}

read_binary_matrix <- function(path) {
  sheet <- readxl::excel_sheets(path)[1]
  raw <- readxl::read_excel(path, sheet = sheet, .name_repair = "minimal")
  ids <- trimws(as.character(raw[[1]]))
  mat <- as.data.frame(raw[, -1, drop = FALSE], check.names = FALSE)
  names(mat) <- clean_gene(names(mat))
  mat[] <- lapply(mat, function(x) {
    y <- suppressWarnings(as.numeric(x))
    y[is.na(y)] <- 0
    as.integer(y > 0)
  })
  keep <- !empty_value(ids)
  mat <- mat[keep, , drop = FALSE]
  ids <- ids[keep]
  if (anyDuplicated(names(mat))) {
    split_cols <- split(seq_along(names(mat)), names(mat))
    mat <- as.data.frame(
      lapply(split_cols, function(idx) as.integer(rowSums(mat[, idx, drop = FALSE]) > 0)),
      check.names = FALSE
    )
  }
  rownames(mat) <- ids
  if (anyDuplicated(rownames(mat))) {
    current_genes <- colnames(mat)
    split_rows <- split(seq_len(nrow(mat)), rownames(mat))
    mat <- as.data.frame(
      t(vapply(split_rows, function(idx) as.integer(colSums(mat[idx, , drop = FALSE]) > 0), integer(ncol(mat)))),
      check.names = FALSE
    )
    colnames(mat) <- current_genes
  }
  mat
}

align_matrices <- function(mats) {
  all_genes <- sort(unique(unlist(lapply(mats, colnames))))
  lapply(mats, function(mat) {
    missing <- setdiff(all_genes, colnames(mat))
    if (length(missing) > 0) {
      zeros <- as.data.frame(matrix(0L, nrow = nrow(mat), ncol = length(missing)), check.names = FALSE)
      colnames(zeros) <- missing
      mat <- cbind(mat, zeros)
    }
    mat[, all_genes, drop = FALSE]
  })
}

read_scz_comparison_cohort_sex_meta <- function(level, cohort) {
  aud <- read.delim(scz_comparison_cohort_audit, check.names = FALSE, stringsAsFactors = FALSE)
  aud <- aud[aud$record_type == level & aud$cohort == cohort, , drop = FALSE]
  data.frame(
    record_id = aud$record_id,
    sample_id = aud$sample_id,
    sex = normalise_sex(aud$sex),
    sex_source = aud$source_sheet,
    stringsAsFactors = FALSE
  )
}

read_public_reference_sex_meta <- function(level) {
  if (level == "sample") {
    x <- read.delim(public_reference_sample_meta_file, check.names = FALSE, stringsAsFactors = FALSE)
    data.frame(
      record_id = x$sample_id,
      sample_id = x$sample_id,
      sex = normalise_sex(x$sex_final),
      sex_source = x$sex_final_source,
      stringsAsFactors = FALSE
    )
  } else {
    x <- read.delim(public_reference_haplotype_meta_file, check.names = FALSE, stringsAsFactors = FALSE)
    data.frame(
      record_id = x$haplotype_id,
      sample_id = x$sample_id,
      sex = normalise_sex(x$sex_final),
      sex_source = x$sex_final_source,
      stringsAsFactors = FALSE
    )
  }
}

make_sex_meta <- function(level) {
  rbind(
    transform(read_scz_comparison_cohort_sex_meta(level, "SCZ_filtered_cases"), cohort = "SCZ"),
    transform(read_scz_comparison_cohort_sex_meta(level, "comparison_cohort_comparison_site"), cohort = "comparison_cohort"),
    transform(read_public_reference_sex_meta(level), cohort = "public_reference")
  )
}

safe_glm <- function(formula, data) {
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      glm(formula, data = data, family = binomial(link = "logit")),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(fit = fit, warnings = unique(warnings))
}

logistic_deviance_result <- function(df, predictor) {
  if (length(unique(df[[predictor]])) < 2 || sum(df[[predictor]] > 0, na.rm = TRUE) == 0) {
    return(list(beta = NA_real_, se = NA_real_, or = NA_real_, p_wald = 1,
                p_dev = 1, deviance_chisq = 0, df = 1, converged = NA,
                warning = "Predictor has no usable variation"))
  }
  df$sex <- factor(df$sex)
  if (nlevels(df$sex) < 2) {
    return(list(beta = NA_real_, se = NA_real_, or = NA_real_, p_wald = 1,
                p_dev = 1, deviance_chisq = NA_real_, df = 1, converged = NA,
                warning = "Sex covariate has <2 levels after filtering"))
  }
  null <- safe_glm(case_status ~ sex + global_rare_CNV_gene_count, df)
  full <- safe_glm(as.formula(paste("case_status ~ sex + global_rare_CNV_gene_count +", predictor)), df)
  if (inherits(null$fit, "error") || inherits(full$fit, "error")) {
    msg <- paste(c(
      if (inherits(null$fit, "error")) paste("null:", conditionMessage(null$fit)),
      if (inherits(full$fit, "error")) paste("full:", conditionMessage(full$fit))
    ), collapse = "; ")
    return(list(beta = NA_real_, se = NA_real_, or = NA_real_, p_wald = 1,
                p_dev = 1, deviance_chisq = NA_real_, df = 1, converged = FALSE,
                warning = msg))
  }
  an <- tryCatch(anova(null$fit, full$fit, test = "Chisq"), error = function(e) e)
  sm <- summary(full$fit)
  coef_table <- sm$coefficients
  beta <- if (predictor %in% rownames(coef_table)) coef_table[predictor, "Estimate"] else NA_real_
  se <- if (predictor %in% rownames(coef_table)) coef_table[predictor, "Std. Error"] else NA_real_
  p_wald <- if (predictor %in% rownames(coef_table)) coef_table[predictor, "Pr(>|z|)"] else NA_real_
  if (inherits(an, "error")) {
    p_dev <- 1
    dev <- NA_real_
    an_df <- 1
    warn <- c(null$warnings, full$warnings, conditionMessage(an))
  } else {
    p_dev <- an$`Pr(>Chi)`[2]
    dev <- an$Deviance[2]
    an_df <- an$Df[2]
    if (is.na(p_dev)) p_dev <- 1
    warn <- c(null$warnings, full$warnings)
  }
  list(beta = beta, se = se, or = ifelse(is.na(beta), NA_real_, exp(beta)),
       p_wald = p_wald, p_dev = p_dev, deviance_chisq = dev, df = an_df,
       converged = isTRUE(full$fit$converged), warning = paste(unique(warn), collapse = " | "))
}

fisher_result <- function(case_binary, ctrl_binary) {
  a <- sum(case_binary > 0)
  b <- length(case_binary) - a
  c <- sum(ctrl_binary > 0)
  d <- length(ctrl_binary) - c
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft_greater <- tryCatch(fisher.test(tab, alternative = "greater"), error = function(e) e)
  ft_two <- tryCatch(fisher.test(tab, alternative = "two.sided"), error = function(e) e)
  list(
    scz_carriers = a,
    ctrl_carriers = c,
    scz_noncarriers = b,
    ctrl_noncarriers = d,
    fisher_or = if (inherits(ft_greater, "error")) NA_real_ else unname(ft_greater$estimate),
    fisher_p_greater = if (inherits(ft_greater, "error")) 1 else ft_greater$p.value,
    fisher_p_two_sided = if (inherits(ft_two, "error")) 1 else ft_two$p.value
  )
}

run_one <- function(level, comparison_name, ctrl_names, mats, sex_meta, gene_sets) {
  case_mat <- mats$SCZ
  ctrl_list <- unname(lapply(ctrl_names, function(x) mats[[x]]))
  ctrl_mat <- do.call(rbind, ctrl_list)
  aligned <- align_matrices(list(case = case_mat, ctrl = ctrl_mat, comparison_cohort = mats$comparison_cohort, public_reference = mats$public_reference, SCZ = mats$SCZ))
  case_mat <- aligned$case
  ctrl_mat <- aligned$ctrl
  comparison_cohort_mat <- aligned$comparison_cohort
  public_reference_mat <- aligned$public_reference
  scz_mat <- aligned$SCZ

  comparison_cohort_count <- colSums(comparison_cohort_mat > 0)
  public_reference_count <- colSums(public_reference_mat > 0)
  scz_count <- colSums(scz_mat > 0)
  comparison_cohort_freq <- comparison_cohort_count / nrow(comparison_cohort_mat)
  public_reference_freq <- public_reference_count / nrow(public_reference_mat)
  rare_genes <- names(scz_count)[comparison_cohort_freq < 0.01 & public_reference_freq < 0.01]
  rare_filter <- data.frame(
    Level = level,
    Gene = names(scz_count),
    Pass_filter = names(scz_count) %in% rare_genes,
    SCZ_count = as.integer(scz_count),
    SCZ_frequency = as.numeric(scz_count / nrow(scz_mat)),
    comparison_cohort_count = as.integer(comparison_cohort_count),
    comparison_cohort_frequency = as.numeric(comparison_cohort_freq),
    public_reference_count = as.integer(public_reference_count),
    public_reference_frequency = as.numeric(public_reference_freq),
    stringsAsFactors = FALSE
  )

  record_ids <- c(rownames(case_mat), rownames(ctrl_mat))
  group <- c(rep("SCZ", nrow(case_mat)), rep("CTRL", nrow(ctrl_mat)))
  ctrl_source <- if (length(ctrl_names) == 1) rep(ctrl_names, nrow(ctrl_mat)) else {
    unlist(lapply(ctrl_names, function(x) rep(x, nrow(mats[[x]]))), use.names = FALSE)
  }
  cohort <- c(rep("SCZ", nrow(case_mat)), ctrl_source)
  pheno <- data.frame(
    record_id = record_ids,
    group = group,
    cohort = cohort,
    case_status = c(rep(1L, nrow(case_mat)), rep(0L, nrow(ctrl_mat))),
    stringsAsFactors = FALSE
  )
  meta <- sex_meta[, c("record_id", "sample_id", "sex", "sex_source", "cohort")]
  colnames(meta)[5] <- "metadata_cohort"
  pheno <- merge(pheno, meta, by = "record_id", all.x = TRUE, sort = FALSE)
  rare_combined <- rbind(case_mat[, rare_genes, drop = FALSE], ctrl_mat[, rare_genes, drop = FALSE])
  pheno$global_rare_CNV_gene_count <- rowSums(rare_combined)
  pheno$usable_for_sex_logistic <- !empty_value(pheno$sex)

  results <- list()
  contributors <- list()
  for (set_name in names(gene_sets)) {
    set_genes <- intersect(clean_gene(gene_sets[[set_name]]), colnames(case_mat))
    rare_set_genes <- intersect(set_genes, rare_genes)
    case_count <- if (length(rare_set_genes) == 0) rep(0L, nrow(case_mat)) else rowSums(case_mat[, rare_set_genes, drop = FALSE])
    ctrl_count <- if (length(rare_set_genes) == 0) rep(0L, nrow(ctrl_mat)) else rowSums(ctrl_mat[, rare_set_genes, drop = FALSE])

    model_df <- pheno
    model_df$gene_set_rare_CNV_gene_count <- c(case_count, ctrl_count)
    model_df$gene_set_rare_CNV_binary <- as.integer(model_df$gene_set_rare_CNV_gene_count > 0)
    fit_df <- model_df[model_df$usable_for_sex_logistic, , drop = FALSE]

    count_res <- logistic_deviance_result(fit_df, "gene_set_rare_CNV_gene_count")
    binary_res <- logistic_deviance_result(fit_df, "gene_set_rare_CNV_binary")
    fish <- fisher_result(c(case_count > 0), c(ctrl_count > 0))

    results[[length(results) + 1]] <- data.frame(
      Level = level,
      Comparison = comparison_name,
      RareDefinition = paste0(level, "_comparison_cohort_frequency<0.01 AND public_reference_frequency<0.01"),
      RareGenesTotal = length(rare_genes),
      GeneSet = set_name,
      GeneSetSize_original_unique_symbols = length(unique(clean_gene(gene_sets[[set_name]]))),
      GeneSetGenes_present_in_CNV_matrix = length(set_genes),
      RareGenesInSet = length(rare_set_genes),
      SCZ_total_records = nrow(case_mat),
      CTRL_total_records = nrow(ctrl_mat),
      Logistic_records_used = nrow(fit_df),
      Logistic_records_dropped_missing_sex = nrow(model_df) - nrow(fit_df),
      SCZ_samples_ge1_gene = fish$scz_carriers,
      CTRL_samples_ge1_gene = fish$ctrl_carriers,
      SCZ_samples_ge2_genes = sum(case_count >= 2),
      CTRL_samples_ge2_genes = sum(ctrl_count >= 2),
      SCZ_mean_gene_set_rare_CNV_gene_count = mean(case_count),
      CTRL_mean_gene_set_rare_CNV_gene_count = mean(ctrl_count),
      Fisher_OR_binary = fish$fisher_or,
      Fisher_P_greater_SCZ_gt_CTRL = fish$fisher_p_greater,
      Fisher_P_two_sided = fish$fisher_p_two_sided,
      LogisticModel_count = "case_status ~ sex + global_rare_CNV_gene_count + gene_set_rare_CNV_gene_count",
      Beta_count = count_res$beta,
      SE_count = count_res$se,
      OR_per_additional_gene_set_gene = count_res$or,
      Deviance_Chisq_count = count_res$deviance_chisq,
      Deviance_df_count = count_res$df,
      Logistic_P_dev_count = count_res$p_dev,
      Logistic_P_wald_count = count_res$p_wald,
      Converged_count = count_res$converged,
      Warning_count = count_res$warning,
      LogisticModel_binary = "case_status ~ sex + global_rare_CNV_gene_count + gene_set_rare_CNV_binary",
      Beta_binary = binary_res$beta,
      SE_binary = binary_res$se,
      OR_binary_ge1_gene = binary_res$or,
      Deviance_Chisq_binary = binary_res$deviance_chisq,
      Deviance_df_binary = binary_res$df,
      Logistic_P_dev_binary = binary_res$p_dev,
      Logistic_P_wald_binary = binary_res$p_wald,
      Converged_binary = binary_res$converged,
      Warning_binary = binary_res$warning,
      stringsAsFactors = FALSE
    )

    if (length(rare_set_genes) > 0) {
      case_carriers <- colSums(case_mat[, rare_set_genes, drop = FALSE] > 0)
      ctrl_carriers <- colSums(ctrl_mat[, rare_set_genes, drop = FALSE] > 0)
      contributors[[length(contributors) + 1]] <- data.frame(
        Level = level,
        Comparison = comparison_name,
        GeneSet = set_name,
        Gene = rare_set_genes,
        SCZ_carrier_records = as.integer(case_carriers[rare_set_genes]),
        CTRL_carrier_records = as.integer(ctrl_carriers[rare_set_genes]),
        comparison_cohort_frequency = as.numeric(comparison_cohort_freq[rare_set_genes]),
        public_reference_frequency = as.numeric(public_reference_freq[rare_set_genes]),
        SCZ_filter_count = as.integer(scz_count[rare_set_genes]),
        stringsAsFactors = FALSE
      )
    }
  }

  res <- do.call(rbind, results)
  res$FDR_BH_Fisher_greater <- p.adjust(res$Fisher_P_greater_SCZ_gt_CTRL, method = "BH")
  res$Bonferroni_Fisher_greater <- p.adjust(res$Fisher_P_greater_SCZ_gt_CTRL, method = "bonferroni")
  res$FDR_BH_Logistic_dev_count <- p.adjust(res$Logistic_P_dev_count, method = "BH")
  res$Bonferroni_Logistic_dev_count <- p.adjust(res$Logistic_P_dev_count, method = "bonferroni")
  res$FDR_BH_Logistic_dev_binary <- p.adjust(res$Logistic_P_dev_binary, method = "BH")
  res$Bonferroni_Logistic_dev_binary <- p.adjust(res$Logistic_P_dev_binary, method = "bonferroni")
  res <- res[order(res$Logistic_P_dev_count, res$Fisher_P_greater_SCZ_gt_CTRL), ]
  contrib <- if (length(contributors) > 0) do.call(rbind, contributors) else data.frame()
  list(results = res, contributors = contrib, pheno = pheno, rare_filter = rare_filter)
}

write_one_workbook <- function(path, result, meta) {
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  writeData(wb, "Results", result$results)
  addWorksheet(wb, "ContributingGenes")
  writeData(wb, "ContributingGenes", result$contributors)
  addWorksheet(wb, "PhenotypeBurden")
  writeData(wb, "PhenotypeBurden", result$pheno)
  addWorksheet(wb, "RareGeneFilter")
  writeData(wb, "RareGeneFilter", result$rare_filter)
  addWorksheet(wb, "Meta")
  writeData(wb, "Meta", meta)
  for (sh in names(wb)) {
    freezePane(wb, sh, firstRow = TRUE)
    setColWidths(wb, sh, cols = 1:30, widths = "auto")
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

gene_sets_raw <- jsonlite::fromJSON(gene_set_json, simplifyVector = FALSE)
gene_sets <- gene_sets_raw$gene_sets
all_results <- list()
all_meta <- list()

for (level in c("sample", "haplotype")) {
  mats <- align_matrices(lapply(inputs[[level]], read_binary_matrix))
  sex_meta <- make_sex_meta(level)
  for (comparison_name in names(comparisons)) {
    ctrl_names <- comparisons[[comparison_name]]
    result <- run_one(level, comparison_name, ctrl_names, mats, sex_meta, gene_sets)
    meta <- data.frame(
      Metric = c(
        "analysis_date", "level", "comparison", "case_input", "control_inputs",
        "gene_set_json", "rare_definition", "primary_logistic_model",
        "fisher_test", "sex_covariate_source", "global_covariate",
        "important_note"
      ),
      Value = c(
        analysis_date, level, comparison_name, inputs[[level]]$SCZ,
        paste(unlist(inputs[[level]][ctrl_names]), collapse = "; "),
        gene_set_json,
        paste0(level, "_comparison_cohort_frequency<0.01 AND public_reference_frequency<0.01"),
        "Deviance test comparing case_status ~ sex + global_rare_CNV_gene_count vs case_status ~ sex + global_rare_CNV_gene_count + gene_set_rare_CNV_gene_count",
        "One-sided Fisher exact test on binary carrier status, alternative SCZ > CTRL; no covariates",
        "SCZ/comparison_cohort from eoscz_ngs_hifi_rm_parents_sample_metadata.xlsx; public_reference from curated public_reference metadata with selected_public_reference_samples local idxstats sex supplement; records without sex excluded from logistic only",
        "global_rare_CNV_gene_count = total number of rare CNV-related genes per record after rare filter",
        "Exploratory control-rare model: comparison_cohort and public_reference frequencies are both <1%; SCZ recurrence is not required"
      ),
      stringsAsFactors = FALSE
    )
    prefix <- paste(level, "comparison_and_public_reference_controlRare_noSCZrec_26-7-17", sep = "_")
    xlsx <- file.path(out_dir, paste0(prefix, ".xlsx"))
    write_one_workbook(xlsx, result, meta)
    write.table(result$results, file.path(out_dir, paste0(prefix, "_results.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
    write.table(result$contributors, file.path(out_dir, paste0(prefix, "_contributing_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
    write.table(result$pheno, file.path(out_dir, paste0(prefix, "_phenotype_burden.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
    all_results[[length(all_results) + 1]] <- result$results
    all_meta[[length(all_meta) + 1]] <- meta
  }
}

all_res <- do.call(rbind, all_results)
write.table(all_res, file.path(out_dir, "all_results_36geneSets_controlRare_noSCZrecurrent_logistic_26-7-17.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

sig_summary <- data.frame(
  Level = character(), Comparison = character(), Test = character(),
  N_nominal_positive = integer(), N_FDR_positive = integer(), TopGeneSet = character(),
  TopP = numeric(), TopOR = numeric(), stringsAsFactors = FALSE
)
for (lv in unique(all_res$Level)) {
  for (cmp in unique(all_res$Comparison)) {
    sub <- all_res[all_res$Level == lv & all_res$Comparison == cmp, ]
    tests <- list(
      Fisher = list(p = "Fisher_P_greater_SCZ_gt_CTRL", fdr = "FDR_BH_Fisher_greater", or = "Fisher_OR_binary"),
      Logistic_count = list(p = "Logistic_P_dev_count", fdr = "FDR_BH_Logistic_dev_count", or = "OR_per_additional_gene_set_gene"),
      Logistic_binary = list(p = "Logistic_P_dev_binary", fdr = "FDR_BH_Logistic_dev_binary", or = "OR_binary_ge1_gene")
    )
    for (tn in names(tests)) {
      pp <- tests[[tn]]$p
      ff <- tests[[tn]]$fdr
      oo <- tests[[tn]]$or
      ord <- order(sub[[pp]], decreasing = FALSE)
      top <- sub[ord[1], ]
      sig_summary <- rbind(sig_summary, data.frame(
        Level = lv,
        Comparison = cmp,
        Test = tn,
        N_nominal_positive = sum(sub[[pp]] < 0.05 & !is.na(sub[[oo]]) & sub[[oo]] > 1, na.rm = TRUE),
        N_FDR_positive = sum(sub[[ff]] < 0.05 & !is.na(sub[[oo]]) & sub[[oo]] > 1, na.rm = TRUE),
        TopGeneSet = top$GeneSet,
        TopP = top[[pp]],
        TopOR = top[[oo]],
        stringsAsFactors = FALSE
      ))
    }
  }
}
write.table(sig_summary, file.path(out_dir, "summary_significance_counts_controlRare_noSCZrecurrent_26-7-17.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

wb <- createWorkbook()
addWorksheet(wb, "AllResults")
writeData(wb, "AllResults", all_res)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", sig_summary)
addWorksheet(wb, "Nominal_Logistic_count")
writeData(wb, "Nominal_Logistic_count", subset(all_res, Logistic_P_dev_count < 0.05 & OR_per_additional_gene_set_gene > 1))
addWorksheet(wb, "Nominal_Fisher")
writeData(wb, "Nominal_Fisher", subset(all_res, Fisher_P_greater_SCZ_gt_CTRL < 0.05 & Fisher_OR_binary > 1))
addWorksheet(wb, "FDR_Logistic_count")
writeData(wb, "FDR_Logistic_count", subset(all_res, FDR_BH_Logistic_dev_count < 0.05 & OR_per_additional_gene_set_gene > 1))
addWorksheet(wb, "FDR_Fisher")
writeData(wb, "FDR_Fisher", subset(all_res, FDR_BH_Fisher_greater < 0.05 & Fisher_OR_binary > 1))
addWorksheet(wb, "Meta")
writeData(wb, "Meta", do.call(rbind, all_meta))
for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:35, widths = "auto")
}
saveWorkbook(wb, file.path(out_dir, "summary_36geneSets_controlRare_noSCZrecurrent_logistic_26-7-17.xlsx"), overwrite = TRUE)

cat("DONE\n")
cat("Output directory:", out_dir, "\n")
print(sig_summary)


