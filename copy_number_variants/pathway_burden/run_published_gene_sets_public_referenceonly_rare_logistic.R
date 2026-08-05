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
out_dir <- file.path(burden_root, "public_referenceOnly_rare_logistic_case_cohort_26-7-17")
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
gene_set_group_file <- file.path(burden_root, "Supplementary_Table_3_gene_sets_summary.tsv")
scz_comparison_cohort_audit <- file.path(project_root, "eoscz_ngs_hifi_rm_parents_sex_age_audit_26-7-16/SCZ_comparison_cohort_CNV_sample_haplotype_sex_age_match_all_records_26-7-16.tsv")
public_reference_sample_meta_file <- file.path(cnv_root, "slurm_scripts_public_reference/public_reference_metadata_26-7-15/public_reference_sample_sex_age_metadata_final_26-7-15.tsv")
public_reference_haplotype_meta_file <- file.path(cnv_root, "slurm_scripts_public_reference/public_reference_metadata_26-7-15/public_reference_haplotype_sex_age_metadata_final_26-7-15.tsv")

inputs <- list(
  sample = list(
    SCZ = sample_binary_path(file.path(cnv_root, "slurm_scripts_case_cohort")),
    public_reference = sample_binary_path(file.path(cnv_root, "slurm_scripts_public_reference"))
  ),
  haplotype = list(
    SCZ = file.path(cnv_root, "slurm_scripts_case_cohort/haplotype_CN.xlsx"),
    public_reference = file.path(cnv_root, "slurm_scripts_public_reference/haplotype_CN.xlsx")
  )
)

empty_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == "" | toupper(trimws(as.character(x))) %in% c(
    "NA", "N/A", "UNKNOWN", "NO DATA", "MISSING", "RESTRICTED ACCESS",
    "NOT PROVIDED", "NOT COLLECTED", "NOT APPLICABLE", "NONE", "AGE", "unknown", "unspecified"
  )
}

clean_gene <- function(x) toupper(trimws(as.character(x)))

normalise_sex <- function(x) {
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[y %in% c("male", "m", "1")] <- "male"
  out[y %in% c("female", "f", "2")] <- "female"
  out
}

read_binary_matrix <- function(path) {
  if (!file.exists(path)) stop("Input not found: ", path)
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

read_scz_sex_meta <- function(level) {
  aud <- read.delim(scz_comparison_cohort_audit, check.names = FALSE, stringsAsFactors = FALSE)
  aud <- aud[aud$record_type == level & aud$cohort == "SCZ_filtered_cases", , drop = FALSE]
  data.frame(
    record_id = aud$record_id,
    sample_id = aud$sample_id,
    sex = normalise_sex(aud$sex),
    sex_source = aud$source_sheet,
    cohort = "SCZ",
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
      cohort = "public_reference",
      stringsAsFactors = FALSE
    )
  } else {
    x <- read.delim(public_reference_haplotype_meta_file, check.names = FALSE, stringsAsFactors = FALSE)
    data.frame(
      record_id = x$haplotype_id,
      sample_id = x$sample_id,
      sex = normalise_sex(x$sex_final),
      sex_source = x$sex_final_source,
      cohort = "public_reference",
      stringsAsFactors = FALSE
    )
  }
}

make_sex_meta <- function(level) {
  rbind(read_scz_sex_meta(level), read_public_reference_sex_meta(level))
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

run_level <- function(level, gene_sets) {
  mats <- align_matrices(lapply(inputs[[level]], read_binary_matrix))
  case_mat <- mats$SCZ
  ctrl_mat <- mats$public_reference
  public_reference_count <- colSums(ctrl_mat > 0)
  scz_count <- colSums(case_mat > 0)
  public_reference_freq <- public_reference_count / nrow(ctrl_mat)
  rare_genes <- names(scz_count)[public_reference_freq < 0.01 & scz_count >= 2]
  rare_filter <- data.frame(
    Level = level,
    Gene = names(scz_count),
    Pass_filter = names(scz_count) %in% rare_genes,
    SCZ_count = as.integer(scz_count),
    SCZ_frequency = as.numeric(scz_count / nrow(case_mat)),
    public_reference_count = as.integer(public_reference_count),
    public_reference_frequency = as.numeric(public_reference_freq),
    stringsAsFactors = FALSE
  )

  pheno <- data.frame(
    record_id = c(rownames(case_mat), rownames(ctrl_mat)),
    group = c(rep("SCZ", nrow(case_mat)), rep("CTRL", nrow(ctrl_mat))),
    cohort = c(rep("SCZ", nrow(case_mat)), rep("public_reference", nrow(ctrl_mat))),
    case_status = c(rep(1L, nrow(case_mat)), rep(0L, nrow(ctrl_mat))),
    stringsAsFactors = FALSE
  )
  sex_meta <- make_sex_meta(level)
  pheno <- merge(pheno, sex_meta[, c("record_id", "sample_id", "sex", "sex_source", "cohort")],
                 by = "record_id", all.x = TRUE, sort = FALSE)
  colnames(pheno)[colnames(pheno) == "cohort.x"] <- "cohort"
  colnames(pheno)[colnames(pheno) == "cohort.y"] <- "metadata_cohort"
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

    results[[length(results) + 1]] <- data.frame(
      Level = level,
      Comparison = "filtered_cases_vs_public_reference",
      RareDefinition = paste0(level, "_public_reference_frequency<0.01 AND SCZ_count>=2"),
      RareGenesTotal = length(rare_genes),
      GeneSet = set_name,
      GeneSetSize_original_unique_symbols = length(unique(clean_gene(gene_sets[[set_name]]))),
      GeneSetGenes_present_in_CNV_matrix = length(set_genes),
      RareGenesInSet = length(rare_set_genes),
      SCZ_total_records = nrow(case_mat),
      CTRL_total_records = nrow(ctrl_mat),
      Logistic_records_used = nrow(fit_df),
      Logistic_records_dropped_missing_sex = nrow(model_df) - nrow(fit_df),
      SCZ_records_ge1_gene = sum(case_count > 0),
      CTRL_records_ge1_gene = sum(ctrl_count > 0),
      SCZ_records_ge2_genes = sum(case_count >= 2),
      CTRL_records_ge2_genes = sum(ctrl_count >= 2),
      SCZ_mean_gene_set_rare_CNV_gene_count = mean(case_count),
      CTRL_mean_gene_set_rare_CNV_gene_count = mean(ctrl_count),
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
      contributors[[length(contributors) + 1]] <- data.frame(
        Level = level,
        Comparison = "filtered_cases_vs_public_reference",
        GeneSet = set_name,
        Gene = rare_set_genes,
        SCZ_carrier_records = as.integer(colSums(case_mat[, rare_set_genes, drop = FALSE] > 0)[rare_set_genes]),
        CTRL_carrier_records = as.integer(colSums(ctrl_mat[, rare_set_genes, drop = FALSE] > 0)[rare_set_genes]),
        public_reference_frequency = as.numeric(public_reference_freq[rare_set_genes]),
        SCZ_filter_count = as.integer(scz_count[rare_set_genes]),
        stringsAsFactors = FALSE
      )
    }
  }

  res <- do.call(rbind, results)
  groups <- read.delim(gene_set_group_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- groups[, c("GeneSet ID (Suppl DataSets)", "GeneSet Group", "Figure Label", "GeneSet FullName")]
  colnames(groups) <- c("GeneSet", "GeneSetGroup", "FigureLabel", "GeneSetFullName")
  res <- merge(res, groups, by = "GeneSet", all.x = TRUE, sort = FALSE)
  res$CNV_type_for_FDR <- "Liftoff_extra_copy_or_CNV_gene_presence"
  res$FDR_BH_Logistic_dev_count_36wide <- p.adjust(res$Logistic_P_dev_count, method = "BH")
  res$Bonferroni_Logistic_dev_count_36wide <- p.adjust(res$Logistic_P_dev_count, method = "bonferroni")
  res$FDR_BH_Logistic_dev_binary_36wide <- p.adjust(res$Logistic_P_dev_binary, method = "BH")
  res$Bonferroni_Logistic_dev_binary_36wide <- p.adjust(res$Logistic_P_dev_binary, method = "bonferroni")
  res$FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype <- NA_real_
  res$FDR_BH_Logistic_dev_binary_by_GeneSetGroup_CNVtype <- NA_real_
  for (grp in unique(res$GeneSetGroup)) {
    idx <- which(res$GeneSetGroup == grp)
    res$FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype[idx] <- p.adjust(res$Logistic_P_dev_count[idx], method = "BH")
    res$FDR_BH_Logistic_dev_binary_by_GeneSetGroup_CNVtype[idx] <- p.adjust(res$Logistic_P_dev_binary[idx], method = "BH")
  }
  res <- res[order(res$Logistic_P_dev_count, res$Logistic_P_dev_binary), ]
  list(
    results = res,
    contributors = if (length(contributors) > 0) do.call(rbind, contributors) else data.frame(),
    pheno = pheno,
    rare_filter = rare_filter
  )
}

write_workbook <- function(path, result, meta) {
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
    setColWidths(wb, sh, cols = 1:40, widths = "auto")
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

gene_sets <- jsonlite::fromJSON(gene_set_json, simplifyVector = FALSE)$gene_sets
all_results <- list()
for (level in c("sample", "haplotype")) {
  result <- run_level(level, gene_sets)
  meta <- data.frame(
    Metric = c("analysis_date", "level", "comparison", "case_input", "control_input",
               "gene_set_json", "rare_definition", "primary_logistic_model",
               "groupwise_fdr_note", "important_note"),
    Value = c(analysis_date, level, "filtered_cases_vs_public_reference", inputs[[level]]$SCZ, inputs[[level]]$public_reference,
              gene_set_json, paste0(level, "_public_reference_frequency<0.01 AND SCZ_count>=2"),
              "Deviance test comparing case_status ~ sex + global_rare_CNV_gene_count vs case_status ~ sex + global_rare_CNV_gene_count + gene_set_rare_CNV_gene_count",
              "BH-FDR calculated both across all 36 gene sets and within original GeneSetGroup as a single Liftoff CNV type",
              "Exploratory SCZ vs public_reference-only analysis; comparison_cohort is excluded from both control group and rare-frequency filter"),
    stringsAsFactors = FALSE
  )
  prefix <- paste(level, "filtered_cases_vs_public_reference", "36geneSets_public_referenceOnly_rare_logistic_26-7-17", sep = "_")
  write_workbook(file.path(out_dir, paste0(prefix, ".xlsx")), result, meta)
  write.table(result$results, file.path(out_dir, paste0(prefix, "_results.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  write.table(result$contributors, file.path(out_dir, paste0(prefix, "_contributing_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  write.table(result$pheno, file.path(out_dir, paste0(prefix, "_phenotype_burden.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  all_results[[length(all_results) + 1]] <- result$results
}

all_res <- do.call(rbind, all_results)
write.table(all_res, file.path(out_dir, "all_results_36geneSets_public_referenceOnly_rare_logistic_26-7-17.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

summary_rows <- list()
for (lv in unique(all_res$Level)) {
  sub <- all_res[all_res$Level == lv, ]
  for (test in c("count", "binary")) {
    p_col <- if (test == "count") "Logistic_P_dev_count" else "Logistic_P_dev_binary"
    or_col <- if (test == "count") "OR_per_additional_gene_set_gene" else "OR_binary_ge1_gene"
    fdr36_col <- if (test == "count") "FDR_BH_Logistic_dev_count_36wide" else "FDR_BH_Logistic_dev_binary_36wide"
    fdrgrp_col <- if (test == "count") "FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype" else "FDR_BH_Logistic_dev_binary_by_GeneSetGroup_CNVtype"
    ord <- order(sub[[p_col]])
    top <- sub[ord[1], ]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Level = lv,
      Comparison = "filtered_cases_vs_public_reference",
      Test = paste0("Logistic_", test),
      RareGenesTotal = unique(sub$RareGenesTotal),
      SCZ_total_records = unique(sub$SCZ_total_records),
      CTRL_total_records = unique(sub$CTRL_total_records),
      Logistic_records_used = unique(sub$Logistic_records_used),
      Logistic_records_dropped_missing_sex = unique(sub$Logistic_records_dropped_missing_sex),
      N_nominal_SCZ_gt_CTRL = sum(sub[[p_col]] < 0.05 & sub[[or_col]] > 1, na.rm = TRUE),
      N_FDR05_36wide_SCZ_gt_CTRL = sum(sub[[fdr36_col]] <= 0.05 & sub[[or_col]] > 1, na.rm = TRUE),
      N_FDR10_36wide_SCZ_gt_CTRL = sum(sub[[fdr36_col]] <= 0.10 & sub[[or_col]] > 1, na.rm = TRUE),
      N_FDR05_groupwise_SCZ_gt_CTRL = sum(sub[[fdrgrp_col]] <= 0.05 & sub[[or_col]] > 1, na.rm = TRUE),
      N_FDR10_groupwise_SCZ_gt_CTRL = sum(sub[[fdrgrp_col]] <= 0.10 & sub[[or_col]] > 1, na.rm = TRUE),
      N_warnings = sum(if (test == "count") sub$Warning_count != "" else sub$Warning_binary != "", na.rm = TRUE),
      TopGeneSet = top$GeneSet,
      TopP = top[[p_col]],
      TopOR = top[[or_col]],
      TopFDR36 = top[[fdr36_col]],
      TopFDRGroupwise = top[[fdrgrp_col]],
      stringsAsFactors = FALSE
    )
  }
}
summary <- do.call(rbind, summary_rows)
write.table(summary, file.path(out_dir, "summary_public_referenceOnly_rare_logistic_26-7-17.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

wb <- createWorkbook()
addWorksheet(wb, "AllResults")
writeData(wb, "AllResults", all_res)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary)
addWorksheet(wb, "FDR10_groupwise_count")
writeData(wb, "FDR10_groupwise_count", subset(all_res, FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype <= 0.10 & OR_per_additional_gene_set_gene > 1))
addWorksheet(wb, "FDR10_groupwise_binary")
writeData(wb, "FDR10_groupwise_binary", subset(all_res, FDR_BH_Logistic_dev_binary_by_GeneSetGroup_CNVtype <= 0.10 & OR_binary_ge1_gene > 1))
for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:45, widths = "auto")
}
saveWorkbook(wb, file.path(out_dir, "summary_public_referenceOnly_rare_logistic_26-7-17.xlsx"), overwrite = TRUE)

cat("DONE\n")
cat("Output directory:", out_dir, "\n")
print(summary)


