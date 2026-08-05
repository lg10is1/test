options(stringsAsFactors = FALSE)

analysis_date <- "2026-07-16"
args_cli <- commandArgs(trailingOnly = TRUE)
get_cli_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args_cli[startsWith(args_cli, prefix)]
  if (length(hit) > 0) return(sub(prefix, "", hit[1], fixed = TRUE))
  default
}
default_project_root <- normalizePath(
  Sys.getenv("EOSCZ_PROJECT_ROOT", getwd()),
  winslash = "/",
  mustWork = FALSE
)
burden_root <- get_cli_arg("burden-root", file.path(default_project_root, "cnv_analysis/pathway_burden"))
in_dir <- file.path(burden_root, "sex_globalRare_logistic_fisher_case_cohort_26-7-16")
results_file <- file.path(in_dir, "all_results_36geneSets_sex_globalRare_logistic_fisher_26-7-16.tsv")
group_file <- file.path(burden_root, "Supplementary_Table_3_gene_sets_summary.tsv")
out_file <- file.path(in_dir, "all_results_36geneSets_with_originalStyle_groupwiseFDR_26-7-16.tsv")
summary_file <- file.path(in_dir, "summary_originalStyle_groupwiseFDR_26-7-16.tsv")

res <- read.delim(results_file, check.names = FALSE)
groups <- read.delim(group_file, check.names = FALSE)
groups <- groups[, c("GeneSet ID (Suppl DataSets)", "GeneSet Group", "Figure Label", "GeneSet FullName")]
colnames(groups) <- c("GeneSet", "GeneSetGroup", "FigureLabel", "GeneSetFullName")
res <- merge(res, groups, by = "GeneSet", all.x = TRUE, sort = FALSE)
res$CNV_type_for_FDR <- "Liftoff_extra_copy_or_CNV_gene_presence"

if (any(is.na(res$GeneSetGroup))) {
  warning("Missing GeneSetGroup for: ", paste(unique(res$GeneSet[is.na(res$GeneSetGroup)]), collapse = ", "))
  res$GeneSetGroup[is.na(res$GeneSetGroup)] <- "Unknown"
}

add_groupwise_fdr <- function(df, p_col, out_col) {
  df[[out_col]] <- NA_real_
  split_key <- interaction(df$Level, df$Comparison, df$CNV_type_for_FDR, df$GeneSetGroup, drop = TRUE, sep = "||")
  for (key in levels(split_key)) {
    idx <- which(split_key == key)
    df[[out_col]][idx] <- p.adjust(df[[p_col]][idx], method = "BH")
  }
  df
}

res <- add_groupwise_fdr(res, "Fisher_P_greater_SCZ_gt_CTRL", "FDR_BH_Fisher_greater_by_GeneSetGroup_CNVtype")
res <- add_groupwise_fdr(res, "Logistic_P_dev_count", "FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype")
res <- add_groupwise_fdr(res, "Logistic_P_dev_binary", "FDR_BH_Logistic_dev_binary_by_GeneSetGroup_CNVtype")

write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

rows <- list()
for (lv in unique(res$Level)) {
  for (cmp in unique(res$Comparison[res$Level == lv])) {
    sub <- res[res$Level == lv & res$Comparison == cmp, ]
    tests <- list(
      Fisher = list(p = "Fisher_P_greater_SCZ_gt_CTRL", fdr36 = "FDR_BH_Fisher_greater", fdrgrp = "FDR_BH_Fisher_greater_by_GeneSetGroup_CNVtype", or = "Fisher_OR_binary"),
      Logistic_count = list(p = "Logistic_P_dev_count", fdr36 = "FDR_BH_Logistic_dev_count", fdrgrp = "FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype", or = "OR_per_additional_gene_set_gene"),
      Logistic_binary = list(p = "Logistic_P_dev_binary", fdr36 = "FDR_BH_Logistic_dev_binary", fdrgrp = "FDR_BH_Logistic_dev_binary_by_GeneSetGroup_CNVtype", or = "OR_binary_ge1_gene")
    )
    for (tn in names(tests)) {
      spec <- tests[[tn]]
      rows[[length(rows) + 1]] <- data.frame(
        Level = lv,
        Comparison = cmp,
        Test = tn,
        N_gene_sets = nrow(sub),
        N_nominal_SCZ_gt_CTRL = sum(sub[[spec$p]] < 0.05 & sub[[spec$or]] > 1, na.rm = TRUE),
        N_FDR05_36wide = sum(sub[[spec$fdr36]] <= 0.05 & sub[[spec$or]] > 1, na.rm = TRUE),
        N_FDR10_36wide = sum(sub[[spec$fdr36]] <= 0.10 & sub[[spec$or]] > 1, na.rm = TRUE),
        N_FDR05_by_GeneSetGroup_CNVtype = sum(sub[[spec$fdrgrp]] <= 0.05 & sub[[spec$or]] > 1, na.rm = TRUE),
        N_FDR10_by_GeneSetGroup_CNVtype = sum(sub[[spec$fdrgrp]] <= 0.10 & sub[[spec$or]] > 1, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}
summary <- do.call(rbind, rows)
write.table(summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

cat("OUTPUT\t", out_file, "\n", sep = "")
cat("SUMMARY\t", summary_file, "\n", sep = "")
cat("\nGeneSetGroup sizes in 36 sets:\n")
print(table(unique(res[, c("GeneSet", "GeneSetGroup")])$GeneSetGroup))
cat("\nSummary:\n")
print(summary, row.names = FALSE)
cat("\nCount logistic FDR<=0.10 by group/CNV type:\n")
print(res[res$FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype <= 0.10 & res$OR_per_additional_gene_set_gene > 1,
          c("Level", "Comparison", "GeneSet", "GeneSetGroup", "RareGenesInSet", "SCZ_samples_ge1_gene", "CTRL_samples_ge1_gene", "OR_per_additional_gene_set_gene", "Logistic_P_dev_count", "FDR_BH_Logistic_dev_count", "FDR_BH_Logistic_dev_count_by_GeneSetGroup_CNVtype", "Warning_count")], row.names = FALSE)

