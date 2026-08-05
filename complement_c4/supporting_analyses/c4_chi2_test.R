#!/usr/bin/env Rscript

# Historical C4 comparison script.
#
# The original file was named as a chi-square test but the active analysis
# used a one-sided Fisher's exact test.  This compatibility script supports
# both tests explicitly; use the main C4 Fisher script for the final analysis.

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
})

parse_named_args <- function(arguments) {
  result <- list()
  for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) {
      stop("Arguments must use the form --name=value: ", argument)
    }
    parts <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    result[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  result
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
required_args <- c("input-xlsx", "sheet", "output-xlsx")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", ")
  )
}

input_xlsx <- args[["input-xlsx"]]
sheet_name <- args[["sheet"]]
output_xlsx <- args[["output-xlsx"]]
test_name <- if ("test" %in% names(args)) args[["test"]] else "fisher"

if (!test_name %in% c("fisher", "chisq")) {
  stop("--test must be either fisher or chisq.")
}
if (!file.exists(input_xlsx)) {
  stop("Input workbook not found: ", input_xlsx)
}

data <- read_excel(input_xlsx, sheet = sheet_name)
required_columns <- c("Type", "SCZ Count", "HC Count")
missing_columns <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0L) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

data[["SCZ Count"]] <- as.numeric(data[["SCZ Count"]])
data[["HC Count"]] <- as.numeric(data[["HC Count"]])
if (anyNA(data[["SCZ Count"]]) || anyNA(data[["HC Count"]])) {
  stop("SCZ Count and HC Count must contain numeric, non-missing values.")
}

total_scz <- if ("scz-total" %in% names(args)) {
  as.numeric(args[["scz-total"]])
} else {
  sum(data[["SCZ Count"]])
}
total_control <- if ("control-total" %in% names(args)) {
  as.numeric(args[["control-total"]])
} else {
  sum(data[["HC Count"]])
}

p_values <- numeric(nrow(data))
for (i in seq_len(nrow(data))) {
  contingency_table <- matrix(
    c(
      data[["SCZ Count"]][i],
      data[["HC Count"]][i],
      total_scz - data[["SCZ Count"]][i],
      total_control - data[["HC Count"]][i]
    ),
    nrow = 2,
    byrow = TRUE
  )
  if (any(contingency_table < 0)) {
    stop("Negative contingency-table cell at row ", i, ".")
  }

  p_values[i] <- if (test_name == "fisher") {
    fisher.test(contingency_table, alternative = "greater")$p.value
  } else {
    suppressWarnings(chisq.test(contingency_table, correct = FALSE)$p.value)
  }
}

data[["p_value"]] <- p_values
data[["p_value_FDR"]] <- p.adjust(p_values, method = "fdr")
data[["p_value_Bonferroni"]] <- p.adjust(p_values, method = "bonferroni")

dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)
write_xlsx(data, output_xlsx)

cat("Test:", test_name, "\n")
cat("SCZ total:", total_scz, "\n")
cat("Control total:", total_control, "\n")
cat("Saved:", output_xlsx, "\n")
