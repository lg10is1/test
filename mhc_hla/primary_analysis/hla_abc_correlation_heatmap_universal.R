suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
})

usage_text <- paste(
  "Usage:",
  "Rscript hla_abc_correlation_heatmap_universal.R",
  "--scz-xlsx=PATH",
  "--row-gene=C",
  "--col-gene=B",
  "--row-result=PATH",
  "--col-result=PATH",
  "--min-scz-count=5",
  "--output-dir=PATH",
  "--output-prefix=NAME",
  sep = " "
)

parse_args <- function(args) {
  defaults <- list(
    scz_xlsx = NA_character_,
    row_gene = "C",
    col_gene = "B",
    row_result = NA_character_,
    col_result = NA_character_,
    row_alleles = "",
    col_alleles = "",
    significance_col = "Significant_FDR_0.05",
    fdr_threshold = 0.05,
    allow_nominal_fallback = FALSE,
    nominal_p_threshold = 0.05,
    min_scz_count = 5,
    drop_null = TRUE,
    output_dir = ".",
    output_prefix = "",
    plot_title = "",
    width = NA_real_,
    height = NA_real_,
    label_digits = 2
  )

  coerce_value <- function(key, value) {
    if (key %in% c("min_scz_count", "label_digits")) {
      return(as.integer(value))
    }
    if (key %in% c("fdr_threshold", "nominal_p_threshold", "width", "height")) {
      return(as.numeric(value))
    }
    if (key %in% c("drop_null", "allow_nominal_fallback")) {
      return(tolower(value) %in% c("true", "1", "yes", "y"))
    }
    value
  }

  opts <- defaults
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    arg_body <- substring(arg, 3)
    pieces <- strsplit(arg_body, "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", pieces[1])
    value <- if (length(pieces) > 1) paste(pieces[-1], collapse = "=") else "TRUE"
    if (!key %in% names(opts)) {
      stop(paste("Unknown argument:", key, "\n", usage_text), call. = FALSE)
    }
    opts[[key]] <- coerce_value(key, value)
  }
  opts
}

normalize_values <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == ""] <- NA_character_
  x[toupper(x) == "NULL"] <- NA_character_
  x[tolower(x) == "unknown"] <- NA_character_
  x
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(character())
  }
  parts <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  parts[parts != ""]
}

load_gene_calls <- function(xlsx_path, gene) {
  sheets <- excel_sheets(xlsx_path)
  public_sheet <- paste0(gene, "_qc_passed")
  legacy_sheet <- paste0(gene, "_validated")
  chosen_sheet <- if (public_sheet %in% sheets) {
    public_sheet
  } else if (legacy_sheet %in% sheets) {
    legacy_sheet
  } else {
    sheets[1]
  }
  df <- suppressMessages(read_excel(xlsx_path, sheet = chosen_sheet, .name_repair = "unique"))
  if (!"Sample_name" %in% names(df)) {
    stop(paste("Sheet", chosen_sheet, "missing column Sample_name"), call. = FALSE)
  }
  if (!gene %in% names(df)) {
    stop(paste("Sheet", chosen_sheet, "missing column", gene), call. = FALSE)
  }
  out <- df[, c("Sample_name", gene)]
  names(out) <- c("Sample_name", "Genotype")
  out$Sample_name <- as.character(out$Sample_name)
  out$Genotype <- normalize_values(out$Genotype)
  out$Gene <- gene
  out$Source_Sheet <- chosen_sheet
  distinct(out)
}

prepare_pair_dataset <- function(xlsx_path, row_gene, col_gene, drop_null = TRUE) {
  row_df <- load_gene_calls(xlsx_path, row_gene) %>%
    transmute(
      Sample_name = Sample_name,
      Row_Gene = Gene,
      Row_Genotype = Genotype,
      Row_Source_Sheet = Source_Sheet
    )

  col_df <- load_gene_calls(xlsx_path, col_gene) %>%
    transmute(
      Sample_name = Sample_name,
      Col_Gene = Gene,
      Col_Genotype = Genotype,
      Col_Source_Sheet = Source_Sheet
    )

  merged <- inner_join(col_df, row_df, by = "Sample_name")
  if (drop_null) {
    merged <- merged %>%
      filter(!is.na(Col_Genotype), !is.na(Row_Genotype))
  }
  distinct(merged)
}

read_result_alleles <- function(
  result_path,
  manual_alleles,
  min_scz_count,
  significance_col,
  fdr_threshold,
  allow_nominal_fallback,
  nominal_p_threshold
) {
  manual_values <- normalize_values(split_csv(manual_alleles))
  manual_values <- unique(manual_values[!is.na(manual_values)])

  ordered_from_result <- character()
  result_table <- data.frame()
  selection_source <- "manual_only"

  if (!is.na(result_path) && nzchar(result_path)) {
    raw_res <- suppressMessages(read_excel(result_path, .name_repair = "unique"))
    res <- raw_res
    genotype_col <- if ("Genotype" %in% names(res)) {
      "Genotype"
    } else if ("Subtype" %in% names(res)) {
      "Subtype"
    } else {
      stop(paste("No Genotype/Subtype column found in", result_path), call. = FALSE)
    }

    if (significance_col %in% names(res)) {
      keep <- res[[significance_col]]
      keep[is.na(keep)] <- FALSE
      res <- res[keep, , drop = FALSE]
    } else if ("FDR_BH" %in% names(res)) {
      res <- res %>% filter(!is.na(FDR_BH), FDR_BH <= fdr_threshold)
    }

    if ("SCZ_Count" %in% names(res)) {
      res <- res %>% filter(!is.na(SCZ_Count), SCZ_Count >= min_scz_count)
    }

    if (nrow(res) == 0 && allow_nominal_fallback && "P-Value" %in% names(raw_res)) {
      nominal_res <- raw_res
      if ("SCZ_Count" %in% names(nominal_res)) {
        nominal_res <- nominal_res %>% filter(!is.na(SCZ_Count), SCZ_Count >= min_scz_count)
      }
      nominal_res <- nominal_res %>% filter(!is.na(`P-Value`), `P-Value` <= nominal_p_threshold)
      if (nrow(nominal_res) > 0) {
        res <- nominal_res
        selection_source <- paste0("nominal_p<=", nominal_p_threshold)
      }
    } else if (nrow(res) > 0) {
      selection_source <- "fdr_or_significance_column"
    }

    if (nrow(res) > 0) {
      if (all(c("SCZ_Count", "FDR_BH") %in% names(res))) {
        res <- res %>% arrange(desc(SCZ_Count), FDR_BH, .data[[genotype_col]])
      } else if ("SCZ_Count" %in% names(res)) {
        res <- res %>% arrange(desc(SCZ_Count), .data[[genotype_col]])
      } else {
        res <- res %>% arrange(.data[[genotype_col]])
      }
      res[[genotype_col]] <- normalize_values(res[[genotype_col]])
      res <- res %>% filter(!is.na(.data[[genotype_col]]))
      ordered_from_result <- unique(res[[genotype_col]])
      result_table <- as.data.frame(res)
    } else {
      result_table <- as.data.frame(res)
    }
  }

  combined <- ordered_from_result
  manual_extras <- setdiff(manual_values, combined)
  if (length(manual_extras) > 0) {
    combined <- c(combined, sort(manual_extras))
  }

  combined <- unique(combined[!is.na(combined)])
  list(
    alleles = combined,
    result_table = result_table,
    manual_alleles = manual_values,
    selection_source = selection_source
  )
}

filter_observed_alleles <- function(selected, observed_vector) {
  observed_vector <- normalize_values(observed_vector)
  observed_vector <- unique(observed_vector[!is.na(observed_vector)])
  selected[selected %in% observed_vector]
}

order_by_observed_count <- function(selected, observed_vector) {
  if (length(selected) == 0) {
    return(selected)
  }
  counts <- vapply(
    selected,
    function(x) as.numeric(sum(observed_vector == x, na.rm = TRUE)),
    numeric(1)
  )
  ordered <- selected[order(-counts, selected)]
  unique(ordered)
}

compute_pairwise_stats <- function(pair_df, col_alleles, row_alleles, col_gene, row_gene) {
  records <- vector("list", length(col_alleles) * length(row_alleles))
  idx <- 1
  total_n <- nrow(pair_df)

  for (col_allele in col_alleles) {
    col_vec <- pair_df$Col_Genotype == col_allele
    col_count <- sum(col_vec, na.rm = TRUE)

    for (row_allele in row_alleles) {
      row_vec <- pair_df$Row_Genotype == row_allele
      row_count <- sum(row_vec, na.rm = TRUE)

      both <- sum(col_vec & row_vec, na.rm = TRUE)
      col_only <- sum(col_vec & !row_vec, na.rm = TRUE)
      row_only <- sum(!col_vec & row_vec, na.rm = TRUE)
      neither <- sum(!col_vec & !row_vec, na.rm = TRUE)

      pearson_correlation <- suppressWarnings(cor(as.numeric(col_vec), as.numeric(row_vec)))
      if (is.nan(pearson_correlation)) {
        pearson_correlation <- NA_real_
      }

      contingency <- matrix(c(both, col_only, row_only, neither), nrow = 2, byrow = TRUE)
      fisher_res <- fisher.test(contingency, alternative = "greater")
      odds_ratio <- if (!is.null(fisher_res$estimate)) unname(fisher_res$estimate[[1]]) else NA_real_
      expected_both <- if (total_n > 0) (col_count * row_count) / total_n else NA_real_
      enrichment <- if (!is.na(expected_both) && expected_both > 0) both / expected_both else NA_real_

      records[[idx]] <- data.frame(
        Col_Gene = col_gene,
        Col_Allele = col_allele,
        Row_Gene = row_gene,
        Row_Allele = row_allele,
        Pair_Total = total_n,
        Col_Count = col_count,
        Row_Count = row_count,
        Both_Count = both,
        Col_Only = col_only,
        Row_Only = row_only,
        Neither = neither,
        Col_Frequency = col_count / total_n,
        Row_Frequency = row_count / total_n,
        Both_Frequency = both / total_n,
        Both_in_Col = if (col_count > 0) both / col_count else NA_real_,
        Both_in_Row = if (row_count > 0) both / row_count else NA_real_,
        Expected_Both = expected_both,
        Enrichment_vs_Expected = enrichment,
        Pearson_Correlation = pearson_correlation,
        Odds_Ratio = odds_ratio,
        P_Value = fisher_res$p.value,
        Fisher_Alternative = "greater",
        Contingency_Table = paste0("[[", both, ", ", col_only, "], [", row_only, ", ", neither, "]]"),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  stats <- bind_rows(records)
  stats$FDR_BH <- p.adjust(stats$P_Value, method = "BH")
  stats$Significant_FDR_0.05 <- stats$FDR_BH <= 0.05
  stats
}

build_selection_sheet <- function(result_info, final_alleles, observed_vector, gene_label) {
  observed_counts <- sapply(final_alleles, function(x) sum(observed_vector == x, na.rm = TRUE))
  result_genotypes <- character()
  if (nrow(result_info$result_table) > 0) {
    genotype_col <- if ("Genotype" %in% names(result_info$result_table)) {
      "Genotype"
    } else if ("Subtype" %in% names(result_info$result_table)) {
      "Subtype"
    } else {
      NULL
    }
    if (!is.null(genotype_col)) {
      result_genotypes <- normalize_values(result_info$result_table[[genotype_col]])
    }
  }
  sheet <- data.frame(
    Gene = gene_label,
    Allele = final_alleles,
    Observed_In_Pair_Data = observed_counts,
    Selected_From_Result = final_alleles %in% result_genotypes,
    Selected_Manually = final_alleles %in% result_info$manual_alleles,
    Selection_Source = result_info$selection_source,
    stringsAsFactors = FALSE
  )

  if (nrow(result_info$result_table) > 0) {
    result_table <- result_info$result_table
    genotype_col <- if ("Genotype" %in% names(result_table)) "Genotype" else if ("Subtype" %in% names(result_table)) "Subtype" else NULL
    if (!is.null(genotype_col)) {
      result_table[[genotype_col]] <- normalize_values(result_table[[genotype_col]])
      result_table <- result_table[result_table[[genotype_col]] %in% final_alleles, , drop = FALSE]
      names(result_table)[names(result_table) == genotype_col] <- "Allele"
      sheet <- left_join(sheet, result_table, by = "Allele")
    }
  }

  sheet
}

write_sheets <- function(sheets, path) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    for (sheet_name in names(sheets)) {
      clean_name <- substr(sheet_name, 1, 31)
      openxlsx::addWorksheet(wb, clean_name)
      openxlsx::writeData(wb, clean_name, sheets[[sheet_name]])
    }
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    return(path)
  }

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path)
    return(path)
  }

  csv_dir <- sub("\\.xlsx$", "_csv", path)
  if (identical(csv_dir, path)) {
    csv_dir <- paste0(path, "_csv")
  }
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  for (sheet_name in names(sheets)) {
    csv_path <- file.path(csv_dir, paste0(sheet_name, ".csv"))
    write.csv(sheets[[sheet_name]], csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  csv_dir
}

plot_heatmap <- function(stats_df, row_alleles, col_alleles, output_pdf, output_png, plot_title, width, height, label_digits) {
  plot_df <- stats_df %>%
    mutate(
      Col_Allele = factor(Col_Allele, levels = col_alleles),
      Row_Allele = factor(Row_Allele, levels = row_alleles)
    )

  value_labels <- ifelse(
    is.na(plot_df$Pearson_Correlation),
    "NA",
    sprintf(paste0("%.", label_digits, "f"), plot_df$Pearson_Correlation)
  )
  plot_df$Label <- ifelse(plot_df$Significant_FDR_0.05, paste0(value_labels, "*"), value_labels)

  label_size <- 8 / 2.845276
  axis_text_size_x <- 9
  axis_text_size_y <- 9

  p <- ggplot(plot_df, aes(x = Col_Allele, y = Row_Allele, fill = Pearson_Correlation)) +
    geom_tile(color = "grey88", linewidth = 0.5) +
    geom_text(aes(label = Label), size = label_size, family = "Arial") +
    scale_fill_gradient2(
      low = "#1f78b4",
      mid = "white",
      high = "#da7271",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90"
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = NULL,
      fill = "Pearson Correlation"
    ) +
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size = axis_text_size_x),
      axis.text.y = element_text(size = axis_text_size_y),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.position = "right"
    )

  ggsave(output_pdf, p, width = width, height = height, device = grDevices::cairo_pdf, family = "Arial")
  ggsave(output_png, p, width = width, height = height, dpi = 300)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.na(args$scz_xlsx) || !nzchar(args$scz_xlsx)) {
  stop(paste("Missing --scz-xlsx\n", usage_text), call. = FALSE)
}

if (!file.exists(args$scz_xlsx)) {
  stop(paste("Input file not found:", args$scz_xlsx), call. = FALSE)
}

if (is.na(args$output_prefix) || !nzchar(args$output_prefix)) {
  args$output_prefix <- paste0("HLA_", args$col_gene, "_", args$row_gene, "_correlation")
}

if (is.na(args$plot_title) || !nzchar(args$plot_title)) {
  args$plot_title <- paste0("HLA-", args$col_gene, " × HLA-", args$row_gene, " correlation heatmap")
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

pair_df <- prepare_pair_dataset(
  xlsx_path = args$scz_xlsx,
  row_gene = args$row_gene,
  col_gene = args$col_gene,
  drop_null = args$drop_null
)

if (nrow(pair_df) == 0) {
  stop("No paired haplotypes remained after merging and filtering.", call. = FALSE)
}

row_info <- read_result_alleles(
  result_path = args$row_result,
  manual_alleles = args$row_alleles,
  min_scz_count = args$min_scz_count,
  significance_col = args$significance_col,
  fdr_threshold = args$fdr_threshold,
  allow_nominal_fallback = args$allow_nominal_fallback,
  nominal_p_threshold = args$nominal_p_threshold
)

col_info <- read_result_alleles(
  result_path = args$col_result,
  manual_alleles = args$col_alleles,
  min_scz_count = args$min_scz_count,
  significance_col = args$significance_col,
  fdr_threshold = args$fdr_threshold,
  allow_nominal_fallback = args$allow_nominal_fallback,
  nominal_p_threshold = args$nominal_p_threshold
)

row_alleles <- filter_observed_alleles(row_info$alleles, pair_df$Row_Genotype)
col_alleles <- filter_observed_alleles(col_info$alleles, pair_df$Col_Genotype)

row_alleles <- order_by_observed_count(row_alleles, pair_df$Row_Genotype)
col_alleles <- order_by_observed_count(col_alleles, pair_df$Col_Genotype)

if (length(row_alleles) == 0) {
  stop("No row-gene alleles were selected. Check --row-result / --row-alleles.", call. = FALSE)
}

if (length(col_alleles) == 0) {
  stop("No col-gene alleles were selected. Check --col-result / --col-alleles.", call. = FALSE)
}

stats_df <- compute_pairwise_stats(
  pair_df = pair_df,
  col_alleles = col_alleles,
  row_alleles = row_alleles,
  col_gene = args$col_gene,
  row_gene = args$row_gene
)

row_sheet <- build_selection_sheet(row_info, row_alleles, pair_df$Row_Genotype, args$row_gene)
col_sheet <- build_selection_sheet(col_info, col_alleles, pair_df$Col_Genotype, args$col_gene)

settings_sheet <- data.frame(
  Parameter = c(
    "scz_xlsx", "row_gene", "col_gene", "row_result", "col_result",
    "min_scz_count", "drop_null", "allow_nominal_fallback", "nominal_p_threshold",
    "pair_total", "row_sheet", "col_sheet", "row_selection_source", "col_selection_source"
  ),
  Value = c(
    args$scz_xlsx, args$row_gene, args$col_gene, args$row_result, args$col_result,
    as.character(args$min_scz_count), as.character(args$drop_null),
    as.character(args$allow_nominal_fallback), as.character(args$nominal_p_threshold),
    as.character(nrow(pair_df)), unique(pair_df$Row_Source_Sheet)[1], unique(pair_df$Col_Source_Sheet)[1],
    row_info$selection_source, col_info$selection_source
  ),
  stringsAsFactors = FALSE
)

base_width <- max(7, 1.6 * length(col_alleles) + 2)
base_height <- max(6, 1.2 * length(row_alleles) + 2.5)
plot_width <- if (is.na(args$width)) base_width else args$width
plot_height <- if (is.na(args$height)) base_height else args$height

pdf_path <- file.path(args$output_dir, paste0(args$output_prefix, "_heatmap.pdf"))
png_path <- file.path(args$output_dir, paste0(args$output_prefix, "_heatmap.png"))
summary_path <- file.path(args$output_dir, paste0(args$output_prefix, "_summary.xlsx"))

plot_heatmap(
  stats_df = stats_df,
  row_alleles = row_alleles,
  col_alleles = col_alleles,
  output_pdf = pdf_path,
  output_png = png_path,
  plot_title = args$plot_title,
  width = plot_width,
  height = plot_height,
  label_digits = args$label_digits
)

actual_summary_path <- write_sheets(
  list(
    settings = settings_sheet,
    selected_cols = col_sheet,
    selected_rows = row_sheet,
    pairwise_stats = stats_df,
    pair_data = pair_df
  ),
  summary_path
)

cat("Heatmap PDF:", pdf_path, "\n")
cat("Heatmap PNG:", png_path, "\n")
cat("Summary output:", actual_summary_path, "\n")
