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
base <- get_cli_arg("burden-root", file.path(default_project_root, "cnv_analysis/pathway_burden"))
summary_path <- file.path(base, "Supplementary_Table_3_gene_sets_summary.tsv")
rdata_path <- file.path(base, "gene_set_reconstruction/cnvGSAdata_src/cnvGSAdata/data/gs_data_example.RData")
gene_map_path <- file.path(base, "gene_set_reconstruction/cnvGSAdata_src/cnvGSAdata/inst/extdata/gene_ID_demo.txt")

out_prefix <- file.path(base, "original_36_gene_sets_from_gs_data_example_26-7-14")
out_entrez_tsv <- paste0(out_prefix, "_membership_entrez.tsv")
out_symbol_tsv <- paste0(out_prefix, "_membership_symbols.tsv")
out_status_tsv <- paste0(out_prefix, "_status.tsv")
out_dup_tsv <- paste0(out_prefix, "_symbol_duplicate_mapping.tsv")
out_symbol_rows_json <- paste0(out_prefix, "_symbols_original_rows.json")
out_symbol_unique_json <- paste0(out_prefix, "_symbols_for_local_burden_unique.json")
out_entrez_json <- paste0(out_prefix, "_entrez_original.json")
# Backward-compatible aliases from the first extraction. _symbols.json is kept as the local-burden unique-symbol cache.
out_symbol_alias_json <- paste0(out_prefix, "_symbols.json")
out_entrez_alias_json <- paste0(out_prefix, "_entrez.json")

summary_df <- read.delim(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
load(rdata_path)  # gs_all.ls, gsid2name.chv

gene_map <- read.delim(gene_map_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
gene_map$geneID <- as.character(gene_map$geneID)
id_to_symbol <- setNames(gene_map$Symbol, gene_map$geneID)

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  x <- gsub("\r", "\\r", x)
  x <- gsub("\n", "\\n", x)
  x <- gsub("\t", "\\t", x)
  x
}
json_string <- function(x) paste0('"', json_escape(x), '"')
json_array <- function(x, indent = "    ") {
  if (length(x) == 0) return("[]")
  paste0("[\n", paste0(indent, json_string(x), collapse = ",\n"), "\n  ]")
}
write_gene_set_json <- function(path, gene_sets, metadata) {
  lines <- c("{")
  lines <- c(lines, "  \"metadata\": {")
  meta_names <- names(metadata)
  for (i in seq_along(meta_names)) {
    nm <- meta_names[i]
    val <- metadata[[nm]]
    comma <- if (i < length(meta_names)) "," else ""
    if (is.numeric(val) || is.integer(val)) lines <- c(lines, paste0("    ", json_string(nm), ": ", val, comma))
    else lines <- c(lines, paste0("    ", json_string(nm), ": ", json_string(as.character(val)), comma))
  }
  lines <- c(lines, "  },")
  lines <- c(lines, "  \"gene_sets\": {")
  set_names <- names(gene_sets)
  for (i in seq_along(set_names)) {
    set_name <- set_names[i]
    comma <- if (i < length(set_names)) "," else ""
    arr_lines <- strsplit(json_array(gene_sets[[set_name]], indent = "      "), "\n", fixed = TRUE)[[1]]
    lines <- c(lines, paste0("    ", json_string(set_name), ": ", arr_lines[1]))
    if (length(arr_lines) > 1) lines <- c(lines, paste0("    ", arr_lines[-1]))
    lines[length(lines)] <- paste0(lines[length(lines)], comma)
  }
  lines <- c(lines, "  }", "}")
  writeLines(lines, con = path, useBytes = TRUE)
}

entrez_rows <- list(); symbol_rows <- list(); status_rows <- list(); duplicate_rows <- list()
gene_sets_entrez <- list(); gene_sets_symbols_rows <- list(); gene_sets_symbols_unique <- list()

for (i in seq_len(nrow(summary_df))) {
  set_id <- summary_df[["GeneSet ID (Suppl DataSets)"]][i]
  paper_count <- as.integer(summary_df[["#Genes in Set"]][i])
  if (!set_id %in% names(gs_all.ls)) stop(paste("Set missing from gs_all.ls:", set_id))
  entrez <- as.character(gs_all.ls[[set_id]])
  entrez <- entrez[!is.na(entrez) & nzchar(entrez)]
  mapped <- ifelse(entrez %in% names(id_to_symbol), id_to_symbol[entrez], paste0("ENTREZ_", entrez))
  mapped_upper <- toupper(mapped)
  mapping_note <- ifelse(entrez %in% names(id_to_symbol), "", "No symbol in gene_ID_demo; kept as ENTREZ_<ID>")
  recovered_count <- length(entrez)
  unique_symbol_count <- length(unique(mapped_upper))
  status <- if (recovered_count == paper_count) "Exact_original_entrez_count" else "Count_mismatch"
  note <- if (unique_symbol_count < recovered_count) {
    paste0("Original Entrez rows are exact; symbol mapping collapses ", recovered_count - unique_symbol_count, " duplicated symbol(s).")
  } else {
    "Original Entrez rows are exact; symbol mapping is one-to-one at unique-symbol level."
  }
  if (set_id == "Kirov_ARC") note <- paste0(note, " Kirov_ARC is recovered directly from gs_data_example.RData, not reconstructed from GMT.")

  gene_sets_entrez[[set_id]] <- entrez
  gene_sets_symbols_rows[[set_id]] <- mapped_upper
  gene_sets_symbols_unique[[set_id]] <- sort(unique(mapped_upper))

  entrez_rows[[length(entrez_rows) + 1]] <- data.frame(
    GeneSet = set_id, GeneSetGroup = summary_df[["GeneSet Group"]][i], FigureLabel = summary_df[["Figure Label"]][i],
    PaperCount = paper_count, EntrezID = entrez, stringsAsFactors = FALSE
  )
  symbol_rows[[length(symbol_rows) + 1]] <- data.frame(
    GeneSet = set_id, GeneSetGroup = summary_df[["GeneSet Group"]][i], FigureLabel = summary_df[["Figure Label"]][i],
    PaperCount = paper_count, EntrezID = entrez, MappedSymbol = mapped_upper, MappingNote = mapping_note, stringsAsFactors = FALSE
  )
  dup_symbols <- names(which(table(mapped_upper) > 1))
  if (length(dup_symbols) > 0) {
    for (sym in dup_symbols) {
      idx <- which(mapped_upper == sym)
      duplicate_rows[[length(duplicate_rows) + 1]] <- data.frame(
        GeneSet = set_id, DuplicatedMappedSymbol = sym, EntrezIDs = paste(entrez[idx], collapse = ";"),
        OriginalRowCount = length(idx), UniqueSymbolCount = 1, stringsAsFactors = FALSE
      )
    }
  }
  status_rows[[length(status_rows) + 1]] <- data.frame(
    GeneSet = set_id, FigureLabel = summary_df[["Figure Label"]][i], FullName = summary_df[["GeneSet FullName"]][i],
    GeneSetGroup = summary_df[["GeneSet Group"]][i], PaperCount = paper_count, OriginalEntrezCount = recovered_count,
    SymbolRowsCount = length(mapped_upper), UniqueMappedSymbolCount = unique_symbol_count, Status = status,
    Source = "cnvGSAdata::gs_data_example.RData / gs_all.ls", Note = note, stringsAsFactors = FALSE
  )
}

entrez_df <- do.call(rbind, entrez_rows)
symbol_df <- do.call(rbind, symbol_rows)
status_df <- do.call(rbind, status_rows)
dup_df <- if (length(duplicate_rows) > 0) do.call(rbind, duplicate_rows) else data.frame(GeneSet=character(), DuplicatedMappedSymbol=character(), EntrezIDs=character(), OriginalRowCount=integer(), UniqueSymbolCount=integer())

write.table(entrez_df, out_entrez_tsv, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
write.table(symbol_df, out_symbol_tsv, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
write.table(status_df, out_status_tsv, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
write.table(dup_df, out_dup_tsv, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

base_metadata <- list(
  name = "Nature Genetics 2017 ng.3725 original 36 gene sets",
  created_date = "2026-07-14",
  source_rdata = rdata_path,
  source_object = "gs_all.ls",
  source_summary = summary_path,
  total_gene_sets = length(gene_sets_entrez),
  total_gene_set_groups = length(unique(status_df$GeneSetGroup)),
  source_note = "Gene-set members are extracted directly from cnvGSAdata gs_data_example.RData / gs_all.ls, including Kirov_ARC and the 3038-row BrainSpan VHM sets.",
  burden_note = "This file changes gene-set definitions only; downstream burden model remains the local haplotype-level Liftoff extra-copy CNV burden test unless the analysis script is modified."
)
metadata_entrez <- base_metadata
metadata_entrez$gene_identifier <- "EntrezID_original_rows"
metadata_entrez$count_note <- "The gene_sets arrays preserve original Entrez rows and match all 36 paper counts exactly."
write_gene_set_json(out_entrez_json, gene_sets_entrez, metadata_entrez)
write_gene_set_json(out_entrez_alias_json, gene_sets_entrez, metadata_entrez)

metadata_symbol_rows <- base_metadata
metadata_symbol_rows$gene_identifier <- "MappedSymbol_original_rows"
metadata_symbol_rows$count_note <- "The gene_sets arrays preserve original row counts after Entrez-to-symbol mapping. Some mapped symbols are duplicated because distinct Entrez IDs map to the same current symbol."
write_gene_set_json(out_symbol_rows_json, gene_sets_symbols_rows, metadata_symbol_rows)

metadata_symbol_unique <- base_metadata
metadata_symbol_unique$gene_identifier <- "MappedSymbol_unique_for_local_burden"
metadata_symbol_unique$count_note <- "The gene_sets arrays are unique uppercase symbols for local CNV matrix matching. BspanVHM_PreNat and BspanVHM_EqlNat contain 3037 unique symbols from 3038 original Entrez rows due duplicate symbol mapping."
write_gene_set_json(out_symbol_unique_json, gene_sets_symbols_unique, metadata_symbol_unique)
write_gene_set_json(out_symbol_alias_json, gene_sets_symbols_unique, metadata_symbol_unique)

cat("WROTE\n")
cat(out_entrez_tsv, "\n", out_symbol_tsv, "\n", out_status_tsv, "\n", out_dup_tsv, "\n", out_entrez_json, "\n", out_symbol_rows_json, "\n", out_symbol_unique_json, "\n", out_symbol_alias_json, "\n", sep = "")
cat("STATUS_COUNTS\n")
print(table(status_df$Status))
cat("DUPLICATE_SYMBOL_MAPPINGS\n")
print(dup_df, row.names = FALSE)
if (any(status_df$OriginalEntrezCount != status_df$PaperCount)) quit(status = 2)

