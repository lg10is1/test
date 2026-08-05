############################################################
## Convert rare_variants.txt to clean sample-level burden table
## Robust version: tolerate absent TR_contraction / TR_expansion / RareSV classes
############################################################

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

############################################################
## 0. Path
############################################################

cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
base_dir <- Sys.getenv("RESULT_ROOT", file.path(dirname(script_dir), "results", "prs"))
data_dir <- file.path(base_dir, "data")

rare_file <- file.path(data_dir, "rare_variants.txt")
out_file <- file.path(data_dir, "rare_variants_clean_burden.tsv")
out_detail_file <- file.path(data_dir, "rare_variants_clean_burden_with_TR_detail.tsv")
out_unmatched_file <- file.path(data_dir, "rare_variants_samples_not_in_case_list.tsv")

############################################################
## 1. Extract all case samples from SCZ PRS file
## Samples beginning with C are controls; all others are cases.
############################################################

scz_score_file <- Sys.getenv("SCZ_PRS_FILE", file.path(base_dir, "SCZ/prscsx_eur.sscore"))

scz_score <- fread(scz_score_file)

case_samples <- scz_score %>%
  transmute(sample = IID) %>%
  filter(!grepl("^C", sample)) %>%
  distinct(sample) %>%
  arrange(sample)

cat("[INFO] Number of case samples:", nrow(case_samples), "\n")

############################################################
## 2. Read rare_variants.txt
## format: variant_id, variant_type, sample_list
############################################################

rv <- fread(
  rare_file,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  col.names = c("variant_id", "variant_type", "sample_list")
)

rv <- rv %>%
  mutate(
    variant_id = str_trim(variant_id),
    variant_type = str_trim(variant_type),
    sample_list = str_trim(sample_list)
  ) %>%
  filter(
    !is.na(variant_id),
    !is.na(variant_type),
    !is.na(sample_list),
    sample_list != ""
  )

cat("[INFO] Number of variant rows:", nrow(rv), "\n")
print(table(rv$variant_type))

############################################################
## 3. Expand sample list
## Repeated samples in column 3 are retained and contribute multiple counts.
############################################################

rv_long <- rv %>%
  separate_rows(sample_list, sep = ",") %>%
  mutate(
    sample = str_trim(sample_list)
  ) %>%
  filter(
    !is.na(sample),
    sample != ""
  ) %>%
  select(variant_id, variant_type, sample)

cat("[INFO] Number of expanded sample-variant records:", nrow(rv_long), "\n")

############################################################
## 4. Helper: add missing count columns after pivot_wider
############################################################

add_missing_count_cols <- function(dat, cols) {
  for (cc in cols) {
    if (!cc %in% colnames(dat)) {
      dat[[cc]] <- 0L
    }
  }
  dat
}

############################################################
## 5. Map variant type
############################################################

rv_long2 <- rv_long %>%
  mutate(
    burden_class = case_when(
      variant_type == "rare SV" ~ "RareSV",
      variant_type %in% c("TR expansion", "TR contraction") ~ "RareTR",
      variant_type == "rare SNV" ~ "RareSNV",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(burden_class))

############################################################
## 6. Count burden per sample
############################################################

count_old <- rv_long2 %>%
  count(sample, burden_class, name = "n") %>%
  pivot_wider(
    names_from = burden_class,
    values_from = n,
    values_fill = 0
  ) %>%
  add_missing_count_cols(c("RareSV", "RareTR", "RareSNV"))

############################################################
## 7. Clean burden table
############################################################

clean_burden <- case_samples %>%
  left_join(count_old, by = "sample") %>%
  mutate(
    RareSV = ifelse(is.na(RareSV), 0L, as.integer(RareSV)),
    RareTR = ifelse(is.na(RareTR), 0L, as.integer(RareTR)),
    RareSNV = ifelse(is.na(RareSNV), 0L, as.integer(RareSNV)),
    sum_rare = RareSV + RareTR + RareSNV,
    sum_all = sum_rare
  ) %>%
  select(
    sample,
    RareSV,
    RareTR,
    RareSNV,
    sum_all,
    sum_rare
  ) %>%
  arrange(
    desc(sum_all),
    desc(RareSV),
    desc(RareTR),
    sample
  )

fwrite(
  clean_burden,
  file = out_file,
  sep = "\t",
  quote = FALSE
)

cat("[INFO] Clean burden table saved to:\n", out_file, "\n")
cat("[INFO] Output sample number:", nrow(clean_burden), "\n")
print(head(clean_burden, 30))

############################################################
## 8. Optional detail file: split TR expansion / contraction
## Do not assume TR contraction exists; fill the missing class with zero.
############################################################

count_detail <- rv_long %>%
  mutate(
    detail_class = case_when(
      variant_type == "rare SV" ~ "RareSV",
      variant_type == "rare SNV" ~ "RareSNV",
      variant_type == "TR expansion" ~ "TR_expansion",
      variant_type == "TR contraction" ~ "TR_contraction",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(detail_class)) %>%
  count(sample, detail_class, name = "n") %>%
  pivot_wider(
    names_from = detail_class,
    values_from = n,
    values_fill = 0
  ) %>%
  add_missing_count_cols(c("RareSV", "RareSNV", "TR_expansion", "TR_contraction"))

clean_burden_detail <- case_samples %>%
  left_join(count_detail, by = "sample") %>%
  mutate(
    RareSV = ifelse(is.na(RareSV), 0L, as.integer(RareSV)),
    RareSNV = ifelse(is.na(RareSNV), 0L, as.integer(RareSNV)),
    TR_expansion = ifelse(is.na(TR_expansion), 0L, as.integer(TR_expansion)),
    TR_contraction = ifelse(is.na(TR_contraction), 0L, as.integer(TR_contraction)),
    RareTR = TR_expansion + TR_contraction,
    sum_rare = RareSV + RareTR + RareSNV,
    sum_all = sum_rare
  ) %>%
  select(
    sample,
    RareSV,
    TR_expansion,
    TR_contraction,
    RareTR,
    RareSNV,
    sum_all,
    sum_rare
  ) %>%
  arrange(
    desc(sum_all),
    desc(RareSV),
    desc(RareTR),
    sample
  )

fwrite(
  clean_burden_detail,
  file = out_detail_file,
  sep = "\t",
  quote = FALSE
)

cat("[INFO] Detail table saved to:\n", out_detail_file, "\n")

############################################################
## 9. Check samples in rare_variants.txt but not in case list
############################################################

samples_not_in_case <- setdiff(unique(rv_long$sample), case_samples$sample)

cat("[INFO] Samples in rare_variants.txt but not in case list:", length(samples_not_in_case), "\n")

if (length(samples_not_in_case) > 0) {
  fwrite(
    data.frame(sample = samples_not_in_case),
    file = out_unmatched_file,
    sep = "\t",
    quote = FALSE
  )
  cat("[INFO] Unmatched samples saved to:\n", out_unmatched_file, "\n")
}
