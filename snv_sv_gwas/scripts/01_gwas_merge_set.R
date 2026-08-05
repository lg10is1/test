library(data.table)
library(qqman)

# ============================================================
# 0.       
# ============================================================

edge_bp <- 1e6
sig_p <- 5e-6
sv_cutoff <- 50

out_dir <- "/path/to/EOSCZ_PROJECT/figure_analysis/01.GWAS_figure.public"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  ## Public example: set00 only. Add set01/set02 here if available.
  set00 = "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/gwas/SCZ.mlm.ngspc.fastGWA"
)

# ============================================================
# 1.          
# ============================================================

chr_lengths <- data.table(
  CHR = as.character(1:22),
  chr_len = c(
    248387328, 242696752, 201105948, 193574945, 182045439,
    172126628, 160567428, 146259331, 150617247, 134758134,
    135127769, 133324548, 113566686, 101161492, 99753195,
    96330374, 84276897, 80542538, 61707364, 66210255,
    45090682, 51324926
  )
)

# ============================================================
# 2.       
# ============================================================

calc_lambda <- function(p) {
  p <- p[!is.na(p) & p > 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  
  chisq_obs <- qchisq(p, df = 1, lower.tail = FALSE)
  lambda <- median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
  
  return(lambda)
}

standardize_alleles <- function(dt) {
  
  dt[, A1 := as.character(A1)]
  dt[, A2 := as.character(A2)]
  
  flip_idx <- which(
    !is.na(dt$A1) &
      !is.na(dt$A2) &
      dt$A1 > dt$A2
  )
  
  if (length(flip_idx) > 0) {
    
    old_A1 <- dt$A1[flip_idx]
    old_A2 <- dt$A2[flip_idx]
    
    set(dt, i = flip_idx, j = "A1", value = old_A2)
    set(dt, i = flip_idx, j = "A2", value = old_A1)
    
    if ("BETA" %in% names(dt)) {
      set(dt, i = flip_idx, j = "BETA", value = -dt$BETA[flip_idx])
    }
    
    if ("T" %in% names(dt)) {
      set(dt, i = flip_idx, j = "T", value = -dt$T[flip_idx])
    }
    
    if ("AF1" %in% names(dt)) {
      set(dt, i = flip_idx, j = "AF1", value = 1 - dt$AF1[flip_idx])
    }
  }
  
  return(dt)
}

process_gwas <- function(file_path, set_name, chr_lengths, edge_bp, out_dir) {
  
  message("Processing: ", set_name)
  
  dt <- fread(file_path)
  
  required_cols <- c("CHR", "POS", "SNP", "A1", "A2", "P")
  missing_cols <- setdiff(required_cols, names(dt))
  
  if (length(missing_cols) > 0) {
    stop(
      "          ? ",
      paste(missing_cols, collapse = ", "),
      "\n   : ",
      file_path
    )
  }
  
  #     CHR    
  dt[, CHR := as.character(CHR)]
  dt[, CHR := sub("^chr", "", CHR, ignore.case = TRUE)]
  n_before_autosome <- nrow(dt)
  dt <- dt[CHR %in% as.character(1:22)]
  message("[AUTOSOME] ", set_name, ": kept ", nrow(dt), " / ", n_before_autosome, " rows on chr1-22")
  
  #           ?  dt[, POS := as.numeric(POS)]
  dt[, P := as.numeric(P)]
  dt[, SNP := as.character(SNP)]
  
  #           ?  dt <- merge(dt, chr_lengths, by = "CHR", all.x = TRUE)
  
  #                      
  dt <- dt[!is.na(chr_len)]
  
  #              ?1Mb
  dt <- dt[POS > edge_bp & POS < (chr_len - edge_bp)]
  
  #        chr_len
  dt[, chr_len := NULL]
  
  #     A1 / A2           ?BETA / T / AF1
  dt <- standardize_alleles(dt)
  
  #     SV    
  dt[, sv_len := pmax(nchar(A2), nchar(A1))]
  
  #     set    
  dt[, set := set_name]
  
  #     ?set        ?       ?  out_file_all <- file.path(
    out_dir,
    paste0(set_name, ".remove_chr_edge_1Mb.standardized.fastGWA.tsv")
  )
  
  fwrite(dt, out_file_all, sep = "\t")
  
  # ==========================================================
  #        ?set        SV  ?SNV/indel
  # ==========================================================
  
  dt_sv <- dt[sv_len >= sv_cutoff]
  dt_snv_indel <- dt[sv_len < sv_cutoff]
  
  out_file_sv <- file.path(
    out_dir,
    paste0(set_name, ".SV_ge50bp.remove_chr_edge_1Mb.standardized.fastGWA.tsv")
  )
  
  out_file_snv_indel <- file.path(
    out_dir,
    paste0(set_name, ".SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.fastGWA.tsv")
  )
  
  fwrite(dt_sv, out_file_sv, sep = "\t")
  fwrite(dt_snv_indel, out_file_snv_indel, sep = "\t")
  
  message("Saved full set: ", out_file_all)
  message("Saved SV set: ", out_file_sv)
  message("Saved SNV/indel set: ", out_file_snv_indel)
  message(set_name, " total variants: ", nrow(dt))
  message(set_name, " SV >= 50bp: ", nrow(dt_sv))
  message(set_name, " SNV/indel < 50bp: ", nrow(dt_snv_indel))
  
  return(dt)
}

prepare_for_merge <- function(dt, suffix) {
  
  dt2 <- copy(dt)
  
  merge_keys <- c("CHR", "POS", "A1", "A2", "sv_len")
  
  #        ?set                 merge      
  if (anyDuplicated(dt2, by = merge_keys)) {
    message("    ?, suffix, "           CHR/POS/A1/A2/sv_len    ?P           ?)
    setorder(dt2, P)
    dt2 <- unique(dt2, by = merge_keys)
  }
  
  cols_to_rename <- setdiff(names(dt2), merge_keys)
  setnames(dt2, cols_to_rename, paste0(cols_to_rename, "_", suffix))
  
  return(dt2)
}

# ============================================================
# 3.           set          set        1Mb     ?# ============================================================

dt_list <- Map(
  f = function(file_path, set_name) {
    process_gwas(
      file_path = file_path,
      set_name = set_name,
      chr_lengths = chr_lengths,
      edge_bp = edge_bp,
      out_dir = out_dir
    )
  },
  file_path = files,
  set_name = names(files)
)

# ============================================================
# 4.              ?set       
# ============================================================

per_set_data <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

fwrite(
  per_set_data,
  file.path(out_dir, "all_sets.long_format.remove_chr_edge_1Mb.standardized.tsv"),
  sep = "\t"
)

# ============================================================
# 5.  ?SV / SNV-indel        ?set pooled    
# ============================================================

sv_sites_by_set <- per_set_data[sv_len >= sv_cutoff]
snv_indel_sites_by_set <- per_set_data[sv_len < sv_cutoff]

fwrite(
  sv_sites_by_set,
  file.path(out_dir, "all_sets.long_format.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv"),
  sep = "\t"
)

fwrite(
  snv_indel_sites_by_set,
  file.path(out_dir, "all_sets.long_format.SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.tsv"),
  sep = "\t"
)

#           ?set        pooled     ?sv_sig_by_set <- sv_sites_by_set[P < sig_p]
snv_indel_sig_by_set <- snv_indel_sites_by_set[P < sig_p]

fwrite(
  sv_sig_by_set,
  file.path(out_dir, "SV.significant.by_set.remove_chr_edge_1Mb.tsv"),
  sep = "\t"
)

fwrite(
  snv_indel_sig_by_set,
  file.path(out_dir, "SNV_INDEL.significant.by_set.remove_chr_edge_1Mb.tsv"),
  sep = "\t"
)

lambda_sv_by_set <- calc_lambda(sv_sites_by_set$P)
lambda_snv_indel_by_set <- calc_lambda(snv_indel_sites_by_set$P)

cat("SV Lambda, by-set pooled P:", round(lambda_sv_by_set, 4), "\n")
cat("SNV/indel Lambda, by-set pooled P:", round(lambda_snv_indel_by_set, 4), "\n")

# ============================================================
# 6.     set          merge
#    key = CHR + POS + A1 + A2 + sv_len
# ============================================================

merge_list <- Map(
  f = prepare_for_merge,
  dt = dt_list,
  suffix = names(files)
)

merge_keys <- c("CHR", "POS", "A1", "A2", "sv_len")

merged_data <- Reduce(
  f = function(x, y) merge(x, y, by = merge_keys, all = TRUE),
  x = merge_list
)

# ============================================================
# 7.     min_P       SNP    
# ============================================================

p_cols <- grep("^P_", names(merged_data), value = TRUE)

merged_data[, min_P := do.call(
  pmin,
  c(.SD, list(na.rm = TRUE))
), .SDcols = p_cols]

merged_data[is.infinite(min_P), min_P := NA_real_]

snp_cols <- grep("^SNP_", names(merged_data), value = TRUE)

merged_data[, unified_SNP := do.call(
  fcoalesce,
  .SD
), .SDcols = snp_cols]

fwrite(
  merged_data,
  file.path(out_dir, "all_sets.merged_minP.remove_chr_edge_1Mb.standardized.tsv"),
  sep = "\t"
)

# ============================================================
# 8.  ?SV / SNV-indel      erged min_P    
# ============================================================

sv_sites_merged <- merged_data[sv_len >= sv_cutoff]
snv_indel_sites_merged <- merged_data[sv_len < sv_cutoff]

#        ?merged     SV  ?SNV/indel       
fwrite(
  sv_sites_merged,
  file.path(out_dir, "all_sets.merged_minP.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv"),
  sep = "\t"
)

fwrite(
  snv_indel_sites_merged,
  file.path(out_dir, "all_sets.merged_minP.SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.tsv"),
  sep = "\t"
)

# merged min_P       
sv_sig_merged <- sv_sites_merged[min_P < sig_p]
snv_indel_sig_merged <- snv_indel_sites_merged[min_P < sig_p]

fwrite(
  sv_sig_merged,
  file.path(out_dir, "SV.significant.merged_minP.remove_chr_edge_1Mb.tsv"),
  sep = "\t"
)

fwrite(
  snv_indel_sig_merged,
  file.path(out_dir, "SNV_INDEL.significant.merged_minP.remove_chr_edge_1Mb.tsv"),
  sep = "\t"
)

lambda_sv_minP <- calc_lambda(sv_sites_merged$min_P)
lambda_snv_indel_minP <- calc_lambda(snv_indel_sites_merged$min_P)

cat("SV Lambda, merged min_P:", round(lambda_sv_minP, 4), "\n")
cat("SNV/indel Lambda, merged min_P:", round(lambda_snv_indel_minP, 4), "\n")

# ============================================================
# 9.     lambda      
# ============================================================

lambda_summary <- data.table(
  category = c(
    "SV_by_set_pooled_P",
    "SNV_INDEL_by_set_pooled_P",
    "SV_merged_min_P",
    "SNV_INDEL_merged_min_P"
  ),
  lambda = c(
    lambda_sv_by_set,
    lambda_snv_indel_by_set,
    lambda_sv_minP,
    lambda_snv_indel_minP
  )
)

fwrite(
  lambda_summary,
  file.path(out_dir, "lambda_summary.remove_chr_edge_1Mb.tsv"),
  sep = "\t"
)

# ============================================================
# 10.           ?# ============================================================

set_count_summary <- rbindlist(lapply(names(dt_list), function(set_name) {
  data.table(dataset = paste0(set_name, "_all"), n = nrow(dt_list[[set_name]]))
}))

count_summary <- rbindlist(list(
  set_count_summary,
  data.table(
    dataset = c(
    "all_sets_long_all",
    "all_sets_long_SV_ge50bp",
    "all_sets_long_SNV_INDEL_lt50bp",
    "merged_all",
    "merged_SV_ge50bp",
    "merged_SNV_INDEL_lt50bp",
    "significant_SV_by_set",
    "significant_SNV_INDEL_by_set",
    "significant_SV_merged_minP",
      "significant_SNV_INDEL_merged_minP"
    ),
    n = c(
    nrow(per_set_data),
    nrow(sv_sites_by_set),
    nrow(snv_indel_sites_by_set),
    nrow(merged_data),
    nrow(sv_sites_merged),
    nrow(snv_indel_sites_merged),
    nrow(sv_sig_by_set),
    nrow(snv_indel_sig_by_set),
    nrow(sv_sig_merged),
      nrow(snv_indel_sig_merged)
    )
  )
))

fwrite(
  count_summary,
  file.path(out_dir, "variant_count_summary.remove_chr_edge_1Mb.tsv"),
  sep = "\t"
)

message("All done.")

