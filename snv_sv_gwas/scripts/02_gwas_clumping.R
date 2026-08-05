library(data.table)
library(qqman)

# ============================================================
# 0. Configuration
# ============================================================

out_dir <- "/path/to/EOSCZ_PROJECT/figure_analysis/01.GWAS_figure.public"

clump_dir <- file.path(out_dir, "clumping_by_set_subtype")
dir.create(clump_dir, recursive = TRUE, showWarnings = FALSE)

# Set this to the PLINK 1.9 executable if it is not available as "plink".
plink_bin <- "plink"

# Clumping parameters: p1 defines lead signals, p2 defines candidate variants,
# and r2/kb define the LD threshold and window size.
clump_p1 <- 5e-6
clump_p2 <- 0.05
clump_r2 <- 0.01
clump_kb <- 1000

# ============================================================
# 1. PLINK binary-file prefix for each call set.
# ============================================================

bfiles <- c(
  ## Public example: set00 only. Add set01/set02 here if available.
  set00 = "/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/NGS.QCsite.QCind"
)

# ============================================================
# 2. GWAS result files for each call set and variant class.
# ============================================================

gwas_files <- list(
  set00 = list(
    SV = file.path(out_dir, "set00.SV_ge50bp.remove_chr_edge_1Mb.standardized.fastGWA.tsv"),
    SNV_INDEL = file.path(out_dir, "set00.SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.fastGWA.tsv")
  )
)

# ============================================================
# 3. Clumping helper.
# ============================================================

run_clumping_one <- function(set_name, subtype, gwas_file, bfile_prefix) {
  
  message("============================================================")
  message("Clumping: ", set_name, " | ", subtype)
  message("GWAS file: ", gwas_file)
  message("BFILE: ", bfile_prefix)
  
  if (!file.exists(gwas_file)) {
    stop("GWAS file does not exist: ", gwas_file)
  }
  
  for (suffix in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(bfile_prefix, suffix))) {
      stop("Missing PLINK bfile component: ", paste0(bfile_prefix, suffix))
    }
  }
  
  dt <- fread(gwas_file)
  
  required_cols <- c("CHR", "SNP", "P")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("GWAS file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dt[, SNP := as.character(SNP)]
  dt[, P := as.numeric(P)]
  dt[, CHR := sub("^chr", "", as.character(CHR), ignore.case = TRUE)]
  dt <- dt[CHR %in% as.character(1:22)]
  
  # PLINK clumping requires valid SNP identifiers and P values.
  clump_input <- dt[
    !is.na(SNP) &
      SNP != "" &
      !is.na(P) &
      P > 0 &
      P <= 1,
    .(SNP, P)
  ]
  
  # Deduplicate SNP identifiers while retaining the smallest P value.
  setorder(clump_input, P)
  clump_input <- unique(clump_input, by = "SNP")
  
  if (nrow(clump_input) == 0) {
    warning(set_name, " ", subtype, ": no valid SNP/P records available for clumping.")
    return(NULL)
  }
  
  prefix <- file.path(
    clump_dir,
    paste0(set_name, ".", subtype, ".clump_p1_", clump_p1, ".r2_", clump_r2, ".kb_", clump_kb)
  )
  
  clump_input_file <- paste0(prefix, ".plink_clump_input.tsv")
  
  fwrite(
    clump_input,
    clump_input_file,
    sep = "\t"
  )
  
  # ==========================================================
  #     PLINK clumping
  # ==========================================================
  
  cmd <- sprintf(
    '%s --bfile %s --threads 1 --clump %s --clump-snp-field SNP --clump-field P --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %g --allow-extra-chr --out %s',
    shQuote(plink_bin),
    shQuote(bfile_prefix),
    shQuote(clump_input_file),
    clump_p1,
    clump_p2,
    clump_r2,
    clump_kb,
    shQuote(prefix)
  )
  
  message("Running command:")
  message(cmd)
  
  status <- system(cmd)
  
  if (status != 0) {
    warning("PLINK clumping       : ", set_name, " ", subtype)
    return(NULL)
  }
  
  clumped_file <- paste0(prefix, ".clumped")
  
  if (!file.exists(clumped_file)) {
    warning("No .clumped file was produced; there may be no variants below clump-p1: ", set_name, " ", subtype)
    return(NULL)
  }
  
  clumped <- tryCatch(
    fread(clumped_file, fill = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(clumped) || nrow(clumped) == 0) {
    warning(".clumped       : ", clumped_file)
    return(NULL)
  }
  
  # PLINK records lead-variant identifiers in the SNP column.
  if (!"SNP" %in% names(clumped)) {
    warning("The .clumped file has no SNP column; lead signals cannot be extracted: ", clumped_file)
    return(NULL)
  }
  
  lead_snps <- unique(clumped$SNP)
  
  independent_leads <- dt[SNP %in% lead_snps]
  
  # Sort independent signals by P value.
  setorder(independent_leads, P)
  
  independent_out <- paste0(prefix, ".independent_lead_signals.tsv")
  
  fwrite(
    independent_leads,
    independent_out,
    sep = "\t"
  )
  
  message("Saved PLINK clumped file: ", clumped_file)
  message("Saved independent lead signals: ", independent_out)
  message("Number of independent lead signals: ", nrow(independent_leads))
  
  return(list(
    set = set_name,
    subtype = subtype,
    input_n = nrow(dt),
    clump_input_n = nrow(clump_input),
    lead_n = nrow(independent_leads),
    clumped_file = clumped_file,
    independent_file = independent_out
  ))
}

# ============================================================
# 4. Run clumping for every call set and variant class.
# ============================================================

clump_results <- list()

for (set_name in names(gwas_files)) {
  
  for (subtype in names(gwas_files[[set_name]])) {
    
    res <- run_clumping_one(
      set_name = set_name,
      subtype = subtype,
      gwas_file = gwas_files[[set_name]][[subtype]],
      bfile_prefix = bfiles[[set_name]]
    )
    
    clump_results[[paste(set_name, subtype, sep = ".")]] <- res
  }
}

# ============================================================
# 5. Write a clumping summary.
# ============================================================

clump_summary <- rbindlist(
  lapply(clump_results, function(x) {
    if (is.null(x)) return(NULL)
    data.table(
      set = x$set,
      subtype = x$subtype,
      input_n = x$input_n,
      clump_input_n = x$clump_input_n,
      lead_n = x$lead_n,
      clumped_file = x$clumped_file,
      independent_file = x$independent_file
    )
  }),
  fill = TRUE
)

fwrite(
  clump_summary,
  file.path(clump_dir, "clumping_summary.by_set_subtype.tsv"),
  sep = "\t"
)

message("All clumping done.")

