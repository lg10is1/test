## Figure 2 SV/SNV GWAS analysis configuration.
## Copy this file to config/figure2_config.R and edit paths before running.

## Output directory. Relative paths are interpreted from the directory where
## run_gwas_workflow.sh is launched.
output_root <- "results/figure2_gwas"

## Software.
plink_bin <- Sys.getenv("PLINK_BIN", "plink")

## Core analysis parameters.
edge_bp <- 1000000L
sv_length_cutoff_bp <- 50L
genomewide_p <- 5e-8
suggestive_p <- 5e-6

## PLINK clumping parameters.
clump_p1 <- 5e-6
clump_p2_pangenie <- 0.05
clump_p2_external <- 1
clump_r2 <- 0.01
clump_kb <- 1000
plink_threads <- 1

## Additional PLINK arguments by source.
## For PanGenie, keep --maf 0.01 if the genotype bfile has not already been
## filtered to the same MAF threshold as the GWAS.
plink_extra_args <- list(
  pangenie = "--maf 0.01",
  deepvariant = "",
  paragraph = ""
)

## PanGenie input. The public example uses set00 only. Add set01/set02 rows if
## you want to process all PanGenie sets.
pangenie_sets <- data.frame(
  dataset_id = c("pangenie_set00"),
  set_id = c("set00"),
  gwas_file = c("/path/to/Pangenie_v3/06.gwas/set00/gwas/SCZ.mlm.ngspc.fastGWA"),
  bfile_prefix = c("/path/to/Pangenie_v3/06.gwas/set00/NGS.QCsite.QCind"),
  stringsAsFactors = FALSE
)

## Optional external GWAS files. Set enabled = FALSE to skip a row.
external_gwas <- data.frame(
  dataset_id = c("deepvariant_snv", "paragraph_sv"),
  source = c("deepvariant", "paragraph"),
  variant_class = c("SNV_INDEL", "SV"),
  gwas_file = c(
    "/path/to/Deepvariant_paragraph/deepvar_gwas/deepvar/03_gwas/SCZ.deepvar.mlm.geno0.1.maf0.01.fastGWA",
    "/path/to/Deepvariant_paragraph/deepvar_gwas/paragraph_test/03_gwas/SCZ.paragraph_test.MODEL.mlm.geno0.1.maf0.01.fastGWA"
  ),
  bfile_prefix = c(
    "/path/to/Deepvariant_paragraph/chr_all2.strict_step2_genimi.common_samples.merged",
    "/path/to/Deepvariant_paragraph/chr_all2.strict_step2_genimi.common_samples.merged"
  ),
  enabled = c(TRUE, TRUE),
  stringsAsFactors = FALSE
)

## Column aliases are auto-detected. Override here only if your files use
## nonstandard names.
column_aliases <- list(
  chr = c("CHR", "#CHROM", "CHROM", "chrom", "chr"),
  pos = c("POS", "BP", "position", "pos"),
  snp = c("SNP", "ID", "variant_id", "MarkerName", "rsid"),
  a1 = c("A1", "ALT", "Allele1", "effect_allele"),
  a2 = c("A2", "REF", "Allele2", "other_allele"),
  p = c("P", "PVAL", "P_VALUE", "p", "pvalue")
)
