# Pipeline step map

This file maps the public scripts to analysis modules. It is intended for GitHub readers who want to inspect the workflow without access to project data.

## GWAS

- `00_run_gcta_fastgwa_binary_example.sh` gives an example GCTA command template.
- `01_gwas_merge_set.R` processes PanGenie set00 GWAS results, filters autosomes/chromosome edges, splits SV and SNV/indel, and writes merged minimum-P tables.
- `02_gwas_clumping.R` runs PLINK clumping.
- `03_gwas_draw_manhattan.R` draws Manhattan plots using PanGenie, DeepVariant, and Paragraph GWAS results.

## LD analyses

- `04_prepare_pangenie_maf01_lists.sh` and `04_prepare_deepvariant_paragraph_maf02_list.sh` prepare MAF-filtered variant lists.
- `04_ld_prepare_ids_and_plink_cmd.R` prepares significant and matched-null index variants and PLINK LD jobs.
- `05_ld_summarise_decay.R` summarizes LD decay and bootstrap AUC.
- `06_ld_friends_scores.R` calculates LD-friends/LD-score style summaries.
- `15_sv_detectable_ld_tests.R` performs formal SV detectable-LD tests against matched null SVs.

## Lead signal tables and annotation

- `06_extract_lead_sig_from_gwas.R` creates source-level lead tables.
- `06_select_canonical_window_leads.py` selects canonical 1 Mb window representatives.
- `07_run_annovar_t2t_refgene.sh` runs ANNOVAR.
- `08_annotate_gwas_catalog_scz_related.py` annotates GWAS Catalog schizophrenia-related overlaps.
- `09_merge_annotations_and_windows.py` merges final annotation and window tables.

## Locus plots

- `10_prepare_locuszoom_ld.R` prepares regional LD files.
- `11_plot_locuszoom_with_nature.R` draws locuszoom-style plots.

## meQTL and SV matching

- `12_match_sig_sv_to_pav.R` matches significant SVs to PAV/truvari SV calls.
- `13_run_sv_meqtl_cis.R` runs cis-meQTL models for matched SVs.
- `13b_qc_sv_meqtl_genotypes.R` and `13c_filter_sv_meqtl_by_mac.R` provide genotype QC and MAC filtering.
- `14_sgv_meqtl_disabled.R` documents the SGV meQTL policy/placeholder in this workflow.

## Carrier samples and summary collection

- `20_collect_results.sh` copies key outputs into summary packages.
- `extract_sv_nonref_samples.py` extracts non-reference carrier sample IDs for SV leads.
- `extract_snv_indel_nonref_samples.py` extracts carrier sample IDs for SNV/indel leads and prepares IGV batch snapshots.

