# MHC / HLA

HLA subtype and haplotype summaries, novel-allele analyses, association tests, and heatmaps (Fig. 4 and supplements). Inputs are Immuannot result workbooks.

## Requirements

- Python 3 with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`
- R with `ggplot2`, `pheatmap`/`ComplexHeatmap`, `readxl`

## Scripts

| Script | Purpose |
|---|---|
| `primary_analysis/haplotype_summary_barplot_universal.py` | Fig. 4A HLA haplotype summary bar plot. |
| `primary_analysis/subtype_summary_barplot_universal.py` | Fig. 4B HLA subtype summary bar plot. |
| `primary_analysis/fig4c_true_novel_typing_digit_barplot_universal.py` | Fig. 4C true novel HLA allele counts by 2-/3-/4-field class. |
| `primary_analysis/abc_fisher_one_sided_unknown_removed_universal.py` | Fig. 4F HLA-A/B/C one-sided frequency tests, case cohort vs comparison cohort. |
| `primary_analysis/novel_typing_barplot_scz_unique_universal.py` | Novel-typing bar plot for case-unique alleles. |
| `primary_analysis/hla_abc_correlation_heatmap_universal.R` | HLA-A/B/C typing correlation heatmap. |
| `primary_analysis/hla_abc_correlation_heatmap_batch_universal.py` | Batch driver for the correlation heatmap. |
| `primary_analysis/hla_abc_significant_typing_correlation_heatmap_universal.py` | Heatmap restricted to significant typings. |
| `supporting_analyses/novel_typing_digit_distribution_universal.py` | Novel-typing counts by digit level. |
| `supporting_analyses/export_novel_candidate_haplotype_map_case_cohort.py` | Exports novel-candidate haplotype maps. |
| `supporting_analyses/export_reverse_comparison_known_and_novel_candidates_case_cohort.py` | Exports known and novel candidates for the reverse comparison. |
| `supporting_analyses/filter_known_hla_reverse_case_cohort.py` | Filters known HLA alleles for the reverse comparison. |
| `supporting_analyses/run_reverse_comparison_hla_case_cohort_comparison_cohort_public_reference_east_asian_subset.py` | Reverse-comparison Fisher tests against the comparison cohort and the public East-Asian reference subset. |

## Usage

```bash
python3 primary_analysis/subtype_summary_barplot_universal.py immuannot_case.xlsx --total-count 420

python3 primary_analysis/abc_fisher_one_sided_unknown_removed_universal.py \
  immuannot_comparison.xlsx immuannot_case.xlsx
```

Run `python3 <script> --help` for the full argument list.
