# Complement C4

C4 haplotype- and sample-level copy-number summaries and case/control association tests.

## Requirements

- Python 3 with `pandas`, `numpy`, `openpyxl`
- R with `readxl`, `writexl`

## Scripts

| Script | Purpose |
|---|---|
| `primary_analysis/c4_haplotype_cnv.py` | Counts `C4AL`/`C4BL`/`C4AS`/`C4BS` structures per haplotype row. |
| `primary_analysis/c4_sample_cnv.py` | Aggregates haplotype-level structures to sample-level counts. Cohort presets: `scz`, `comparison`, `public`; override with `--sample-regex`. |
| `primary_analysis/c4_fisher_test.R` | One-sided Fisher's exact test of C4 structure enrichment in SCZ vs a control group. |
| `strict_filtering/annotate_haplotype_calls_region.py` | Annotates C4 haplotype calls with region metadata. |
| `strict_filtering/verify_region_annotation.py` | Verifies region-annotation columns. |
| `supporting_analyses/c4_chi2_test.R` | Chi-square test of C4 structure counts. |

## Usage

```bash
python3 primary_analysis/c4_haplotype_cnv.py input.xlsx -o haplotype_counts.xlsx

python3 primary_analysis/c4_sample_cnv.py input.xlsx -o sample_counts.xlsx \
  --cohort-pattern public

Rscript primary_analysis/c4_fisher_test.R \
  --input=sample_counts.xlsx --sheet=Sheet1 --output=fisher_results.xlsx
```

Run `python3 <script> --help` or inspect script headers for the full argument list.
