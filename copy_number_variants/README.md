# Copy number variants

CNV copy-number matrices, case/control association tests, pathway burden analyses, PCA, and figure workflows (Fig. 5 and supplements).

## Requirements

- Python 3 with `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `upsetplot`, `openpyxl`
- R with `data.table`, `ggplot2`, `readxl`, `writexl`, `openxlsx`, `karyoploteR`

## Layout

| Directory | Contents |
|---|---|
| `primary_analysis/` | Sample-level CNV matrix generation, primary Fisher tests, and Fig. 5a-5e plots. |
| `pathway_burden/` | Published gene-set extraction and rare-CNV pathway burden models (logistic/Fisher, with FDR control). |
| `copy_number_range_plots/` | Fig. 5d/5e copy-number range plots (0-100 crop, 1-100 crop with legend, full range). |
| `supporting_analyses/` | Sensitivity Fisher tests (common/rare, control-only, East-Asian subset, reverse comparison), CNV PCA plots, haplotype-level plots, and QC scripts. |

## Usage

Scripts read input workbooks/tables from paths given on the command line or at the top of each script (`/path/to/EOSCZ_PROJECT` placeholders must be replaced locally). Examples:

```bash
python3 primary_analysis/generate_sample_cn_filtered_cases.py --help

Rscript pathway_burden/run_published_gene_sets_36genesets_sex_globalrare_logistic_fisher.R \
  --burden-root=/path/to/cnv_analysis/pathway_burden
```

The two `run_published_gene_sets_controlrare_*` scripts are model variants: `..._with_sc.R` adds a cohort-size scaling covariate, `...nosczrecurrent_....R` does not.
