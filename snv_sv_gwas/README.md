# SNV/SV GWAS

GWAS result standardization, clumping, Manhattan/QQ plots, LD analyses, canonical lead selection, ANNOVAR and GWAS Catalog annotation, locus plots, SV-to-PAV matching, cis-meQTL analysis, and result collection for the PanGenie, DeepVariant, and Paragraph call sets.

## Requirements

- Bash, Python 3, R
- PLINK 1.9/2, GCTA (fastGWA), bcftools, ANNOVAR
- Local GWAS summary statistics and reference panels (controlled inputs, not included)

## Layout

| Path | Contents |
|---|---|
| `run_gwas_workflow.sh` | Module entry point for the numbered workflow. |
| `scripts/` | Numbered analysis scripts; see [docs/pipeline_steps.md](docs/pipeline_steps.md) for the step map. |
| `scripts/individual_level_tools/` | Carrier-sample extraction tools; controlled-data environment only, see [scripts/individual_level_tools/README.md](scripts/individual_level_tools/README.md). |
| `configs/gwas_config.example.R` | Path configuration example; see [configs/path_placeholders.md](configs/path_placeholders.md). |
| `docs/methods_gwas.md` | GWAS methods summary. |
| `docs/gwas_catalog_annotation.md` | GWAS Catalog annotation notes. |

## Configuration

Copy `configs/gwas_config.example.R` and replace every `/path/to/EOSCZ_PROJECT` placeholder with local paths. Several R scripts also accept `--burden-root`-style command-line overrides.

## Usage

```bash
bash run_gwas_workflow.sh [--start-at N] [--stop-after N] [--force]
```

The runner delegates to `scripts/run_pipeline.sh`, which executes numbered steps with checkpoints. Numbered scripts can also be run individually in the order given in [docs/pipeline_steps.md](docs/pipeline_steps.md). `scripts/14_sgv_meqtl_disabled.R` is an intentional placeholder: SGV cis-meQTL is disabled in this workflow (SV-only meQTL policy).

This directory retains its original [license](license).
