# Configuration

## General template

Copy the template and edit the ignored local file:

```bash
cp configs/config.example.yaml configs/config.yaml
```

Not every historical script reads this YAML directly. It is a centralized checklist for paths and resources; retained shell scripts commonly use environment variables instead.

## Common environment variables

| Variable | Used for |
|---|---|
| `EOSCZ_PROJECT_ROOT` | Local project root for scripts that require a shared root |
| `EOSCZ_REFERENCE_DIR` | Reference and annotation directory |
| `EOSCZ_REFERENCE_FASTA` | Reference FASTA path |
| `BULK_RNA_WORK_DIR` | Bulk RNA-seq working directory |
| `BULK_RNA_INPUT_DIR` | Authorized bulk RNA-seq FASTQ directory |
| `FASTQ_DIR`, `OUT_ROOT` | Genome assembly input and output roots |
| `EOSCZ_SV_CONFIG` | Alternate structural-variant shell configuration path |
| `EOSCZ_SCRNA_INPUT_DIR` | Single-cell Matrix Market input directory |
| `EOSCZ_SNAPSEED_MARKERS` | Marker YAML path |
| `GWAS_CATALOG_ASSOC_ZIP` | Existing local GWAS Catalog associations zip |

Scripts use fail-fast checks for required variables where it could be done without changing the statistical logic. Other source-derived Python and R scripts still define placeholder paths near the top; edit only local copies or add an author-reviewed launcher.

## Slurm

Partition, account, module names, memory, and thread limits vary by site. Set them locally and do not commit internal hostnames, node exclusions, account names, or queue names.

## Secrets and privacy

No retained script requires a public token. If a local system requires credentials, inject them through an approved secret manager. Do not store credentials in YAML, `.env`, command history, scheduler logs, or Git.
