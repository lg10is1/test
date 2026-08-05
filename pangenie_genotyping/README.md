# PanGenie processing and GWAS workflow

This directory contains the PanGenie processing, genotype QC, downstream
variant/sample QC, PCA, GRM, and GWAS scripts used in this project. Personal
HPC account paths and batch-generation code are not included.

## Workflow

Run the scripts in this order:

```text
01_prepare_graph_vcf.sh
02_run_pangenie.slurm
03_qc_single_sample.slurm
04_run_bfile_qc.sh
05_run_gwas.sh
```

Genotype consistency is calculated inside `04_run_bfile_qc.sh`; there is no
separate consistency script or additional consistency step.

## Requirements

- Bash
- `bcftools` with the `+setGT` plugin and `tabix`
- Snakemake and the configured PanGenie workflow
- PLINK 1.9, PLINK 2, and GCTA
- Python 3
- `pandas` and `plotly` for PCA HTML plots
- A Slurm cluster for the `.slurm` scripts

The consistency calculation calls `genotype_consistency.py`. By default, this
program must be available at:

```text
Pangenie_v3/06.gwas/genotype_consistency.py
```

Set `PYTHON_SCRIPT=/path/to/genotype_consistency.py` if it is stored elsewhere.

## Expected project layout

```text
Pangenie_v3/
|-- 01.split_vcf/
|   |-- filtered_part_00.vcf.gz
|   `-- sample_map.txt
|-- 02.pangenome/
|-- 03.pangenie_results/
|   `-- set00/
|       `-- SAMPLE001_genotyping.vcf.gz
|-- 04.pangenie_qc/
|-- 05.merge/
|   `-- set00/
|       `-- step1.split.vcf.gz
`-- 06.gwas/
    |-- genotype_consistency.py
    |-- SCZ_pheno.txt
    |-- SCZ_batch.txt
    |-- tgs_sample.txt
    |-- bad.sample
    |-- eas.sample
    `-- exclude_complex
```

## 1. Prepare the graph VCF

```bash
bash 01_prepare_graph_vcf.sh /path/to/Pangenie_v3
```

This creates `02.pangenome/set00/graph.vcf` after changing fully missing
genotypes (`./.`) to phased reference genotypes (`0|0`).

## 2. Run PanGenie

```bash
sbatch 02_run_pangenie.slurm /path/to/pangenie/pipeline
```

The default Conda/Mamba environment is named `pangenie`. Override it when
needed:

```bash
PANGENIE_ENV=/path/to/conda/envs/pangenie \
    sbatch --export=ALL 02_run_pangenie.slurm /path/to/pangenie/pipeline
```

The PanGenie/Snakemake configuration must point to the prepared graph VCF and
the intended sample BAM/CRAM.

## 3. Apply GQ60 genotype QC

Run this job once for each sample:

```bash
sbatch 03_qc_single_sample.slurm SAMPLE001 /path/to/Pangenie_v3
```

The output files are:

```text
04.pangenie_qc/set00/SAMPLE001.gq60.vcf.gz
04.pangenie_qc/set00/SAMPLE001.gq60.vcf.gz.tbi
```

Genotypes with `FORMAT/GQ < 60` are changed to missing. The `GQ`, `GL`, and
`KC` FORMAT fields are removed, and records with 2-5 alleles are retained.

## 4. Run BFILE QC, consistency, PCA, and GRM

```bash
bash 04_run_bfile_qc.sh /path/to/Pangenie_v3
```

This script performs the following operations in one pipeline:

1. Converts `05.merge/set00/step1.split.vcf.gz` to PLINK BED files.
2. Calculates chromosome-level genotype consistency and merges the results.
3. Retains variants with consistency-report column 18 `R2 > 0.4` and column 13
   `concordance > 0.7`, then applies missingness and complex-region filters.
4. Performs sample QC and PCA.
5. Creates `pca.ngs.eigenvec` containing 20 NGS principal components.
6. Creates sparse and full GCTA GRMs.

The chromosome consistency calculation runs up to 12 jobs concurrently. Use
`MAX_CHROMOSOME_JOBS` to change this value:

```bash
MAX_CHROMOSOME_JOBS=6 bash 04_run_bfile_qc.sh /path/to/Pangenie_v3
```

## 5. Run GWAS adjusted for 20 PCs

```bash
bash 05_run_gwas.sh /path/to/Pangenie_v3
```

Only the following GWAS model is run:

```text
SCZ.mlm.ngspc.fastGWA
```

The model is adjusted for population structure using 20 NGS-derived principal
components (PC1-PC20). All 20 PCs from
`06.gwas/set00/pca.ngs.eigenvec` are supplied through GCTA's `--qcovar` option.
The model also uses the sparse GRM and `SCZ_batch.txt`. The script stops if the
PC file does not contain exactly 20 PC columns.
