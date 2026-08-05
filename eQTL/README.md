# cis-eQTL analysis and fine-mapping

## Requirements

- R with `data.table` and `peer`
- `rnanorm`, `bgzip`, and `tabix`
- QTLtools
- PLINK 1.9 and SuSiEx

## 1. RNA-seq QC and phenotype preparation

```bash
Rscript 01_rnaseq_qc.R \
    gene_counts.csv \
    sample_names.txt \
    gene_coordinates.tsv \
    base_covariates.tsv \
    pca_projection.sscore \
    output_dir
```

The script retains genes with counts greater than 6 in more than 80% of
samples, performs TMM normalization and inverse-normal transformation, creates
15 PEER factors, and writes:

```text
rna_tmm_inv.bed.gz
rna_tmm_inv.bed.gz.tbi
covariates_qtltools.tsv
```

`gene_coordinates.tsv` must contain `gene`, `chr`, `start`, `end`, and `strand`
columns. The base covariate and PCA files must contain an `IID` column.

## 2. Nominal cis-eQTL mapping

```bash
bash 02_run_cis_eqtl.sh \
    genotypes.vcf.gz \
    output_dir/rna_tmm_inv.bed.gz \
    output_dir/covariates_qtltools.tsv \
    nominal_cis_eqtl.txt
```

The nominal threshold is 0.001. The script does not set `--window`, so QTLtools
uses its default 1 Mb cis window. Change the reporting threshold with
`NOMINAL_P`. See the
[QTLtools cis documentation](https://qtltools.github.io/qtltools/pages/QTLtools-cis.1.html).

## 3. Fine-mapping

```bash
Rscript 03_finemapping.R \
    nominal_cis_eqtl.txt \
    /path/to/reference_bfile \
    finemapping_results \
    139 \
    16
```

The last two arguments are the GWAS sample size and thread count. Fine-mapping
is run for non-X phenotypes whose top cis-eQTL has `P < 5e-3`.
