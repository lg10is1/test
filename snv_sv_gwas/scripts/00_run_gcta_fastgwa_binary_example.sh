#!/usr/bin/env bash
set -euo pipefail

## Example only: this script is not called by run_gwas_workflow.sh.
## Edit the paths below to match your cohort and genotype source.

GCTA_BIN="${GCTA_BIN:-gcta64}"

## PLINK genotype prefix for one variant call set, for example:
## - PanGenie set00
## - DeepVariant SNV/indel
## - Paragraph SV
BFILE_PREFIX="/path/to/genotypes"

## Phenotype file accepted by GCTA. For binary traits, cases/controls are
## usually coded as 1/0 or 2/1 depending on the project convention.
PHENO_FILE="/path/to/phenotype.txt"

## Fixed-effect covariates, e.g. batch/sex if included.
COVAR_FILE="/path/to/covariates.txt"

## Quantitative covariates: genetic PC1-PC20. The file should contain 20 PC columns.
QCOVAR_FILE="/path/to/genetic_PC1_to_PC20.qcovar"

## Sparse genetic relationship matrix generated from a matched genotype source.
SPARSE_GRM_PREFIX="/path/to/sparse_grm"

OUT_PREFIX="results/gcta_fastgwa/example.mlm.geno0.1.maf0.01"

mkdir -p "$(dirname "${OUT_PREFIX}")"

"${GCTA_BIN}" \
  --bfile "${BFILE_PREFIX}" \
  --grm-sparse "${SPARSE_GRM_PREFIX}" \
  --fastGWA-mlm-binary \
  --pheno "${PHENO_FILE}" \
  --covar "${COVAR_FILE}" \
  --qcovar "${QCOVAR_FILE}" \
  --geno 0.1 \
  --maf 0.01 \
  --thread-num 8 \
  --out "${OUT_PREFIX}"



