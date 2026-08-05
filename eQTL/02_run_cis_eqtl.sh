#!/usr/bin/env bash

# Usage: bash 02_run_cis_eqtl.sh GENOTYPES.vcf.gz PHENOTYPES.bed.gz COVARIATES.tsv OUTPUT.txt

set -euo pipefail

VCF="$1"
BED="$2"
COVARIATES="$3"
OUTPUT="$4"
NOMINAL_P="${NOMINAL_P:-0.001}"

QTLtools cis \
    --vcf "${VCF}" \
    --bed "${BED}" \
    --cov "${COVARIATES}" \
    --nominal "${NOMINAL_P}" \
    --out "${OUTPUT}"

