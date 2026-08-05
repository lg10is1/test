#!/usr/bin/env bash

# Run the single retained GCTA GWAS model:
#   SCZ.mlm.ngspc.fastGWA
#
# The model uses the sparse GRM, the batch covariate, and all 20 NGS PCs from
# pca.ngs.eigenvec.
#
# Usage:
#   bash 05_run_gwas.sh PROJECT_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PROJECT_DIR" >&2
    exit 2
fi

PROJECT_DIR="$(cd "$1" && pwd)"
BASE_WORK_DIR="${PROJECT_DIR}/06.gwas"
DATA_DIR="set00"
THREADS="${THREADS:-6}"
PHENO_FILE="${PHENO_FILE:-${BASE_WORK_DIR}/SCZ_pheno.txt}"
BATCH_FILE="${BATCH_FILE:-${BASE_WORK_DIR}/SCZ_batch.txt}"

if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer." >&2
    exit 2
fi

for tool in gcta64 awk; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: required command is not available in PATH: ${tool}" >&2
        exit 1
    }
done

for file in "${PHENO_FILE}" "${BATCH_FILE}"; do
    [[ -f "${file}" ]] || {
        echo "ERROR: required GWAS input not found: ${file}" >&2
        exit 1
    }
done

run_gwas() {
    local set_dir="${BASE_WORK_DIR}/${DATA_DIR}"
    local gwas_dir="${set_dir}/gwas"
    local bfile="${set_dir}/NGS.QCsite.QCind"
    local sparse_grm="${set_dir}/sp_grm"
    local pc_ngs="${set_dir}/pca.ngs.eigenvec"
    local output_prefix="${gwas_dir}/SCZ.mlm.ngspc"
    local ngs_pc_count

    echo "[$(date +'%H:%M:%S')] [START] SCZ.mlm.ngspc"

    for suffix in bed bim fam; do
        [[ -f "${bfile}.${suffix}" ]] || {
            echo "ERROR: missing PLINK input: ${bfile}.${suffix}" >&2
            return 1
        }
    done

    for file in \
        "${sparse_grm}.grm.sp" \
        "${sparse_grm}.grm.id" \
        "${pc_ngs}"; do
        [[ -s "${file}" ]] || {
            echo "ERROR: required set-level input not found: ${file}" >&2
            return 1
        }
    done

    # PLINK --pca 20 writes FID, IID, followed by exactly 20 PC columns.
    ngs_pc_count="$(awk 'NF >= 3 { print NF - 2; exit }' "${pc_ngs}")"
    if [[ "${ngs_pc_count}" != "20" ]]; then
        echo "ERROR: ${pc_ngs} contains ${ngs_pc_count:-0} PCs; expected 20." >&2
        return 1
    fi
    echo "[CHECK] Using all 20 PCs from ${pc_ngs}"

    mkdir -p "${gwas_dir}"

    if [[ -s "${output_prefix}.fastGWA" ]]; then
        echo "[SKIP] Existing result: ${output_prefix}.fastGWA"
        return
    fi

    gcta64 \
        --bfile "${bfile}" \
        --pheno "${PHENO_FILE}" \
        --thread-num "${THREADS}" \
        --geno 0.1 \
        --maf 0.01 \
        --grm-sparse "${sparse_grm}" \
        --fastGWA-mlm-binary \
        --qcovar "${pc_ngs}" \
        --covar "${BATCH_FILE}" \
        --out "${output_prefix}"

    [[ -s "${output_prefix}.fastGWA" ]] || {
        echo "ERROR: expected output was not created: ${output_prefix}.fastGWA" >&2
        return 1
    }

    echo "[$(date +'%H:%M:%S')] [DONE] ${output_prefix}.fastGWA"
}

mkdir -p "${BASE_WORK_DIR}"
run_gwas
