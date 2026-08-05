#!/usr/bin/env bash
#SBATCH --job-name=pangenie_qc
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=pangenie_qc.%j.out
#SBATCH --error=pangenie_qc.%j.err

# Apply the final PanGenie genotype QC to one sample.
# Genotypes with FORMAT/GQ < 60 are changed to missing. FORMAT/GQ, FORMAT/GL,
# and FORMAT/KC are then removed, and sites with 2-5 alleles are retained.
#
# Usage:
#   sbatch 03_qc_single_sample.slurm SAMPLE_ID /path/to/Pangenie_v3
#
# Optional environment variables:
#   PANGENIE_ENV     Conda/Mamba environment name or path (default: pangenie)
#   CONDA_INIT_FILE  Shell initialization file (default: $HOME/.bashrc)

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 SAMPLE_ID PROJECT_DIR" >&2
    exit 2
fi

SAMPLE_ID="$1"
PROJECT_DIR="$(cd "$2" && pwd)"
PANGENIE_ENV="${PANGENIE_ENV:-pangenie}"
CONDA_INIT_FILE="${CONDA_INIT_FILE:-${HOME}/.bashrc}"

INPUT_VCF="${PROJECT_DIR}/03.pangenie_results/set00/${SAMPLE_ID}_genotyping.vcf.gz"
OUTPUT_DIR="${PROJECT_DIR}/04.pangenie_qc/set00"
OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE_ID}.gq60.vcf.gz"
TMP_VCF="${OUTPUT_DIR}/.${SAMPLE_ID}.gq60.tmp.$$.vcf.gz"

[[ -f "${INPUT_VCF}" ]] || {
    echo "ERROR: input VCF not found: ${INPUT_VCF}" >&2
    exit 1
}

if [[ -f "${CONDA_INIT_FILE}" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_INIT_FILE}"
fi

if command -v mamba >/dev/null 2>&1; then
    mamba activate "${PANGENIE_ENV}"
elif command -v conda >/dev/null 2>&1; then
    conda activate "${PANGENIE_ENV}"
else
    echo "ERROR: neither mamba nor conda is available after initialization." >&2
    exit 1
fi

for tool in bcftools tabix; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: ${tool} is not available in the activated environment." >&2
        exit 1
    }
done

mkdir -p "${OUTPUT_DIR}"
trap 'rm -f "${TMP_VCF}" "${TMP_VCF}.tbi"' EXIT

echo "[START] ${SAMPLE_ID}"

bcftools +setGT "${INPUT_VCF}" -Ou -- -t q -n . -i 'FMT/GQ<60' | \
    bcftools annotate -Ou -x FORMAT/GQ,FORMAT/GL,FORMAT/KC | \
    bcftools view -m2 -M5 -Oz -o "${TMP_VCF}"

tabix -f -p vcf "${TMP_VCF}"
mv -f "${TMP_VCF}" "${OUTPUT_VCF}"
mv -f "${TMP_VCF}.tbi" "${OUTPUT_VCF}.tbi"

trap - EXIT
echo "[DONE] ${OUTPUT_VCF}"
