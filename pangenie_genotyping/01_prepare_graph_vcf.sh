#!/usr/bin/env bash

# Prepare the graph VCF used by the PanGenie workflow.
#
# Usage:
#   bash 01_prepare_graph_vcf.sh /path/to/Pangenie_v3

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PROJECT_DIR" >&2
    exit 2
fi

PROJECT_DIR="$(cd "$1" && pwd)"
SET_ID="00"
INPUT_VCF="${PROJECT_DIR}/01.split_vcf/filtered_part_${SET_ID}.vcf.gz"
WORK_DIR="${PROJECT_DIR}/02.pangenome/set${SET_ID}"
LINKED_VCF="${WORK_DIR}/filtered_part_${SET_ID}.vcf.gz"
OUTPUT_VCF="${WORK_DIR}/graph.vcf"
TMP_VCF="${OUTPUT_VCF}.tmp.$$"

command -v bcftools >/dev/null 2>&1 || {
    echo "ERROR: bcftools is not available in PATH." >&2
    exit 1
}

[[ -f "${INPUT_VCF}" ]] || {
    echo "ERROR: input VCF not found: ${INPUT_VCF}" >&2
    exit 1
}

mkdir -p "${WORK_DIR}"
ln -sfn "${INPUT_VCF}" "${LINKED_VCF}"

# Link the available index as well. Either TBI or CSI is accepted.
if [[ -f "${INPUT_VCF}.tbi" ]]; then
    ln -sfn "${INPUT_VCF}.tbi" "${LINKED_VCF}.tbi"
elif [[ -f "${INPUT_VCF}.csi" ]]; then
    ln -sfn "${INPUT_VCF}.csi" "${LINKED_VCF}.csi"
else
    echo "WARNING: no .tbi or .csi index found for ${INPUT_VCF}." >&2
fi

trap 'rm -f "${TMP_VCF}"' EXIT

echo "[START] Preparing graph VCF: ${INPUT_VCF}"

# Replace fully missing genotypes (./.) with phased reference genotypes (0|0).
# -Ov is intentional: graph.vcf is written as a plain-text VCF, not as BCF.
bcftools +setGT "${LINKED_VCF}" -Ov -o "${TMP_VCF}" -- -t './.' -n 0p
mv -f "${TMP_VCF}" "${OUTPUT_VCF}"

trap - EXIT
echo "[DONE] ${OUTPUT_VCF}"
