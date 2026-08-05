#!/usr/bin/env bash
#SBATCH --job-name=pangenie
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=pangenie.%j.out
#SBATCH --error=pangenie.%j.err

# Run the configured PanGenie Snakemake workflow.
#
# Usage:
#   sbatch 02_run_pangenie.slurm /path/to/pangenie/pipeline
#
# Optional environment variables:
#   PANGENIE_ENV     Conda/Mamba environment name or path (default: pangenie)
#   CONDA_INIT_FILE  Shell initialization file (default: $HOME/.bashrc)

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 PANGENIE_PIPELINE_DIR" >&2
    exit 2
fi

PIPELINE_DIR="$1"
PANGENIE_ENV="${PANGENIE_ENV:-pangenie}"
CONDA_INIT_FILE="${CONDA_INIT_FILE:-${HOME}/.bashrc}"

[[ -d "${PIPELINE_DIR}" ]] || {
    echo "ERROR: pipeline directory not found: ${PIPELINE_DIR}" >&2
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

command -v snakemake >/dev/null 2>&1 || {
    echo "ERROR: snakemake is not available in the activated environment." >&2
    exit 1
}

cd "${PIPELINE_DIR}"

echo "[START] PanGenie workflow in ${PIPELINE_DIR}"
snakemake \
    --cores "${SLURM_CPUS_PER_TASK:-20}" \
    --printshellcmds \
    --rerun-incomplete
echo "[DONE] PanGenie workflow"
