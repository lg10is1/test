#!/usr/bin/env bash
#SBATCH --job-name=call_5mc
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --output=call_5mc.%j.out
#SBATCH --error=call_5mc.%j.err

# Usage: sbatch 01_call_5mc.slurm INPUT.bam OUTPUT_PREFIX

set -euo pipefail

BAM="$1"
OUTPUT_PREFIX="$2"
THREADS="${SLURM_CPUS_PER_TASK:-40}"
TMP_ROOT="${TMP_ROOT:-${TMPDIR:-/tmp}}"
TASK_TMP="$(mktemp -d "${TMP_ROOT%/}/pb_cpg.XXXXXX")"

trap 'rm -rf "${TASK_TMP}"' EXIT
export TMPDIR="${TASK_TMP}"

mkdir -p "$(dirname "${OUTPUT_PREFIX}")"

aligned_bam_to_cpg_scores \
    --bam "${BAM}" \
    --output-prefix "${OUTPUT_PREFIX}" \
    --threads "${THREADS}"

