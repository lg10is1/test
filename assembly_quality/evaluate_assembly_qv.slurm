#!/usr/bin/env bash
#SBATCH --job-name=assembly_qv
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${MERYL_INPUT_DIR:?Set MERYL_INPUT_DIR}"
: "${MERYL_OUTPUT_DIR:?Set MERYL_OUTPUT_DIR}"
: "${ASSEMBLY_FASTA_DIR:?Set ASSEMBLY_FASTA_DIR}"
: "${MERQURY_SH:?Set MERQURY_SH to merqury.sh}"
SAMPLE_ID="${1:?Usage: sbatch evaluate_assembly_qv.slurm SAMPLE_ID [ASSEMBLY_ID]}"
ASSEMBLY_ID="${2:-$SAMPLE_ID}"
THREADS="${THREADS:-20}"
MEMORY_GB="${MEMORY_GB:-80}"
KMER_SIZE="${KMER_SIZE:-21}"
ANALYSIS_DIR="${MERYL_ANALYSIS_DIR:-${MERYL_OUTPUT_DIR}/analysis}"

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-merqury}"
mkdir -p "$MERYL_OUTPUT_DIR" "$ANALYSIS_DIR/$ASSEMBLY_ID" task_file
cd "$MERYL_OUTPUT_DIR"
meryl k="$KMER_SIZE" threads="$THREADS" memory="$MEMORY_GB" count \
  "$MERYL_INPUT_DIR/${SAMPLE_ID}_R1.fastq.gz" output "${SAMPLE_ID}_R1.meryl"
meryl k="$KMER_SIZE" threads="$THREADS" memory="$MEMORY_GB" count \
  "$MERYL_INPUT_DIR/${SAMPLE_ID}_R2.fastq.gz" output "${SAMPLE_ID}_R2.meryl"
meryl threads="$THREADS" memory="$MEMORY_GB" union-sum \
  output "${SAMPLE_ID}.merged.meryl" "${SAMPLE_ID}_R1.meryl" "${SAMPLE_ID}_R2.meryl"
cd "$ANALYSIS_DIR/$ASSEMBLY_ID"
"$MERQURY_SH" "$MERYL_OUTPUT_DIR/${SAMPLE_ID}.merged.meryl" \
  "$ASSEMBLY_FASTA_DIR/${ASSEMBLY_ID}.1.fa" \
  "$ASSEMBLY_FASTA_DIR/${ASSEMBLY_ID}.2.fa" "$ASSEMBLY_ID"
