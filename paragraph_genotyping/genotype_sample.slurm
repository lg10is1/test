#!/usr/bin/env bash
#SBATCH --job-name=paragraph_genotype
#SBATCH -N 1
#SBATCH --ntasks-per-node=28
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${PARAGRAPH_MULTIGRMPY:?Set PARAGRAPH_MULTIGRMPY to multigrmpy.py}"
: "${INPUT_VCF:?Set INPUT_VCF}"
: "${SAMPLE_MANIFEST:?Set SAMPLE_MANIFEST}"
: "${EOSCZ_REFERENCE_FASTA:?Set EOSCZ_REFERENCE_FASTA}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR}"
THREADS="${THREADS:-28}"
TMPDIR="${TMPDIR:-${OUTPUT_DIR}/tmp}"
export TMPDIR

module load miniconda3
source activate "${CONDA_ENV:-paragraph}"
mkdir -p "$OUTPUT_DIR" "$TMPDIR" task_file
python3 "$PARAGRAPH_MULTIGRMPY" -i "$INPUT_VCF" -m "$SAMPLE_MANIFEST" \
  -r "$EOSCZ_REFERENCE_FASTA" -o "$OUTPUT_DIR" -t "$THREADS"
