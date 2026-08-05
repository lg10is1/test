#!/usr/bin/env bash
#SBATCH --job-name=bulk_rna_quantify
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${BULK_RNA_WORK_DIR:?Set BULK_RNA_WORK_DIR}"
: "${SAMPLE_NAME:?Set SAMPLE_NAME}"
THREADS="${THREADS:-8}"

module load miniconda3
source activate "${CONDA_ENV:-rnaseq}"
cd "$BULK_RNA_WORK_DIR"
mkdir -p "ballgown/$SAMPLE_NAME" task_file
stringtie -e -B -p "$THREADS" -G merged_gtf/stringtie_merge.gtf \
  -o "ballgown/$SAMPLE_NAME/$SAMPLE_NAME.gtf" "bam/$SAMPLE_NAME.bam"
