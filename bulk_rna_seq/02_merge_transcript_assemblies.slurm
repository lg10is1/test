#!/usr/bin/env bash
#SBATCH --job-name=bulk_rna_merge
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${BULK_RNA_WORK_DIR:?Set BULK_RNA_WORK_DIR}"
: "${EOSCZ_REFERENCE_DIR:?Set EOSCZ_REFERENCE_DIR}"
THREADS="${THREADS:-20}"

module load miniconda3
source activate "${CONDA_ENV:-rnaseq}"
cd "$BULK_RNA_WORK_DIR"
mkdir -p merged_gtf task_file
find gtf -maxdepth 1 -type f -name '*.gtf' -print | sort > merged_gtf/stringtie_mergelist.txt
[[ -s merged_gtf/stringtie_mergelist.txt ]] || { echo "No GTF files found in $BULK_RNA_WORK_DIR/gtf" >&2; exit 1; }
stringtie --merge -p "$THREADS" -G "$EOSCZ_REFERENCE_DIR/T2T-CHM13V2.0.gff" \
  -o merged_gtf/stringtie_merge.gtf merged_gtf/stringtie_mergelist.txt
