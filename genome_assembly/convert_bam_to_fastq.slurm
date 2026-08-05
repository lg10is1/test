#!/usr/bin/env bash
#SBATCH --job-name=bam_to_fastq
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${SOURCE_DIR:?Set SOURCE_DIR to the input BAM directory}"
DEST_DIR="${DEST_DIR:-${SOURCE_DIR}/fastq}"

module load miniconda3
source activate "${CONDA_ENV:-pacbio_tools}"
mkdir -p "$DEST_DIR" task_file

find "$SOURCE_DIR" -maxdepth 1 -type f -name '*.bam' -print0 | while IFS= read -r -d '' bam_file; do
  filename=$(basename "$bam_file" .bam)
  bam2fastq -o "$DEST_DIR/$filename" "$bam_file"
done

echo "Conversion complete. FASTQ files were written to $DEST_DIR."
