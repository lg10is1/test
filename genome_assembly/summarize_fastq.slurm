#!/usr/bin/env bash
#SBATCH --job-name=fastq_summary
#SBATCH -n 4
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${FASTQ_DIR:?Set FASTQ_DIR}"
OUTPUT_FILE="${FASTQ_STATS_OUTPUT:-fastq_stats.txt}"
THREADS="${THREADS:-4}"
mkdir -p task_file
seqkit stats --threads "$THREADS" "$FASTQ_DIR"/*.fastq.gz > "$OUTPUT_FILE"
