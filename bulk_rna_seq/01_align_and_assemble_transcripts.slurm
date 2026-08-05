#!/usr/bin/env bash
#SBATCH --job-name=bulk_rna_align
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${BULK_RNA_WORK_DIR:?Set BULK_RNA_WORK_DIR}"
: "${BULK_RNA_INPUT_DIR:?Set BULK_RNA_INPUT_DIR}"
: "${EOSCZ_REFERENCE_DIR:?Set EOSCZ_REFERENCE_DIR}"
: "${SAMPLE_NAME:?Set SAMPLE_NAME to the FASTQ basename without _R1/_R2 suffixes}"
THREADS="${THREADS:-20}"

module load miniconda3
source activate "${CONDA_ENV:-rnaseq}"
cd "$BULK_RNA_WORK_DIR"
mkdir -p fastqc sam bam gtf merged_gtf ballgown task_file

hisat2 -p "$THREADS" --dta \
  -x "$EOSCZ_REFERENCE_DIR/GCF_009914755.1_T2T-CHM13v2.0_genomic" \
  -1 "$BULK_RNA_INPUT_DIR/${SAMPLE_NAME}_R1_001.fastq.gz" \
  # The source used the R1 file for both mates. This behavior is retained pending author review.
  -2 "$BULK_RNA_INPUT_DIR/${SAMPLE_NAME}_R1_001.fastq.gz" \
  -S "sam/${SAMPLE_NAME}.sam"
samtools sort -@ "$THREADS" -o "bam/${SAMPLE_NAME}.bam" "sam/${SAMPLE_NAME}.sam"
stringtie -p "$THREADS" -G "$EOSCZ_REFERENCE_DIR/T2T-CHM13V2.0.gff" \
  -o "gtf/${SAMPLE_NAME}.gtf" "bam/${SAMPLE_NAME}.bam"

echo "Finished sample: $SAMPLE_NAME"
