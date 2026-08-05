# Bulk RNA-seq

HISAT2 alignment, StringTie transcript assembly and quantification, and gene/transcript count matrices.

## Requirements

- HISAT2, samtools, StringTie (Conda environment, default name `rnaseq`)
- Python 3
- Slurm

## Workflow

```bash
# 1. Align reads and assemble transcripts per sample
BULK_RNA_WORK_DIR=/path/to/work \
BULK_RNA_INPUT_DIR=/path/to/fastq \
EOSCZ_REFERENCE_DIR=/path/to/hisat2_and_annotation \
SAMPLE_NAME=<sample> \
sbatch 01_align_and_assemble_transcripts.slurm

# 2. Merge per-sample assemblies
BULK_RNA_WORK_DIR=/path/to/work \
sbatch 02_merge_transcript_assemblies.slurm

# 3. Quantify against the merged assembly
BULK_RNA_WORK_DIR=/path/to/work \
SAMPLE_NAME=<sample> \
sbatch 03_quantify_transcripts.slurm

# 4. Build count matrices
BULK_RNA_WORK_DIR=/path/to/work \
sbatch 04_prepare_count_matrices.slurm
```

`prepare_count_matrices.py` is the StringTie `prepDE.py` utility used by step 4 and can also be run directly:

```bash
python3 prepare_count_matrices.py \
  -i ../examples/synthetic/gtf_files.tsv \
  -g gene_count_matrix.csv \
  -t transcript_count_matrix.csv
```

Optional variables: `THREADS` (20), `CONDA_ENV` (`rnaseq`).
