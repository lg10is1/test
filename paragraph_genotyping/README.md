# Paragraph genotyping

Multi-sample SV genotyping with Paragraph `multigrmpy.py`.

## Requirements

- Paragraph (Conda environment, default name `paragraph`)
- Slurm

## Usage

```bash
PARAGRAPH_MULTIGRMPY=/path/to/multigrmpy.py \
INPUT_VCF=/path/to/svs.vcf \
SAMPLE_MANIFEST=/path/to/manifest.tsv \
EOSCZ_REFERENCE_FASTA=/path/to/reference.fasta \
OUTPUT_DIR=/path/to/output \
sbatch genotype_sample.slurm
```

Optional variables: `THREADS` (28), `TMPDIR`, `CONDA_ENV` (`paragraph`).
