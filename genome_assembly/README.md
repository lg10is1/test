# Genome assembly

BAM-to-FASTQ conversion, FASTQ statistics, and hifiasm assembly.

## Requirements

- `bam2fastq`, `seqkit`, `hifiasm`, `gfatools` (Conda environments)
- Slurm

## Scripts

| Script | Purpose |
|---|---|
| `convert_bam_to_fastq.slurm` | Converts all BAM files in `SOURCE_DIR` to FASTQ with bam2fastq. |
| `summarize_fastq.slurm` | Writes seqkit statistics for all FASTQ files in `FASTQ_DIR`. |
| `assemble_hifiasm_array.sh` | Submits hifiasm assembly jobs for samples selected by sorted position (`-n`, `-r`, or `-l` to list). |

## Usage

```bash
SOURCE_DIR=/path/to/bam sbatch convert_bam_to_fastq.slurm

FASTQ_DIR=/path/to/fastq sbatch summarize_fastq.slurm

FASTQ_DIR=/path/to/fastq OUT_ROOT=/path/to/assemblies SLURM_PARTITION=<partition> \
bash assemble_hifiasm_array.sh -l
```

Optional variables: `CONDA_ENV` (`hifiasm`), `GFATOOLS_ENV` (`gfatools`), `THREADS` (40), `DEST_DIR`, `FASTQ_STATS_OUTPUT`.
