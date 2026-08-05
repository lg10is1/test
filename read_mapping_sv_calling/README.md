# Read mapping and SV calling

PacBio read mapping and pbsv/Sniffles2 structural-variant calling as a Slurm array job.

## Requirements

- `pbmm2`, `pbsv`, `sniffles`, `samtools` (Conda environments, defaults `pacbio_tools` and `sniffles`)
- Reference FASTA, minimap2 index (`.mmi`), and TRF BED
- Slurm

## Configuration

```bash
cp config.example.sh config.sh
# Edit config.sh: TMP_DIR, RAW_BAM_DIR, ALIGN_DIR, FASTQ_DIR, PBSV_DIR,
# SNF_DIR, SNF_TMP, DISCOVER_DIR, REF_MMI, REF_FASTA, TRF_BED, THREADS,
# SAMPLE_LIST (one de-identified sample ID per line).
```

`config.sh` must remain untracked. The script also accepts `EOSCZ_SV_CONFIG=/path/to/config.sh`.

## Usage

```bash
sbatch --array=1-N call_structural_variants.sh
```

The array index selects the sample from `SAMPLE_LIST`. Per sample, the script converts BAM to FASTQ when needed, maps reads with pbmm2, runs pbsv discover/call, and runs Sniffles2.
