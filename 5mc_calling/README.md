# 5mC calling and merge

## Requirements

- pb-CpG-tools 3.x
- Python 3 with `numpy` and `pandas`
- Slurm for the calling script

## 1. Call 5mC scores

Run one job per indexed aligned BAM:

```bash
sbatch 01_call_5mc.slurm sample.bam /path/to/calls/sample
```

The script uses the default pb-CpG-tools calling modes and writes:

```text
sample.combined.bed.gz
sample.combined.bed.gz.tbi
sample.combined.bw
```

The current pb-CpG-tools release includes the model internally, so no `--model`
argument is used. See the
[official pb-CpG-tools documentation](https://github.com/PacificBiosciences/pb-CpG-tools).

## 2. Merge samples by chromosome

```bash
python 02_merge_5mc.py \
    --input-dir /path/to/calls \
    --output-dir /path/to/merged \
    --chrom chr1
```

The script reads all `*.combined.bed.gz` files and creates position, methylation
score, coverage, and discretized-score matrices. Run it once per chromosome.

## 3. Filter CpG sites

```bash
python 03_filter_5mc.py \
    --input-dir /path/to/merged \
    --output-dir /path/to/qc \
    --chrom chr1
```

Sites missing in more than 30% of samples are removed. Change the threshold with
`--max-missing`.
