# Reference-population analysis

VCF utilities, single-sample inversion clustering, public 1000 Genomes Project (IGSR) ONT manifest handling, and a Slurm-array template for remapping CRAM reads to CHM13. This module is a set of utilities, not an end-to-end workflow.

## Contents

| Path | Purpose |
|---|---|
| `filter_vcf_by_positions.py` | Streams records matching exact CHROM/POS pairs to a new VCF. |
| `list_vcf_samples.sh` | Extracts sample names with `bcftools query -l`. |
| `combine_large_inversions.py` | Clusters single-sample inversion calls and calculates case/control carrier statistics. |
| `1kgp/ont/download_igsr_ont.sh` | Validates and optionally downloads the public IGSR ONT manifest. |
| `1kgp/ont/file_list.txt` | Public IGSR HTTPS URLs (2,038 entries). |
| `1kgp/ont/md5_checksums.txt` | Expected IGSR MD5 checksums matched to the URL basenames. |
| `1kgp/ngs/remap_cram_to_chm13.sh` | Resumable Slurm-array CRAM-to-FASTQ-to-CHM13 remapping template. |

## Exact-site VCF filtering

The positions file is tab-separated with `CHROM` and one-based `POS` in the first two columns. Blank lines and rows beginning with `#` are ignored.

```bash
python3 reference_population_analysis/filter_vcf_by_positions.py \
  --position-file /path/to/positions.tsv \
  --vcf-file /path/to/input.vcf.gz \
  --output-vcf /path/to/filtered.vcf.gz
```

## Inversion clustering and carrier comparison

Input VCFs must be named `<sample>_filtered.vcf[.gz]` with one sample genotype column. Case and control manifests contain one de-identified sample ID per line and must not overlap.

```bash
python3 reference_population_analysis/combine_large_inversions.py \
  --input-dir /path/to/single_sample_vcfs \
  --case-file /path/to/cases.list \
  --control-file /path/to/controls.list \
  --dist 1000 \
  --overlap 0.5 \
  --output /path/to/INV_Stats_Full.txt
```

Retained assumptions: clusters require both breakpoint distances within `--dist` and reciprocal overlap of at least `--overlap`; absent samples are treated as `0/0`; only ALT allele `1` contributes to dosage; the p-value is a one-sided Fisher exact test for case enrichment with Benjamini-Hochberg FDR; a Haldane-Anscombe correction of 0.5 is used for the odds ratio.

## Public IGSR ONT manifest

The downloader validates locally and prints its plan by default (no network request):

```bash
bash reference_population_analysis/1kgp/ont/download_igsr_ont.sh \
  --workdir /path/to/approved/large_storage
```

Add `--execute` to run the downloads. Only URLs on the listed IGSR HTTPS host are accepted, and a downloaded file is removed only when its checksum fails. Downloaded files, logs, and state files must remain untracked.

## NGS remapping template

Requires a Slurm array and the environment variables `SOURCE_REFERENCE`, `TARGET_REFERENCE`, `INPUT_LIST`, and `OUT_DIR`. The input list contains one authorized local CRAM path per line. The array index is one-based; set `OFFSET=1` for a zero-based array.

```bash
SOURCE_REFERENCE=/path/to/GRCh38.fa \
TARGET_REFERENCE=/path/to/CHM13v2.0.fa \
INPUT_LIST=/path/to/deidentified_crams.list \
OUT_DIR=/path/to/output \
sbatch --array=1-N --export=ALL \
  reference_population_analysis/1kgp/ngs/remap_cram_to_chm13.sh
```

Requires `samtools`, `bwa-mem2`, and `gzip`. This is a template; confirm reference checksums, CRAM reference compatibility, and read-group policy before submission.
