#!/usr/bin/env bash
#SBATCH --job-name=pacbio_sv
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG_FILE="${EOSCZ_SV_CONFIG:-${SCRIPT_DIR}/config.sh}"
[[ -f "$CONFIG_FILE" ]] || { echo "Configuration not found: $CONFIG_FILE" >&2; exit 2; }
# shellcheck source=/dev/null
source "$CONFIG_FILE"

for name in TMP_DIR RAW_BAM_DIR ALIGN_DIR FASTQ_DIR PBSV_DIR SNF_DIR SNF_TMP DISCOVER_DIR REF_MMI REF_FASTA TRF_BED THREADS; do
  [[ -n "${!name:-}" ]] || { echo "Missing required configuration variable: $name" >&2; exit 2; }
done
SAMPLE_LIST="${SAMPLE_LIST:-${SCRIPT_DIR}/samples.list}"
[[ -f "$SAMPLE_LIST" ]] || { echo "Sample list not found: $SAMPLE_LIST" >&2; exit 2; }
: "${SLURM_ARRAY_TASK_ID:?Run this script as a Slurm array job}"
export TMPDIR="$TMP_DIR"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
[[ -n "$SAMPLE" ]] || { echo "No sample for array index $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
[[ "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe sample identifier" >&2; exit 2; }

RAW_SAMPLE_DIR="$RAW_BAM_DIR/$SAMPLE"
ALIGN_BAM="$ALIGN_DIR/$SAMPLE.align.hg38.bam"
FASTQ="$FASTQ_DIR/$SAMPLE.fastq.gz"
SNF_VCF="$SNF_DIR/$SAMPLE.hg38.snf.vcf"
PBSV_VCF="$PBSV_DIR/$SAMPLE.pbsv.hg38.vcf"

[[ ! -f "$SNF_VCF" ]] || { echo "Sample already finished: $SAMPLE"; exit 0; }
mkdir -p "$ALIGN_DIR" "$FASTQ_DIR" "$PBSV_DIR" "$SNF_DIR" "$SNF_TMP" "$DISCOVER_DIR" logs
module load miniconda3
source activate "${PACBIO_CONDA_ENV:-pacbio_tools}"

if [[ ! -f "$FASTQ" ]]; then
  bam2fastq "$RAW_SAMPLE_DIR/$SAMPLE.fofn" -o "$FASTQ_DIR/$SAMPLE" -j "$THREADS"
fi
if [[ ! -f "$ALIGN_BAM" ]]; then
  pbmm2 align --preset CCS --sort -j "$THREADS" -J 4 -m 1G --sample "$SAMPLE" \
    "$REF_MMI" "$RAW_SAMPLE_DIR/$SAMPLE.fofn" "$ALIGN_BAM"
fi

while IFS= read -r chromosome; do
  [[ -n "$chromosome" ]] || continue
  out="$DISCOVER_DIR/$SAMPLE.$chromosome.svsig.gz"
  [[ -f "$out" ]] || pbsv discover --region "$chromosome" --tandem-repeats "$TRF_BED" "$ALIGN_BAM" "$out"
done < <(samtools view -H "$ALIGN_BAM" | awk -F'\t' '$1=="@SQ"{split($2,a,":");print a[2]}')

if [[ ! -f "$PBSV_VCF" ]]; then
  mapfile -t signatures < <(find "$DISCOVER_DIR" -maxdepth 1 -type f -name "$SAMPLE.*.svsig.gz" -print | sort)
  (( ${#signatures[@]} > 0 )) || { echo "No pbsv signature files found" >&2; exit 1; }
  pbsv call -j "$THREADS" --ccs --min-sv-length 50 "$REF_FASTA" "${signatures[@]}" "$PBSV_VCF"
fi

conda activate "${SNIFFLES_CONDA_ENV:-sniffles}"
if [[ ! -f "$SNF_VCF" ]]; then
  sniffles --input "$ALIGN_BAM" --reference "$REF_FASTA" --tandem-repeats "$TRF_BED" \
    --threads "$THREADS" --snf "$SNF_TMP/$SAMPLE.snf"
  sniffles --input "$ALIGN_BAM" --reference "$REF_FASTA" --tandem-repeats "$TRF_BED" \
    --threads "$THREADS" --vcf "$SNF_VCF"
fi
echo "Finished sample: $SAMPLE"
