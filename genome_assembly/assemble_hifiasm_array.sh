#!/usr/bin/env bash
set -euo pipefail

: "${FASTQ_DIR:?Set FASTQ_DIR}"
: "${OUT_ROOT:?Set OUT_ROOT}"
: "${SLURM_PARTITION:?Set SLURM_PARTITION for your cluster}"
CONDA_ENV="${CONDA_ENV:-hifiasm}"
GFATOOLS_ENV="${GFATOOLS_ENV:-gfatools}"
THREADS="${THREADS:-40}"

usage() {
  echo "Usage: $0 [-n NUM | -r START-END] [-l]"
  echo "  -n NUM        Submit one sample by sorted position."
  echo "  -r START-END  Submit an inclusive range of sorted positions."
  echo "  -l            List samples without submitting jobs."
  echo "  -h            Show this help."
}

LIST_ONLY=false
RANGE=""
while getopts ":n:r:lh" opt; do
  case "$opt" in
    n) RANGE="$OPTARG-$OPTARG" ;;
    r) RANGE="$OPTARG" ;;
    l) LIST_ONLY=true ;;
    h) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
done

mapfile -t FQS < <(find "$FASTQ_DIR" -type f -name '*.fastq.gz' -print | sort)
TOTAL=${#FQS[@]}
(( TOTAL > 0 )) || { echo "No *.fastq.gz files found under $FASTQ_DIR" >&2; exit 1; }

if [[ -n "$RANGE" ]]; then
  [[ "$RANGE" =~ ^[0-9]+-[0-9]+$ ]] || { echo "Invalid range: $RANGE" >&2; exit 2; }
  START=${RANGE%-*}
  END=${RANGE#*-}
  (( START > 0 && END <= TOTAL && START <= END )) || { echo "Range must be within 1-$TOTAL" >&2; exit 2; }
else
  START=1
  END=$TOTAL
fi

if "$LIST_ONLY"; then
  for ((i=1; i<=TOTAL; i++)); do
    fq="${FQS[$((i-1))]}"
    printf '%s  %s\n' "$i" "$(basename "$fq" .fastq.gz)"
  done
  exit 0
fi

for ((i=START; i<=END; i++)); do
  fq="${FQS[$((i-1))]}"
  sample=$(basename "$fq" .fastq.gz)
  outdir="$OUT_ROOT/$sample"
  flagfa="$outdir/$sample.asm.bp.hap2.p_ctg.gfa"
  [[ ! -f "$flagfa" ]] || { echo "[SKIP] $i/$TOTAL $sample"; continue; }
  mkdir -p "$outdir"
  jobscript=$(mktemp)
  trap 'rm -f -- "$jobscript"' EXIT
  cat > "$jobscript" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_${sample}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH -n ${THREADS}
#SBATCH --ntasks-per-node=${THREADS}
#SBATCH --output=${outdir}/%x_%j.out
#SBATCH --error=${outdir}/%x_%j.err
set -euo pipefail
module purge
module load miniconda3
source activate "${CONDA_ENV}"
hifiasm -o "${outdir}/${sample}.asm" -t "${THREADS}" "${fq}"
source activate "${GFATOOLS_ENV}"
gfatools gfa2fa "${outdir}/${sample}.asm.bp.hap1.p_ctg.gfa" > "${outdir}/${sample}.1.fa"
gfatools gfa2fa "${outdir}/${sample}.asm.bp.hap2.p_ctg.gfa" > "${outdir}/${sample}.2.fa"
EOF
  sbatch "$jobscript"
  echo "[SUBMITTED] $i/$TOTAL $sample"
  rm -f -- "$jobscript"
  trap - EXIT
done
