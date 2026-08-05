#!/usr/bin/env bash
set -euo pipefail

PANGENIE_BASE="${PANGENIE_BASE:-/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas}"
LD_OUTDIR="${LD_OUTDIR:-/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public}"
PLINK_BIN="${PLINK_BIN:-plink}"
MAF="${LD_MAF:-0.01}"
OUTDIR="${LD_OUTDIR}/maf01"
MAX_PLINK_JOBS=2

mkdir -p "$OUTDIR"

pids=()
labels=()
status=0

wait_batch() {
  local i
  for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
      echo "[DONE MAF] ${labels[$i]}"
    else
      echo "[ERROR MAF] ${labels[$i]}" >&2
      status=1
    fi
  done
  pids=()
  labels=()
}

for set_name in set00; do
  bfile="${PANGENIE_BASE}/${set_name}/NGS.QCsite.QCind"
  out="${OUTDIR}/${set_name}.maf0.01"
  echo "[LAUNCH MAF] ${set_name}: ${bfile}"
  (
    "$PLINK_BIN" \
      --bfile "$bfile" \
      --threads 1 \
      --autosome \
      --maf "$MAF" \
      --write-snplist \
      --out "$out"
  ) &
  pids+=("$!")
  labels+=("$set_name")
  if (( ${#pids[@]} >= MAX_PLINK_JOBS )); then
    wait_batch
  fi
done

wait_batch

if (( status != 0 )); then
  exit "$status"
fi

for set_name in set00; do
  file="${OUTDIR}/${set_name}.maf0.01.snplist"
  [[ -s "$file" ]] || { echo "[ERROR] Missing or empty MAF list: $file" >&2; exit 1; }
done

echo "[DONE] Pangenie MAF >= ${MAF} variant lists"
