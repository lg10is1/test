#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/EOSCZ_PROJECT"
LD_OUTDIR="${LD_OUTDIR:-${ROOT}/figure_analysis/SV_SNV_LD/LD_decay_public}"
MAF="${DEEPVARIANT_PARAGRAPH_NULL_MAF:-0.02}"
PLINK_BIN="${PLINK_BIN:-plink}"
OUTDIR="${LD_OUTDIR}/maf02"
BFILE="${ROOT}/GWAS/Deepvariant_paragraph/chr_all2.strict_step2_genimi.common_samples.merged"
OUT_PREFIX="${OUTDIR}/deepvariant_paragraph.maf0.02"

mkdir -p "$OUTDIR"

if [[ -s "${OUT_PREFIX}.snplist" && "${OUT_PREFIX}.snplist" -nt "${BFILE}.bed" ]]; then
  echo "[SKIP MAF] deepvariant_paragraph MAF >= ${MAF}: ${OUT_PREFIX}.snplist"
else
  echo "[RUN MAF] deepvariant_paragraph MAF >= ${MAF}: ${BFILE}"
  "$PLINK_BIN" \
    --bfile "$BFILE" \
    --threads 1 \
    --autosome \
    --maf "$MAF" \
    --write-snplist \
    --out "$OUT_PREFIX"
fi

[[ -s "${OUT_PREFIX}.snplist" ]] || {
  echo "[ERROR] Missing or empty MAF list: ${OUT_PREFIX}.snplist" >&2
  exit 1
}

echo "[DONE] deepvariant_paragraph MAF >= ${MAF} variant list: ${OUT_PREFIX}.snplist"
