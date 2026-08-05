#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "${script_dir}/.." && pwd)

if [[ $# -gt 0 ]]; then
    echo "Usage: bash run_prs_scores.sh" >&2
    exit 1
fi

RESULT_ROOT="${RESULT_ROOT:-${repo_dir}/results/prs}"

: "${BFILE:?Set BFILE to the PLINK genotype prefix used for PRS scoring}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
PRS_THREADS="${PRS_THREADS:-40}"
PRS_MEMORY_MB="${PRS_MEMORY_MB:-100000}"
DEDUP_MODE="${DEDUP_MODE:-first}"
FORCE="${FORCE:-0}"

declare -A LOADING_DIRS=(
    [BIP]="${BIP_LOADING_DIR:?Set BIP_LOADING_DIR}"
    [SCZ]="${SCZ_LOADING_DIR:?Set SCZ_LOADING_DIR}"
    [ADHD]="${ADHD_LOADING_DIR:?Set ADHD_LOADING_DIR}"
    [ASD]="${ASD_LOADING_DIR:?Set ASD_LOADING_DIR}"
    [MDD]="${MDD_LOADING_DIR:?Set MDD_LOADING_DIR}"
)

traits=(BIP SCZ ADHD ASD MDD)

mkdir -p "${RESULT_ROOT}"

for ext in bed bim fam; do
    if [[ ! -s "${BFILE}.${ext}" ]]; then
        echo "ERROR: missing PLINK input file: ${BFILE}.${ext}" >&2
        exit 1
    fi
done

if ! command -v "${PLINK2_BIN}" >/dev/null 2>&1; then
    echo "ERROR: plink2 not found: ${PLINK2_BIN}" >&2
    exit 1
fi

dedup_prefix="${RESULT_ROOT}/EOSCZ_PROJECT_hg19_hapmap3_rename.dedup"
duplicate_ids="${RESULT_ROOT}/EOSCZ_PROJECT_hg19_hapmap3_rename.duplicate_ids.txt"
dedup_rmdup_list="${RESULT_ROOT}/EOSCZ_PROJECT_hg19_hapmap3_rename.dedup.rmdup.list"

if [[ "${DEDUP_MODE}" != "first" ]]; then
    echo "ERROR: unsupported DEDUP_MODE=${DEDUP_MODE}; currently only DEDUP_MODE=first is implemented." >&2
    exit 1
fi

awk '{count[$2]++} END {for (id in count) if (count[id] > 1) print id}' \
    "${BFILE}.bim" | sort > "${duplicate_ids}"

duplicate_count=$(wc -l < "${duplicate_ids}")
echo "[INFO] Duplicate variant IDs in input BIM: ${duplicate_count}"

if [[ "${duplicate_count}" -gt 0 ]]; then
    regenerate_dedup=0
    if [[ ! -s "${dedup_prefix}.bed" ||
          ! -s "${dedup_prefix}.bim" ||
          ! -s "${dedup_prefix}.fam" ||
          "${FORCE}" -eq 1 ]]; then
        regenerate_dedup=1
    else
        dedup_duplicate_count=$(awk '{count[$2]++} END {n=0; for (id in count) if (count[id] > 1) n++; print n}' "${dedup_prefix}.bim")
        if [[ "${dedup_duplicate_count}" -gt 0 ]]; then
            echo "[INFO] Existing deduplicated PLINK still has ${dedup_duplicate_count} duplicated IDs; regenerating."
            regenerate_dedup=1
        else
            echo "[SKIP] Existing deduplicated PLINK files: ${dedup_prefix}.{bed,bim,fam}"
        fi
    fi

    if [[ "${regenerate_dedup}" -eq 1 ]]; then
        echo "[INFO] Creating deduplicated PLINK files with plink2 --rm-dup force-first."
        "${PLINK2_BIN}" \
            --bfile "${BFILE}" \
            --rm-dup force-first list \
            --make-bed \
            --out "${dedup_prefix}" \
            --threads "${PRS_THREADS}" \
            --memory "${PRS_MEMORY_MB}"
    fi

    BFILE_SCORE="${dedup_prefix}"
else
    : > "${dedup_rmdup_list}"
    BFILE_SCORE="${BFILE}"
fi

for ext in bed bim fam; do
    if [[ ! -s "${BFILE_SCORE}.${ext}" ]]; then
        echo "ERROR: missing scoring PLINK file after deduplication: ${BFILE_SCORE}.${ext}" >&2
        exit 1
    fi
done

score_duplicate_count=$(awk '{count[$2]++} END {n=0; for (id in count) if (count[id] > 1) n++; print n}' "${BFILE_SCORE}.bim")
if [[ "${score_duplicate_count}" -gt 0 ]]; then
    echo "ERROR: scoring BIM still contains ${score_duplicate_count} duplicated variant IDs: ${BFILE_SCORE}.bim" >&2
    echo "       Remove stale dedup files or rerun with FORCE=1 after updating this script." >&2
    exit 1
fi

summary="${RESULT_ROOT}/prs_score_outputs.tsv"
printf 'trait\tloading\tout_prefix\tsscore\tbfile\n' > "${summary}"

echo "[INFO] PLINK input: ${BFILE}"
echo "[INFO] PLINK scoring input: ${BFILE_SCORE}"
echo "[INFO] Result root: ${RESULT_ROOT}"
echo "[INFO] Threads: ${PRS_THREADS}, memory MB: ${PRS_MEMORY_MB}"

for trait in "${traits[@]}"; do
    loading="${LOADING_DIRS[${trait}]}/loadings"
    out_dir="${RESULT_ROOT}/${trait}"
    out_prefix="${out_dir}/prscsx_eur"
    sscore="${out_prefix}.sscore"

    mkdir -p "${out_dir}"

    if [[ ! -s "${loading}" ]]; then
        echo "ERROR: missing loading file for ${trait}: ${loading}" >&2
        exit 1
    fi

    printf '%s\n' "${loading}" > "${out_dir}/source_loading.txt"
    printf '%s\t%s\t%s\t%s\t%s\n' "${trait}" "${loading}" "${out_prefix}" "${sscore}" "${BFILE_SCORE}" >> "${summary}"

    if [[ -s "${sscore}" && "${FORCE}" -eq 0 ]]; then
        echo "[SKIP] ${trait}: existing ${sscore}"
        continue
    fi

    echo "[INFO] Scoring ${trait}"
    "${PLINK2_BIN}" \
        --bfile "${BFILE_SCORE}" \
        --out "${out_prefix}" \
        --score "${loading}" 2 4 6 \
        --threads "${PRS_THREADS}" \
        --memory "${PRS_MEMORY_MB}"

    if [[ ! -s "${sscore}" ]]; then
        echo "ERROR: plink2 did not create ${sscore}" >&2
        exit 1
    fi

    echo "[DONE] ${trait}: ${sscore}"
done

echo
echo "[DONE] PRS scoring finished."
echo "Summary: ${summary}"
