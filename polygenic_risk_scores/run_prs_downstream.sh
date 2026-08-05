#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ $# -gt 0 ]]; then
  echo "Usage: bash run_prs_downstream.sh" >&2
  exit 1
fi

RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
RESULT_ROOT="${RESULT_ROOT:-${REPO_DIR}/results/prs}"

export RESULT_ROOT
export N_PCS="${N_PCS:-20}"
export N_QUANTILE="${N_QUANTILE:-5}"
export RARE_THRESHOLDS="${RARE_THRESHOLDS:-1,2,3}"

echo "[INFO] Script dir: ${SCRIPT_DIR}"
echo "[INFO] Result root: ${RESULT_ROOT}"
echo "[INFO] N_PCS=${N_PCS}; N_QUANTILE=${N_QUANTILE}; RARE_THRESHOLDS=${RARE_THRESHOLDS}"
echo "[INFO] Step 03 rare burden is not run here; existing data/rare_variants_clean_burden.tsv is reused."
echo "[INFO] LRS exclude samples: ${LRS_EXCLUDE_FILE:-NOT_SET}"
echo "[INFO] Version5 case/control violin annotations use logistic-regression P values from glm(y ~ PRS)."
echo "[INFO] Version5 quantile plotting keeps SCZ for the main figure and BIP/ADHD/MDD/ASD for supplementary figure."

mkdir -p "${RESULT_ROOT}/figure"

echo "[RUN] LRS previous PRS with 20 PCs and sample exclusion"
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/02_lrs_prs_pc20.R"

echo "[RUN] NGS previous PRS with DeepVariant 20 PCs + batch"
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/02_ngs_prs_deepvar_pc20_batch.R"

echo "[RUN] LRS SCZ PRS rare-threshold analysis with sample exclusion"
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/04_lrs_prs_rare_threshold_pc20.R"

echo "[RUN] Combined Nagelkerke + PRS violin panel"
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/05_plot_combined_prs_panels.R"

echo "[RUN] SCZ main-text quantile figure and non-SCZ supplementary quantile figure"
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/06_plot_quantile_main_supplement.R"

echo "[DONE] Downstream PRS version5 analyses finished."
echo "[DONE] Figures: ${RESULT_ROOT}/figure"
