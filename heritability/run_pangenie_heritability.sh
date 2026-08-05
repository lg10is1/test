#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PROJECT_BASE="${PROJECT_BASE:-${REPO_DIR}/data}"
PANGENIE_BASE="${PANGENIE_BASE:-${PROJECT_BASE}/TGS_callset/Pangenie_v3}"
GWAS_BASE="${GWAS_BASE:-${PANGENIE_BASE}/06.gwas}"
FIGURE7_BASE="${FIGURE7_BASE:-${PROJECT_BASE}/Figure7_heritability}"
SOURCE_NAME="${SOURCE_NAME:-pangenie}"
MERGE_SOURCES="${MERGE_SOURCES:-set00,set01,set02}"
OUT_BASE="${OUT_BASE:-${REPO_DIR}/results/heritability/pangenie}"

TR_GRM_PREFIX="${TR_GRM_PREFIX:-${FIGURE7_BASE}/tr_grm/ngstr}"
REMOVE_FILE="${REMOVE_FILE:-${FIGURE7_BASE}/tr_grm/remove_sample}"
PHENO_FILE="${PHENO_FILE:-${GWAS_BASE}/SCZ_pheno.txt}"
COVAR_FILE="${COVAR_FILE:-${GWAS_BASE}/SCZ_batch.txt}"

GCTA_BIN="${GCTA_BIN:-gcta64}"
PLINK_BIN="${PLINK_BIN:-plink}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
THREADS="${THREADS:-16}"
GCTA_THREADS="${GCTA_THREADS:-$THREADS}"
PLINK_THREADS="${PLINK_THREADS:-$THREADS}"
REML_JOBS="${REML_JOBS:-1}"
MEMORY_MB="${MEMORY_MB:-64000}"
PREVALENCE="0.007"
MAF="${MAF:-0.02}"
REML_MAXIT="${REML_MAXIT:-1000}"
PCA_COUNT="${PCA_COUNT:-20}"
KING_CUTOFF="${KING_CUTOFF:-0.0884}"
PRUNE_WINDOW="${PRUNE_WINDOW:-50}"
PRUNE_STEP="${PRUNE_STEP:-5}"
PRUNE_R2="${PRUNE_R2:-0.2}"
SMALL_VARIANT_MODE="${SMALL_VARIANT_MODE:-SNV_INDEL_LT50}"
MODE="${MODE:-all}"
FORCE="${FORCE:-0}"
DRY_RUN="${DRY_RUN:-0}"

usage() {
  cat <<'EOF'
Usage:
  bash run_pangenie_heritability.sh [options]

Modes:
  --mode all, prepare, pca, plot, grm, inputs, reml, summarise

Common options:
  --out-base DIR
  --extract-base DIR       Existing source extract root, default Figure7_heritability/res_remove_sample
  --gcta-threads N
  --plink-threads N
  --reml-jobs N
  --memory-mb N
  --force

Example:
  bash run_pangenie_heritability.sh --mode all \
    --plink-threads 8 --gcta-threads 2 --reml-jobs 20 --memory-mb 24000
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --out-base) OUT_BASE="$2"; shift 2 ;;
    --project-base) PROJECT_BASE="$2"; PANGENIE_BASE="${PROJECT_BASE}/TGS_callset/Pangenie_v3"; GWAS_BASE="${PANGENIE_BASE}/06.gwas"; FIGURE7_BASE="${PROJECT_BASE}/Figure7_heritability"; shift 2 ;;
    --gwas-base) GWAS_BASE="$2"; shift 2 ;;
    --source-name) SOURCE_NAME="$2"; shift 2 ;;
    --merge-sources) MERGE_SOURCES="$2"; shift 2 ;;
    --tr-grm-prefix) TR_GRM_PREFIX="$2"; shift 2 ;;
    --remove-samples) REMOVE_FILE="$2"; shift 2 ;;
    --pheno) PHENO_FILE="$2"; shift 2 ;;
    --covar) COVAR_FILE="$2"; shift 2 ;;
    --maf) MAF="$2"; shift 2 ;;
    --reml-maxit) REML_MAXIT="$2"; shift 2 ;;
    --pca-count) PCA_COUNT="$2"; shift 2 ;;
    --king-cutoff) KING_CUTOFF="$2"; shift 2 ;;
    --small-variant-mode) SMALL_VARIANT_MODE="$2"; shift 2 ;;
    --threads) THREADS="$2"; GCTA_THREADS="$2"; PLINK_THREADS="$2"; shift 2 ;;
    --gcta-threads) GCTA_THREADS="$2"; shift 2 ;;
    --plink-threads) PLINK_THREADS="$2"; shift 2 ;;
    --reml-jobs) REML_JOBS="$2"; shift 2 ;;
    --memory-mb) MEMORY_MB="$2"; shift 2 ;;
    --gcta) GCTA_BIN="$2"; shift 2 ;;
    --plink) PLINK_BIN="$2"; shift 2 ;;
    --plink2) PLINK2_BIN="$2"; shift 2 ;;
    --rscript) RSCRIPT_BIN="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$MODE" in
  all|prepare|pca|plot|grm|inputs|reml|summarise) ;;
  *) echo "[ERROR] Invalid --mode: ${MODE}" >&2; exit 2 ;;
esac
if ! [[ "$REML_JOBS" =~ ^[0-9]+$ ]] || (( REML_JOBS < 1 )); then
  echo "[ERROR] --reml-jobs must be an integer >=1" >&2
  exit 2
fi
if ! [[ "$PCA_COUNT" =~ ^[0-9]+$ ]] || (( PCA_COUNT < 20 )); then
  echo "[ERROR] --pca-count must be an integer >=20" >&2
  exit 2
fi

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log_info() { echo "[$(timestamp)] $*"; }

print_cmd() { printf '  '; printf '%q ' "$@"; printf '\n'; }

require_file() {
  local file="$1"
  if [[ "$DRY_RUN" == "1" ]]; then return 0; fi
  [[ -s "$file" ]] || { echo "[ERROR] Missing/empty file: $file" >&2; exit 1; }
}

require_bfile() {
  local prefix="$1"
  require_file "${prefix}.bed"
  require_file "${prefix}.bim"
  require_file "${prefix}.fam"
}

grm_complete() {
  local prefix="$1"
  [[ -s "${prefix}.grm.bin" && -s "${prefix}.grm.N.bin" && -s "${prefix}.grm.id" ]]
}

require_grm() {
  local prefix="$1"
  if [[ "$DRY_RUN" == "1" ]]; then return 0; fi
  grm_complete "$prefix" || { echo "[ERROR] Incomplete GRM: $prefix" >&2; exit 1; }
}

run_logged() {
  local log_file="$1"; shift
  if [[ "$DRY_RUN" == "1" ]]; then print_cmd "$@"; return 0; fi
  mkdir -p "$(dirname "$log_file")"
  { printf '[CMD] '; printf '%q ' "$@"; printf '\n'; } > "$log_file"
  "$@" >> "$log_file" 2>&1
}

try_logged() {
  local log_file="$1"; shift
  if [[ "$DRY_RUN" == "1" ]]; then print_cmd "$@"; return 0; fi
  mkdir -p "$(dirname "$log_file")"
  { printf '[CMD] '; printf '%q ' "$@"; printf '\n'; } > "$log_file"
  set +e
  "$@" >> "$log_file" 2>&1
  local rc=$?
  set -e
  return "$rc"
}

run_rscript() {
  if [[ "$DRY_RUN" == "1" ]]; then print_cmd "$RSCRIPT_BIN" "$@"; else "$RSCRIPT_BIN" "$@"; fi
}

maf_args=()
case "$MAF" in
  ""|none|NONE|NA|na) maf_args=() ;;
  *) maf_args=(--maf "$MAF") ;;
esac

mkdir -p "$OUT_BASE" "$OUT_BASE/shared" "$OUT_BASE/plots" "$OUT_BASE/summary"
PANGENIE_BFILE_DIR="${OUT_BASE}/pangenie_bfiles"
FILTER_REMOVE_FILE="${OUT_BASE}/shared/remove_sample.matched.tsv"
FILTER_TR_KEEP="${OUT_BASE}/shared/TR.filtered.keep"
remove_args=()
IFS=',' read -r -a MERGE_SOURCE_ARRAY <<< "$MERGE_SOURCES"

prepare_global_remove() {
  require_file "$REMOVE_FILE"
  require_file "${TR_GRM_PREFIX}.grm.id"
  : > "$FILTER_REMOVE_FILE"
  awk -v matched="$FILTER_REMOVE_FILE" '
    BEGIN { FS="[ ,\t]+"; OFS="\t" }
    NR==FNR {
      gsub(/\r/, "")
      if ($0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/) next
      u1=toupper($1); u2=toupper($2)
      if (u1=="FID" || u1=="IID" || u2=="IID") next
      if (NF>=2) pair[$1 SUBSEP $2]=1; else iid[$1]=1
      next
    }
    {
      if (($1 SUBSEP $2) in pair || $2 in iid) print $1, $2 > matched
      else print $1, $2
    }
  ' "$REMOVE_FILE" "${TR_GRM_PREFIX}.grm.id" > "$FILTER_TR_KEEP"
  local removed_n keep_n
  removed_n="$(awk 'NF>=2 {n++} END {print n+0}' "$FILTER_REMOVE_FILE")"
  keep_n="$(awk 'NF>=2 {n++} END {print n+0}' "$FILTER_TR_KEEP")"
  if (( removed_n == 0 )); then
    echo "[ERROR] No remove_sample entry matched ${TR_GRM_PREFIX}.grm.id" >&2
    exit 1
  fi
  log_info "Global sample exclusion: removed=${removed_n}, TR keep=${keep_n}, list=${REMOVE_FILE}"
  remove_args=(--remove "$FILTER_REMOVE_FILE")
}

if [[ "$MODE" == "all" || "$MODE" == "pca" || "$MODE" == "inputs" || "$MODE" == "reml" ]]; then
  prepare_global_remove
fi

cat <<EOF
[CONFIG] SCRIPT_DIR=${SCRIPT_DIR}
[CONFIG] OUT_BASE=${OUT_BASE}
[CONFIG] SOURCE_NAME=${SOURCE_NAME}
[CONFIG] MERGE_SOURCES=${MERGE_SOURCES}
[CONFIG] MODE=${MODE}
[CONFIG] TR_GRM_PREFIX=${TR_GRM_PREFIX}
[CONFIG] REMOVE_FILE=${REMOVE_FILE}
[CONFIG] MAF=${MAF} PREVALENCE=${PREVALENCE} REML_MAXIT=${REML_MAXIT}
[CONFIG] PCA_COUNT=${PCA_COUNT} KING_CUTOFF=${KING_CUTOFF}
[CONFIG] GCTA_THREADS=${GCTA_THREADS} PLINK_THREADS=${PLINK_THREADS} REML_JOBS=${REML_JOBS}
[CONFIG] DRY_RUN=${DRY_RUN} FORCE=${FORCE}
EOF

prepare_pangenie_bfiles() {
  log_info "Preparing Pangenie SNV/SV bfiles with source-suffixed variant IDs"
  mkdir -p "$PANGENIE_BFILE_DIR/logs"
  local source bfile
  for source in "${MERGE_SOURCE_ARRAY[@]}"; do
    bfile="${GWAS_BASE}/${source}/NGS.QCsite.QCind"
    require_bfile "$bfile"
  done

  run_rscript "${SCRIPT_DIR}/prepare_pangenie_bfiles.R" \
    --sources "$MERGE_SOURCES" --gwas-base "$GWAS_BASE" \
    --out-dir "$PANGENIE_BFILE_DIR" --small-variant-mode "$SMALL_VARIANT_MODE"

  local component extract rename raw_prefix merge_prefix list_file first_prefix n_prefixes
  for component in SNV_INDEL SV; do
    merge_prefix="${PANGENIE_BFILE_DIR}/${SOURCE_NAME}.${component}"
    if [[ "$FORCE" != "1" && -s "${merge_prefix}.bed" && -s "${merge_prefix}.bim" && -s "${merge_prefix}.fam" ]]; then
      log_info "Reusing merged ${component} bfile: ${merge_prefix}"
      continue
    fi
    list_file="${PANGENIE_BFILE_DIR}/${component}/${component}.merge_list.txt"
    mkdir -p "${PANGENIE_BFILE_DIR}/${component}"
    : > "$list_file"
    first_prefix=""
    n_prefixes=0
    for source in "${MERGE_SOURCE_ARRAY[@]}"; do
      bfile="${GWAS_BASE}/${source}/NGS.QCsite.QCind"
      extract="${PANGENIE_BFILE_DIR}/${component}/${source}.${component}.old_ids.extract"
      rename="${PANGENIE_BFILE_DIR}/${component}/${source}.${component}.rename.tsv"
      require_file "$extract"
      require_file "$rename"
      raw_prefix="${PANGENIE_BFILE_DIR}/${component}/${source}.${component}.renamed"
      if [[ "$FORCE" == "1" || ! -s "${raw_prefix}.bed" || ! -s "${raw_prefix}.bim" || ! -s "${raw_prefix}.fam" ]]; then
        run_logged "${PANGENIE_BFILE_DIR}/logs/${source}.${component}.make_filtered_bfile.log" \
          "$PLINK2_BIN" --bfile "$bfile" --extract "$extract" --make-bed \
          --threads "$PLINK_THREADS" --memory "$MEMORY_MB" --out "$raw_prefix"
        if [[ "$DRY_RUN" != "1" ]]; then
          awk 'BEGIN{FS=OFS="\t"} NR==FNR {rename[$1]=$2; next} \
            {if ($2 in rename) $2=rename[$2]; print}' \
            "$rename" "${raw_prefix}.bim" > "${raw_prefix}.bim.tmp"
          mv "${raw_prefix}.bim.tmp" "${raw_prefix}.bim"
        fi
      fi
      if [[ -z "$first_prefix" ]]; then
        first_prefix="$raw_prefix"
      else
        printf '%s\t%s\t%s\n' "${raw_prefix}.bed" "${raw_prefix}.bim" "${raw_prefix}.fam" >> "$list_file"
      fi
      n_prefixes=$((n_prefixes + 1))
    done
    if (( n_prefixes == 0 )); then
      echo "[ERROR] No ${component} variants selected for merged bfile" >&2
      exit 1
    elif (( n_prefixes == 1 )); then
      cp "${first_prefix}.bed" "${merge_prefix}.bed"
      cp "${first_prefix}.bim" "${merge_prefix}.bim"
      cp "${first_prefix}.fam" "${merge_prefix}.fam"
    else
      run_logged "${PANGENIE_BFILE_DIR}/logs/${component}.plink_merge.log" \
        "$PLINK_BIN" --bfile "$first_prefix" --merge-list "$list_file" --make-bed --out "$merge_prefix"
    fi
    require_bfile "$merge_prefix"
    awk '{print $2}' "${merge_prefix}.bim" > "${merge_prefix}.all_variants.extract"
  done
}

build_projected_pca() {
  local component="$1" bfile="$2" extract="$3" reference_keep="$4" pca_method="$5"
  local pc_dir="${OUT_BASE}/${SOURCE_NAME}/pc"
  local log_dir="${OUT_BASE}/${SOURCE_NAME}/logs"
  local stem="${pc_dir}/${SOURCE_NAME}.${component}"
  local prune="${stem}.pruned"
  local ref="${stem}.unrelated_pca"
  local proj="${stem}.all.projected_pca"
  local method_file="${stem}.pca_method.txt"
  local recompute="$FORCE"
  local previous_method=""
  mkdir -p "$pc_dir" "$log_dir"
  if [[ "$DRY_RUN" != "1" && -s "$method_file" ]]; then
    previous_method="$(head -n 1 "$method_file")"
    if [[ "$previous_method" != "$pca_method" ]]; then recompute=1; fi
  fi
  if [[ "$recompute" == "1" || ! -s "${prune}.prune.in" ]]; then
    run_logged "${log_dir}/${SOURCE_NAME}.${component}.prune.log" \
      "$PLINK2_BIN" --bfile "$bfile" --extract "$extract" "${maf_args[@]}" "${remove_args[@]}" \
      --autosome --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_R2" \
      --threads "$PLINK_THREADS" --memory "$MEMORY_MB" --out "$prune"
  fi
  require_file "${prune}.prune.in"

  local keep_args=()
  if [[ -n "$reference_keep" ]]; then keep_args=(--keep "$reference_keep"); fi
  if [[ "$recompute" == "1" || ! -s "${ref}.eigenvec.allele" || ! -s "${ref}.acount" ]]; then
    if ! try_logged "${log_dir}/${SOURCE_NAME}.${component}.reference_pca.log" \
      "$PLINK2_BIN" --bfile "$bfile" --extract "${prune}.prune.in" "${keep_args[@]}" "${remove_args[@]}" \
      --freq counts --pca allele-wts "$PCA_COUNT" vcols=chrom,ref,alt \
      --threads "$PLINK_THREADS" --memory "$MEMORY_MB" --out "$ref"; then
      log_info "${SOURCE_NAME} ${component}: unrelated-reference PCA failed; retrying with all samples"
      reference_keep=""
      keep_args=()
      pca_method="full_sample_projected_fallback"
      run_logged "${log_dir}/${SOURCE_NAME}.${component}.reference_pca.full_sample_fallback.log" \
        "$PLINK2_BIN" --bfile "$bfile" --extract "${prune}.prune.in" "${remove_args[@]}" \
        --freq counts --pca allele-wts "$PCA_COUNT" vcols=chrom,ref,alt \
        --threads "$PLINK_THREADS" --memory "$MEMORY_MB" --out "$ref"
    fi
  fi
  require_file "${ref}.eigenvec.allele"
  require_file "${ref}.eigenval"
  require_file "${ref}.acount"

  local score_end=$((5 + PCA_COUNT))
  if [[ "$recompute" == "1" || ! -s "${proj}.sscore" ]]; then
    run_logged "${log_dir}/${SOURCE_NAME}.${component}.projection.log" \
      "$PLINK2_BIN" --bfile "$bfile" --extract "${prune}.prune.in" "${remove_args[@]}" \
      --read-freq "${ref}.acount" \
      --score "${ref}.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
      --score-col-nums "6-${score_end}" \
      --threads "$PLINK_THREADS" --memory "$MEMORY_MB" --out "$proj"
  fi
  require_file "${proj}.sscore"

  local unrelated_args=()
  if [[ -n "$reference_keep" ]]; then unrelated_args=(--unrelated "$reference_keep"); fi
  run_rscript "${SCRIPT_DIR}/prepare_pc_outputs.R" \
    --source "$SOURCE_NAME" --component "$component" --method "$pca_method" \
    --scores "${proj}.sscore" --eigenval "${ref}.eigenval" "${unrelated_args[@]}" \
    --pheno "$PHENO_FILE" --batch "$COVAR_FILE" --out-prefix "$stem"
  if [[ "$DRY_RUN" != "1" ]]; then printf '%s\n' "$pca_method" > "$method_file"; fi
}

run_pca() {
  require_file "$PHENO_FILE"
  require_file "$COVAR_FILE"
  require_grm "$TR_GRM_PREFIX"
  local pc_dir="${OUT_BASE}/${SOURCE_NAME}/pc"
  local log_dir="${OUT_BASE}/${SOURCE_NAME}/logs"
  mkdir -p "$pc_dir" "$log_dir"
  local snv_bfile="${PANGENIE_BFILE_DIR}/${SOURCE_NAME}.SNV_INDEL"
  local sv_bfile="${PANGENIE_BFILE_DIR}/${SOURCE_NAME}.SV"
  require_bfile "$snv_bfile"
  require_bfile "$sv_bfile"

  local shared_tr_prefix="${OUT_BASE}/shared/TR.full_grm_pca"
  if [[ "$FORCE" == "1" || ! -s "${shared_tr_prefix}.eigenvec" ]]; then
    run_logged "${OUT_BASE}/shared/TR.full_grm_pca.log" \
      "$GCTA_BIN" --grm "$TR_GRM_PREFIX" --keep "$FILTER_TR_KEEP" --pca "$PCA_COUNT" \
      --thread-num "$GCTA_THREADS" --out "$shared_tr_prefix"
  fi
  require_file "${shared_tr_prefix}.eigenvec"
  require_file "${shared_tr_prefix}.eigenval"

  local snv_prune="${pc_dir}/${SOURCE_NAME}.SNV_INDEL.king_pruned"
  local unrelated="${pc_dir}/${SOURCE_NAME}.SNV_INDEL.king.king.cutoff.in.id"
  if [[ "$FORCE" == "1" || ! -s "${snv_prune}.prune.in" ]]; then
    run_logged "${log_dir}/${SOURCE_NAME}.SNV_INDEL.king_prune.log" \
      "$PLINK2_BIN" --bfile "$snv_bfile" --extract "${snv_bfile}.all_variants.extract" "${maf_args[@]}" "${remove_args[@]}" \
      --autosome --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_R2" \
      --threads "$PLINK_THREADS" --memory "$MEMORY_MB" --out "$snv_prune"
  fi
  require_file "${snv_prune}.prune.in"
  local pca_reference_keep="$unrelated"
  local pca_method="unrelated_projected"
  if [[ "$FORCE" == "1" || ! -s "$unrelated" ]]; then
    if ! try_logged "${log_dir}/${SOURCE_NAME}.SNV_INDEL.king.log" \
      "$PLINK2_BIN" --bfile "$snv_bfile" --extract "${snv_prune}.prune.in" "${remove_args[@]}" \
      --king-cutoff "$KING_CUTOFF" --threads "$PLINK_THREADS" --memory "$MEMORY_MB" \
      --out "${pc_dir}/${SOURCE_NAME}.SNV_INDEL.king"; then
      log_info "${SOURCE_NAME}: KING failed; using full-sample PC spaces"
      pca_reference_keep=""
      pca_method="full_sample_projected_fallback"
    fi
  fi
  if [[ -n "$pca_reference_keep" && "$DRY_RUN" != "1" ]]; then
    local n_unrelated
    n_unrelated="$(awk 'NF >= 2 {n++} END {print n+0}' "$pca_reference_keep")"
    if (( n_unrelated <= PCA_COUNT )); then
      log_info "${SOURCE_NAME}: only ${n_unrelated} unrelated samples; using full-sample PC spaces"
      pca_reference_keep=""
      pca_method="full_sample_projected_fallback"
    fi
  fi

  build_projected_pca SNV_INDEL "$snv_bfile" "${snv_bfile}.all_variants.extract" "$pca_reference_keep" "$pca_method"
  build_projected_pca SV "$sv_bfile" "${sv_bfile}.all_variants.extract" "$pca_reference_keep" "$pca_method"
  run_rscript "${SCRIPT_DIR}/prepare_pc_outputs.R" \
    --source "$SOURCE_NAME" --component TR --method full_grm_fallback \
    --scores "${shared_tr_prefix}.eigenvec" --eigenval "${shared_tr_prefix}.eigenval" \
    --pheno "$PHENO_FILE" --batch "$COVAR_FILE" \
    --out-prefix "${pc_dir}/${SOURCE_NAME}.TR"
}

run_plot() {
  local component prefix
  for component in SNV_INDEL SV TR; do
    prefix="${OUT_BASE}/${SOURCE_NAME}/pc/${SOURCE_NAME}.${component}"
    require_file "${prefix}.plot.tsv"
    require_file "${prefix}.eigenval.tsv"
    run_rscript "${SCRIPT_DIR}/plot_pca.R" \
      --plot-data "${prefix}.plot.tsv" --eigenval "${prefix}.eigenval.tsv" \
      --out-prefix "${OUT_BASE}/plots/${SOURCE_NAME}.${component}"
  done
}

run_grm() {
  local grm_dir="${OUT_BASE}/${SOURCE_NAME}/grm"
  local log_dir="${OUT_BASE}/${SOURCE_NAME}/logs"
  mkdir -p "$grm_dir" "$log_dir"
  local component bfile grm
  for component in SNV_INDEL SV; do
    bfile="${PANGENIE_BFILE_DIR}/${SOURCE_NAME}.${component}"
    grm="${grm_dir}/${SOURCE_NAME}.${component}"
    require_bfile "$bfile"
    if [[ "$FORCE" == "1" ]] || ! grm_complete "$grm"; then
      run_logged "${log_dir}/${SOURCE_NAME}.${component}.make_grm.log" \
        "$GCTA_BIN" --bfile "$bfile" "${maf_args[@]}" --make-grm \
        --thread-num "$GCTA_THREADS" --out "$grm"
    else
      log_info "Reusing complete GRM: ${grm}"
    fi
    require_grm "$grm"
  done
}

build_inputs() {
  local pc_dir="${OUT_BASE}/${SOURCE_NAME}/pc"
  local grm_dir="${OUT_BASE}/${SOURCE_NAME}/grm"
  local input_dir="${OUT_BASE}/${SOURCE_NAME}/inputs"
  mkdir -p "$input_dir"
  require_grm "${grm_dir}/${SOURCE_NAME}.SNV_INDEL"
  require_grm "${grm_dir}/${SOURCE_NAME}.SV"
  require_grm "$TR_GRM_PREFIX"
  run_rscript "${SCRIPT_DIR}/build_analysis_inputs.R" \
    --source "$SOURCE_NAME" --pheno "$PHENO_FILE" --batch "$COVAR_FILE" \
    --snv-pc10 "${pc_dir}/${SOURCE_NAME}.SNV_INDEL.pc10.qcovar.tsv" \
    --snv-pc20 "${pc_dir}/${SOURCE_NAME}.SNV_INDEL.pc20.qcovar.tsv" \
    --sv-pc10 "${pc_dir}/${SOURCE_NAME}.SV.pc10.qcovar.tsv" \
    --sv-pc20 "${pc_dir}/${SOURCE_NAME}.SV.pc20.qcovar.tsv" \
    --tr-pc10 "${pc_dir}/${SOURCE_NAME}.TR.pc10.qcovar.tsv" \
    --tr-pc20 "${pc_dir}/${SOURCE_NAME}.TR.pc20.qcovar.tsv" \
    --snv-grm-id "${grm_dir}/${SOURCE_NAME}.SNV_INDEL.grm.id" \
    --sv-grm-id "${grm_dir}/${SOURCE_NAME}.SV.grm.id" \
    --tr-grm-id "${TR_GRM_PREFIX}.grm.id" --remove-file "$FILTER_REMOVE_FILE" --out-dir "$input_dir"

  local snv_grm="${grm_dir}/${SOURCE_NAME}.SNV_INDEL"
  local sv_grm="${grm_dir}/${SOURCE_NAME}.SV"
  printf '%s\n%s\n' "$snv_grm" "$sv_grm" > "${input_dir}/${SOURCE_NAME}.SNV_INDEL_SV.mgrm"
  printf '%s\n%s\n' "$snv_grm" "$TR_GRM_PREFIX" > "${input_dir}/${SOURCE_NAME}.SNV_INDEL_TR.mgrm"
  printf '%s\n%s\n' "$sv_grm" "$TR_GRM_PREFIX" > "${input_dir}/${SOURCE_NAME}.SV_TR.mgrm"
  printf '%s\n%s\n%s\n' "$snv_grm" "$sv_grm" "$TR_GRM_PREFIX" > "${input_dir}/${SOURCE_NAME}.SNV_INDEL_SV_TR.mgrm"

  run_rscript "${SCRIPT_DIR}/make_reml_tasks.R" \
    --source "$SOURCE_NAME" --input-dir "$input_dir" --reml-dir "${OUT_BASE}/${SOURCE_NAME}/reml" \
    --log-dir "${OUT_BASE}/${SOURCE_NAME}/logs/reml" --snv-grm "$snv_grm" --sv-grm "$sv_grm" \
    --tr-grm "$TR_GRM_PREFIX" --output "${input_dir}/${SOURCE_NAME}.reml_tasks.tsv"
}

clean_reason() {
  local log_file="$1" reason=""
  [[ -s "$log_file" ]] || { echo "log file missing"; return; }
  reason="$(grep -Eim1 'not invertible|non-invertible|singular matrix|matrix.*singular|not positive definite|collinear|rank deficient|not converg|failed to converge|iteration limit|no individual|no common individual|too few|cannot open|failed|error' "$log_file" || true)"
  if [[ -z "$reason" ]]; then reason="$(tail -n 1 "$log_file" 2>/dev/null || true)"; fi
  printf '%s' "$reason" | tr '\t\r\n' '   '
}

append_status() {
  local file="$1" source="$2" task_id="$3" model="$4" pc_n="$5" adjustment="$6" status="$7" rc="$8" reason="$9" log_file="${10}" hsq="${11}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$source" "$task_id" "$model" "$pc_n" "$adjustment" "$status" "$rc" "$reason" "$log_file" "$hsq" "$(timestamp)" >> "$file"
}

run_reml_task() {
  local status_file="$1" source_name="$2" task_id="$3" model="$4" pc_n="$5" adjustment="$6" input_type="$7" input_path="$8" qcovar="$9" out_prefix="${10}" log_file="${11}"
  local hsq="${out_prefix}.hsq" rc=0 status="success" reason=""
  mkdir -p "$(dirname "$out_prefix")" "$(dirname "$log_file")"
  if [[ "$DRY_RUN" == "1" ]]; then
    if [[ "$input_type" == "grm" ]]; then
      print_cmd "$GCTA_BIN" --grm "$input_path" --keep "${OUT_BASE}/${source_name}/inputs/${source_name}.analysis.keep" --pheno "$PHENO_FILE" --covar "$COVAR_FILE" --qcovar "$qcovar" --reml --reml-no-constrain --reml-maxit "$REML_MAXIT" --prevalence "$PREVALENCE" --thread-num "$GCTA_THREADS" --out "$out_prefix"
    else
      print_cmd "$GCTA_BIN" --mgrm "$input_path" --keep "${OUT_BASE}/${source_name}/inputs/${source_name}.analysis.keep" --pheno "$PHENO_FILE" --covar "$COVAR_FILE" --qcovar "$qcovar" --reml --reml-no-constrain --reml-maxit "$REML_MAXIT" --prevalence "$PREVALENCE" --thread-num "$GCTA_THREADS" --out "$out_prefix"
    fi
    append_status "$status_file" "$source_name" "$task_id" "$model" "$pc_n" "$adjustment" dry_run 0 "" "$log_file" "$hsq"
    return 0
  fi
  if [[ "$FORCE" != "1" && -s "$hsq" ]]; then
    append_status "$status_file" "$source_name" "$task_id" "$model" "$pc_n" "$adjustment" reused 0 "existing hsq reused" "$log_file" "$hsq"
    return 0
  fi
  rm -f "$hsq"
  local cmd=("$GCTA_BIN")
  if [[ "$input_type" == "grm" ]]; then cmd+=(--grm "$input_path"); else cmd+=(--mgrm "$input_path"); fi
  cmd+=(--keep "${OUT_BASE}/${source_name}/inputs/${source_name}.analysis.keep"
    --pheno "$PHENO_FILE" --covar "$COVAR_FILE" --qcovar "$qcovar"
    --reml --reml-no-constrain --reml-maxit "$REML_MAXIT"
    --prevalence "$PREVALENCE" --thread-num "$GCTA_THREADS" --out "$out_prefix")
  { printf '[CMD] '; printf '%q ' "${cmd[@]}"; printf '\n'; } > "$log_file"
  set +e
  "${cmd[@]}" >> "$log_file" 2>&1
  rc=$?
  set -e
  if (( rc != 0 )); then
    status="failed"; reason="$(clean_reason "$log_file")"
  elif [[ ! -s "$hsq" ]]; then
    status="failed"; reason="GCTA returned zero but no non-empty .hsq was produced"
  elif grep -Eiq 'not converg|failed to converge|iteration limit' "$log_file"; then
    status="warning"; reason="$(clean_reason "$log_file")"
  fi
  append_status "$status_file" "$source_name" "$task_id" "$model" "$pc_n" "$adjustment" "$status" "$rc" "$reason" "$log_file" "$hsq"
}

run_reml() {
  local task_file="${OUT_BASE}/${SOURCE_NAME}/inputs/${SOURCE_NAME}.reml_tasks.tsv"
  local status_dir="${OUT_BASE}/${SOURCE_NAME}/status"
  local status_file="${status_dir}/${SOURCE_NAME}.reml_status.tsv"
  local task_status_dir="${status_dir}/task_status"
  local rc=0
  local -a pids=()
  local -a status_files=()
  local task_source task_id model pc_n adjustment input_type input_path qcovar out_prefix log_file components task_status_file
  require_file "$task_file"
  mkdir -p "$status_dir" "$task_status_dir"
  printf 'source\ttask_id\tmodel\tpc_n\tadjustment\tstatus\texit_code\terror_reason\tlog_file\thsq_file\tfinished_at\n' > "$status_file"
  log_info "Running REML tasks: reml_jobs=${REML_JOBS}, gcta_threads=${GCTA_THREADS}"
  while IFS=$'\t' read -r task_source task_id model pc_n adjustment input_type input_path qcovar out_prefix log_file components; do
    [[ "$task_source" == "source" ]] && continue
    task_status_file="${task_status_dir}/${task_id}.status.tsv"
    rm -f "$task_status_file"
    status_files+=("$task_status_file")
    if (( REML_JOBS <= 1 )); then
      run_reml_task "$task_status_file" "$task_source" "$task_id" "$model" "$pc_n" "$adjustment" "$input_type" "$input_path" "$qcovar" "$out_prefix" "$log_file"
    else
      run_reml_task "$task_status_file" "$task_source" "$task_id" "$model" "$pc_n" "$adjustment" "$input_type" "$input_path" "$qcovar" "$out_prefix" "$log_file" &
      pids+=("$!")
      if (( ${#pids[@]} >= REML_JOBS )); then
        if ! wait "${pids[0]}"; then rc=1; fi
        pids=("${pids[@]:1}")
      fi
    fi
  done < "$task_file"
  while (( ${#pids[@]} > 0 )); do
    if ! wait "${pids[0]}"; then rc=1; fi
    pids=("${pids[@]:1}")
  done
  for task_status_file in "${status_files[@]}"; do
    if [[ -s "$task_status_file" ]]; then
      cat "$task_status_file" >> "$status_file"
    else
      echo "[WARN] Missing per-task status file: ${task_status_file}" >&2
      rc=1
    fi
  done
  return "$rc"
}

run_summarise() {
  run_rscript "${SCRIPT_DIR}/collect_results.R" \
    --sources "$SOURCE_NAME" --out-base "$OUT_BASE" --summary-dir "${OUT_BASE}/summary" \
    --plot-dir "${OUT_BASE}/plots"
}

if [[ "$MODE" == "all" || "$MODE" == "prepare" ]]; then prepare_pangenie_bfiles; fi
if [[ "$MODE" == "all" || "$MODE" == "pca" ]]; then run_pca; fi
if [[ "$MODE" == "all" || "$MODE" == "plot" ]]; then run_plot; fi
if [[ "$MODE" == "all" || "$MODE" == "grm" ]]; then run_grm; fi
if [[ "$MODE" == "all" || "$MODE" == "inputs" || "$MODE" == "reml" ]]; then build_inputs; fi
if [[ "$MODE" == "all" || "$MODE" == "reml" ]]; then run_reml; fi
if [[ "$MODE" == "all" || "$MODE" == "summarise" ]]; then run_summarise; fi

log_info "Merged set00-set02 heritability workflow finished (mode=${MODE})"
