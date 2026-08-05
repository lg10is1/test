#!/usr/bin/env bash
set -euo pipefail

START_AT=1
STOP_AFTER=20
RESUME=1
RESET_CHECKPOINTS=0
SEED=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-at)
      START_AT="$2"
      shift 2
      ;;
    --stop-after)
      STOP_AFTER="$2"
      shift 2
      ;;
    --force)
      RESUME=0
      shift
      ;;
    --reset-checkpoints)
      RESET_CHECKPOINTS=1
      shift
      ;;
    --seed)
      SEED="$2"
      shift 2
      ;;
    -h|--help)
      cat <<'EOF'
Usage:
  bash run_pipeline.sh [--start-at N] [--stop-after N] [--force] [--reset-checkpoints] [--seed N]

Steps:
  1  Validate the final DeepVariant and selected Paragraph GWAS/clumping files
  2  Seed public GWAS/clumping workspace from v20 and remove manual-excluded lead
  3  Validate reused Pangenie clumping files
  4  Rscript 03_gwas_draw_manhattan.R
  5  Select canonical leads, then prepare LD IDs and PLINK commands
  6  Run generated PLINK LD command script
  7  Rscript 05_ld_summarise_decay.R
  8  Rscript 06_ld_friends_scores.R
  9  Refresh source-level and canonical lead tables
  10 ANNOVAR standard + t2t_noloc
  11 GWAS Catalog schizophrenia-related annotation
  12 Final annotation merge and 1000 kb windows
  13 Prepare SGV/SV lead-specific locuszoom LD jobs
  14 Run locuszoom PLINK LD jobs
  15 Draw SGV/SV locuszoom plots with Nature comparison
  16 Match final SV leads to the PAV SV VCF
  17 Run SV cis-meQTL and save SV tables/plots
  18 Record SV-only meQTL policy (SGV meQTL intentionally disabled)
  19 Test detectable LD and per-SV maximum r2 for lead versus null SVs
  20 Collect public full and main-text summary packages

Resume behavior:
  Completed steps are skipped using .pipeline_checkpoints/step_XX.done.
  --force              Run selected steps even if checkpoint files exist.
  --reset-checkpoints  Remove checkpoint files before running.
  --seed N             Integer random seed for null sampling, LD-decay bootstrap, and permutations (default: 1).

Examples:
  bash run_pipeline.sh --force
  bash run_pipeline.sh --seed 20260713 --start-at 5 --force
  bash run_pipeline.sh --start-at 1 --stop-after 1 --force
  bash run_pipeline.sh --start-at 2 --force
  bash run_pipeline.sh --start-at 13 --stop-after 15 --force
  bash run_pipeline.sh --start-at 16 --stop-after 18 --force
  bash run_pipeline.sh --start-at 19 --stop-after 19 --force
  bash run_pipeline.sh --start-at 20 --stop-after 20 --force
EOF
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if ! [[ "$SEED" =~ ^-?[0-9]+$ ]]; then
  echo "[ERROR] --seed must be an integer: ${SEED}" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DEEPVARIANT_PARAGRAPH_DIR="/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph"
DEEPVAR_ROOT="${DEEPVARIANT_PARAGRAPH_DIR}/deepvar_gwas/deepvar"
PARAGRAPH_MODEL="pcsrc_deepvar_pc20_grm_deepvar_with_batch"
PARAGRAPH_ROOT="${DEEPVARIANT_PARAGRAPH_DIR}/deepvar_gwas/paragraph_test"
DEEPVAR_GWAS="${DEEPVAR_ROOT}/03_gwas/SCZ.deepvar.mlm.geno0.1.maf0.01.fastGWA"
PARAGRAPH_GWAS="${PARAGRAPH_ROOT}/03_gwas/SCZ.paragraph_test.${PARAGRAPH_MODEL}.mlm.geno0.1.maf0.01.fastGWA"
PARAGRAPH_VALID_GWAS="${PARAGRAPH_ROOT}/03_gwas/SCZ.paragraph_test.${PARAGRAPH_MODEL}.mlm.geno0.1.maf0.01.valid_statistics.fastGWA.tsv"
DEEPVAR_CLUMP="${DEEPVAR_ROOT}/04_clumping/SCZ.deepvar.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped"
PARAGRAPH_CLUMP_PREFIX="${PARAGRAPH_ROOT}/04_clumping/SCZ.paragraph_test.${PARAGRAPH_MODEL}.mlm.clump_p1_5e-6.r2_0.01.kb_1000"
PARAGRAPH_CLUMP="${PARAGRAPH_CLUMP_PREFIX}.clumped"
PARAGRAPH_CLUMP_VALIDATION="${PARAGRAPH_CLUMP_PREFIX}.clump_validation.audit.tsv"
PARAGRAPH_CLUMP_COMPLETE="${PARAGRAPH_CLUMP_PREFIX}.clump.complete"
PARAGRAPH_MODEL_SUMMARY="${PARAGRAPH_ROOT}/gwas_model_summary.tsv"

GWAS_FIGURE_DIR="/path/to/EOSCZ_PROJECT/figure_analysis/01.GWAS_figure.public"
MANHATTAN_SV="/path/to/EOSCZ_PROJECT/figure_analysis/sv_gwas.public.tiff"
MANHATTAN_COMBINED="/path/to/EOSCZ_PROJECT/figure_analysis/sv_snv_gwas.public.tiff"
CLUMP_DIR="${GWAS_FIGURE_DIR}/clumping_by_set_subtype"
SOURCE_GWAS_FIGURE_DIR="/path/to/EOSCZ_PROJECT/Figure3/01.GWAS_figure.version20"
MANUAL_EXCLUSION_FILE="${SCRIPT_DIR}/manual_excluded_leads.tsv"
LD_OUTDIR="/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public"
PLINK_LD_SCRIPT="${LD_OUTDIR}/cmd/run_plink_ld.sh"
ANNOVAR_ANNOTATED="${LD_OUTDIR}/tables/lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.tsv"
ANNOVAR_ANNOTATED_NOLOC="${LD_OUTDIR}/tables/lead_sig_from_gwas.all.annovar_noloc.refGene.annotated.fixed.tsv"
GWAS_CATALOG_ZIP="gwas_catalog_v1.0.2-associations_full.zip"
GWAS_CATALOG_ANNOTATED="${LD_OUTDIR}/tables/lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.gwas_scz_related.tsv"
GWAS_CATALOG_ANNOTATED_NOLOC="${LD_OUTDIR}/tables/lead_sig_from_gwas.all.annovar_noloc.refGene.annotated.fixed.gwas_scz_related.tsv"
FINAL_MERGED="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.tsv"
FINAL_WINDOW="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.window_1000kb.tsv"
FINAL_WINDOW_SNV_INDEL="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.window_1000kb.snv_indel.tsv"
FINAL_WINDOW_SV="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.window_1000kb.sv.tsv"
FINAL_MERGED_NOLOC="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.noloc.tsv"
FINAL_WINDOW_NOLOC="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.noloc.window_1000kb.tsv"
FINAL_WINDOW_NOLOC_SNV_INDEL="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.noloc.window_1000kb.snv_indel.tsv"
FINAL_WINDOW_NOLOC_SV="${LD_OUTDIR}/tables/lead_sig_from_gwas.final_merged.cleaned.noloc.window_1000kb.sv.tsv"
LD_FRIENDS_TABLE="${LD_OUTDIR}/ld_friends_scores/ld_friends_scores.wilcoxon.tsv"
LD_FRIENDS_FIGURE="${LD_OUTDIR}/figures/ld_friends_scores_sig_vs_null.all_partners.pdf"
LD_FRIENDS_TABLE_250KB="${LD_OUTDIR}/ld_friends_scores/ld_friends_scores.wilcoxon.250kb.tsv"
LD_FRIENDS_FIGURE_250KB="${LD_OUTDIR}/figures/ld_friends_scores_sig_vs_null.all_partners.250kb.pdf"
SV_DETECTABLE_PREFIX="${LD_OUTDIR}/rdata/sv_sig_vs_null.proximal_snv_indel"
SV_DETECTABLE_MAXR2="${SV_DETECTABLE_PREFIX}.maxR2.tsv"
SV_DETECTABLE_SUMMARY="${SV_DETECTABLE_PREFIX}.summary.tsv"
SV_DETECTABLE_CONTINGENCY="${SV_DETECTABLE_PREFIX}.contingency.tsv"
SV_DETECTABLE_TESTS="${SV_DETECTABLE_PREFIX}.formal_tests.tsv"
SV_DETECTABLE_CONFIG="${SV_DETECTABLE_PREFIX}.run_config.tsv"
SV_DETECTABLE_THRESHOLDS="${SV_DETECTABLE_PREFIX}.threshold_sensitivity.tsv"
SV_DETECTABLE_DISTRIBUTIONS="${SV_DETECTABLE_PREFIX}.distribution_tests.tsv"
SV_DETECTABLE_STRATA="${SV_DETECTABLE_PREFIX}.permutation_strata.tsv"
LEAD_TABLE_ALL="${LD_OUTDIR}/tables/lead_sig_from_gwas.all.tsv"
LEAD_TABLE="${LD_OUTDIR}/tables/lead_sig_from_gwas.canonical_1000kb.tsv"
LEAD_MAPPING="${LD_OUTDIR}/tables/lead_sig_from_gwas.canonical_1000kb.mapping.tsv"
LEAD_CANONICAL_SUMMARY="${LD_OUTDIR}/tables/lead_sig_from_gwas.canonical_1000kb.summary.tsv"
LEAD_MANUAL_EXACT_AUDIT="${LD_OUTDIR}/tables/lead_sig_from_gwas.manual_excluded_exact_representatives.tsv"
LOCUSZOOM_ROOT="${LD_OUTDIR}/locuszoom"
LOCUSZOOM_LD_SCRIPT="${LOCUSZOOM_ROOT}/run_locuszoom_ld.sh"
LOCUSZOOM_SUMMARY="${LOCUSZOOM_ROOT}/locuszoom_plot_summary.tsv"
MEQTL_ROOT="/path/to/EOSCZ_PROJECT/figure_analysis/02.meQTL/public"
MEQTL_SV_DIR="${MEQTL_ROOT}/SV"
MEQTL_SGV_DIR="${MEQTL_ROOT}/SGV"
SV_MATCH_DIR="${MEQTL_SV_DIR}/matching"
SV_MATCH_FILE="${SV_MATCH_DIR}/sig_sv_to_pav_sv_len50.best.tsv"
SV_MEQTL_DEBUG_PREFIX="${MEQTL_SV_DIR}/tables/sv_meqtl.genotype_debug"
SV_MEQTL_MAC_SUMMARY="${MEQTL_SV_DIR}/tables/meqtl_cis.mac_filter_summary.tsv"
PAV_SV_VCF="/path/to/EOSCZ_PROJECT/TGS_SV_merge_SCZ/truvari_single_sample/truvari_merged_sort_pP0.5.sv_len_gt50.sorted.vcf.gz"
BCFTOOLS_BIN="${HOME}/mambaforge/envs/truvari5/bin/bcftools"
CHECKPOINT_DIR="${SCRIPT_DIR}/.pipeline_checkpoints"

mkdir -p "$CHECKPOINT_DIR"
if (( RESET_CHECKPOINTS == 1 )); then
  rm -f "${CHECKPOINT_DIR}"/step_*.done
  echo "[INFO] Reset checkpoints in ${CHECKPOINT_DIR}"
fi

checkpoint_file() {
  printf "%s/step_%02d.done" "$CHECKPOINT_DIR" "$1"
}

run_step() {
  local step="$1"
  local label="$2"
  shift 2
  local marker
  marker="$(checkpoint_file "$step")"

  if (( step < START_AT || step > STOP_AFTER )); then
    echo "[SKIP] Step ${step}: ${label}"
    return 0
  fi
  if (( RESUME == 1 )) && [[ -s "$marker" ]]; then
    if grep -q '^seed=' "$marker"; then
      local marker_seed
      marker_seed="$(awk -F= '$1 == "seed" {print $2; exit}' "$marker")"
      if [[ "$marker_seed" == "$SEED" ]]; then
        echo "[RESUME] Step ${step} already completed: ${label}"
        echo "         marker: ${marker}"
        return 0
      fi
      echo "[RERUN] Step ${step}: checkpoint seed=${marker_seed}, requested seed=${SEED}"
    elif [[ "$SEED" == "1" ]]; then
      echo "[RESUME] Step ${step} already completed: ${label}"
      echo "         marker: ${marker}"
      return 0
    else
      echo "[RERUN] Step ${step}: checkpoint lacks seed, requested seed=${SEED}"
    fi
  fi

  echo
  echo "============================================================"
  echo "[RUN] Step ${step}: ${label}"
  echo "============================================================"
  "$@"
  {
    echo "step=${step}"
    echo "label=${label}"
    echo "completed_at=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "script_dir=${SCRIPT_DIR}"
    echo "seed=${SEED}"
  } > "$marker"
  echo "[CHECKPOINT] Wrote ${marker}"
}

check_file() {
  [[ -s "$1" ]] || { echo "[ERROR] Required file not found or empty: $1" >&2; exit 1; }
}

check_autosome_tsv() {
  local file="$1"
  local column="$2"
  check_file "$file"
  awk -F '\t' -v target="$column" '
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == target) idx = i
      if (!idx) {
        printf("[ERROR] Column %s not found in %s\n", target, FILENAME) > "/dev/stderr"
        exit 2
      }
      next
    }
    {
      value = $idx
      sub(/^[Cc][Hh][Rr]/, "", value)
      if (value !~ /^([1-9]|1[0-9]|2[0-2])$/) {
        printf("[ERROR] Non-autosomal value %s in column %s at %s:%d\n", $idx, target, FILENAME, NR) > "/dev/stderr"
        exit 1
      }
    }
  ' "$file"
  echo "[AUTOSOME CHECK OK] ${file} | column=${column}"
}

verify_selected_source_results() {
  check_file "$DEEPVAR_GWAS"
  check_file "$DEEPVAR_CLUMP"
  check_file "$PARAGRAPH_MODEL_SUMMARY"
  check_file "$PARAGRAPH_GWAS"
  check_file "$PARAGRAPH_VALID_GWAS"
  check_file "$PARAGRAPH_CLUMP"
  check_file "$PARAGRAPH_CLUMP_VALIDATION"
  check_file "$PARAGRAPH_CLUMP_COMPLETE"

  grep -Fq "$PARAGRAPH_MODEL" "$PARAGRAPH_MODEL_SUMMARY" || {
    echo "[ERROR] Selected Paragraph model is absent from ${PARAGRAPH_MODEL_SUMMARY}: ${PARAGRAPH_MODEL}" >&2
    exit 1
  }
  awk -F '\t' 'NR == 2 && $1 ~ /^PASS/ {ok=1} END {exit !ok}' "$PARAGRAPH_CLUMP_VALIDATION" || {
    echo "[ERROR] Selected Paragraph clump validation did not PASS: ${PARAGRAPH_CLUMP_VALIDATION}" >&2
    exit 1
  }
  echo "[SOURCE CHECK OK] DeepVariant GWAS/clump and selected Paragraph model are complete."
}

check_pangenie_clumping_outputs() {
  local set_name subtype
  for set_name in set00; do
    for subtype in SNV_INDEL SV; do
      check_file "${CLUMP_DIR}/${set_name}.${subtype}.clump_p1_5e-06.r2_0.01.kb_1000.clumped"
    done
  done
}

apply_manual_exclusions() {
  check_file "$MANUAL_EXCLUSION_FILE"
  check_file "${GWAS_FIGURE_DIR}/all_sets.merged_minP.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv"
  python3 00_apply_manual_exclusions.py \
    --root "$GWAS_FIGURE_DIR" \
    --exclusions "$MANUAL_EXCLUSION_FILE" \
    --audit "${GWAS_FIGURE_DIR}/manual_exclusion_filter.audit.tsv"
}

seed_gwas_workspace_from_v20() {
  check_file "${SOURCE_GWAS_FIGURE_DIR}/all_sets.merged_minP.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv"
  if [[ ! -d "$GWAS_FIGURE_DIR" ]]; then
    mkdir -p "$(dirname "$GWAS_FIGURE_DIR")"
    cp -a "$SOURCE_GWAS_FIGURE_DIR" "$GWAS_FIGURE_DIR"
    echo "[SEED] Copied v20 GWAS/clumping workspace to ${GWAS_FIGURE_DIR}"
  else
    echo "[SEED] Reusing existing GWAS/clumping workspace: ${GWAS_FIGURE_DIR}"
  fi
  apply_manual_exclusions
}

draw_manhattan() {
  apply_manual_exclusions
  Rscript 03_gwas_draw_manhattan.R
}

check_commands() {
  local cmd
  for cmd in Rscript python3 plink plink2 gcta64 xargs gzip; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "[ERROR] Command not found: $cmd" >&2; exit 1; }
  done
  Rscript -e 'suppressPackageStartupMessages({library(data.table); library(dplyr); library(ggplot2); library(cowplot); library(ggsci); library(qqman)})' >/dev/null
  python3 - <<'PY' >/dev/null
import pandas
PY
}

check_detectable_ld_commands() {
  command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Command not found: Rscript" >&2; exit 1; }
  Rscript -e 'suppressPackageStartupMessages({library(data.table); library(ggplot2); library(cowplot)})' >/dev/null
}

run_gwas_catalog_annotations() {
  python3 08_annotate_gwas_catalog_scz_related.py --input "$ANNOVAR_ANNOTATED" --gwas-zip "$GWAS_CATALOG_ZIP"
  python3 08_annotate_gwas_catalog_scz_related.py --input "$ANNOVAR_ANNOTATED_NOLOC" --gwas-zip "$GWAS_CATALOG_ZIP"
}

run_final_merges() {
  python3 09_merge_annotations_and_windows.py \
    --table-dir "${LD_OUTDIR}/tables" --window-kb 1000 \
    --lead-table "$LEAD_TABLE"
  python3 09_merge_annotations_and_windows.py \
    --table-dir "${LD_OUTDIR}/tables" --window-kb 1000 \
    --lead-table "$LEAD_TABLE" \
    --catalog-table "$GWAS_CATALOG_ANNOTATED_NOLOC" \
    --out-final "$FINAL_MERGED_NOLOC" --out-window "$FINAL_WINDOW_NOLOC" \
    --out-window-snv-indel "$FINAL_WINDOW_NOLOC_SNV_INDEL" \
    --out-window-sv "$FINAL_WINDOW_NOLOC_SV"
}
prepare_ld_ids_and_commands() {
  prepare_canonical_leads
  LD_OUTDIR="$LD_OUTDIR" bash 04_prepare_pangenie_maf01_lists.sh
  LD_OUTDIR="$LD_OUTDIR" bash 04_prepare_deepvariant_paragraph_maf02_list.sh
  Rscript 04_ld_prepare_ids_and_plink_cmd.R --seed "$SEED"
}

prepare_canonical_leads() {
  Rscript 06_extract_lead_sig_from_gwas.R
  python3 06_select_canonical_window_leads.py \
    --table-dir "${LD_OUTDIR}/tables" \
    --window-kb 1000 \
    --manual-exclusions "$MANUAL_EXCLUSION_FILE"
  check_file "$LEAD_TABLE_ALL"
  check_file "$LEAD_TABLE"
  check_file "$LEAD_MAPPING"
  check_file "$LEAD_CANONICAL_SUMMARY"
  check_file "$LEAD_MANUAL_EXACT_AUDIT"
  check_autosome_tsv "$LEAD_TABLE" lead_chr
}

check_bcftools() {
  [[ -x "$BCFTOOLS_BIN" ]] || {
    echo "[ERROR] Fixed bcftools executable not found: $BCFTOOLS_BIN" >&2
    exit 1
  }
}

run_sv_match() {
  check_bcftools
  Rscript 12_match_sig_sv_to_pav.R \
    --sig-file "$LEAD_TABLE" \
    --pav-sv-vcf "$PAV_SV_VCF" \
    --bcftools "$BCFTOOLS_BIN" \
    --out-dir "$SV_MATCH_DIR" \
    --window-bp 1000
}

run_sv_meqtl() {
  check_bcftools
  Rscript 13_run_sv_meqtl_cis.R \
    --matches "$SV_MATCH_FILE" \
    --pav-sv-vcf "$PAV_SV_VCF" \
    --bcftools "$BCFTOOLS_BIN" \
    --out-dir "$MEQTL_SV_DIR" \
    --window-bp 1000000 \
    --match-types exact_ref_alt,fuzzy_ref_same
  Rscript 13b_qc_sv_meqtl_genotypes.R \
    --matches "$SV_MATCH_FILE" \
    --pav-sv-vcf "$PAV_SV_VCF" \
    --bcftools "$BCFTOOLS_BIN" \
    --out-prefix "$SV_MEQTL_DEBUG_PREFIX" \
    --min-n 50 \
    --match-types exact_ref_alt,fuzzy_ref_same
  Rscript 13c_filter_sv_meqtl_by_mac.R \
    --input "${MEQTL_SV_DIR}/tables/meqtl_cis.all_results.tsv.gz" \
    --out-dir "${MEQTL_SV_DIR}/tables" \
    --mac-thresholds 2,5 \
    --fdr-threshold 0.05 \
    --p-threshold 1e-8
}

record_sgv_meqtl_disabled() {
  mkdir -p "$MEQTL_SGV_DIR"
  printf 'analysis\tstatus\treason\nSGV_cis_meQTL\tDISABLED\tSV-only meQTL policy; do not read NGS or SGV genotype files\n' \
    > "${MEQTL_SGV_DIR}/sgv_meqtl.disabled.tsv"
  echo "[INFO] SGV meQTL disabled by SV-only policy; no NGS/SGV genotype file was read."
}

if (( START_AT == 19 && STOP_AFTER == 19 )); then
  check_detectable_ld_commands
elif (( START_AT == 20 && STOP_AFTER == 20 )); then
  command -v python3 >/dev/null 2>&1 || { echo "[ERROR] Command not found: python3" >&2; exit 1; }
  python3 - <<'PY' >/dev/null
import pandas
PY
else
  check_commands
fi
echo "[INFO] Pipeline code version: public"
echo "[INFO] Random seed: ${SEED}"
echo "[INFO] Output paths: public-compatible paths are intentionally reused and may be overwritten."
echo "[INFO] Chromosome policy: autosomes only (chr1-22); chr23/X and all non-autosomes are excluded."
echo "[INFO] Step 1 validates existing final source results; it does not rerun GWAS."
echo "[INFO] DeepVariant final GWAS: ${DEEPVAR_GWAS}"
echo "[INFO] Paragraph final model: ${PARAGRAPH_MODEL}"
echo "[INFO] Paragraph model settings: DeepVariant PC20 + DeepVariant sparse GRM + SCZ_batch"
echo "[INFO] Paragraph final GWAS: ${PARAGRAPH_GWAS}"
echo "[INFO] Paragraph final clump: ${PARAGRAPH_CLUMP}"
echo "[INFO] Source clumping: invalid or missing P values are removed before PLINK"
echo "[INFO] LD genotype: common-sample merged DeepVariant + Paragraph bfile"
echo "[INFO] New downstream output: ${LD_OUTDIR}"
echo "[INFO] Locuszoom: SGV/SV separate; own GWAS colored by lead-specific r2; Nature comparison below"
echo "[INFO] meQTL: SV cis-meQTL tables, per-lead results, and plots saved under ${MEQTL_ROOT}"
echo "[INFO] Fixed bcftools: ${BCFTOOLS_BIN}"
echo "[INFO] Step 19: Fisher test of detectable LD (per-SV max r2 > 0.1) and Wilcoxon test of per-SV max r2"
echo "[INFO] LD null sampling: 100 null variants per significant index variant"
echo "[INFO] deepvariant_paragraph null sampling: null candidates require MAF >= 0.02; significant index IDs are not filtered by this threshold"
echo "[INFO] Canonical lead counts are data-driven and recorded in ${LEAD_CANONICAL_SUMMARY}"
echo "[INFO] Manual lead exclusion: exact source_set + variant_type + lead_id rows are removed; if an excluded row was the original canonical representative, no replacement representative is selected"
echo "[INFO] Paragraph SV rule: retain every already-clumped Paragraph SV independently unless that exact Paragraph lead row is manually excluded"
echo "[INFO] meQTL SV matching used for modeling: exact_ref_alt and fuzzy_ref_same"
echo "[INFO] meQTL scope: SV only; SGV cis-meQTL is intentionally disabled"

run_step 1 "validate final DeepVariant and selected Paragraph GWAS/clumping" verify_selected_source_results

run_step 2 "seed public GWAS/clumping workspace from v20 and apply manual exclusion" seed_gwas_workspace_from_v20
if (( START_AT <= 2 && STOP_AFTER >= 2 )); then
  check_file "${GWAS_FIGURE_DIR}/all_sets.merged_minP.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv"
  check_file "${GWAS_FIGURE_DIR}/all_sets.merged_minP.SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.tsv"
  check_autosome_tsv "${GWAS_FIGURE_DIR}/all_sets.merged_minP.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv" CHR
  check_autosome_tsv "${GWAS_FIGURE_DIR}/all_sets.merged_minP.SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.tsv" CHR
fi

run_step 3 "validate reused Pangenie clumping files" check_pangenie_clumping_outputs
if (( START_AT <= 3 && STOP_AFTER >= 3 )); then check_pangenie_clumping_outputs; fi

run_step 4 "draw Manhattan plots from manually filtered GWAS results" draw_manhattan
if (( START_AT <= 4 && STOP_AFTER >= 4 )); then
  check_file "$MANHATTAN_SV"
  check_file "$MANHATTAN_COMBINED"
fi
run_step 5 "select canonical leads and prepare MAF-filtered LD-decay IDs/PLINK commands" prepare_ld_ids_and_commands
if (( START_AT <= 5 && STOP_AFTER >= 5 )); then
  check_file "${LD_OUTDIR}/maf01/set00.maf0.01.snplist"
  check_file "${LD_OUTDIR}/maf02/deepvariant_paragraph.maf0.02.snplist"
  check_file "${LD_OUTDIR}/ld_jobs.metadata.tsv"
  check_file "$PLINK_LD_SCRIPT"
fi

run_step 6 "run PLINK LD jobs" bash "$PLINK_LD_SCRIPT"
run_step 7 "summarise LD decay and draw LD plots" Rscript 05_ld_summarise_decay.R --seed "$SEED"
run_step 8 "calculate LD friends and LD scores" Rscript 06_ld_friends_scores.R
if (( START_AT <= 8 && STOP_AFTER >= 8 )); then
  check_file "$LD_FRIENDS_TABLE"; check_file "$LD_FRIENDS_FIGURE"
  check_file "$LD_FRIENDS_TABLE_250KB"; check_file "$LD_FRIENDS_FIGURE_250KB"
fi

run_step 9 "refresh source-level and canonical lead tables" prepare_canonical_leads
if (( START_AT <= 9 && STOP_AFTER >= 9 )); then
  check_file "$LEAD_TABLE_ALL"
  check_file "$LEAD_TABLE"
  check_file "$LEAD_MANUAL_EXACT_AUDIT"
  check_autosome_tsv "$LEAD_TABLE" lead_chr
fi

run_step 10 "ANNOVAR t2t/refGene and t2t_noloc/refGene" bash 07_run_annovar_t2t_refgene.sh
if (( START_AT <= 10 && STOP_AFTER >= 10 )); then
  check_file "$ANNOVAR_ANNOTATED"; check_file "$ANNOVAR_ANNOTATED_NOLOC"
fi

run_step 11 "GWAS Catalog schizophrenia-related annotation" run_gwas_catalog_annotations
if (( START_AT <= 11 && STOP_AFTER >= 11 )); then
  check_file "$GWAS_CATALOG_ANNOTATED"; check_file "$GWAS_CATALOG_ANNOTATED_NOLOC"
fi

run_step 12 "merge final annotations and 1000 kb windows" run_final_merges
if (( START_AT <= 12 && STOP_AFTER >= 12 )); then
  for file in "$FINAL_MERGED" "$FINAL_WINDOW" "$FINAL_WINDOW_SNV_INDEL" "$FINAL_WINDOW_SV" \
    "$FINAL_MERGED_NOLOC" "$FINAL_WINDOW_NOLOC" "$FINAL_WINDOW_NOLOC_SNV_INDEL" "$FINAL_WINDOW_NOLOC_SV"; do
    check_file "$file"
  done
fi

run_step 13 "prepare SGV/SV lead-specific locuszoom LD jobs" Rscript 10_prepare_locuszoom_ld.R
if (( START_AT <= 13 && STOP_AFTER >= 13 )); then
  check_file "${LOCUSZOOM_ROOT}/locuszoom_jobs.metadata.tsv"
  check_autosome_tsv "${LOCUSZOOM_ROOT}/locuszoom_jobs.metadata.tsv" lead_chr
  check_file "$LOCUSZOOM_LD_SCRIPT"
fi

run_step 14 "run lead-specific PLINK LD for locuszoom" bash "$LOCUSZOOM_LD_SCRIPT"
run_step 15 "draw SGV/SV locuszoom plots with Nature comparison" Rscript 11_plot_locuszoom_with_nature.R
if (( START_AT <= 15 && STOP_AFTER >= 15 )); then
  check_file "$LOCUSZOOM_SUMMARY"
  check_file "${LOCUSZOOM_ROOT}/SGV/locuszoom_plot_summary.tsv"
  check_file "${LOCUSZOOM_ROOT}/SV/locuszoom_plot_summary.tsv"
fi

run_step 16 "match final SV leads to PAV SV VCF" run_sv_match
if (( START_AT <= 16 && STOP_AFTER >= 16 )); then check_file "$SV_MATCH_FILE"; fi

run_step 17 "run SV cis-meQTL" run_sv_meqtl
if (( START_AT <= 17 && STOP_AFTER >= 17 )); then
  check_file "${MEQTL_SV_DIR}/tables/meqtl_cis.summary.tsv"
  check_file "${SV_MEQTL_DEBUG_PREFIX}.tsv"
  check_file "${SV_MEQTL_DEBUG_PREFIX}.reason_summary.tsv"
  check_file "$SV_MEQTL_MAC_SUMMARY"
  check_file "${MEQTL_SV_DIR}/tables/meqtl_cis.per_sv_mac.tsv"
fi

run_step 18 "record SV-only meQTL policy; SGV meQTL disabled" record_sgv_meqtl_disabled
if (( START_AT <= 18 && STOP_AFTER >= 18 )); then
  check_file "${MEQTL_SGV_DIR}/sgv_meqtl.disabled.tsv"
fi

run_step 19 "test SV detectable LD and per-SV maximum r2" \
  Rscript 15_sv_detectable_ld_tests.R \
    --outdir "$LD_OUTDIR" \
    --threshold 0.1 \
    --max-dist-bp 1000000 \
    --seed "$SEED"
if (( START_AT <= 19 && STOP_AFTER >= 19 )); then
  check_file "$SV_DETECTABLE_MAXR2"
  check_file "$SV_DETECTABLE_SUMMARY"
  check_file "$SV_DETECTABLE_CONTINGENCY"
  check_file "$SV_DETECTABLE_TESTS"
  check_file "$SV_DETECTABLE_CONFIG"
  check_file "$SV_DETECTABLE_THRESHOLDS"
  check_file "$SV_DETECTABLE_DISTRIBUTIONS"
  check_file "$SV_DETECTABLE_STRATA"
fi

run_step 20 "collect public summary and main-text result packages" bash 20_collect_results.sh

echo
echo "============================================================"
echo "[DONE] public pipeline finished."
echo "============================================================"
