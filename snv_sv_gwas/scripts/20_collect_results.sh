#!/usr/bin/env bash
set -euo pipefail

ROOT="${EOSCZ_PROJECT_ROOT:?Set EOSCZ_PROJECT_ROOT}"
FIGURE3="${ROOT}/figure_analysis"
SCRIPT_DIR="${FIGURE3}/scripts"
SUMMARY="${FIGURE3}/summary_public"
SUMMARY_MAIN="${FIGURE3}/summary_maintext_public"
GWAS_DIR="${FIGURE3}/01.GWAS_figure.public"
LD_DIR="${FIGURE3}/SV_SNV_LD/LD_decay_public"
MEQTL_DIR="${FIGURE3}/02.meQTL/public"
SV_MEQTL_ALL="${MEQTL_DIR}/SV/tables/meqtl_cis.all_results.tsv.gz"
SV_MEQTL_MAC_SUMMARY="${MEQTL_DIR}/SV/tables/meqtl_cis.mac_filter_summary.tsv"
SV_MEQTL_MAC_SCRIPT="${SCRIPT_DIR}/13c_filter_sv_meqtl_by_mac.R"
DP_DIR="${ROOT}/GWAS/Deepvariant_paragraph"
PARAGRAPH_MODEL="pcsrc_deepvar_pc20_grm_deepvar_with_batch"
PARAGRAPH_DIR="${DP_DIR}/deepvar_gwas/paragraph_test"
DEEPVAR_DIR="${DP_DIR}/deepvar_gwas/deepvar"

# Keep derived MAC summaries synchronized with the completed SV meQTL table.
# This makes a standalone Step 20 sufficient for both MAC post-processing and collection.
if [[ -s "${SV_MEQTL_ALL}" ]]; then
  if [[ ! -s "${SV_MEQTL_MAC_SCRIPT}" ]]; then
    printf '[ERROR] Missing MAC post-processing script: %s\n' "${SV_MEQTL_MAC_SCRIPT}" >&2
    exit 1
  fi
  if [[ ! -s "${SV_MEQTL_MAC_SUMMARY}" || "${SV_MEQTL_ALL}" -nt "${SV_MEQTL_MAC_SUMMARY}" ]]; then
    printf '[RUN] Refreshing MAC>=2/5 SV meQTL summaries before collection\n'
    Rscript "${SV_MEQTL_MAC_SCRIPT}" \
      --input "${SV_MEQTL_ALL}" \
      --out-dir "${MEQTL_DIR}/SV/tables" \
      --mac-thresholds 2,5 \
      --fdr-threshold 0.05 \
      --p-threshold 1e-8
  else
    printf '[INFO] MAC-filtered SV meQTL summaries are up to date\n'
  fi
else
  printf '[WARN] Cannot run MAC filtering; missing: %s\n' "${SV_MEQTL_ALL}" >&2
fi

# This is a generated summary location. Remove any SGV meQTL files copied by an
# earlier run so the SV-only package cannot silently retain stale SGV results.
[[ -n "${SUMMARY:-}" && "$SUMMARY" != "/" ]] || { echo "Refusing unsafe SUMMARY path" >&2; exit 2; }
rm -rf -- "${SUMMARY}/09_meQTL_SGV"

mkdir -p \
  "${SUMMARY}/01_GWAS" \
  "${SUMMARY}/02_leads_and_annotations" \
  "${SUMMARY}/03_LD_decay" \
  "${SUMMARY}/04_SV_detectable_LD" \
  "${SUMMARY}/05_LD_friends_optional" \
  "${SUMMARY}/06_locuszoom/summaries/SGV" \
  "${SUMMARY}/06_locuszoom/summaries/SV" \
  "${SUMMARY}/06_locuszoom/plots" \
  "${SUMMARY}/07_SV_PAV_matching" \
  "${SUMMARY}/08_meQTL_SV/tables" \
  "${SUMMARY}/08_meQTL_SV/plots" \
  "${SUMMARY}/09_meQTL_SGV/tables" \
  "${SUMMARY}/09_meQTL_SGV/plots" \
  "${SUMMARY}/99_metadata" \
  "${SUMMARY_MAIN}/01_GWAS_leads" \
  "${SUMMARY_MAIN}/02_LD" \
  "${SUMMARY_MAIN}/03_SV_meQTL" \
  "${SUMMARY_MAIN}/03_SV_meQTL/matching" \
  "${SUMMARY_MAIN}/99_metadata"

MANIFEST="${SUMMARY}/99_metadata/copy_manifest.tsv"
printf 'status\tcategory\tsource\tdestination\n' > "${MANIFEST}"

copy_one() {
  local category="$1" source="$2" destination_dir="$3"
  local destination="${destination_dir}/$(basename "${source}")"
  if [[ -f "${source}" ]]; then
    cp -f -- "${source}" "${destination}"
    printf 'COPIED\t%s\t%s\t%s\n' "${category}" "${source}" "${destination}" >> "${MANIFEST}"
  else
    printf 'MISSING\t%s\t%s\t%s\n' "${category}" "${source}" "${destination}" >> "${MANIFEST}"
    printf '[WARN] Missing: %s\n' "${source}" >&2
  fi
}

copy_top_glob() {
  local category="$1" source_dir="$2" destination_dir="$3" pattern="$4"
  local found=0 source
  while IFS= read -r -d '' source; do
    found=1
    copy_one "${category}" "${source}" "${destination_dir}"
  done < <(find "${source_dir}" -maxdepth 1 -type f -name "${pattern}" -print0 2>/dev/null || true)
  if (( found == 0 )); then
    printf 'NO_MATCH\t%s\t%s/%s\t%s\n' "${category}" "${source_dir}" "${pattern}" "${destination_dir}" >> "${MANIFEST}"
  fi
}

copy_plot_tree() {
  local category="$1" source_root="$2" destination_root="$3"
  local source relative destination
  [[ -d "${source_root}" ]] || {
    printf 'MISSING_DIR\t%s\t%s\t%s\n' "${category}" "${source_root}" "${destination_root}" >> "${MANIFEST}"
    return 0
  }
  while IFS= read -r -d '' source; do
    relative="${source#${source_root}/}"
    destination="${destination_root}/${relative}"
    mkdir -p "$(dirname "${destination}")"
    cp -f -- "${source}" "${destination}"
    printf 'COPIED\t%s\t%s\t%s\n' "${category}" "${source}" "${destination}" >> "${MANIFEST}"
  done < <(find "${source_root}" -type f \( -name '*.pdf' -o -name '*.png' -o -name '*.tiff' -o -name '*.svg' \) -print0)
}

# 01: GWAS quality summaries, significant loci, clumping summaries, and Manhattan figures.
for file in \
  lambda_summary.remove_chr_edge_1Mb.tsv \
  variant_count_summary.remove_chr_edge_1Mb.tsv \
  SV.significant.by_set.remove_chr_edge_1Mb.tsv \
  SNV_INDEL.significant.by_set.remove_chr_edge_1Mb.tsv \
  SV.significant.merged_minP.remove_chr_edge_1Mb.tsv \
  SNV_INDEL.significant.merged_minP.remove_chr_edge_1Mb.tsv; do
  copy_one GWAS "${GWAS_DIR}/${file}" "${SUMMARY}/01_GWAS"
done
copy_one GWAS "${GWAS_DIR}/clumping_by_set_subtype/clumping_summary.by_set_subtype.tsv" "${SUMMARY}/01_GWAS"
copy_one GWAS "${PARAGRAPH_DIR}/gwas_model_summary.tsv" "${SUMMARY}/01_GWAS"
copy_one GWAS "${PARAGRAPH_DIR}/04_clumping/SCZ.paragraph_test.${PARAGRAPH_MODEL}.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped" "${SUMMARY}/01_GWAS"
copy_one GWAS "${PARAGRAPH_DIR}/04_clumping/SCZ.paragraph_test.${PARAGRAPH_MODEL}.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clump_validation.audit.tsv" "${SUMMARY}/01_GWAS"
copy_one GWAS "${DEEPVAR_DIR}/04_clumping/SCZ.deepvar.mlm.clump_p1_5e-6.r2_0.01.kb_1000.clumped" "${SUMMARY}/01_GWAS"
copy_one GWAS "${FIGURE3}/sv_gwas.public.tiff" "${SUMMARY}/01_GWAS"
copy_one GWAS "${FIGURE3}/sv_snv_gwas.public.tiff" "${SUMMARY}/01_GWAS"

# 02: Lead variants and all compact annotation/window tables.
copy_top_glob LEADS_ANNOTATIONS "${LD_DIR}/tables" "${SUMMARY}/02_leads_and_annotations" 'lead_sig_from_gwas*.tsv'
copy_top_glob LEADS_ANNOTATIONS "${LD_DIR}/tables" "${SUMMARY}/02_leads_and_annotations" 'lead_sig_from_gwas*.csv'

# Create a compact locus-level table from the already merged 1 Mb-window table.
# The window is used for merging, but all window columns are omitted from the
# final compact output. The original wide files above remain untouched.
FULL_LEAD_TABLE="${LD_DIR}/tables/lead_sig_from_gwas.final_merged.cleaned.noloc.window_1000kb.tsv"
COMPACT_LEAD_TABLE="${SUMMARY}/02_leads_and_annotations/lead_sig_from_gwas.final_merged.cleaned.noloc.compact.tsv"
COMPACT_AUDIT="${SUMMARY}/02_leads_and_annotations/lead_sig_from_gwas.final_merged.cleaned.noloc.compact.audit.tsv"
MAIN_GWAS_COUNTS="${SUMMARY_MAIN}/01_GWAS_leads/maintext_gwas_lead_counts.tsv"
if [[ -f "${FULL_LEAD_TABLE}" ]]; then
  python3 - "${FULL_LEAD_TABLE}" "${COMPACT_LEAD_TABLE}" "${COMPACT_AUDIT}" "${MAIN_GWAS_COUNTS}" <<'PY'
import sys
import pandas as pd

source, compact_out, audit_out, counts_out = sys.argv[1:]
df = pd.read_csv(source, sep="\t", dtype=str)

compact_columns = [
    "source_set", "variant_type", "lead_id",
    "SNP", "CHR", "POS", "A1", "A2",
    "BETA", "SE", "T", "OR", "Z_STAT",
    "A1_FREQ", "AF1", "N", "OBS_CT", "P", "P_noSPA",
    "Func_refGene", "Gene_refGene", "GeneDetail_refGene",
    "ExonicFunc_refGene", "AAChange_refGene",
    "GWASCatalog_SCZ_related", "GWASCatalog_SCZ_related_genes",
    "GWASCatalog_SCZ_related_traits",
    "GWASCatalog_SCZ_related_association_count",
    "GWASCatalog_SCZ_related_study_accessions",
]
compact_columns = [column for column in compact_columns if column in df.columns]
compact = df.loc[:, compact_columns].copy()
compact.to_csv(compact_out, sep="\t", index=False)

duplicate_coordinate_rows = 0
coordinate_keys = [column for column in ["variant_type", "CHR", "POS"] if column in compact.columns]
if len(coordinate_keys) == 3:
    duplicate_coordinate_rows = int(compact.duplicated(coordinate_keys, keep=False).sum())
audit = pd.DataFrame([
    {"metric": "source_file", "value": source},
    {"metric": "source_data_rows", "value": len(df)},
    {"metric": "compact_data_rows", "value": len(compact)},
    {"metric": "source_unique_window_ids", "value": df["window_id"].nunique() if "window_id" in df.columns else "NA"},
    {"metric": "compact_duplicate_variant_type_chr_pos_rows", "value": duplicate_coordinate_rows},
])
audit.to_csv(audit_out, sep="\t", index=False)
if len(compact) != len(df):
    raise RuntimeError("Compact output row count differs from merged-window source")

p = pd.to_numeric(df.get("final_p"), errors="coerce")
count_rows = []
for variant_type, sub in df.assign(_p=p).groupby("variant_type", dropna=False):
    valid = sub[sub["_p"].notna()]
    genomewide = valid[valid["_p"] < 5e-8]
    suggestive_only = valid[(valid["_p"] >= 5e-8) & (valid["_p"] < 5e-6)]
    count_rows.append({
        "variant_type": variant_type,
        "n_rows_with_p": len(valid),
        "n_unique_leads_with_p": valid["lead_id"].nunique(),
        "n_genomewide_rows_p_lt_5e_8": len(genomewide),
        "n_genomewide_unique_leads_p_lt_5e_8": genomewide["lead_id"].nunique(),
        "n_suggestive_only_rows_5e_8_to_5e_6": len(suggestive_only),
        "n_suggestive_only_unique_leads_5e_8_to_5e_6": suggestive_only["lead_id"].nunique(),
    })
pd.DataFrame(count_rows).to_csv(counts_out, sep="\t", index=False)
PY
  printf 'GENERATED\tLEADS_COMPACT\t%s\t%s\n' "${FULL_LEAD_TABLE}" "${COMPACT_LEAD_TABLE}" >> "${MANIFEST}"
  copy_one MAIN_TEXT_GWAS "${COMPACT_LEAD_TABLE}" "${SUMMARY_MAIN}/01_GWAS_leads"
  copy_one MAIN_TEXT_GWAS "${COMPACT_AUDIT}" "${SUMMARY_MAIN}/01_GWAS_leads"
else
  printf 'MISSING\tLEADS_COMPACT\t%s\t%s\n' "${FULL_LEAD_TABLE}" "${COMPACT_LEAD_TABLE}" >> "${MANIFEST}"
  printf '[WARN] Cannot generate compact lead table; missing: %s\n' "${FULL_LEAD_TABLE}" >&2
fi

copy_one MAIN_TEXT_GWAS "${LD_DIR}/tables/lead_sig_from_gwas.summary.tsv" "${SUMMARY_MAIN}/01_GWAS_leads"
copy_one MAIN_TEXT_GWAS "${LD_DIR}/tables/lead_sig_from_gwas.canonical_1000kb.tsv" "${SUMMARY_MAIN}/01_GWAS_leads"
copy_one MAIN_TEXT_GWAS "${LD_DIR}/tables/lead_sig_from_gwas.canonical_1000kb.mapping.tsv" "${SUMMARY_MAIN}/01_GWAS_leads"
copy_one MAIN_TEXT_GWAS "${LD_DIR}/tables/lead_sig_from_gwas.canonical_1000kb.summary.tsv" "${SUMMARY_MAIN}/01_GWAS_leads"
copy_one MAIN_TEXT_GWAS "${LD_DIR}/tables/lead_sig_from_gwas.manual_excluded_exact_representatives.tsv" "${SUMMARY_MAIN}/01_GWAS_leads"
copy_one MAIN_TEXT_GWAS "${GWAS_DIR}/variant_count_summary.remove_chr_edge_1Mb.tsv" "${SUMMARY_MAIN}/01_GWAS_leads"

# 03: LD-decay summaries, inference tables, and publication figures.
for file in null_matching_config.tsv ld_jobs.metadata.tsv id_summary.tsv; do
  copy_one LD_DECAY "${LD_DIR}/${file}" "${SUMMARY}/03_LD_decay"
done
for file in observed_decay.tsv bootstrap_ci.tsv plot_data.tsv ld_decay_auc_observed.tsv ld_decay_auc_pairwise_diff.tsv sv_sig_vs_null.maxR2.tsv; do
  copy_one LD_DECAY "${LD_DIR}/rdata/${file}" "${SUMMARY}/03_LD_decay"
done
for file in \
  f3b_2_v2.svg f3b_2_v2.pdf f3b_2_v2.tiff \
  sv_maxR2_hist_v2.svg sv_maxR2_hist_v2.pdf sv_maxR2_hist_v2.tiff \
  sv_maxR2_pair_observed_hist_v2.svg sv_maxR2_pair_observed_hist_v2.pdf sv_maxR2_pair_observed_hist_v2.tiff \
  sv_maxR2_pair_observed_density_v2.svg sv_maxR2_pair_observed_density_v2.pdf; do
  copy_one LD_DECAY "${LD_DIR}/figures/${file}" "${SUMMARY}/03_LD_decay"
done

# 04: Prespecified detectable-LD test and all sensitivity analyses.
copy_top_glob SV_DETECTABLE_LD "${LD_DIR}/rdata" "${SUMMARY}/04_SV_detectable_LD" 'sv_sig_vs_null.proximal_snv_indel.*.tsv'
for file in \
  sv_sig_vs_null.proximal_snv_indel.summary.tsv \
  sv_sig_vs_null.proximal_snv_indel.formal_tests.tsv \
  sv_sig_vs_null.proximal_snv_indel.threshold_sensitivity.tsv \
  sv_sig_vs_null.proximal_snv_indel.distribution_tests.tsv \
  sv_sig_vs_null.proximal_snv_indel.maxR2.tsv; do
  copy_one MAIN_TEXT_LD "${LD_DIR}/rdata/${file}" "${SUMMARY_MAIN}/02_LD"
done
copy_one MAIN_TEXT_LD "${LD_DIR}/rdata/observed_decay.tsv" "${SUMMARY_MAIN}/02_LD"
copy_one MAIN_TEXT_LD "${LD_DIR}/rdata/ld_decay_auc_pairwise_diff.tsv" "${SUMMARY_MAIN}/02_LD"
copy_one MAIN_TEXT_LD "${LD_DIR}/null_matching_config.tsv" "${SUMMARY_MAIN}/02_LD"
copy_one MAIN_TEXT_LD "${LD_DIR}/figures/f3b_2_v2.pdf" "${SUMMARY_MAIN}/02_LD"
copy_one MAIN_TEXT_LD "${LD_DIR}/figures/sv_maxR2_hist_v2.pdf" "${SUMMARY_MAIN}/02_LD"

# 05: Optional LD-friends/LD-score analysis (pipeline step 8).
copy_top_glob LD_FRIENDS "${LD_DIR}/ld_friends_scores" "${SUMMARY}/05_LD_friends_optional" 'ld_friends_scores.*.tsv'
copy_top_glob LD_FRIENDS "${LD_DIR}/figures" "${SUMMARY}/05_LD_friends_optional" 'ld_friends_scores_sig_vs_null.*.pdf'
copy_top_glob LD_FRIENDS "${LD_DIR}/figures" "${SUMMARY}/05_LD_friends_optional" 'ld_friends_scores_sig_vs_null.*.tiff'
copy_top_glob MAIN_TEXT_LD_FRIENDS "${LD_DIR}/ld_friends_scores" "${SUMMARY_MAIN}/02_LD" 'ld_friends_scores.*.tsv'
copy_top_glob MAIN_TEXT_LD_FRIENDS "${LD_DIR}/figures" "${SUMMARY_MAIN}/02_LD" 'ld_friends_scores_sig_vs_null.*.pdf'
copy_top_glob MAIN_TEXT_LD_FRIENDS "${LD_DIR}/figures" "${SUMMARY_MAIN}/02_LD" 'ld_friends_scores_sig_vs_null.*.tiff'

# 06: Locuszoom job/plot summaries and final regional plots; raw LD is not duplicated.
copy_top_glob LOCUSZOOM "${LD_DIR}/locuszoom" "${SUMMARY}/06_locuszoom/summaries" 'locuszoom*.tsv'
copy_one LOCUSZOOM "${LD_DIR}/locuszoom/SGV/locuszoom_plot_summary.tsv" "${SUMMARY}/06_locuszoom/summaries/SGV"
copy_one LOCUSZOOM "${LD_DIR}/locuszoom/SV/locuszoom_plot_summary.tsv" "${SUMMARY}/06_locuszoom/summaries/SV"
copy_plot_tree LOCUSZOOM "${LD_DIR}/locuszoom/SGV/plots" "${SUMMARY}/06_locuszoom/plots/SGV"
copy_plot_tree LOCUSZOOM "${LD_DIR}/locuszoom/SV/plots" "${SUMMARY}/06_locuszoom/plots/SV"

# 07: Exact matching of significant SV leads to the PAV VCF.
copy_top_glob SV_PAV_MATCH "${MEQTL_DIR}/SV/matching" "${SUMMARY}/07_SV_PAV_matching" 'sig_sv_to_pav_sv_len50.*.tsv'
copy_top_glob MAIN_TEXT_SV_PAV_MATCH "${MEQTL_DIR}/SV/matching" "${SUMMARY_MAIN}/03_SV_meQTL/matching" 'sig_sv_to_pav_sv_len50.*.tsv'

# 08: SV-only meQTL result tables and plots; all_results/per_lead stay in place.
for file in meqtl_cis.summary.tsv meqtl_cis.error_summary.tsv meqtl_cis.top_per_lead.tsv meqtl_cis.global_fdr05.tsv.gz meqtl_selected_matches.tsv meqtl_cis.genotype_policy.tsv meqtl_cis.mac_filter_summary.tsv meqtl_cis.per_sv_mac.tsv; do
  [[ -f "${MEQTL_DIR}/SV/tables/${file}" ]] && copy_one MEQTL_SV "${MEQTL_DIR}/SV/tables/${file}" "${SUMMARY}/08_meQTL_SV/tables"
done
copy_top_glob MEQTL_SV_DEBUG "${MEQTL_DIR}/SV/tables" "${SUMMARY}/08_meQTL_SV/tables" 'sv_meqtl.genotype_debug*.tsv*'
copy_top_glob MEQTL_SV_MAC_HITS "${MEQTL_DIR}/SV/tables" "${SUMMARY}/08_meQTL_SV/tables" 'meqtl_cis.mac_ge*.fdr_lt_0.05.tsv.gz'
copy_top_glob MEQTL_SV_MAC_HITS "${MEQTL_DIR}/SV/tables" "${SUMMARY}/08_meQTL_SV/tables" 'meqtl_cis.mac_ge*.p_lt_1e-8.tsv.gz'
copy_top_glob MEQTL_SV_MAC_FULL "${MEQTL_DIR}/SV/tables" "${SUMMARY}/08_meQTL_SV/tables" 'meqtl_cis.mac_ge*.all_results.tsv.gz'
copy_plot_tree MEQTL_SV "${MEQTL_DIR}/SV/plots" "${SUMMARY}/08_meQTL_SV/plots"
copy_one MEQTL_POLICY "${MEQTL_DIR}/SGV/sgv_meqtl.disabled.tsv" "${SUMMARY}/09_meQTL_SGV"

for file in meqtl_cis.summary.tsv meqtl_cis.top_per_lead.tsv meqtl_cis.global_fdr05.tsv.gz meqtl_selected_matches.tsv meqtl_cis.genotype_policy.tsv meqtl_cis.mac_filter_summary.tsv meqtl_cis.per_sv_mac.tsv; do
  copy_one MAIN_TEXT_SV_MEQTL "${MEQTL_DIR}/SV/tables/${file}" "${SUMMARY_MAIN}/03_SV_meQTL"
done
copy_top_glob MAIN_TEXT_SV_MEQTL_DEBUG "${MEQTL_DIR}/SV/tables" "${SUMMARY_MAIN}/03_SV_meQTL" 'sv_meqtl.genotype_debug*.tsv*'
copy_top_glob MAIN_TEXT_SV_MEQTL_MAC_HITS "${MEQTL_DIR}/SV/tables" "${SUMMARY_MAIN}/03_SV_meQTL" 'meqtl_cis.mac_ge*.fdr_lt_0.05.tsv.gz'
copy_top_glob MAIN_TEXT_SV_MEQTL_MAC_HITS "${MEQTL_DIR}/SV/tables" "${SUMMARY_MAIN}/03_SV_meQTL" 'meqtl_cis.mac_ge*.p_lt_1e-8.tsv.gz'

# Documentation and reproducibility metadata.
copy_one DOCUMENTATION "${SCRIPT_DIR}/SUMMARY_readme.md" "${SUMMARY}"
if [[ -f "${SUMMARY}/SUMMARY_readme.md" ]]; then
  mv -f "${SUMMARY}/SUMMARY_readme.md" "${SUMMARY}/readme.md"
  sed -i "s#${SUMMARY}/SUMMARY_readme.md#${SUMMARY}/readme.md#g" "${MANIFEST}" || true
fi
copy_one DOCUMENTATION "${SCRIPT_DIR}/run_pipeline.sh" "${SUMMARY}/99_metadata"
copy_one DOCUMENTATION "${SCRIPT_DIR}/readme.md" "${SUMMARY}/99_metadata"
copy_one DOCUMENTATION "${SCRIPT_DIR}/RUNNING.md" "${SUMMARY}/99_metadata"
copy_one DOCUMENTATION "${SCRIPT_DIR}/SUMMARY_MAIN_TEXT_readme.md" "${SUMMARY_MAIN}"
if [[ -f "${SUMMARY_MAIN}/SUMMARY_MAIN_TEXT_readme.md" ]]; then
  mv -f "${SUMMARY_MAIN}/SUMMARY_MAIN_TEXT_readme.md" "${SUMMARY_MAIN}/readme.md"
fi

INVENTORY="${SUMMARY}/99_metadata/summary_inventory.tsv"
printf 'relative_path\tsize_bytes\tsha256\n' > "${INVENTORY}"
while IFS= read -r -d '' file; do
  relative="${file#${SUMMARY}/}"
  size="$(stat -c '%s' "${file}")"
  checksum="$(sha256sum "${file}" | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "${relative}" "${size}" "${checksum}" >> "${INVENTORY}"
done < <(find "${SUMMARY}" -type f ! -path "${INVENTORY}" -print0 | sort -z)

MAIN_INVENTORY="${SUMMARY_MAIN}/99_metadata/summary_maintext_inventory.tsv"
printf 'relative_path\tsize_bytes\tsha256\n' > "${MAIN_INVENTORY}"
while IFS= read -r -d '' file; do
  relative="${file#${SUMMARY_MAIN}/}"
  size="$(stat -c '%s' "${file}")"
  checksum="$(sha256sum "${file}" | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "${relative}" "${size}" "${checksum}" >> "${MAIN_INVENTORY}"
done < <(find "${SUMMARY_MAIN}" -type f ! -path "${MAIN_INVENTORY}" -print0 | sort -z)

copied="$(awk -F '\t' '$1=="COPIED"{n++} END{print n+0}' "${MANIFEST}")"
missing="$(awk -F '\t' '$1=="MISSING" || $1=="MISSING_DIR" || $1=="NO_MATCH"{n++} END{print n+0}' "${MANIFEST}")"
printf '[DONE] Summary directory: %s\n' "${SUMMARY}"
printf '[DONE] Main-text package: %s\n' "${SUMMARY_MAIN}"
printf '[DONE] Copied entries: %s; missing/no-match entries: %s\n' "${copied}" "${missing}"
printf '[INFO] Inspect: %s and %s\n' "${MANIFEST}" "${INVENTORY}"
