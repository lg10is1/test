# -*- coding: utf-8 -*-

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

OUTPUT_DIR = Path(os.environ.get("EOSCZ_CNV_SUPPORT_OUTPUT_DIR", Path(__file__).resolve().parent / "outputs"))
HELPER_SCRIPT = Path(
    os.environ.get(
        "EOSCZ_CNV_SAMPLE_HELPER_SCRIPT",
        Path(__file__).resolve().with_name("sample_cn_vs_haplotype_cn_qc_and_sample_level_common_fisher.py"),
    )
)
OUTPUT_PREFIX = "true_case_cohort_comparison_cohort_public_reference_sample_rare_comparison_cohortlt1pct_public_referencelt1pct_SCZrecurrent2"
RARE_DEFINITION = "Frequency_comparison_cohort_sample < 0.01 and Frequency_public_reference_sample < 0.01 and Count_Scz_sample >= 2"


def load_helper_module():
    spec = importlib.util.spec_from_file_location("sample_qc_helper", HELPER_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import helper script: {HELPER_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_sample_level_gene_table(rebuilt_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    scz = rebuilt_tables["SCZ_filtered_cases"]
    comparison_cohort = rebuilt_tables["comparison_cohort"]
    public_reference = rebuilt_tables["public_reference"]

    genes = sorted(set(scz.columns) | set(comparison_cohort.columns) | set(public_reference.columns))
    scz = scz.reindex(columns=genes, fill_value=0)
    comparison_cohort = comparison_cohort.reindex(columns=genes, fill_value=0)
    public_reference = public_reference.reindex(columns=genes, fill_value=0)
    ctrl = pd.concat([comparison_cohort, public_reference], axis=0)

    scz_binary = (scz > 0).astype(int)
    comparison_cohort_binary = (comparison_cohort > 0).astype(int)
    public_reference_binary = (public_reference > 0).astype(int)
    ctrl_binary = (ctrl > 0).astype(int)

    result = pd.DataFrame({"Gene": genes})
    result["Count_Scz_sample"] = [int(scz_binary[gene].sum()) for gene in genes]
    result["Frequency_Scz_sample"] = result["Count_Scz_sample"] / scz_binary.shape[0]
    result["Count_comparison_cohort_sample"] = [int(comparison_cohort_binary[gene].sum()) for gene in genes]
    result["Frequency_comparison_cohort_sample"] = result["Count_comparison_cohort_sample"] / comparison_cohort_binary.shape[0]
    result["Count_public_reference_sample"] = [int(public_reference_binary[gene].sum()) for gene in genes]
    result["Frequency_public_reference_sample"] = result["Count_public_reference_sample"] / public_reference_binary.shape[0]
    result["Count_comparison_cohort+public_reference_sample"] = [int(ctrl_binary[gene].sum()) for gene in genes]
    result["Frequency_comparison_cohort+public_reference_sample"] = result["Count_comparison_cohort+public_reference_sample"] / ctrl_binary.shape[0]
    result["SCZ_other_sample"] = scz_binary.shape[0] - result["Count_Scz_sample"]
    result["comparison_cohort_other_sample"] = comparison_cohort_binary.shape[0] - result["Count_comparison_cohort_sample"]
    result["public_reference_other_sample"] = public_reference_binary.shape[0] - result["Count_public_reference_sample"]
    result["comparison_cohort+public_reference_other_sample"] = ctrl_binary.shape[0] - result["Count_comparison_cohort+public_reference_sample"]
    result.attrs["SCZ_samples"] = scz_binary.shape[0]
    result.attrs["comparison_cohort_samples"] = comparison_cohort_binary.shape[0]
    result.attrs["public_reference_samples"] = public_reference_binary.shape[0]
    result.attrs["comparison_cohort_public_reference_samples"] = ctrl_binary.shape[0]
    return result


def fisher_p_and_or(row: pd.Series) -> tuple[float, float]:
    contingency = np.array(
        [
            [int(row["Count_Scz_sample"]), int(row["SCZ_other_sample"])],
            [int(row["Count_comparison_cohort+public_reference_sample"]), int(row["comparison_cohort+public_reference_other_sample"])],
        ]
    )
    odds_ratio, p_value = fisher_exact(contingency, alternative="greater")
    return float(p_value), float(odds_ratio)


def run_sample_rare_fisher(rebuilt_tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_genes = build_sample_level_gene_table(rebuilt_tables)
    rare = all_genes[
        (all_genes["Frequency_comparison_cohort_sample"] < 0.01)
        & (all_genes["Frequency_public_reference_sample"] < 0.01)
        & (all_genes["Count_Scz_sample"] >= 2)
    ].copy()

    if len(rare):
        p_or = rare.apply(fisher_p_and_or, axis=1, result_type="expand")
        rare["P_value"] = p_or[0]
        rare["Odds Ratio"] = p_or[1]
        rare["FDR"] = multipletests(rare["P_value"].to_numpy(dtype=float), method="fdr_bh")[1]
        rare["Bonferroni"] = multipletests(rare["P_value"].to_numpy(dtype=float), method="bonferroni")[1]
    else:
        rare["P_value"] = []
        rare["Odds Ratio"] = []
        rare["FDR"] = []
        rare["Bonferroni"] = []

    rare["Rare_definition"] = RARE_DEFINITION
    rare = rare.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

    raw_sig = rare[rare["P_value"] < 0.05].copy()
    fdr_sig = rare[rare["FDR"] < 0.05].copy()
    bonf_sig = rare[rare["Bonferroni"] < 0.05].copy()

    summary = pd.DataFrame(
        [
            {"Metric": "SCZ_samples", "Value": all_genes.attrs["SCZ_samples"]},
            {"Metric": "comparison_cohort_samples", "Value": all_genes.attrs["comparison_cohort_samples"]},
            {"Metric": "public_reference_samples", "Value": all_genes.attrs["public_reference_samples"]},
            {"Metric": "comparison_cohort_public_reference_samples", "Value": all_genes.attrs["comparison_cohort_public_reference_samples"]},
            {"Metric": "all_tested_genes_union", "Value": len(all_genes)},
            {"Metric": "rare_genes", "Value": len(rare)},
            {"Metric": "rare_definition", "Value": RARE_DEFINITION},
            {"Metric": "raw_p_lt_0.05", "Value": len(raw_sig)},
            {"Metric": "FDR_lt_0.05", "Value": len(fdr_sig)},
            {"Metric": "Bonferroni_lt_0.05", "Value": len(bonf_sig)},
            {"Metric": "SCZ_higher_raw_p_lt_0.05", "Value": int((raw_sig["Frequency_Scz_sample"] > raw_sig["Frequency_comparison_cohort+public_reference_sample"]).sum())},
            {"Metric": "min_p_value", "Value": rare["P_value"].min() if len(rare) else np.nan},
            {"Metric": "min_FDR", "Value": rare["FDR"].min() if len(rare) else np.nan},
            {"Metric": "min_Bonferroni", "Value": rare["Bonferroni"].min() if len(rare) else np.nan},
        ]
    )
    return rare, raw_sig, fdr_sig, bonf_sig, summary


def main() -> None:
    helper = load_helper_module()
    rebuilt_tables, qc_summary, hap_parse_summary, mismatch_examples, hap_group_counts = helper.cohort_qc()
    rare, raw_sig, fdr_sig, bonf_sig, summary = run_sample_rare_fisher(rebuilt_tables)

    merge_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_tested_genes.xlsx"
    raw_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_nominal_p_lt_0.05_all.xlsx"
    scz_higher_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_nominal_p_lt_0.05_scz_higher_frequency.xlsx"
    fdr_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_fdr_lt_0.05.xlsx"
    bonf_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_bonferroni_lt_0.05.xlsx"
    summary_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_summary.xlsx"

    scz_higher = raw_sig[raw_sig["Frequency_Scz_sample"] > raw_sig["Frequency_comparison_cohort+public_reference_sample"]].copy()
    rare.to_excel(merge_out, index=False)
    raw_sig.to_excel(raw_out, index=False)
    scz_higher.to_excel(scz_higher_out, index=False)
    fdr_sig.to_excel(fdr_out, index=False)
    bonf_sig.to_excel(bonf_out, index=False)
    with pd.ExcelWriter(summary_out, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        rare.head(100).to_excel(writer, sheet_name="Top100_by_P", index=False)
        raw_sig.to_excel(writer, sheet_name="RawP_lt_0.05", index=False)
        fdr_sig.to_excel(writer, sheet_name="FDR_lt_0.05", index=False)
        bonf_sig.to_excel(writer, sheet_name="Bonferroni_lt_0.05", index=False)
        qc_summary.to_excel(writer, sheet_name="QC_sample_CN", index=False)
        hap_parse_summary.to_excel(writer, sheet_name="QC_haplotype_parse", index=False)
        if not mismatch_examples.empty:
            mismatch_examples.to_excel(writer, sheet_name="QC_mismatch_examples", index=False)

    print("Sample-level rare Fisher summary:")
    print(summary.to_string(index=False))
    print("\nTop 30 sample-level rare Fisher results:")
    cols = [
        "Gene",
        "Count_Scz_sample",
        "Frequency_Scz_sample",
        "Count_comparison_cohort+public_reference_sample",
        "Frequency_comparison_cohort+public_reference_sample",
        "Count_comparison_cohort_sample",
        "Frequency_comparison_cohort_sample",
        "Count_public_reference_sample",
        "Frequency_public_reference_sample",
        "P_value",
        "FDR",
        "Bonferroni",
        "Odds Ratio",
    ]
    print(rare[cols].head(30).to_string(index=False))
    if len(fdr_sig):
        print("\nFDR significant genes:")
        print(fdr_sig[cols].to_string(index=False))
    else:
        print("\nNo FDR significant genes in the sample-level restricted rare-gene Fisher test.")

    print("\nQC sample_CN summary:")
    print(qc_summary.to_string(index=False))
    for path in [merge_out, raw_out, scz_higher_out, fdr_out, bonf_out, summary_out]:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()

