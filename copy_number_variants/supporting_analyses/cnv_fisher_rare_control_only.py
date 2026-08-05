# -*- coding: utf-8 -*-

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


OUTPUT_DIR = Path(
    os.environ.get(
        "EOSCZ_CNV_SUPPORT_OUTPUT_DIR",
        Path(__file__).resolve().parent / "outputs",
    )
)
HAPLOTYPE_INPUT = OUTPUT_DIR / "true_case_cohort_comparison_cohort_public_reference_haplotype_pretest_merge.xlsx"
LEGACY_HAPLOTYPE_INPUT = OUTPUT_DIR / "true_case_cohort_comparison_cohort_public_reference_haplotype_merge_pre_test.xlsx"
SAMPLE_HELPER = Path(
    os.environ.get(
        "EOSCZ_CNV_SAMPLE_HELPER_SCRIPT",
        Path(__file__).resolve().with_name("sample_cn_vs_haplotype_cn_qc_and_sample_level_common_fisher.py"),
    )
)

HAPLOTYPE_PREFIX = "true_case_cohort_comparison_cohort_public_reference_haplotype_rare_comparison_cohortlt1pct_public_referencelt1pct_noSCZrecurrent"
SAMPLE_PREFIX = "true_case_cohort_comparison_cohort_public_reference_sample_rare_comparison_cohortlt1pct_public_referencelt1pct_noSCZrecurrent"


def load_sample_helper():
    spec = importlib.util.spec_from_file_location("sample_helper", SAMPLE_HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import helper script: {SAMPLE_HELPER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_existing_path(primary: Path, legacy: Path) -> Path:
    if primary.exists():
        return primary
    if legacy.exists():
        return legacy
    raise FileNotFoundError(f"Input table not found. Tried: {primary}; {legacy}")


def fisher_p_and_or(case_count: int, case_other: int, ctrl_count: int, ctrl_other: int) -> tuple[float, float]:
    odds_ratio, p_value = fisher_exact(
        np.array([[case_count, case_other], [ctrl_count, ctrl_other]]),
        alternative="greater",
    )
    return float(p_value), float(odds_ratio)


def add_adjustments(table: pd.DataFrame) -> pd.DataFrame:
    if table.empty:
        table["FDR"] = []
        table["Bonferroni"] = []
        return table
    p_values = table["P_value"].to_numpy(dtype=float)
    table["FDR"] = multipletests(p_values, method="fdr_bh")[1]
    table["Bonferroni"] = multipletests(p_values, method="bonferroni")[1]
    return table


def write_result(prefix: str, result: pd.DataFrame, summary: pd.DataFrame) -> None:
    raw_sig = result[result["P_value"] < 0.05].copy()
    fdr_sig = result[result["FDR"] < 0.05].copy()
    bonf_sig = result[result["Bonferroni"] < 0.05].copy()
    scz_higher = raw_sig[raw_sig["Frequency_Scz"] > raw_sig["Frequency_comparison_cohort+public_reference"]].copy()

    merge_out = OUTPUT_DIR / f"{prefix}_tested_genes.xlsx"
    raw_out = OUTPUT_DIR / f"{prefix}_nominal_p_lt_0.05_all.xlsx"
    scz_higher_out = OUTPUT_DIR / f"{prefix}_nominal_p_lt_0.05_scz_higher_frequency.xlsx"
    fdr_out = OUTPUT_DIR / f"{prefix}_fdr_lt_0.05.xlsx"
    bonf_out = OUTPUT_DIR / f"{prefix}_bonferroni_lt_0.05.xlsx"
    summary_out = OUTPUT_DIR / f"{prefix}_summary.xlsx"

    result.to_excel(merge_out, index=False)
    raw_sig.to_excel(raw_out, index=False)
    scz_higher.to_excel(scz_higher_out, index=False)
    fdr_sig.to_excel(fdr_out, index=False)
    bonf_sig.to_excel(bonf_out, index=False)
    with pd.ExcelWriter(summary_out, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        result.head(100).to_excel(writer, sheet_name="Top100_by_P", index=False)
        raw_sig.to_excel(writer, sheet_name="RawP_lt_0.05", index=False)
        fdr_sig.to_excel(writer, sheet_name="FDR_lt_0.05", index=False)
        bonf_sig.to_excel(writer, sheet_name="Bonferroni_lt_0.05", index=False)

    print(f"\n[{prefix}]")
    print(summary.to_string(index=False))
    print("\nTop 30:")
    print(result.head(30).to_string(index=False))
    if fdr_sig.empty:
        print("No FDR significant genes.")
    else:
        print("\nFDR significant genes:")
        print(fdr_sig.to_string(index=False))
    for path in [merge_out, raw_out, scz_higher_out, fdr_out, bonf_out, summary_out]:
        print(f"Saved: {path}")


def run_haplotype_level() -> pd.DataFrame:
    table = pd.read_excel(resolve_existing_path(HAPLOTYPE_INPUT, LEGACY_HAPLOTYPE_INPUT))
    table = table.rename(
        columns={
            "P_value": "P_value_all_tested_genes",
            "FDR": "FDR_all_tested_genes",
            "Bonferroni": "Bonferroni_all_tested_genes",
            "Odds Ratio": "Odds_Ratio_all_tested_genes",
        }
    )

    rare = table[(table["Frequency_comparison_cohort"] < 0.01) & (table["Frequency_public_reference"] < 0.01)].copy()
    rare[["P_value", "Odds Ratio"]] = rare.apply(
        lambda row: pd.Series(
            fisher_p_and_or(
                int(row["Count_Scz"]),
                int(row["SCZ_other"]),
                int(row["Count_comparison_cohort+public_reference"]),
                int(row["comparison_cohort+public_reference_other"]),
            )
        ),
        axis=1,
    )
    rare = add_adjustments(rare)
    rare["Rare_definition"] = "Frequency_comparison_cohort < 0.01 and Frequency_public_reference < 0.01; no SCZ recurrent filter"
    rare = rare.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

    summary = pd.DataFrame(
        [
            {"Metric": "level", "Value": "haplotype"},
            {"Metric": "input_all_tested_genes", "Value": len(table)},
            {"Metric": "rare_genes", "Value": len(rare)},
            {"Metric": "SCZ_total_haplotypes", "Value": int((table["Count_Scz"] + table["SCZ_other"]).iloc[0])},
            {"Metric": "comparison_cohort_total_haplotypes", "Value": 304},
            {"Metric": "public_reference_total_haplotypes", "Value": 727},
            {"Metric": "comparison_cohort_public_reference_total_haplotypes", "Value": int((table["Count_comparison_cohort+public_reference"] + table["comparison_cohort+public_reference_other"]).iloc[0])},
            {"Metric": "rare_definition", "Value": "Frequency_comparison_cohort < 0.01 and Frequency_public_reference < 0.01; no SCZ recurrent filter"},
            {"Metric": "rare_genes_with_SCZ_count_0", "Value": int((rare["Count_Scz"] == 0).sum())},
            {"Metric": "rare_genes_with_SCZ_count_1", "Value": int((rare["Count_Scz"] == 1).sum())},
            {"Metric": "rare_genes_with_SCZ_count_ge2", "Value": int((rare["Count_Scz"] >= 2).sum())},
            {"Metric": "raw_p_lt_0.05", "Value": int((rare["P_value"] < 0.05).sum())},
            {"Metric": "FDR_lt_0.05", "Value": int((rare["FDR"] < 0.05).sum())},
            {"Metric": "Bonferroni_lt_0.05", "Value": int((rare["Bonferroni"] < 0.05).sum())},
            {"Metric": "min_p_value", "Value": rare["P_value"].min()},
            {"Metric": "min_FDR", "Value": rare["FDR"].min()},
        ]
    )

    columns = [
        "Gene",
        "Count_Scz",
        "Frequency_Scz",
        "Count_comparison_cohort+public_reference",
        "Frequency_comparison_cohort+public_reference",
        "Count_comparison_cohort",
        "Frequency_comparison_cohort",
        "Count_public_reference",
        "Frequency_public_reference",
        "SCZ_other",
        "comparison_cohort+public_reference_other",
        "P_value",
        "FDR",
        "Bonferroni",
        "Odds Ratio",
        "P_value_all_tested_genes",
        "FDR_all_tested_genes",
        "Bonferroni_all_tested_genes",
        "Rare_definition",
    ]
    result = rare[columns]
    write_result(HAPLOTYPE_PREFIX, result, summary)
    return result


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
    result["Count_Scz"] = [int(scz_binary[gene].sum()) for gene in genes]
    result["Frequency_Scz"] = result["Count_Scz"] / scz_binary.shape[0]
    result["Count_comparison_cohort"] = [int(comparison_cohort_binary[gene].sum()) for gene in genes]
    result["Frequency_comparison_cohort"] = result["Count_comparison_cohort"] / comparison_cohort_binary.shape[0]
    result["Count_public_reference"] = [int(public_reference_binary[gene].sum()) for gene in genes]
    result["Frequency_public_reference"] = result["Count_public_reference"] / public_reference_binary.shape[0]
    result["Count_comparison_cohort+public_reference"] = [int(ctrl_binary[gene].sum()) for gene in genes]
    result["Frequency_comparison_cohort+public_reference"] = result["Count_comparison_cohort+public_reference"] / ctrl_binary.shape[0]
    result["SCZ_other"] = scz_binary.shape[0] - result["Count_Scz"]
    result["comparison_cohort+public_reference_other"] = ctrl_binary.shape[0] - result["Count_comparison_cohort+public_reference"]
    result.attrs["SCZ_samples"] = scz_binary.shape[0]
    result.attrs["comparison_cohort_samples"] = comparison_cohort_binary.shape[0]
    result.attrs["public_reference_samples"] = public_reference_binary.shape[0]
    result.attrs["comparison_cohort_public_reference_samples"] = ctrl_binary.shape[0]
    return result


def run_sample_level() -> pd.DataFrame:
    helper = load_sample_helper()
    rebuilt_tables, qc_summary, hap_parse_summary, mismatch_examples, _ = helper.cohort_qc()
    table = build_sample_level_gene_table(rebuilt_tables)

    rare = table[(table["Frequency_comparison_cohort"] < 0.01) & (table["Frequency_public_reference"] < 0.01)].copy()
    rare[["P_value", "Odds Ratio"]] = rare.apply(
        lambda row: pd.Series(
            fisher_p_and_or(
                int(row["Count_Scz"]),
                int(row["SCZ_other"]),
                int(row["Count_comparison_cohort+public_reference"]),
                int(row["comparison_cohort+public_reference_other"]),
            )
        ),
        axis=1,
    )
    rare = add_adjustments(rare)
    rare["Rare_definition"] = "Frequency_comparison_cohort_sample < 0.01 and Frequency_public_reference_sample < 0.01; no SCZ recurrent filter"
    rare = rare.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

    summary = pd.DataFrame(
        [
            {"Metric": "level", "Value": "sample"},
            {"Metric": "input_all_tested_genes", "Value": len(table)},
            {"Metric": "rare_genes", "Value": len(rare)},
            {"Metric": "SCZ_samples", "Value": table.attrs["SCZ_samples"]},
            {"Metric": "comparison_cohort_samples", "Value": table.attrs["comparison_cohort_samples"]},
            {"Metric": "public_reference_samples", "Value": table.attrs["public_reference_samples"]},
            {"Metric": "comparison_cohort_public_reference_samples", "Value": table.attrs["comparison_cohort_public_reference_samples"]},
            {"Metric": "rare_definition", "Value": "Frequency_comparison_cohort_sample < 0.01 and Frequency_public_reference_sample < 0.01; no SCZ recurrent filter"},
            {"Metric": "rare_genes_with_SCZ_count_0", "Value": int((rare["Count_Scz"] == 0).sum())},
            {"Metric": "rare_genes_with_SCZ_count_1", "Value": int((rare["Count_Scz"] == 1).sum())},
            {"Metric": "rare_genes_with_SCZ_count_ge2", "Value": int((rare["Count_Scz"] >= 2).sum())},
            {"Metric": "raw_p_lt_0.05", "Value": int((rare["P_value"] < 0.05).sum())},
            {"Metric": "FDR_lt_0.05", "Value": int((rare["FDR"] < 0.05).sum())},
            {"Metric": "Bonferroni_lt_0.05", "Value": int((rare["Bonferroni"] < 0.05).sum())},
            {"Metric": "min_p_value", "Value": rare["P_value"].min()},
            {"Metric": "min_FDR", "Value": rare["FDR"].min()},
        ]
    )

    columns = [
        "Gene",
        "Count_Scz",
        "Frequency_Scz",
        "Count_comparison_cohort+public_reference",
        "Frequency_comparison_cohort+public_reference",
        "Count_comparison_cohort",
        "Frequency_comparison_cohort",
        "Count_public_reference",
        "Frequency_public_reference",
        "SCZ_other",
        "comparison_cohort+public_reference_other",
        "P_value",
        "FDR",
        "Bonferroni",
        "Odds Ratio",
        "Rare_definition",
    ]
    result = rare[columns]

    summary_out = OUTPUT_DIR / f"{SAMPLE_PREFIX}_summary.xlsx"
    write_result(SAMPLE_PREFIX, result, summary)
    with pd.ExcelWriter(summary_out, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        qc_summary.to_excel(writer, sheet_name="QC_sample_CN", index=False)
        hap_parse_summary.to_excel(writer, sheet_name="QC_haplotype_parse", index=False)
        if not mismatch_examples.empty:
            mismatch_examples.to_excel(writer, sheet_name="QC_mismatch_examples", index=False)
    return result


def verify_fdr(label: str, result: pd.DataFrame) -> None:
    if result.empty:
        print(f"{label}: empty result")
        return
    recalculated = multipletests(result["P_value"].to_numpy(dtype=float), method="fdr_bh")[1]
    max_diff = float(np.max(np.abs(result["FDR"].to_numpy(dtype=float) - recalculated)))
    print(f"{label}: FDR max abs diff vs statsmodels = {max_diff}")


def main() -> None:
    hap_result = run_haplotype_level()
    sample_result = run_sample_level()
    verify_fdr("haplotype", hap_result)
    verify_fdr("sample", sample_result)


if __name__ == "__main__":
    main()

