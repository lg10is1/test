# -*- coding: utf-8 -*-

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

OUTPUT_DIR = Path(os.environ.get("EOSCZ_CNV_SUPPORT_OUTPUT_DIR", Path(__file__).resolve().parent / "outputs"))
INPUT_MERGE = OUTPUT_DIR / "true_case_cohort_comparison_cohort_public_reference_haplotype_pretest_merge.xlsx"
LEGACY_INPUT_MERGE = OUTPUT_DIR / "true_case_cohort_comparison_cohort_public_reference_haplotype_merge_pre_test.xlsx"
OUTPUT_PREFIX = "true_case_cohort_comparison_cohort_public_reference_haplotype_rare_comparison_cohortlt1pct_public_referencelt1pct_SCZrecurrent2"
RARE_DEFINITION = "Frequency_comparison_cohort < 0.01 and Frequency_public_reference < 0.01 and Count_Scz >= 2"


def resolve_existing_path(primary: Path, legacy: Path) -> Path:
    if primary.exists():
        return primary
    if legacy.exists():
        return legacy
    raise FileNotFoundError(f"Input table not found. Tried: {primary}; {legacy}")


def fisher_p_and_or(row: pd.Series) -> tuple[float, float]:
    contingency = np.array(
        [
            [int(row["Count_Scz"]), int(row["SCZ_other"])],
            [int(row["Count_comparison_cohort+public_reference"]), int(row["comparison_cohort+public_reference_other"])],
        ]
    )
    odds_ratio, p_value = fisher_exact(contingency, alternative="greater")
    return float(p_value), float(odds_ratio)


def main() -> None:
    merge = pd.read_excel(resolve_existing_path(INPUT_MERGE, LEGACY_INPUT_MERGE))
    merge = merge.rename(
        columns={
            "P_value": "P_value_all_tested_genes",
            "FDR": "FDR_all_tested_genes",
            "Bonferroni": "Bonferroni_all_tested_genes",
            "Odds Ratio": "Odds_Ratio_all_tested_genes",
        }
    )

    scz_totals = sorted((merge["Count_Scz"] + merge["SCZ_other"]).dropna().astype(int).unique().tolist())
    ctrl_totals = sorted((merge["Count_comparison_cohort+public_reference"] + merge["comparison_cohort+public_reference_other"]).dropna().astype(int).unique().tolist())
    if scz_totals != [686]:
        raise ValueError(f"Unexpected SCZ totals: {scz_totals}")
    if ctrl_totals != [1031]:
        raise ValueError(f"Unexpected comparison_cohort+public_reference totals: {ctrl_totals}")

    for column in ["Frequency_comparison_cohort", "Frequency_public_reference", "Count_Scz"]:
        merge[column] = pd.to_numeric(merge[column], errors="coerce").fillna(0)

    rare = merge[
        (merge["Frequency_comparison_cohort"] < 0.01)
        & (merge["Frequency_public_reference"] < 0.01)
        & (merge["Count_Scz"] >= 2)
    ].copy()

    p_or = rare.apply(fisher_p_and_or, axis=1, result_type="expand")
    rare["P_value"] = p_or[0]
    rare["Odds Ratio"] = p_or[1]
    rare["FDR"] = multipletests(rare["P_value"].to_numpy(dtype=float), method="fdr_bh")[1]
    rare["Bonferroni"] = multipletests(rare["P_value"].to_numpy(dtype=float), method="bonferroni")[1]
    rare["Rare_definition"] = RARE_DEFINITION
    rare = rare.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

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
        "Odds_Ratio_all_tested_genes",
        "Rare_definition",
    ]
    rare = rare[columns]

    raw_sig = rare[rare["P_value"] < 0.05].copy()
    fdr_sig = rare[rare["FDR"] < 0.05].copy()
    bonf_sig = rare[rare["Bonferroni"] < 0.05].copy()
    scz_higher = raw_sig[raw_sig["Frequency_Scz"] > raw_sig["Frequency_comparison_cohort+public_reference"]].copy()

    summary = pd.DataFrame(
        [
            {"Metric": "input_all_tested_genes", "Value": len(merge)},
            {"Metric": "rare_genes", "Value": len(rare)},
            {"Metric": "SCZ_total_haplotypes", "Value": 686},
            {"Metric": "comparison_cohort_total_haplotypes", "Value": 304},
            {"Metric": "public_reference_total_haplotypes", "Value": 727},
            {"Metric": "comparison_cohort_public_reference_total_haplotypes", "Value": 1031},
            {"Metric": "rare_definition", "Value": RARE_DEFINITION},
            {"Metric": "raw_p_lt_0.05", "Value": len(raw_sig)},
            {"Metric": "FDR_lt_0.05", "Value": len(fdr_sig)},
            {"Metric": "Bonferroni_lt_0.05", "Value": len(bonf_sig)},
            {"Metric": "SCZ_higher_raw_p_lt_0.05", "Value": len(scz_higher)},
            {"Metric": "min_p_value", "Value": rare["P_value"].min() if len(rare) else np.nan},
            {"Metric": "min_FDR", "Value": rare["FDR"].min() if len(rare) else np.nan},
            {"Metric": "min_Bonferroni", "Value": rare["Bonferroni"].min() if len(rare) else np.nan},
        ]
    )

    merge_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_tested_genes.xlsx"
    raw_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_nominal_p_lt_0.05_all.xlsx"
    scz_higher_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_nominal_p_lt_0.05_scz_higher_frequency.xlsx"
    fdr_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_fdr_lt_0.05.xlsx"
    bonf_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_bonferroni_lt_0.05.xlsx"
    summary_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_summary.xlsx"

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

    print(summary.to_string(index=False))
    print("\nTop 30 rare Fisher results:")
    print(
        rare[
            [
                "Gene",
                "Count_Scz",
                "Frequency_Scz",
                "Count_comparison_cohort+public_reference",
                "Frequency_comparison_cohort+public_reference",
                "Count_comparison_cohort",
                "Frequency_comparison_cohort",
                "Count_public_reference",
                "Frequency_public_reference",
                "P_value",
                "FDR",
                "Bonferroni",
                "Odds Ratio",
            ]
        ].head(30).to_string(index=False)
    )
    if len(fdr_sig):
        print("\nFDR significant genes:")
        print(
            fdr_sig[
                [
                    "Gene",
                    "Count_Scz",
                    "Frequency_Scz",
                    "Count_comparison_cohort+public_reference",
                    "Frequency_comparison_cohort+public_reference",
                    "Count_comparison_cohort",
                    "Frequency_comparison_cohort",
                    "Count_public_reference",
                    "Frequency_public_reference",
                    "P_value",
                    "FDR",
                    "Bonferroni",
                    "Odds Ratio",
                ]
            ].to_string(index=False)
        )
    else:
        print("\nNo FDR significant genes in the restricted rare-gene Fisher test.")

    for path in [merge_out, raw_out, scz_higher_out, fdr_out, bonf_out, summary_out]:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()

