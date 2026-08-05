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
OUTPUT_PREFIX = "true_case_cohort_comparison_cohort_public_reference_haplotype_common_CTRLgt1pct_SCZrecurrent2"


def resolve_existing_path(primary: Path, legacy: Path) -> Path:
    if primary.exists():
        return primary
    if legacy.exists():
        return legacy
    raise FileNotFoundError(f"Input table not found. Tried: {primary}; {legacy}")


def fisher_p(row: pd.Series) -> float:
    contingency = np.array(
        [
            [int(row["Count_Scz"]), int(row["SCZ_other"])],
            [int(row["Count_comparison_cohort+public_reference"]), int(row["comparison_cohort+public_reference_other"])],
        ]
    )
    _, p_value = fisher_exact(contingency, alternative="greater")
    return float(p_value)


def fisher_or(row: pd.Series) -> float:
    contingency = np.array(
        [
            [int(row["Count_Scz"]), int(row["SCZ_other"])],
            [int(row["Count_comparison_cohort+public_reference"]), int(row["comparison_cohort+public_reference_other"])],
        ]
    )
    odds_ratio, _ = fisher_exact(contingency, alternative="greater")
    return float(odds_ratio)


def main() -> None:
    merge = pd.read_excel(resolve_existing_path(INPUT_MERGE, LEGACY_INPUT_MERGE))
    merge = merge.rename(columns={
        "FDR": "FDR_all_tested_genes",
        "Bonferroni": "Bonferroni_all_tested_genes",
        "P_value": "P_value_all_tested_genes",
        "Odds Ratio": "Odds_Ratio_all_tested_genes",
    })

    scz_total_unique = sorted((merge["Count_Scz"] + merge["SCZ_other"]).unique().tolist())
    ctrl_total_unique = sorted((merge["Count_comparison_cohort+public_reference"] + merge["comparison_cohort+public_reference_other"]).unique().tolist())
    if scz_total_unique != [686]:
        raise ValueError(f"Unexpected SCZ totals: {scz_total_unique}")
    if ctrl_total_unique != [1031]:
        raise ValueError(f"Unexpected comparison_cohort+public_reference totals: {ctrl_total_unique}")

    common = merge[
        (merge["Frequency_comparison_cohort+public_reference"] > 0.01)
        & (merge["Count_Scz"] >= 2)
    ].copy()

    common["Common_definition"] = "Frequency_comparison_cohort+public_reference > 0.01 and Count_Scz >= 2"
    common["P_value"] = common.apply(fisher_p, axis=1)
    common["Odds Ratio"] = common.apply(fisher_or, axis=1)
    common["FDR"] = multipletests(common["P_value"].to_numpy(dtype=float), method="fdr_bh")[1]
    common["Bonferroni"] = multipletests(common["P_value"].to_numpy(dtype=float), method="bonferroni")[1]
    common = common.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

    primary_columns = [
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
        "Common_definition",
    ]
    common = common[primary_columns]

    raw_sig = common[common["P_value"] < 0.05].copy()
    fdr_sig = common[common["FDR"] < 0.05].copy()
    bonf_sig = common[common["Bonferroni"] < 0.05].copy()
    scz_higher = raw_sig[raw_sig["Frequency_Scz"] > raw_sig["Frequency_comparison_cohort+public_reference"]].copy()

    summary = pd.DataFrame(
        [
            {"Metric": "input_all_tested_genes", "Value": len(merge)},
            {"Metric": "common_genes", "Value": len(common)},
            {"Metric": "SCZ_total", "Value": 686},
            {"Metric": "comparison_cohort_public_reference_total", "Value": 1031},
            {"Metric": "common_definition", "Value": "Frequency_comparison_cohort+public_reference > 0.01 and Count_Scz >= 2"},
            {"Metric": "raw_p_lt_0.05", "Value": len(raw_sig)},
            {"Metric": "FDR_lt_0.05", "Value": len(fdr_sig)},
            {"Metric": "Bonferroni_lt_0.05", "Value": len(bonf_sig)},
            {"Metric": "SCZ_higher_raw_p_lt_0.05", "Value": len(scz_higher)},
            {"Metric": "min_p_value", "Value": common["P_value"].min() if len(common) else np.nan},
            {"Metric": "min_FDR", "Value": common["FDR"].min() if len(common) else np.nan},
        ]
    )

    merge_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_tested_genes.xlsx"
    raw_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_nominal_p_lt_0.05_all.xlsx"
    scz_higher_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_nominal_p_lt_0.05_scz_higher_frequency.xlsx"
    fdr_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_fdr_lt_0.05.xlsx"
    summary_out = OUTPUT_DIR / f"{OUTPUT_PREFIX}_summary.xlsx"

    common.to_excel(merge_out, index=False)
    raw_sig.to_excel(raw_out, index=False)
    scz_higher.to_excel(scz_higher_out, index=False)
    fdr_sig.to_excel(fdr_out, index=False)
    with pd.ExcelWriter(summary_out, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        common.head(50).to_excel(writer, sheet_name="Top50_by_P", index=False)
        fdr_sig.to_excel(writer, sheet_name="FDR_lt_0.05", index=False)
        raw_sig.to_excel(writer, sheet_name="RawP_lt_0.05", index=False)

    print(summary.to_string(index=False))
    print(f"Saved: {merge_out}")
    print(f"Saved: {raw_out}")
    print(f"Saved: {scz_higher_out}")
    print(f"Saved: {fdr_out}")
    print(f"Saved: {summary_out}")
    if len(fdr_sig):
        print("FDR significant genes:")
        print(fdr_sig[["Gene", "Count_Scz", "Frequency_Scz", "Count_comparison_cohort+public_reference", "Frequency_comparison_cohort+public_reference", "P_value", "FDR", "Bonferroni", "Odds Ratio"]].to_string(index=False))
    else:
        print("No FDR significant genes in the restricted common-gene Fisher test.")


if __name__ == "__main__":
    main()

