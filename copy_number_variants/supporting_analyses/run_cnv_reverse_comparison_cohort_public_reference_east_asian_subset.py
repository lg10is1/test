# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


DATE_TAG = "26-7-10"
SCRIPT_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Dataset:
    key: str
    label: str
    directory_name: str
    file_name: str


DATASETS = (
    Dataset(
        key="SCZ",
        label="case_cohort",
        directory_name="slurm_scripts_case_cohort",
        file_name="haplotype_CN.xlsx",
    ),
    Dataset(
        key="comparison_cohort",
        label="comparison_cohort",
        directory_name="slurm_scripts_comparison_site",
        file_name="haplotype_CN.xlsx",
    ),
    Dataset(
        key="public_reference_east_asian_subset",
        label="public_reference_east_asian_subset",
        directory_name="slurm_scripts_public_reference_east_asian_subset",
        file_name="haplotype_CN.east_asian_subset.xlsx",
    ),
)


def find_cnv_base() -> Path:
    for parent in SCRIPT_DIR.parents:
        cnv_root = parent / "CNV"
        if not cnv_root.is_dir():
            continue
        for candidate in cnv_root.glob("*CNV"):
            if not candidate.is_dir():
                continue
            matches = [path for path in candidate.glob("*protein coding genes*") if path.is_dir()]
            if matches:
                return sorted(matches)[0]
    raise FileNotFoundError("Cannot locate filtered CNV directory")


CNV_BASE = find_cnv_base()


def input_path(dataset: Dataset) -> Path:
    return CNV_BASE / dataset.directory_name / dataset.file_name


def load_presence(dataset: Dataset) -> pd.DataFrame:
    path = input_path(dataset)
    matrix = pd.read_excel(
        path,
        sheet_name="Sheet1",
        index_col=0,
        engine="openpyxl",
    )
    matrix = matrix.dropna(how="all")
    if "[xxx]" in matrix.columns:
        matrix = matrix.drop(columns=["[xxx]"])
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    return (matrix > 0).astype(int)


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    p_values = values.to_numpy(dtype=float)
    n_values = len(p_values)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted_ranked = np.minimum.accumulate(
        (ranked * n_values / np.arange(1, n_values + 1))[::-1]
    )[::-1]
    adjusted = np.empty(n_values, dtype=float)
    adjusted[order] = np.clip(adjusted_ranked, 0, 1)
    return pd.Series(adjusted, index=values.index)


def main() -> None:
    matrices = {dataset.key: load_presence(dataset) for dataset in DATASETS}
    totals = {key: matrix.shape[0] for key, matrix in matrices.items()}
    control_total = totals["comparison_cohort"] + totals["public_reference_east_asian_subset"]

    all_genes = sorted(
        set().union(*(matrix.columns.astype(str) for matrix in matrices.values()))
    )
    counts = {
        key: matrix.sum(axis=0).astype(int).reindex(all_genes, fill_value=0).astype(int)
        for key, matrix in matrices.items()
    }

    rows = []
    for gene in all_genes:
        scz_count = int(counts["SCZ"].loc[gene])
        comparison_cohort_count = int(counts["comparison_cohort"].loc[gene])
        eas_count = int(counts["public_reference_east_asian_subset"].loc[gene])
        control_count = comparison_cohort_count + eas_count
        contingency = np.array(
            [
                [control_count, control_total - control_count],
                [scz_count, totals["SCZ"] - scz_count],
            ],
            dtype=int,
        )
        odds_ratio, p_value = fisher_exact(contingency, alternative="greater")
        rows.append(
            {
                "Gene": gene,
                "Count_Scz": scz_count,
                "Frequency_Scz": scz_count / totals["SCZ"],
                "Count_comparison_cohort": comparison_cohort_count,
                "Frequency_comparison_cohort": comparison_cohort_count / totals["comparison_cohort"],
                "Count_public_reference_east_asian_subset": eas_count,
                "Frequency_public_reference_east_asian_subset": eas_count / totals["public_reference_east_asian_subset"],
                "Count_comparison_cohort+public_reference_east_asian_subset": control_count,
                "Frequency_comparison_cohort+public_reference_east_asian_subset": control_count / control_total,
                "SCZ_other": totals["SCZ"] - scz_count,
                "comparison_cohort+public_reference_east_asian_subset_other": control_total - control_count,
                "Odds_Ratio_Control_vs_SCZ": float(odds_ratio),
                "P_value": float(p_value),
                "Contingency_Table": repr(contingency.tolist()),
            }
        )

    result = pd.DataFrame(rows)
    result["FDR"] = benjamini_hochberg(result["P_value"])
    result["Bonferroni"] = np.minimum(result["P_value"] * len(result), 1.0)
    result["Reverse_Candidate"] = (
        (result["Count_comparison_cohort+public_reference_east_asian_subset"] >= 2)
        & (result["Count_Scz"] == 0)
    )
    result["RawP_lt_0.05"] = result["P_value"] < 0.05
    result["FDR_lt_0.05"] = result["FDR"] < 0.05
    result["Bonferroni_lt_0.05"] = result["Bonferroni"] < 0.05
    result = result.sort_values(
        ["P_value", "Gene"],
        kind="mergesort",
    ).reset_index(drop=True)

    candidates = result[result["Reverse_Candidate"]].copy()
    candidates_raw = candidates[candidates["RawP_lt_0.05"]].copy()
    candidates_fdr = candidates[candidates["FDR_lt_0.05"]].copy()
    candidates_bonf = candidates[candidates["Bonferroni_lt_0.05"]].copy()
    significant_raw = result[result["RawP_lt_0.05"]].copy()

    manifest = pd.DataFrame(
        [
            {
                "Cohort_Key": dataset.key,
                "Label": dataset.label,
                "Input_File": str(input_path(dataset)),
                "Sheet": "Sheet1",
                "N_Haplotypes": totals[dataset.key],
                "N_Gene_Columns": matrices[dataset.key].shape[1],
            }
            for dataset in DATASETS
        ]
    )
    summary = pd.DataFrame(
        [
            {"Metric": "SCZ_haplotypes", "Value": totals["SCZ"]},
            {"Metric": "comparison_cohort_haplotypes", "Value": totals["comparison_cohort"]},
            {"Metric": "public_reference_east_asian_subset_haplotypes", "Value": totals["public_reference_east_asian_subset"]},
            {
                "Metric": "comparison_cohort+public_reference_east_asian_subset_haplotypes",
                "Value": control_total,
            },
            {"Metric": "Tested_gene_union", "Value": len(result)},
            {"Metric": "All_genes_rawP_lt_0.05", "Value": len(significant_raw)},
            {
                "Metric": "All_genes_FDR_lt_0.05",
                "Value": int(result["FDR_lt_0.05"].sum()),
            },
            {
                "Metric": "All_genes_Bonferroni_lt_0.05",
                "Value": int(result["Bonferroni_lt_0.05"].sum()),
            },
            {
                "Metric": "Reverse_candidates_controlGE2_SCZ0",
                "Value": len(candidates),
            },
            {
                "Metric": "Reverse_candidates_rawP_lt_0.05",
                "Value": len(candidates_raw),
            },
            {
                "Metric": "Reverse_candidates_FDR_lt_0.05",
                "Value": len(candidates_fdr),
            },
            {
                "Metric": "Reverse_candidates_Bonferroni_lt_0.05",
                "Value": len(candidates_bonf),
            },
        ]
    )

    prefix = SCRIPT_DIR / (
        "CNV_reverse_comparison_cohort_public_reference_east_asian_subset_vs_"
        f"case_cohort_{DATE_TAG}"
    )
    output_xlsx = prefix.with_suffix(".xlsx")
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        result.to_excel(writer, index=False, sheet_name="all_tested_genes")
        candidates.to_excel(writer, index=False, sheet_name="reverse_candidates")
        candidates_raw.to_excel(writer, index=False, sheet_name="reverse_rawP")
        candidates_fdr.to_excel(writer, index=False, sheet_name="reverse_FDR")
        candidates_bonf.to_excel(writer, index=False, sheet_name="reverse_Bonferroni")
        significant_raw.to_excel(writer, index=False, sheet_name="all_significant_raw")
        summary.to_excel(writer, index=False, sheet_name="summary")
        manifest.to_excel(writer, index=False, sheet_name="input_manifest")

    result.to_csv(
        prefix.with_name(prefix.name + "_all_tested_genes.tsv"),
        sep="\t",
        index=False,
    )
    candidates.to_csv(
        prefix.with_name(prefix.name + "_reverse_candidates.tsv"),
        sep="\t",
        index=False,
    )
    candidates_raw.to_csv(
        prefix.with_name(prefix.name + "_reverse_candidates_rawP_lt_0.05.tsv"),
        sep="\t",
        index=False,
    )
    candidates_raw.loc[:, ["Gene"]].to_csv(
        prefix.with_name(prefix.name + "_reverse_gene_list_rawP_lt_0.05.tsv"),
        sep="\t",
        index=False,
    )
    summary.to_csv(
        prefix.with_name(prefix.name + "_summary.tsv"),
        sep="\t",
        index=False,
    )
    manifest.to_csv(
        prefix.with_name(prefix.name + "_input_manifest.tsv"),
        sep="\t",
        index=False,
    )

    print(summary.to_string(index=False))
    print()
    print("Top reverse candidates with raw P < 0.05:")
    print(
        candidates_raw.loc[
            :,
            [
                "Gene",
                "Count_comparison_cohort+public_reference_east_asian_subset",
                "Count_Scz",
                "P_value",
                "FDR",
                "Bonferroni",
            ],
        ].head(20).to_string(index=False)
    )
    print(f"Saved: {output_xlsx}")


if __name__ == "__main__":
    main()
