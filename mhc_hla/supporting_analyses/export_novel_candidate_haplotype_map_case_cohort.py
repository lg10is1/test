# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


DATE_TAG = "26-6-18"
OUTPUT_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = OUTPUT_DIR.parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

CANDIDATE_TSV = (
    OUTPUT_DIR
    / "reverse_comparison_known_and_new_candidates_case_cohort_comparison_cohort_public_reference_26-6-18_new_candidate_list.tsv"
)
PUBLIC_IMMUANNOT_DIR = Path(os.environ.get("EOSCZ_IMMUANNOT_DIR", PROJECT_ROOT / "immuannot_results"))
LEGACY_IMMUANNOT_DIR = PROJECT_ROOT / "immuannot_results"


def first_existing_path(primary: Path, legacy: Path) -> Path:
    if primary.exists():
        return primary
    if legacy.exists():
        return legacy
    return primary


comparison_cohort_XLSX = first_existing_path(
    PUBLIC_IMMUANNOT_DIR / "Immuannot_comparison_site_results.xlsx",
    LEGACY_IMMUANNOT_DIR / "immuannot_results.xlsx",
)
ALLPUB_XLSX = first_existing_path(
    PUBLIC_IMMUANNOT_DIR / "Immuannot_public_reference_merged.xlsx",
    LEGACY_IMMUANNOT_DIR / "Immuannot_public_reference_merged.xlsx",
)


def normalize(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def load_source_rows(path: Path, sheet_name: str, cohort: str) -> pd.DataFrame:
    table = pd.read_excel(path, sheet_name=sheet_name, dtype=object)
    table = table.copy()
    table["Source_Cohort"] = cohort
    table["Source_Workbook"] = str(path)
    table["Source_Sheet"] = sheet_name
    table["Sample_name"] = table["Sample_name"].map(normalize)
    for gene in ("A", "B", "C"):
        table[gene] = table[gene].map(normalize)
    return table


def build_haplotype_map(candidates: pd.DataFrame) -> pd.DataFrame:
    source_tables = [
        load_source_rows(comparison_cohort_XLSX, "Sheet1", "comparison_cohort"),
        load_source_rows(ALLPUB_XLSX, "public_reference_result", "public_reference"),
    ]
    source = pd.concat(source_tables, ignore_index=True)

    rows = []
    for _, candidate in candidates.iterrows():
        gene = candidate["Gene"]
        genotype = candidate["Genotype"]
        matched = source[source[gene] == genotype].copy()
        matched = matched.sort_values(["Source_Cohort", "Sample_name"], kind="mergesort")
        for order, row in enumerate(matched.itertuples(index=False), start=1):
            rows.append(
                {
                    "Gene": gene,
                    "Genotype": genotype,
                    "Candidate_HC_Count": int(candidate["HC_Count"]),
                    "Candidate_comparison_cohort_Count": int(candidate["comparison_cohort_Count"]),
                    "Candidate_public_reference_Count": int(candidate["public_reference_Count"]),
                    "Candidate_SCZ_Count": int(candidate["SCZ_Count"]),
                    "P-Value": float(candidate["P-Value"]),
                    "FDR_BH": float(candidate["FDR_BH"]),
                    "Significant_FDR_0.05": bool(candidate["Significant_FDR_0_05"]),
                    "Manual_Check": candidate["Manual_Check"],
                    "Haplotype_Order_Within_Genotype": order,
                    "Source_Cohort": row.Source_Cohort,
                    "Source_Sheet": row.Source_Sheet,
                    "Sample_name": row.Sample_name,
                    "A": row.A,
                    "B": row.B,
                    "C": row.C,
                    "Source_Workbook": row.Source_Workbook,
                    "Sequence_Identity_Check": "pending",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    candidates = pd.read_csv(CANDIDATE_TSV, sep="\t")
    # Make column names attribute-safe for itertuples while preserving original output names later.
    candidates = candidates.rename(columns={"Significant_FDR_0.05": "Significant_FDR_0_05"})
    haplotype_map = build_haplotype_map(candidates)
    candidates_for_output = candidates.rename(columns={"Significant_FDR_0_05": "Significant_FDR_0.05"})
    haplotype_map = haplotype_map.rename(columns={"Significant_FDR_0_05": "Significant_FDR_0.05"})

    summary = (
        haplotype_map.groupby(["Gene", "Genotype"], sort=True)
        .agg(
            N_Haplotypes=("Sample_name", "size"),
            comparison_cohort_Haplotypes=("Source_Cohort", lambda values: int((values == "comparison_cohort").sum())),
            public_reference_Haplotypes=("Source_Cohort", lambda values: int((values == "public_reference").sum())),
            P_Value=("P-Value", "first"),
            FDR_BH=("FDR_BH", "first"),
            Significant_FDR_0_05=("Significant_FDR_0.05", "first"),
            Sequence_Identity_Check=("Sequence_Identity_Check", "first"),
        )
        .reset_index()
        .sort_values(["Gene", "P_Value", "N_Haplotypes"], ascending=[True, True, False], kind="mergesort")
    )
    summary = summary.rename(columns={"Significant_FDR_0_05": "Significant_FDR_0.05"})

    output_prefix = OUTPUT_DIR / f"new_candidate_haplotype_map_case_cohort_comparison_cohort_public_reference_{DATE_TAG}"
    output_xlsx = output_prefix.with_suffix(".xlsx")
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, index=False, sheet_name="summary_by_genotype")
        haplotype_map.to_excel(writer, index=False, sheet_name="haplotype_map")
        candidates_for_output.to_excel(writer, index=False, sheet_name="candidate_new_labels")

    summary_tsv = output_prefix.with_name(output_prefix.name + "_summary_by_genotype.tsv")
    map_tsv = output_prefix.with_name(output_prefix.name + "_haplotype_map.tsv")
    summary.to_csv(summary_tsv, sep="\t", index=False)
    haplotype_map.to_csv(map_tsv, sep="\t", index=False)

    print(summary.to_string(index=False))
    print(f"Total haplotype rows: {len(haplotype_map)}")
    print(f"Saved: {output_xlsx}")
    print(f"Saved: {summary_tsv}")
    print(f"Saved: {map_tsv}")


if __name__ == "__main__":
    main()
