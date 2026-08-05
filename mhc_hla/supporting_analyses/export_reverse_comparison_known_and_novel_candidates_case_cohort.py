# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


DATE_TAG = "26-6-18"
RESULT_TAG = "known_reverse_case_cohort_vs_comparison_cohort_public_reference_actual_totals_26-6-18"
GENES = ("A", "B", "C")

SCRIPT_DIR = Path(__file__).resolve().parent
FISHER_DIR = SCRIPT_DIR / "comparison_known_reverse_case_cohort_comparison_cohort_public_reference_26-6-18"
OUTPUT_DIR = SCRIPT_DIR / "reverse_comparison_control_specific_lists_case_cohort_comparison_cohort_public_reference_26-6-18"

PROJECT_ROOT = SCRIPT_DIR.parents[2]
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

INVALID_TOKENS = {"", "NULL", "UNKNOWN", "(blank)"}


def normalize(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def is_valid(value: str) -> bool:
    return bool(value) and value.upper() not in INVALID_TOKENS


def count_gene_values(path: Path, sheet_name: str, gene: str) -> pd.Series:
    table = pd.read_excel(path, sheet_name=sheet_name, dtype=object)
    if gene not in table.columns:
        raise ValueError(f"{path} / {sheet_name} lacks column {gene}")
    values = table[gene].map(normalize)
    values = values[values.map(is_valid)]
    return values.value_counts().astype(int)


def add_control_split_counts(table: pd.DataFrame, gene: str) -> pd.DataFrame:
    comparison_cohort_counts = count_gene_values(comparison_cohort_XLSX, "Sheet1", gene)
    public_reference_counts = count_gene_values(ALLPUB_XLSX, "public_reference_result", gene)
    enriched = table.copy()
    enriched["comparison_cohort_Count"] = enriched["Genotype"].map(comparison_cohort_counts).fillna(0).astype(int)
    enriched["public_reference_Count"] = enriched["Genotype"].map(public_reference_counts).fillna(0).astype(int)
    enriched["Check_HC_Count_Equals_comparison_cohort_plus_public_reference"] = (
        enriched["HC_Count"] == enriched["comparison_cohort_Count"] + enriched["public_reference_Count"]
    )
    return enriched


def load_reverse_results() -> pd.DataFrame:
    frames = []
    for gene in GENES:
        path = FISHER_DIR / f"{gene}_fisher_results_{RESULT_TAG}.xlsx"
        table = pd.read_excel(path, sheet_name="Sheet1")
        table = add_control_split_counts(table, gene)
        frames.append(table)
    return pd.concat(frames, ignore_index=True)


def is_new(genotype: object) -> bool:
    return "new" in str(genotype).lower()


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = load_reverse_results()
    recurrent_scz0 = all_results[
        (all_results["HC_Count"] >= 2)
        & (all_results["SCZ_Count"] == 0)
    ].copy()
    recurrent_scz0["Control_Specific_Recurrent_Criteria"] = "comparison_cohort+public_reference_count>=2; SCZ_count=0"

    known = recurrent_scz0[~recurrent_scz0["Genotype"].map(is_new)].copy()
    known["Category"] = "known"
    known["Manual_Check"] = "not required for known labels"

    new_candidates = recurrent_scz0[recurrent_scz0["Genotype"].map(is_new)].copy()
    new_candidates["Category"] = "new_candidate_label_level"
    new_candidates["Manual_Check"] = (
        "sequence identity not checked; same new label may represent distinct sequences"
    )

    sort_columns = ["Gene", "P-Value", "HC_Count", "Genotype"]
    known = known.sort_values(sort_columns, ascending=[True, True, False, True]).reset_index(drop=True)
    new_candidates = new_candidates.sort_values(sort_columns, ascending=[True, True, False, True]).reset_index(drop=True)

    summary = pd.DataFrame(
        [
            {
                "Category": "known_control_specific_recurrent_SCZ0",
                "N": len(known),
                "RawP_lt_0.05": int((known["P-Value"] < 0.05).sum()),
                "FDR_lt_0.05": int((known["FDR_BH"] < 0.05).sum()),
                "Notes": "known; comparison_cohort+public_reference_count>=2; SCZ_count=0",
            },
            {
                "Category": "new_candidate_control_specific_recurrent_SCZ0_label_level",
                "N": len(new_candidates),
                "RawP_lt_0.05": int((new_candidates["P-Value"] < 0.05).sum()),
                "FDR_lt_0.05": int((new_candidates["FDR_BH"] < 0.05).sum()),
                "Notes": "new label candidates only; require sequence-identity check",
            },
        ]
    )

    by_gene = (
        pd.concat([known, new_candidates], ignore_index=True)
        .groupby(["Category", "Gene"], dropna=False)
        .agg(
            N=("Genotype", "size"),
            RawP_lt_0_05=("P-Value", lambda values: int((values < 0.05).sum())),
            FDR_lt_0_05=("FDR_BH", lambda values: int((values < 0.05).sum())),
        )
        .reset_index()
    )

    output_prefix = OUTPUT_DIR / f"reverse_comparison_known_and_new_candidates_case_cohort_comparison_cohort_public_reference_{DATE_TAG}"
    output_xlsx = output_prefix.with_suffix(".xlsx")
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, index=False, sheet_name="summary")
        by_gene.to_excel(writer, index=False, sheet_name="summary_by_gene")
        known.to_excel(writer, index=False, sheet_name="known_SCZ0_list")
        new_candidates.to_excel(writer, index=False, sheet_name="new_candidate_list")
        recurrent_scz0.to_excel(writer, index=False, sheet_name="all_SCZ0_recurrent")
        all_results.to_excel(writer, index=False, sheet_name="all_reverse_results")

    known_tsv = output_prefix.with_name(output_prefix.name + "_known_SCZ0_list.tsv")
    new_tsv = output_prefix.with_name(output_prefix.name + "_new_candidate_list.tsv")
    summary_tsv = output_prefix.with_name(output_prefix.name + "_summary.tsv")
    known.to_csv(known_tsv, sep="\t", index=False)
    new_candidates.to_csv(new_tsv, sep="\t", index=False)
    summary.to_csv(summary_tsv, sep="\t", index=False)

    print(summary.to_string(index=False))
    print()
    print(by_gene.to_string(index=False))
    print(f"Saved: {output_xlsx}")
    print(f"Saved: {known_tsv}")
    print(f"Saved: {new_tsv}")
    print(f"Saved: {summary_tsv}")


if __name__ == "__main__":
    main()
