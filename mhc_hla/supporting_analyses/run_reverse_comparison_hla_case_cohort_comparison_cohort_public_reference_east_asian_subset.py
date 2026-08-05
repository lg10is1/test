# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


DATE_TAG = "26-7-10"
GENES = ("A", "B", "C")
SCRIPT_DIR = Path(__file__).resolve().parent


def find_project_root() -> Path:
    for parent in SCRIPT_DIR.parents:
        immun_dirs = [path for path in parent.glob("Immannot*") if path.is_dir()]
        if immun_dirs and (parent / "CNV").is_dir():
            return parent
    raise FileNotFoundError("Cannot locate project root")


PROJECT_ROOT = find_project_root()
IMMUN_DIR = next(path for path in PROJECT_ROOT.glob("Immannot*") if path.is_dir())
SCZ_XLSX = IMMUN_DIR / "Immuannot_case_cohort_forFig4_26-6-14.xlsx"
comparison_cohort_XLSX = next(path for path in IMMUN_DIR.glob("Immuannot_comparison_site*.xlsx") if path.is_file())
EAS_XLSX = IMMUN_DIR / "Immuannot_public_reference_east_asian_subset.xlsx"

INVALID_TOKENS = {"", "NULL", "UNKNOWN", "(\u7a7a\u767d)", "NAN"}


def normalize(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def valid_mask(values: pd.Series) -> pd.Series:
    normalized = values.map(normalize)
    return normalized.map(lambda value: value.upper() not in INVALID_TOKENS)


def find_scz_gene_sheet(gene: str) -> str:
    workbook = pd.ExcelFile(SCZ_XLSX, engine="openpyxl")
    candidates = []
    for sheet_name in workbook.sheet_names:
        header = pd.read_excel(
            SCZ_XLSX,
            sheet_name=sheet_name,
            nrows=0,
            engine="openpyxl",
        )
        if list(header.columns) == ["Sample_name", gene] and sheet_name.startswith(gene):
            candidates.append(sheet_name)
    if len(candidates) != 1:
        raise ValueError(f"Expected one manually curated SCZ sheet for {gene}, found {candidates}")
    return candidates[0]


def load_gene_haplotypes(cohort: str, gene: str) -> tuple[pd.DataFrame, str, Path]:
    if cohort == "SCZ":
        path = SCZ_XLSX
        sheet_name = find_scz_gene_sheet(gene)
    elif cohort == "comparison_cohort":
        path = comparison_cohort_XLSX
        sheet_name = "Sheet1"
    elif cohort == "public_reference_east_asian_subset":
        path = EAS_XLSX
        sheet_name = "public_reference_result"
    else:
        raise ValueError(cohort)

    table = pd.read_excel(path, sheet_name=sheet_name, engine="openpyxl", dtype=object)
    required = {"Sample_name", gene}
    if not required.issubset(table.columns):
        raise ValueError(f"{path} / {sheet_name} lacks columns {sorted(required)}")

    subset = table.loc[:, ["Sample_name", gene]].copy()
    subset["Sample_name"] = subset["Sample_name"].map(normalize)
    subset["Genotype"] = subset[gene].map(normalize)
    subset = subset.loc[valid_mask(subset["Genotype"]), ["Sample_name", "Genotype"]]
    subset.insert(0, "Gene", gene)
    subset.insert(1, "Cohort", cohort)
    subset = subset.reset_index(drop=True)
    return subset, sheet_name, path


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


def build_gene_results(
    gene: str,
    cohort_tables: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    totals = {cohort: len(table) for cohort, table in cohort_tables.items()}
    counts = {
        cohort: table["Genotype"].value_counts().astype(int)
        for cohort, table in cohort_tables.items()
    }
    genotypes = sorted(set().union(*(series.index.astype(str) for series in counts.values())))

    rows = []
    control_total = totals["comparison_cohort"] + totals["public_reference_east_asian_subset"]
    for genotype in genotypes:
        scz_count = int(counts["SCZ"].get(genotype, 0))
        comparison_cohort_count = int(counts["comparison_cohort"].get(genotype, 0))
        eas_count = int(counts["public_reference_east_asian_subset"].get(genotype, 0))
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
                "Genotype": genotype,
                "Category": "novel_label" if "new" in genotype.lower() else "known",
                "SCZ_Count": scz_count,
                "SCZ_Total": totals["SCZ"],
                "SCZ_Frequency": scz_count / totals["SCZ"],
                "comparison_cohort_Count": comparison_cohort_count,
                "comparison_cohort_Total": totals["comparison_cohort"],
                "comparison_cohort_Frequency": comparison_cohort_count / totals["comparison_cohort"],
                "public_reference_east_asian_subset_Count": eas_count,
                "public_reference_east_asian_subset_Total": totals["public_reference_east_asian_subset"],
                "public_reference_east_asian_subset_Frequency": eas_count / totals["public_reference_east_asian_subset"],
                "Control_Count": control_count,
                "Control_Total": control_total,
                "Control_Frequency": control_count / control_total,
                "Odds_Ratio_Control_vs_SCZ": float(odds_ratio),
                "P_Value_Control_Greater": float(p_value),
                "Contingency_Table": repr(contingency.tolist()),
            }
        )

    result = pd.DataFrame(rows)
    result["FDR_BH"] = benjamini_hochberg(result["P_Value_Control_Greater"])
    result["Bonferroni"] = np.minimum(
        result["P_Value_Control_Greater"] * len(result),
        1.0,
    )
    result["Control_Specific_Recurrent"] = (
        (result["Control_Count"] >= 2) & (result["SCZ_Count"] == 0)
    )
    result["RawP_lt_0.05"] = result["P_Value_Control_Greater"] < 0.05
    result["FDR_lt_0.05"] = result["FDR_BH"] < 0.05
    result["Bonferroni_lt_0.05"] = result["Bonferroni"] < 0.05
    return result.sort_values(
        ["P_Value_Control_Greater", "Control_Count", "Genotype"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def summarize_candidates(table: pd.DataFrame, category: str) -> pd.DataFrame:
    filtered = table[table["Category"] == category]
    rows = []
    for gene in (*GENES, "ALL"):
        gene_table = filtered if gene == "ALL" else filtered[filtered["Gene"] == gene]
        candidates = gene_table[gene_table["Control_Specific_Recurrent"]]
        rows.append(
            {
                "Category": category,
                "Gene": gene,
                "Tested_Subtypes": len(gene_table),
                "ControlSpecific_Recurrent_SCZ0": len(candidates),
                "RawP_lt_0.05": int(candidates["RawP_lt_0.05"].sum()),
                "FDR_lt_0.05": int(candidates["FDR_lt_0.05"].sum()),
                "Bonferroni_lt_0.05": int(candidates["Bonferroni_lt_0.05"].sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    cohort_gene_tables: dict[tuple[str, str], pd.DataFrame] = {}
    manifest_rows = []

    for gene in GENES:
        for cohort in ("SCZ", "comparison_cohort", "public_reference_east_asian_subset"):
            table, sheet_name, path = load_gene_haplotypes(cohort, gene)
            cohort_gene_tables[(cohort, gene)] = table
            manifest_rows.append(
                {
                    "Gene": gene,
                    "Cohort": cohort,
                    "Input_File": str(path),
                    "Sheet": sheet_name,
                    "Valid_Haplotypes": len(table),
                    "Unique_Subtype_Labels": table["Genotype"].nunique(),
                }
            )

    result_frames = []
    for gene in GENES:
        result_frames.append(
            build_gene_results(
                gene,
                {
                    cohort: cohort_gene_tables[(cohort, gene)]
                    for cohort in ("SCZ", "comparison_cohort", "public_reference_east_asian_subset")
                },
            )
        )
    all_results = pd.concat(result_frames, ignore_index=True)

    recurrent = all_results[all_results["Control_Specific_Recurrent"]].copy()
    known = recurrent[recurrent["Category"] == "known"].copy()
    novel = recurrent[recurrent["Category"] == "novel_label"].copy()
    known_raw = known[known["RawP_lt_0.05"]].copy()
    novel_raw = novel[novel["RawP_lt_0.05"]].copy()

    control_haplotypes = pd.concat(
        [
            cohort_gene_tables[(cohort, gene)]
            for gene in GENES
            for cohort in ("comparison_cohort", "public_reference_east_asian_subset")
        ],
        ignore_index=True,
    )
    novel_keys = novel.loc[:, ["Gene", "Genotype", "Control_Count"]].rename(
        columns={"Control_Count": "Candidate_Control_Count"}
    )
    novel_haplotypes = control_haplotypes.merge(
        novel_keys,
        on=["Gene", "Genotype"],
        how="inner",
        validate="many_to_one",
    )
    novel_haplotypes["Manual_Sequence_Check"] = (
        "required: identical subtype labels may represent distinct sequences"
    )
    novel_haplotypes = novel_haplotypes.sort_values(
        ["Gene", "Genotype", "Cohort", "Sample_name"],
        kind="mergesort",
    )

    summary = pd.concat(
        [
            summarize_candidates(all_results, "known"),
            summarize_candidates(all_results, "novel_label"),
        ],
        ignore_index=True,
    )
    manifest = pd.DataFrame(manifest_rows)

    prefix = SCRIPT_DIR / (
        "reverse_comparison_HLA_case_cohort_"
        f"comparison_cohort_public_reference_east_asian_subset_{DATE_TAG}"
    )
    output_xlsx = prefix.with_suffix(".xlsx")
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, index=False, sheet_name="summary")
        manifest.to_excel(writer, index=False, sheet_name="input_manifest")
        known.to_excel(writer, index=False, sheet_name="known_SCZ0")
        known_raw.to_excel(writer, index=False, sheet_name="known_rawP")
        novel.to_excel(writer, index=False, sheet_name="novel_label_SCZ0")
        novel_raw.to_excel(writer, index=False, sheet_name="novel_label_rawP")
        novel_haplotypes.to_excel(writer, index=False, sheet_name="novel_haplotypes")
        recurrent.to_excel(writer, index=False, sheet_name="all_SCZ0_recurrent")
        all_results.to_excel(writer, index=False, sheet_name="all_reverse_results")

    summary.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)
    known.to_csv(prefix.with_name(prefix.name + "_known_SCZ0.tsv"), sep="\t", index=False)
    known_raw.to_csv(prefix.with_name(prefix.name + "_known_SCZ0_rawP_lt_0.05.tsv"), sep="\t", index=False)
    novel.to_csv(prefix.with_name(prefix.name + "_novel_label_SCZ0.tsv"), sep="\t", index=False)
    novel_raw.to_csv(prefix.with_name(prefix.name + "_novel_label_SCZ0_rawP_lt_0.05.tsv"), sep="\t", index=False)
    novel_haplotypes.to_csv(
        prefix.with_name(prefix.name + "_novel_candidate_haplotypes.tsv"),
        sep="\t",
        index=False,
    )

    print(summary.to_string(index=False))
    print(f"Known recurrent SCZ0: {len(known)}; raw P < 0.05: {len(known_raw)}")
    print(f"Novel label recurrent SCZ0: {len(novel)}; raw P < 0.05: {len(novel_raw)}")
    print(f"Saved: {output_xlsx}")


if __name__ == "__main__":
    main()
