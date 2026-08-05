# -*- coding: utf-8 -*-

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests


SCRIPT_DIR = Path(__file__).resolve().parent
CNV_BASE_DIR = next(path for path in SCRIPT_DIR.parent.glob("*protein coding genes*") if path.is_dir())

CASE_FILE = CNV_BASE_DIR / "slurm_scripts_case_cohort" / "haplotype_CN.xlsx"
comparison_cohort_FILE = CNV_BASE_DIR / "slurm_scripts_comparison_site" / "haplotype_CN.xlsx"
ALLPUB_FILE = CNV_BASE_DIR / "slurm_scripts_public_reference" / "haplotype_CN.xlsx"

GENE_SET_CACHE = SCRIPT_DIR / "pathway_gene_sets_case_cohort_vs_comparison_cohort_public_reference_26-4-27.json"
RESULT_XLSX = SCRIPT_DIR / "case_cohort_vs_comparison_cohort_public_reference_pathway_burden_results_26-6-14.xlsx"
CONTRIBUTOR_XLSX = SCRIPT_DIR / "case_cohort_vs_comparison_cohort_public_reference_pathway_contributing_genes_26-6-14.xlsx"
PRIMARY_CONTRIBUTOR_XLSX = (
    SCRIPT_DIR / "case_cohort_vs_comparison_cohort_public_reference_Neurof_GoNervTransm_contributing_genes_26-6-14.xlsx"
)


@dataclass(frozen=True)
class GroupMeta:
    name: str
    total_haplotypes: int
    rare_gene_count: int


@dataclass(frozen=True)
class RareGeneContext:
    case_freq: pd.Series
    ctrl_freq: pd.Series
    case_rare_genes: set[str]
    ctrl_rare_genes: set[str]


def load_cached_gene_sets() -> tuple[dict[str, list[str]], dict[str, object]]:
    cache = json.loads(GENE_SET_CACHE.read_text(encoding="utf-8"))
    gene_sets = {
        str(set_name): sorted({str(gene).strip().upper() for gene in genes if str(gene).strip()})
        for set_name, genes in cache["gene_sets"].items()
    }
    return gene_sets, cache.get("metadata", {})


def load_haplotype_matrix(file_path: Path) -> pd.DataFrame:
    df = pd.read_excel(file_path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    df.columns = [str(column).strip().upper() for column in df.columns]
    if df.columns.duplicated().any():
        df = df.T.groupby(level=0).sum().T
    return df.astype(int)


def align_gene_tables(*tables: pd.DataFrame) -> list[pd.DataFrame]:
    union_genes = sorted({gene for table in tables for gene in table.columns})
    aligned_tables = []
    for table in tables:
        aligned = table.reindex(columns=union_genes, fill_value=0)
        aligned_tables.append(aligned.astype(int))
    return aligned_tables


def calc_freq(matrix: pd.DataFrame) -> pd.Series:
    return (matrix > 0).sum(axis=0) / matrix.shape[0]


def calc_burden(matrix: pd.DataFrame, rare_genes: set[str], gene_set: list[str]) -> pd.Series:
    genes = [gene for gene in gene_set if gene in rare_genes and gene in matrix.columns]
    if not genes:
        return pd.Series(0, index=matrix.index, dtype=int)
    return matrix[genes].sum(axis=1)


def build_burden_binary(
    case: pd.DataFrame,
    ctrl: pd.DataFrame,
    gene_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, GroupMeta], pd.DataFrame, RareGeneContext]:
    case_freq = calc_freq(case)
    ctrl_freq = calc_freq(ctrl)
    case_rare_genes = {gene for gene in case_freq.index if case_freq[gene] < 0.05}
    ctrl_rare_genes = {gene for gene in ctrl_freq.index if ctrl_freq[gene] < 0.05}

    burden_rows = []
    burden_meta = []
    for group_name, matrix, rare_genes in [
        ("SCZ", case, case_rare_genes),
        ("CTRL", ctrl, ctrl_rare_genes),
    ]:
        burden_df = pd.DataFrame(index=matrix.index)
        for set_name, genes in gene_sets.items():
            burden_df[set_name] = calc_burden(matrix, rare_genes, genes)
        burden_df["group"] = group_name
        burden_rows.append(burden_df)

        for set_name, genes in gene_sets.items():
            rare_in_set = sorted(set(genes).intersection(rare_genes))
            burden_meta.append(
                {
                    "Group": group_name,
                    "GeneSet": set_name,
                    "GeneSetSize": len(genes),
                    "RareGenesInGroup": len(rare_in_set),
                    "RareGeneExamples": ", ".join(rare_in_set[:20]),
                }
            )

    burden_all = pd.concat(burden_rows, axis=0)
    burden_binary = burden_all.copy()
    burden_binary.loc[:, burden_binary.columns != "group"] = (
        burden_binary.loc[:, burden_binary.columns != "group"] > 0
    ).astype(int)

    group_meta = {
        "SCZ": GroupMeta("SCZ", int(case.shape[0]), len(case_rare_genes)),
        "CTRL": GroupMeta("CTRL", int(ctrl.shape[0]), len(ctrl_rare_genes)),
    }
    rare_context = RareGeneContext(
        case_freq=case_freq,
        ctrl_freq=ctrl_freq,
        case_rare_genes=case_rare_genes,
        ctrl_rare_genes=ctrl_rare_genes,
    )
    return burden_binary, group_meta, pd.DataFrame(burden_meta), rare_context


def add_multiple_testing(result_df: pd.DataFrame) -> pd.DataFrame:
    p_values = result_df["p_value"].to_numpy(dtype=float)
    result_df["FDR"] = multipletests(p_values, method="fdr_bh")[1]
    result_df["Bonferroni"] = multipletests(p_values, method="bonferroni")[1]
    return result_df.sort_values(["FDR", "p_value", "GeneSet"], kind="mergesort").reset_index(drop=True)


def fisher_burden_test(burden_binary: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    results = []
    for set_name in gene_sets:
        x = burden_binary[set_name].astype(int)
        scz_with = int(((burden_binary["group"] == "SCZ") & (x == 1)).sum())
        scz_without = int(((burden_binary["group"] == "SCZ") & (x == 0)).sum())
        ctrl_with = int(((burden_binary["group"] == "CTRL") & (x == 1)).sum())
        ctrl_without = int(((burden_binary["group"] == "CTRL") & (x == 0)).sum())
        contingency = np.array([[ctrl_without, ctrl_with], [scz_without, scz_with]])

        if int(contingency[:, 1].sum()) == 0:
            odds_ratio = np.nan
            p_value = 1.0
        else:
            odds_ratio, p_value = fisher_exact(contingency, alternative="greater")

        results.append(
            {
                "GeneSet": set_name,
                "Test": "Fisher_one_sided",
                "SCZ_with_burden": scz_with,
                "SCZ_total": scz_with + scz_without,
                "CTRL_with_burden": ctrl_with,
                "CTRL_total": ctrl_with + ctrl_without,
                "SCZ_rate": scz_with / (scz_with + scz_without),
                "CTRL_rate": ctrl_with / (ctrl_with + ctrl_without),
                "OddsRatio": odds_ratio,
                "p_value": p_value,
            }
        )
    return add_multiple_testing(pd.DataFrame(results))


def chisq_burden_test(burden_binary: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    results = []
    for set_name in gene_sets:
        x = burden_binary[set_name].astype(int)
        scz_with = int(((burden_binary["group"] == "SCZ") & (x == 1)).sum())
        scz_without = int(((burden_binary["group"] == "SCZ") & (x == 0)).sum())
        ctrl_with = int(((burden_binary["group"] == "CTRL") & (x == 1)).sum())
        ctrl_without = int(((burden_binary["group"] == "CTRL") & (x == 0)).sum())
        contingency = np.array([[ctrl_without, ctrl_with], [scz_without, scz_with]])

        if int(contingency[:, 1].sum()) == 0:
            chi_square = np.nan
            p_value = 1.0
        else:
            chi_square, p_value, _, _ = chi2_contingency(contingency, correction=False)

        results.append(
            {
                "GeneSet": set_name,
                "Test": "Chi_square",
                "SCZ_with_burden": scz_with,
                "SCZ_total": scz_with + scz_without,
                "CTRL_with_burden": ctrl_with,
                "CTRL_total": ctrl_with + ctrl_without,
                "SCZ_rate": scz_with / (scz_with + scz_without),
                "CTRL_rate": ctrl_with / (ctrl_with + ctrl_without),
                "Chi_square": chi_square,
                "p_value": p_value,
            }
        )
    return add_multiple_testing(pd.DataFrame(results))


def logistic_burden_test(burden_binary: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    results = []
    for set_name in gene_sets:
        df = pd.DataFrame(
            {
                "group": (burden_binary["group"] == "SCZ").astype(int),
                "burden": burden_binary[set_name].astype(int),
            }
        )
        if int(df["burden"].sum()) == 0:
            results.append({"GeneSet": set_name, "Beta": np.nan, "SE": np.nan, "OR": np.nan, "p_value": 1.0})
            continue

        exog = sm.add_constant(df["burden"], has_constant="add")
        try:
            fit = sm.Logit(df["group"], exog).fit(disp=False)
            beta = float(fit.params["burden"])
            se = float(fit.bse["burden"])
            p_value = float(fit.pvalues["burden"])
            odds_ratio = float(np.exp(beta))
        except Exception:
            beta = np.nan
            se = np.nan
            odds_ratio = np.nan
            p_value = 1.0

        results.append({"GeneSet": set_name, "Beta": beta, "SE": se, "OR": odds_ratio, "p_value": p_value})
    return add_multiple_testing(pd.DataFrame(results))


def build_contributor_table(
    pathway_name: str,
    case: pd.DataFrame,
    ctrl: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    rare_context: RareGeneContext,
) -> pd.DataFrame:
    rows = []
    for gene in sorted(set(gene_sets[pathway_name]).intersection(case.columns)):
        scz_carriers = int((case[gene] > 0).sum())
        ctrl_carriers = int((ctrl[gene] > 0).sum())
        scz_rare = gene in rare_context.case_rare_genes
        ctrl_rare = gene in rare_context.ctrl_rare_genes
        contributes_to_scz = bool(scz_rare and scz_carriers > 0)
        contributes_to_ctrl = bool(ctrl_rare and ctrl_carriers > 0)
        if not contributes_to_scz and not contributes_to_ctrl:
            continue

        if contributes_to_scz and contributes_to_ctrl:
            status = "Shared"
        elif contributes_to_scz:
            status = "SCZ_only"
        else:
            status = "CTRL_only"

        rows.append(
            {
                "Pathway": pathway_name,
                "Gene": gene,
                "SCZ_carriers": scz_carriers,
                "SCZ_total_CN": int(case[gene].sum()),
                "SCZ_freq": float(rare_context.case_freq[gene]),
                "SCZ_rare": scz_rare,
                "CTRL_carriers": ctrl_carriers,
                "CTRL_total_CN": int(ctrl[gene].sum()),
                "CTRL_freq": float(rare_context.ctrl_freq[gene]),
                "CTRL_rare": ctrl_rare,
                "Contributes_to_SCZ_burden": contributes_to_scz,
                "Contributes_to_CTRL_burden": contributes_to_ctrl,
                "Status": status,
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "Pathway",
                "Gene",
                "SCZ_carriers",
                "SCZ_total_CN",
                "SCZ_freq",
                "SCZ_rare",
                "CTRL_carriers",
                "CTRL_total_CN",
                "CTRL_freq",
                "CTRL_rare",
                "Contributes_to_SCZ_burden",
                "Contributes_to_CTRL_burden",
                "Status",
            ]
        )

    return pd.DataFrame(rows).sort_values(
        ["Status", "SCZ_carriers", "CTRL_carriers", "Gene"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )


def write_primary_contributor_workbook(
    pathway_name: str,
    contributor_table: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    group_meta: dict[str, GroupMeta],
) -> None:
    scz_contributors = contributor_table.loc[contributor_table["Contributes_to_SCZ_burden"]].drop(
        columns=["Pathway"]
    )
    ctrl_contributors = contributor_table.loc[contributor_table["Contributes_to_CTRL_burden"]].drop(
        columns=["Pathway"]
    )
    shared = contributor_table.loc[contributor_table["Status"] == "Shared"].drop(columns=["Pathway"])
    all_contributors = contributor_table.drop(columns=["Pathway"])
    summary = pd.DataFrame(
        [
            {"Metric": "PathwayName", "Value": pathway_name},
            {"Metric": "PathwayGeneCount", "Value": len(gene_sets[pathway_name])},
            {"Metric": "SCZ_haplotypes", "Value": group_meta["SCZ"].total_haplotypes},
            {"Metric": "CTRL_haplotypes", "Value": group_meta["CTRL"].total_haplotypes},
            {"Metric": "SCZ_contributing_genes", "Value": int(len(scz_contributors))},
            {"Metric": "CTRL_contributing_genes", "Value": int(len(ctrl_contributors))},
            {"Metric": "Shared_contributing_genes", "Value": int(len(shared))},
            {"Metric": "All_contributing_genes", "Value": int(len(all_contributors))},
        ]
    )
    with pd.ExcelWriter(PRIMARY_CONTRIBUTOR_XLSX, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        scz_contributors.to_excel(writer, sheet_name="SCZ_Contributors", index=False)
        ctrl_contributors.to_excel(writer, sheet_name="CTRL_Contributors", index=False)
        shared.to_excel(writer, sheet_name="Shared", index=False)
        all_contributors.to_excel(writer, sheet_name="All_Contributors", index=False)


def write_results(
    fisher_results: pd.DataFrame,
    chisq_results: pd.DataFrame,
    logistic_results: pd.DataFrame,
    group_meta: dict[str, GroupMeta],
    gene_sets: dict[str, list[str]],
    gene_set_metadata: dict[str, object],
    burden_meta: pd.DataFrame,
    all_contributors: pd.DataFrame,
) -> None:
    meta_rows = [
        {"Metric": "SCZ_haplotypes", "Value": group_meta["SCZ"].total_haplotypes},
        {"Metric": "CTRL_haplotypes", "Value": group_meta["CTRL"].total_haplotypes},
        {"Metric": "SCZ_rare_genes", "Value": group_meta["SCZ"].rare_gene_count},
        {"Metric": "CTRL_rare_genes", "Value": group_meta["CTRL"].rare_gene_count},
        {"Metric": "Gene_sets_tested", "Value": len(gene_sets)},
        {"Metric": "Case_file", "Value": str(CASE_FILE)},
        {"Metric": "comparison_cohort_file", "Value": str(comparison_cohort_FILE)},
        {"Metric": "public_reference_file", "Value": str(ALLPUB_FILE)},
    ]
    gene_set_rows = [
        {"GeneSet": set_name, "GeneSetSize": len(genes), "FirstGenes": ", ".join(genes[:25])}
        for set_name, genes in gene_sets.items()
    ]

    with pd.ExcelWriter(RESULT_XLSX, engine="openpyxl") as writer:
        fisher_results.to_excel(writer, sheet_name="Fisher", index=False)
        chisq_results.to_excel(writer, sheet_name="ChiSquare", index=False)
        logistic_results.to_excel(writer, sheet_name="Logistic", index=False)
        pd.DataFrame(meta_rows).to_excel(writer, sheet_name="Meta", index=False)
        pd.DataFrame(gene_set_rows).to_excel(writer, sheet_name="GeneSets", index=False)
        burden_meta.to_excel(writer, sheet_name="RareGeneMeta", index=False)
        pd.DataFrame(gene_set_metadata.get("reactome_axon_pathways", [])).to_excel(
            writer, sheet_name="ReactomeAxonPathways", index=False
        )
        all_contributors.to_excel(writer, sheet_name="ContributingGenes", index=False)

    summary = (
        all_contributors.groupby("Pathway", dropna=False)
        .agg(
            SCZ_contributing_genes=("Contributes_to_SCZ_burden", "sum"),
            CTRL_contributing_genes=("Contributes_to_CTRL_burden", "sum"),
            Shared_contributing_genes=("Status", lambda values: int((values == "Shared").sum())),
            All_contributing_genes=("Gene", "nunique"),
        )
        .reset_index()
    )
    with pd.ExcelWriter(CONTRIBUTOR_XLSX, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        all_contributors.to_excel(writer, sheet_name="All_Contributors", index=False)


def main() -> None:
    gene_sets, gene_set_metadata = load_cached_gene_sets()

    case = load_haplotype_matrix(CASE_FILE)
    comparison_cohort = load_haplotype_matrix(comparison_cohort_FILE)
    public_reference = load_haplotype_matrix(ALLPUB_FILE)
    case, comparison_cohort, public_reference = align_gene_tables(case, comparison_cohort, public_reference)
    ctrl = pd.concat([comparison_cohort, public_reference], axis=0)

    burden_binary, group_meta, burden_meta, rare_context = build_burden_binary(case, ctrl, gene_sets)
    fisher_results = fisher_burden_test(burden_binary, gene_sets)
    chisq_results = chisq_burden_test(burden_binary, gene_sets)
    logistic_results = logistic_burden_test(burden_binary, gene_sets)

    contributor_tables = [
        build_contributor_table(set_name, case, ctrl, gene_sets, rare_context) for set_name in gene_sets
    ]
    all_contributors = pd.concat(contributor_tables, axis=0, ignore_index=True)
    if "Neurof_GoNervTransm" in gene_sets:
        primary_contributors = build_contributor_table("Neurof_GoNervTransm", case, ctrl, gene_sets, rare_context)
        write_primary_contributor_workbook("Neurof_GoNervTransm", primary_contributors, gene_sets, group_meta)

    write_results(
        fisher_results=fisher_results,
        chisq_results=chisq_results,
        logistic_results=logistic_results,
        group_meta=group_meta,
        gene_sets=gene_sets,
        gene_set_metadata=gene_set_metadata,
        burden_meta=burden_meta,
        all_contributors=all_contributors,
    )

    print(f"SCZ haplotypes: {group_meta['SCZ'].total_haplotypes}")
    print(f"CTRL haplotypes: {group_meta['CTRL'].total_haplotypes}")
    print(f"SCZ rare genes: {group_meta['SCZ'].rare_gene_count}")
    print(f"CTRL rare genes: {group_meta['CTRL'].rare_gene_count}")
    print(f"Gene sets tested: {len(gene_sets)}")
    print("Top Fisher results:")
    print(fisher_results.to_string(index=False))
    print(f"Saved: {RESULT_XLSX}")
    print(f"Saved: {CONTRIBUTOR_XLSX}")
    if PRIMARY_CONTRIBUTOR_XLSX.exists():
        print(f"Saved: {PRIMARY_CONTRIBUTOR_XLSX}")


if __name__ == "__main__":
    main()
