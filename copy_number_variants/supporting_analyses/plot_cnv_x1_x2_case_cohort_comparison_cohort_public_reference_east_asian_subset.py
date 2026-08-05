# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


DATE_TAG = "26-7-10"
SCRIPT_DIR = Path(__file__).resolve().parent
COHORTS = ("SCZ", "comparison_cohort", "PUB-east_asian_subset")
COLORS = {
    "SCZ": "#da7271",
    "comparison_cohort": "#1f78b4",
    "PUB-east_asian_subset": "#d9a627",
}


@dataclass(frozen=True)
class Dataset:
    cohort: str
    directory_name: str
    sample_file: str
    haplotype_file: str


DATASETS = (
    Dataset(
        cohort="SCZ",
        directory_name="slurm_scripts_case_cohort",
        sample_file="sample_CN.xlsx",
        haplotype_file="haplotype_CN.xlsx",
    ),
    Dataset(
        cohort="comparison_cohort",
        directory_name="slurm_scripts_comparison_site",
        sample_file="sample_CN.xlsx",
        haplotype_file="haplotype_CN.xlsx",
    ),
    Dataset(
        cohort="PUB-east_asian_subset",
        directory_name="slurm_scripts_public_reference_east_asian_subset",
        sample_file="sample_CN.east_asian_subset.xlsx",
        haplotype_file="haplotype_CN.east_asian_subset.xlsx",
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


def dataset_dir(dataset: Dataset) -> Path:
    return CNV_BASE / dataset.directory_name


def clean_numeric_matrix(table: pd.DataFrame) -> pd.DataFrame:
    table = table.dropna(how="all")
    if "[xxx]" in table.columns:
        table = table.drop(columns=["[xxx]"])
    return table.apply(pd.to_numeric, errors="coerce").fillna(0)


def load_x1_gene_counts(dataset: Dataset) -> tuple[pd.DataFrame, dict[str, object]]:
    path = dataset_dir(dataset) / dataset.sample_file
    if dataset.cohort in {"SCZ", "comparison_cohort"}:
        table = pd.read_excel(path, sheet_name="Sum_ECN", engine="openpyxl")
        gene_counts = table.loc[:, ["Gene", "Sum_ECN_true"]].copy()
        gene_counts.columns = ["Gene", "Extra_Copy_Number"]
        input_sheet = "Sum_ECN"
    else:
        matrix = pd.read_excel(
            path,
            sheet_name="Sheet1",
            index_col=0,
            engine="openpyxl",
        )
        matrix = clean_numeric_matrix(matrix)
        summed = matrix.sum(axis=0)
        gene_counts = pd.DataFrame(
            {
                "Gene": summed.index.astype(str),
                "Extra_Copy_Number": summed.to_numpy(),
            }
        )
        input_sheet = "Sheet1; summed across 65 samples"

    gene_counts["Extra_Copy_Number"] = pd.to_numeric(
        gene_counts["Extra_Copy_Number"],
        errors="coerce",
    ).fillna(0).astype(int)
    gene_counts = gene_counts[gene_counts["Extra_Copy_Number"] > 0].copy()
    gene_counts.insert(0, "Cohort", dataset.cohort)
    gene_counts = gene_counts.sort_values(
        ["Extra_Copy_Number", "Gene"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    manifest = {
        "Analysis": "X1",
        "Cohort": dataset.cohort,
        "Input_File": str(path),
        "Sheet": input_sheet,
        "N_CNV_Related_Genes": len(gene_counts),
    }
    return gene_counts, manifest


def group_extra_copy(values: pd.Series) -> pd.Series:
    labels = [str(value) for value in range(1, 10)] + ["10-99", "100+"]

    def assign(value: int) -> str:
        if 1 <= value <= 9:
            return str(value)
        if 10 <= value <= 99:
            return "10-99"
        return "100+"

    return values.map(assign).value_counts().reindex(labels, fill_value=0) / len(values)


def plot_x1(distribution: pd.DataFrame, output_prefix: Path) -> None:
    labels = distribution["Bin"].tolist()
    x = np.arange(len(labels))
    width = 0.24

    fig, ax = plt.subplots(figsize=(10.5, 6))
    for index, cohort in enumerate(COHORTS):
        offsets = x + (index - 1) * width
        ax.bar(
            offsets,
            distribution[cohort].to_numpy(),
            width=width,
            color=COLORS[cohort],
            edgecolor="black",
            linewidth=0.4,
            label=cohort,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlabel("Extra Copy Number", fontsize=14)
    ax.set_ylabel("Proportion of CNV-related Protein-coding Genes", fontsize=13)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def load_x2_assembly_counts(dataset: Dataset) -> tuple[pd.DataFrame, dict[str, object]]:
    path = dataset_dir(dataset) / dataset.haplotype_file
    matrix = pd.read_excel(
        path,
        sheet_name="Sheet1",
        index_col=0,
        engine="openpyxl",
    )
    matrix = clean_numeric_matrix(matrix)
    presence = (matrix > 0).astype(int)
    counts = presence.sum(axis=1).astype(int)
    result = pd.DataFrame(
        {
            "Cohort": dataset.cohort,
            "Assembly": counts.index.astype(str),
            "CNV_Related_ProteinCoding_Genes": counts.to_numpy(),
        }
    )
    manifest = {
        "Analysis": "X2",
        "Cohort": dataset.cohort,
        "Input_File": str(path),
        "Sheet": "Sheet1",
        "N_Haplotype_Assemblies": matrix.shape[0],
        "N_Gene_Columns": matrix.shape[1],
        "Unit": "haplotype assembly",
    }
    return result, manifest


def summarize_x2(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort in COHORTS:
        values = table.loc[
            table["Cohort"] == cohort,
            "CNV_Related_ProteinCoding_Genes",
        ]
        rows.append(
            {
                "Cohort": cohort,
                "N_Haplotype_Assemblies": len(values),
                "Mean": float(values.mean()),
                "Median": float(values.median()),
                "SD": float(values.std(ddof=1)),
                "Min": int(values.min()),
                "Q1": float(values.quantile(0.25)),
                "Q3": float(values.quantile(0.75)),
                "Max": int(values.max()),
            }
        )
    return pd.DataFrame(rows)


def plot_x2(table: pd.DataFrame, output_prefix: Path) -> None:
    values_by_cohort = [
        table.loc[
            table["Cohort"] == cohort,
            "CNV_Related_ProteinCoding_Genes",
        ].to_numpy()
        for cohort in COHORTS
    ]
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    box = ax.boxplot(
        values_by_cohort,
        positions=np.arange(1, len(COHORTS) + 1),
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
        boxprops={"linewidth": 0.8},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
    )
    for patch, cohort in zip(box["boxes"], COHORTS):
        patch.set_facecolor(COLORS[cohort])
        patch.set_alpha(0.48)
        patch.set_edgecolor("black")

    rng = np.random.default_rng(20260710)
    for position, cohort, values in zip(
        range(1, len(COHORTS) + 1),
        COHORTS,
        values_by_cohort,
    ):
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        ax.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=8,
            alpha=0.55,
            color=COLORS[cohort],
            edgecolor="white",
            linewidth=0.2,
        )
    ax.set_xticks(np.arange(1, len(COHORTS) + 1))
    ax.set_xticklabels(COHORTS, fontsize=10)
    ax.set_ylabel("CNV-related protein-coding genes per haplotype assembly", fontsize=10.5)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    x1_frames = []
    x2_frames = []
    manifest_rows = []
    for dataset in DATASETS:
        x1_table, x1_manifest = load_x1_gene_counts(dataset)
        x2_table, x2_manifest = load_x2_assembly_counts(dataset)
        x1_frames.append(x1_table)
        x2_frames.append(x2_table)
        manifest_rows.extend([x1_manifest, x2_manifest])

    x1_gene_counts = pd.concat(x1_frames, ignore_index=True)
    distributions = {
        cohort: group_extra_copy(
            x1_gene_counts.loc[
                x1_gene_counts["Cohort"] == cohort,
                "Extra_Copy_Number",
            ]
        )
        for cohort in COHORTS
    }
    x1_distribution = pd.DataFrame(distributions)
    x1_distribution.index.name = "Bin"
    x1_distribution = x1_distribution.reset_index()

    x2_counts = pd.concat(x2_frames, ignore_index=True)
    x2_summary = summarize_x2(x2_counts)
    manifest = pd.DataFrame(manifest_rows)

    x1_prefix = SCRIPT_DIR / (
        "CNV_X1_extra_copy_distribution_case_cohort_"
        f"comparison_cohort_PUB_east_asian_subset_{DATE_TAG}"
    )
    x2_prefix = SCRIPT_DIR / (
        "CNV_X2_per_haplotype_assembly_gene_counts_case_cohort_"
        f"comparison_cohort_PUB_east_asian_subset_{DATE_TAG}"
    )
    output_xlsx = SCRIPT_DIR / (
        "CNV_X1_X2_case_cohort_"
        f"comparison_cohort_PUB_east_asian_subset_{DATE_TAG}.xlsx"
    )

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        x1_gene_counts.to_excel(writer, index=False, sheet_name="X1_gene_counts")
        x1_distribution.to_excel(writer, index=False, sheet_name="X1_distribution")
        x2_counts.to_excel(writer, index=False, sheet_name="X2_per_assembly")
        x2_summary.to_excel(writer, index=False, sheet_name="X2_summary")
        manifest.to_excel(writer, index=False, sheet_name="input_manifest")

    x1_gene_counts.to_csv(
        x1_prefix.with_name(x1_prefix.name + "_gene_counts.tsv"),
        sep="\t",
        index=False,
    )
    x1_distribution.to_csv(
        x1_prefix.with_name(x1_prefix.name + "_distribution.tsv"),
        sep="\t",
        index=False,
    )
    x2_counts.to_csv(
        x2_prefix.with_name(x2_prefix.name + "_per_assembly.tsv"),
        sep="\t",
        index=False,
    )
    x2_summary.to_csv(
        x2_prefix.with_name(x2_prefix.name + "_summary.tsv"),
        sep="\t",
        index=False,
    )
    manifest.to_csv(
        SCRIPT_DIR / f"CNV_X1_X2_input_manifest_{DATE_TAG}.tsv",
        sep="\t",
        index=False,
    )

    plot_x1(x1_distribution, x1_prefix)
    plot_x2(x2_counts, x2_prefix)

    print("X1 CNV-related genes:")
    print(x1_gene_counts.groupby("Cohort").size().to_string())
    print()
    print("X1 distribution:")
    print(x1_distribution.to_string(index=False))
    print()
    print("X2 summary:")
    print(x2_summary.to_string(index=False))
    print(f"Saved: {output_xlsx}")


if __name__ == "__main__":
    main()
