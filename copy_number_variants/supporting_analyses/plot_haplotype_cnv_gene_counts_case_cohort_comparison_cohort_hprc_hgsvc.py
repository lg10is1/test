# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATE_TAG = "26-6-18"
OUTPUT_LABEL = "case_cohort_comparison_cohort_HPRC_HGSVC"

COLORS = {
    "SCZ": "#da7271",
    "comparison_cohort": "#1f78b4",
    "HPRC_HGSVC": "#d9a627",
}


@dataclass(frozen=True)
class Dataset:
    label: str
    input_label: str
    relative_file: Path
    sheet_name: str = "Sheet1"


SCRIPT_DIR = Path(__file__).resolve().parent
CNV_BASE = Path(os.environ.get("EOSCZ_CNV_BASE_DIR", SCRIPT_DIR / "cnv_base"))
OUTPUT_DIR = SCRIPT_DIR / f"output_{OUTPUT_LABEL}_{DATE_TAG}"

DATASETS = (
    Dataset(
        label="SCZ",
        input_label="case_cohort",
        relative_file=Path("slurm_scripts_case_cohort") / "haplotype_CN.xlsx",
    ),
    Dataset(
        label="comparison_cohort",
        input_label="comparison_cohort",
        relative_file=Path("slurm_scripts_comparison_site") / "haplotype_CN.xlsx",
    ),
    Dataset(
        label="HPRC_HGSVC",
        input_label="public_reference",
        relative_file=Path("slurm_scripts_public_reference") / "haplotype_CN.xlsx",
    ),
)


def load_haplotype_matrix(dataset: Dataset) -> pd.DataFrame:
    path = CNV_BASE / dataset.relative_file
    matrix = pd.read_excel(path, sheet_name=dataset.sheet_name, index_col=0)
    matrix = matrix.dropna(how="all")
    if "[xxx]" in matrix.columns:
        matrix = matrix.drop(columns=["[xxx]"])
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return matrix


def count_cnv_genes_per_haplotype(dataset: Dataset) -> tuple[pd.DataFrame, dict[str, object]]:
    matrix = load_haplotype_matrix(dataset)
    counts = (matrix > 0).sum(axis=1).astype(int)
    table = pd.DataFrame(
        {
            "Cohort": dataset.label,
            "Input_Label": dataset.input_label,
            "Haplotype": counts.index.astype(str),
            "CNV_Gene_Count": counts.to_numpy(),
        }
    ).sort_values("CNV_Gene_Count", kind="mergesort")
    manifest = {
        "Cohort": dataset.label,
        "Input_Label": dataset.input_label,
        "Input_File": str(CNV_BASE / dataset.relative_file),
        "Sheet": dataset.sheet_name,
        "N_Haplotypes": matrix.shape[0],
        "N_Genes_In_Matrix": matrix.shape[1],
    }
    return table.reset_index(drop=True), manifest


def summarize(per_haplotype: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cohort, input_label), group in per_haplotype.groupby(["Cohort", "Input_Label"], sort=False):
        values = group["CNV_Gene_Count"]
        rows.append(
            {
                "Cohort": cohort,
                "Input_Label": input_label,
                "N_Haplotypes": int(values.size),
                "Total_CNV_Gene_Counts": int(values.sum()),
                "Mean": float(values.mean()),
                "Median": float(values.median()),
                "SD": float(values.std(ddof=1)) if values.size > 1 else 0.0,
                "Min": int(values.min()),
                "Q1": float(values.quantile(0.25)),
                "Q3": float(values.quantile(0.75)),
                "Max": int(values.max()),
            }
        )
    return pd.DataFrame(rows)


def setup_plot_style() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"


def plot_combined_sorted_bars(per_haplotype: pd.DataFrame, output_prefix: Path) -> None:
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(14.5, 5.2))

    x_offset = 0
    max_count = int(per_haplotype["CNV_Gene_Count"].max())
    centers = []
    labels = []

    for dataset in DATASETS:
        cohort = dataset.label
        group = per_haplotype[per_haplotype["Cohort"] == cohort].sort_values(
            "CNV_Gene_Count", kind="mergesort"
        )
        x = np.arange(len(group)) + x_offset
        ax.bar(
            x,
            group["CNV_Gene_Count"],
            width=0.85,
            color=COLORS[cohort],
            alpha=0.86,
            linewidth=0,
            label=cohort,
        )
        centers.append(x_offset + (len(group) - 1) / 2)
        labels.append(f"{cohort}\n(n={len(group)})")
        x_offset += len(group)
        if cohort != DATASETS[-1].label:
            ax.axvline(x_offset - 0.5, color="#777777", linestyle="--", linewidth=0.8, alpha=0.5)

    ax.set_xticks(centers)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlim(-1, x_offset)
    ax.set_ylim(0, max_count * 1.08)
    ax.set_ylabel("CNV-related protein-coding genes per haplotype", fontsize=11)
    ax.set_xlabel("Haplotypes sorted by CNV gene count within each cohort", fontsize=11)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_boxplot(per_haplotype: pd.DataFrame, output_prefix: Path) -> None:
    setup_plot_style()
    cohorts = [dataset.label for dataset in DATASETS]
    values = [
        per_haplotype.loc[per_haplotype["Cohort"] == cohort, "CNV_Gene_Count"].to_numpy()
        for cohort in cohorts
    ]

    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    box = ax.boxplot(
        values,
        positions=np.arange(1, len(cohorts) + 1),
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.1},
        boxprops={"linewidth": 0.8},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
    )
    for patch, cohort in zip(box["boxes"], cohorts):
        patch.set_facecolor(COLORS[cohort])
        patch.set_alpha(0.5)
        patch.set_edgecolor("black")

    rng = np.random.default_rng(20260618)
    for x, cohort, cohort_values in zip(range(1, len(cohorts) + 1), cohorts, values):
        jitter = rng.uniform(-0.15, 0.15, size=len(cohort_values))
        ax.scatter(
            np.full(len(cohort_values), x) + jitter,
            cohort_values,
            s=8,
            alpha=0.45,
            color=COLORS[cohort],
            edgecolor="white",
            linewidth=0.2,
        )

    ax.set_xticks(np.arange(1, len(cohorts) + 1))
    ax.set_xticklabels(cohorts, fontsize=10)
    ax.set_ylabel("CNV-related protein-coding genes per haplotype", fontsize=11)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    per_haplotype_tables = []
    manifests = []

    for dataset in DATASETS:
        table, manifest = count_cnv_genes_per_haplotype(dataset)
        per_haplotype_tables.append(table)
        manifests.append(manifest)

    per_haplotype = pd.concat(per_haplotype_tables, ignore_index=True)
    summary = summarize(per_haplotype)
    manifest = pd.DataFrame(manifests)

    prefix = OUTPUT_DIR / f"haplotype_CNV_gene_counts_{OUTPUT_LABEL}_{DATE_TAG}"
    with pd.ExcelWriter(prefix.with_suffix(".xlsx"), engine="openpyxl") as writer:
        per_haplotype.to_excel(writer, index=False, sheet_name="per_haplotype_counts")
        summary.to_excel(writer, index=False, sheet_name="summary")
        manifest.to_excel(writer, index=False, sheet_name="input_manifest")

    per_haplotype.to_csv(prefix.with_name(prefix.name + "_per_haplotype_counts.tsv"), sep="\t", index=False)
    summary.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)
    manifest.to_csv(prefix.with_name(prefix.name + "_input_manifest.tsv"), sep="\t", index=False)

    plot_combined_sorted_bars(per_haplotype, prefix.with_name(prefix.name + "_sorted_bar"))
    plot_boxplot(per_haplotype, prefix.with_name(prefix.name + "_boxplot"))

    print(summary.to_string(index=False))
    print(f"Saved output dir: {OUTPUT_DIR}")
    print(f"Saved: {prefix.with_suffix('.xlsx')}")
    print(f"Saved: {prefix.with_name(prefix.name + '_sorted_bar.pdf')}")
    print(f"Saved: {prefix.with_name(prefix.name + '_boxplot.pdf')}")


if __name__ == "__main__":
    main()
