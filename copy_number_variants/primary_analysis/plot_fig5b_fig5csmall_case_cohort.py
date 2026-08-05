# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
CASE_LABEL = "case_cohort"

ANNOTATION_SOURCE = SCRIPT_DIR / "gene_frequencies_with_chr_length_case_cohort_ragtag.txt"
OUTPUT_TABLE = SCRIPT_DIR / "gene_frequencies_with_chr_length_case_cohort.txt"
OUTPUT_BED = SCRIPT_DIR / "gene_info_case_cohort.bed"
OUTPUT_FIG5B_PDF = SCRIPT_DIR / "Fig5B_gene_length_distribution_case_cohort_capped.pdf"
OUTPUT_FIG5B_PNG = SCRIPT_DIR / "Fig5B_gene_length_distribution_case_cohort_capped.png"
OUTPUT_FIG5C_SMALL_PDF = SCRIPT_DIR / "Fig5Csmall_gene_distribution_by_chromosome_case_cohort.pdf"
OUTPUT_FIG5C_SMALL_PNG = SCRIPT_DIR / "Fig5Csmall_gene_distribution_by_chromosome_case_cohort.png"


def find_cnv_base_dir() -> Path:
    for parent in SCRIPT_DIR.parents:
        cnv_root = parent / "CNV"
        if not cnv_root.is_dir():
            continue
        for candidate in cnv_root.glob("*CNV"):
            if not candidate.is_dir():
                continue
            matches = sorted(path for path in candidate.glob("*protein coding genes*") if path.is_dir())
            if matches:
                return matches[0]
    raise FileNotFoundError("Cannot locate filtered CNV directory")


CNV_BASE_DIR = find_cnv_base_dir()
CASE_FREQUENCY_FILE = (
    CNV_BASE_DIR / "slurm_scripts_case_cohort" / "gene_frequencies_filtered.txt"
)


def build_annotation_table() -> pd.DataFrame:
    frequency = pd.read_table(CASE_FREQUENCY_FILE, header=None, names=["Gene", "Col2", "Col3"])
    frequency["Gene"] = frequency["Gene"].astype(str)
    frequency["Col2"] = pd.to_numeric(frequency["Col2"], errors="coerce").fillna(0).astype(int)
    frequency["Col3"] = pd.to_numeric(frequency["Col3"], errors="coerce").fillna(0.0)

    annotation = pd.read_table(ANNOTATION_SOURCE)
    annotation["Gene"] = annotation["Gene"].astype(str)
    annotation = annotation[["Gene", "Chromosome", "Start", "End", "Length"]].drop_duplicates("Gene")

    merged = frequency.merge(annotation, on="Gene", how="left")
    missing = merged.loc[merged["Chromosome"].isna(), "Gene"].tolist()
    if missing:
        raise ValueError(f"Missing annotation for {len(missing)} genes: {missing[:20]}")

    merged["Start"] = pd.to_numeric(merged["Start"], errors="raise").astype(int)
    merged["End"] = pd.to_numeric(merged["End"], errors="raise").astype(int)
    merged["Length"] = pd.to_numeric(merged["Length"], errors="raise").astype(int)
    merged = merged[["Gene", "Col2", "Col3", "Chromosome", "Start", "End", "Length"]]

    merged.to_csv(OUTPUT_TABLE, sep="\t", index=False)
    merged[["Chromosome", "Start", "End", "Gene"]].to_csv(
        OUTPUT_BED, sep="\t", header=False, index=False
    )
    return merged


def plot_gene_length(table: pd.DataFrame) -> None:
    lengths = table["Length"].clip(upper=40000)
    bins = np.linspace(0, 40000, 41)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.hist(lengths, bins=bins, color="#da7271", edgecolor="none", linewidth=0, alpha=0.85)

    try:
        from scipy.stats import gaussian_kde

        x_grid = np.linspace(0, 40000, 400)
        kde = gaussian_kde(lengths)
        bin_width = bins[1] - bins[0]
        ax.plot(x_grid, kde(x_grid) * len(lengths) * bin_width, color="#9d1c2c", linewidth=2)
    except Exception:
        pass

    ax.set_xlabel("Gene Length", fontsize=14)
    ax.set_ylabel("Count of CNV-related Genes", fontsize=14)
    ax.set_xlim(0, 40000)
    ax.set_xticks([0, 10000, 20000, 30000, 40000])
    ax.set_xticklabels(["0", "10000", "20000", "30000", "40000+"])
    ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG5B_PDF, dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_FIG5B_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)


def chromosome_sort_key(chromosome: str) -> int:
    value = str(chromosome).replace("chr", "")
    if value == "X":
        return 23
    if value == "Y":
        return 24
    return int(value)


def plot_chromosome_gene_count(table: pd.DataFrame) -> None:
    autosomes = {f"chr{chromosome}" for chromosome in range(1, 23)}
    table = table[table["Chromosome"].astype(str).isin(autosomes)].copy()
    counts = (
        table.assign(Chromosome=table["Chromosome"].astype(str))
        .groupby("Chromosome", sort=False)
        .size()
        .rename("Count")
        .reset_index()
    )
    counts["SortKey"] = counts["Chromosome"].map(chromosome_sort_key)
    counts = counts.sort_values("SortKey")

    fig, ax = plt.subplots(figsize=(4.2, 5.2))
    y = np.arange(len(counts))
    ax.barh(y, counts["Count"], color="#cd7775", edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels([str(chromosome).replace("chr", "") for chromosome in counts["Chromosome"]], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Count", fontsize=11)
    ax.set_ylabel("Chromosome", fontsize=11)
    ax.set_title("Gene Distribution\non Chromosome", fontsize=11)
    ax.grid(axis="x", alpha=0.25, linestyle="--", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG5C_SMALL_PDF, dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_FIG5C_SMALL_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    table = build_annotation_table()
    plot_gene_length(table)
    plot_chromosome_gene_count(table)

    print(f"CNV_BASE_DIR={CNV_BASE_DIR}")
    print(f"{CASE_LABEL}: {len(table)} annotated genes")
    print(f"Saved: {OUTPUT_TABLE}")
    print(f"Saved: {OUTPUT_BED}")
    print(f"Saved: {OUTPUT_FIG5B_PDF}")
    print(f"Saved: {OUTPUT_FIG5B_PNG}")
    print(f"Saved: {OUTPUT_FIG5C_SMALL_PDF}")
    print(f"Saved: {OUTPUT_FIG5C_SMALL_PNG}")


if __name__ == "__main__":
    main()
