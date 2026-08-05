# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_LABEL = "case_cohort_comparison_cohort_HPRC_HGSVC"
DATE_TAG = "26-6-18"
OUTPUT_DIR = SCRIPT_DIR / f"output_{OUTPUT_LABEL}_{DATE_TAG}"
PREFIX = OUTPUT_DIR / f"haplotype_CNV_gene_counts_{OUTPUT_LABEL}_{DATE_TAG}"

INPUT_TSV = PREFIX.with_name(PREFIX.name + "_per_haplotype_counts.tsv")
OUTPUT_PREFIX = PREFIX.with_name(PREFIX.name + "_sorted_bar")

COHORT_ORDER = ("SCZ", "comparison_cohort", "HPRC_HGSVC")
COLORS = {
    "SCZ": "#da7271",
    "comparison_cohort": "#1f78b4",
    "HPRC_HGSVC": "#d9a627",
}
DISPLAY_LABELS = {
    "SCZ": "SCZ",
    "comparison_cohort": "Control",
    "HPRC_HGSVC": "Public",
}


def setup_plot_style() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"


def plot_combined_sorted_bars(per_haplotype: pd.DataFrame, output_prefix: Path) -> None:
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(14.5, 5.2))

    x_offset = 0
    max_count = int(per_haplotype["CNV_Gene_Count"].max())
    centers: list[float] = []
    labels: list[str] = []

    for cohort in COHORT_ORDER:
        group = per_haplotype[per_haplotype["Cohort"] == cohort].sort_values(
            "CNV_Gene_Count",
            kind="mergesort",
        )
        x = np.arange(len(group)) + x_offset
        ax.bar(
            x,
            group["CNV_Gene_Count"],
            width=0.85,
            color=COLORS[cohort],
            alpha=0.86,
            linewidth=0,
            label=DISPLAY_LABELS[cohort],
        )
        centers.append(x_offset + (len(group) - 1) / 2)
        labels.append(f"{DISPLAY_LABELS[cohort]}\n(n={len(group)})")
        x_offset += len(group)
        if cohort != COHORT_ORDER[-1]:
            ax.axvline(
                x_offset - 0.5,
                color="#777777",
                linestyle="--",
                linewidth=0.8,
                alpha=0.5,
            )

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


def main() -> None:
    per_haplotype = pd.read_csv(INPUT_TSV, sep="\t")
    plot_combined_sorted_bars(per_haplotype, OUTPUT_PREFIX)
    print(f"Read data: {INPUT_TSV}")
    print(f"Updated: {OUTPUT_PREFIX.with_suffix('.pdf')}")
    print(f"Updated: {OUTPUT_PREFIX.with_suffix('.png')}")


if __name__ == "__main__":
    main()
