# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, PercentFormatter


SCRIPT_DIR = Path(__file__).resolve().parent
COMPARISON_NAME = "case_cohort_vs_comparison_cohort_vs_public_reference"
OUTPUT_DIR = SCRIPT_DIR / f"output_Fig5D_Fig5E_{COMPARISON_NAME}_26-6-16"

COLORS = {
    "case_cohort": "#da7271",
    "comparison_cohort": "#1f78b4",
    "public_reference": "#d9a627",
}

INPUTS = {
    "Fig5D": OUTPUT_DIR
    / f"Fig5D_Extra_Copy_Number_{COMPARISON_NAME}_1_10_direct_crop_full_denominator_full_data.tsv",
    "Fig5E": OUTPUT_DIR
    / f"Fig5E_Sample_Count_{COMPARISON_NAME}_1_10_direct_crop_full_denominator_full_data.tsv",
}

OUTPUTS = {
    "Fig5D": OUTPUT_DIR / f"Fig5D_Extra_Copy_Number_{COMPARISON_NAME}_exact_full_range_no_grouping",
    "Fig5E": OUTPUT_DIR / f"Fig5E_Sample_Count_{COMPARISON_NAME}_exact_full_range_no_grouping",
}

XLABELS = {
    "Fig5D": "Extra Copy Number",
    "Fig5E": "Sample Count",
}


def read_exact_distribution(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    table["Bin"] = pd.to_numeric(table["Bin"], errors="raise").astype(int)
    table = table.sort_values("Bin").reset_index(drop=True)
    group_columns = [column for column in table.columns if column != "Bin"]
    table[group_columns] = table[group_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return table


def figure_width(max_bin: int) -> float:
    return min(72.0, max(18.0, max_bin / 80.0))


def plot_exact_full_range(table: pd.DataFrame, output_prefix: Path, xlabel: str) -> None:
    group_columns = [column for column in table.columns if column != "Bin"]
    x = table["Bin"].to_numpy(dtype=float)
    max_bin = int(table["Bin"].max())
    width = min(0.8 / len(group_columns), 0.28)

    fig, ax = plt.subplots(figsize=(figure_width(max_bin), 6.0))
    offset_start = -width * (len(group_columns) - 1) / 2
    for group_index, group_label in enumerate(group_columns):
        ax.bar(
            x + offset_start + group_index * width,
            table[group_label].to_numpy(dtype=float),
            width=width,
            label=group_label,
            color=COLORS.get(group_label, "#777777"),
            edgecolor="none",
            linewidth=0,
        )

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Proportion of CNV-related Genes", fontsize=14)
    ax.set_xlim(0.5, max_bin + 0.5)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=14, integer=True))
    ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.5)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_note(summary_rows: list[dict[str, object]]) -> None:
    lines = [
        "Fig.5D/Fig.5E exact full-range no-grouping note",
        f"Comparison: {COMPARISON_NAME}",
        "These plots use exact integer bins from the full-data TSV files.",
        "No x-axis grouping is applied.",
        "No x-axis range truncation is applied.",
        "The x-axis spans from 1 to the observed maximum bin for each panel.",
        "Proportions are calculated using the full positive CNV-related gene denominator for each group.",
        "",
        "Panel\tRows\tMinBin\tMaxBin\tColumnSums",
    ]
    for row in summary_rows:
        lines.append(
            f"{row['Panel']}\t{row['Rows']}\t{row['MinBin']}\t{row['MaxBin']}\t{row['ColumnSums']}"
        )
    (OUTPUT_DIR / "Fig5D_Fig5E_exact_full_range_no_grouping_NOTE.txt").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    summary_rows: list[dict[str, object]] = []
    for panel, input_path in INPUTS.items():
        table = read_exact_distribution(input_path)
        output_prefix = OUTPUTS[panel]
        plot_exact_full_range(table, output_prefix, XLABELS[panel])
        table.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)

        group_columns = [column for column in table.columns if column != "Bin"]
        column_sums = "; ".join(f"{column}={table[column].sum():.6f}" for column in group_columns)
        summary_rows.append(
            {
                "Panel": panel,
                "Rows": len(table),
                "MinBin": int(table["Bin"].min()),
                "MaxBin": int(table["Bin"].max()),
                "ColumnSums": column_sums,
            }
        )
        print(f"Saved: {output_prefix.with_suffix('.pdf')}")
        print(f"Saved: {output_prefix.with_suffix('.png')}")
        print(f"Saved: {output_prefix.with_suffix('.tsv')}")

    write_note(summary_rows)
    print(f"Saved: {OUTPUT_DIR / 'Fig5D_Fig5E_exact_full_range_no_grouping_NOTE.txt'}")


if __name__ == "__main__":
    main()
