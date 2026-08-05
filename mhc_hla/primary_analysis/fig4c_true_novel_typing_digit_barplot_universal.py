# -*- coding: utf-8 -*-
"""
Create the new Fig.4C: true novel HLA allele species counts by 2-/3-/4-field class.

Input workbook logic:
1. SCZ_all_candidate_new_alleles_1 contains all Immuannot candidate novel allele
   haplotype rows.
2. SCZ_rec_candidate_new_alleles_2 contains recurrent candidate novel allele rows
   after sequence-level manual check.
3. For recurrent rows, the true allele label must use:
   "SCZ Recurrent Candidate New alleles After Sequence-level Manual Check".
4. For non-recurrent rows, the Immuannot candidate label from
   SCZ_all_candidate_new_alleles_1 is retained.
5. The plotted counts are unique true allele labels, not haplotype-row counts.

Important manual-check note:
This script does not infer sequence identity. It trusts the manual-check labels in
SCZ_rec_candidate_new_alleles_2. If the manual labels are revised, rerun the script.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


DEFAULT_ALL_SHEET = "SCZ_all_candidate_new_alleles_1"
DEFAULT_REC_SHEET = "SCZ_rec_candidate_new_alleles_2"
DEFAULT_ALL_LABEL_COL = "All Immauannot Candidate New alleles"
DEFAULT_REC_LABEL_COL = "SCZ Recurrent Candidate New alleles After Sequence-level Manual Check"
DEFAULT_FIG4B_SUMMARY = (
    Path(__file__).resolve().parent
    / "Subtype_Stat_case_cohort_summary_26-6-14.xlsx"
)

DIGIT_CLASSES = ["2-digit", "3-digit", "4-digit"]
GENE_ORDER = ["A", "B", "C"]
GENE_DIGIT_COLORS = {
    "A": {
        "2-digit": "#47a1a2",
        "3-digit": "#6fbfc0",
        "4-digit": "#9edcdd",
    },
    "B": {
        "2-digit": "#da7271",
        "3-digit": "#e99594",
        "4-digit": "#f2bbbb",
    },
    "C": {
        "2-digit": "#1f78b4",
        "3-digit": "#319fd0",
        "4-digit": "#4DBBD5",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count true novel HLA allele species by 2-/3-/4-field class and draw "
            "the new Fig.4C."
        )
    )
    parser.add_argument("input_xlsx", help="Input new_MHC_subtypes xlsx path.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output figure path. Default: script folder/Fig4C_true_new_allele_digit_counts_<input stem>.pdf",
    )
    parser.add_argument(
        "--png-output",
        help="Optional PNG output path. Default: same basename as PDF.",
    )
    parser.add_argument(
        "--summary-output",
        help="Output summary xlsx path. Default: same basename as PDF plus _summary.xlsx.",
    )
    parser.add_argument("--all-sheet", default=DEFAULT_ALL_SHEET)
    parser.add_argument("--rec-sheet", default=DEFAULT_REC_SHEET)
    parser.add_argument("--all-label-col", default=DEFAULT_ALL_LABEL_COL)
    parser.add_argument("--rec-label-col", default=DEFAULT_REC_LABEL_COL)
    parser.add_argument(
        "--fig4b-summary",
        default=str(DEFAULT_FIG4B_SUMMARY),
        help="Optional Fig.4B summary xlsx for consistency check.",
    )
    parser.add_argument("--title", default="")
    parser.add_argument("--xlabel", default="HLA genes")
    parser.add_argument("--ylabel", default="Allele counts")
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--figsize", nargs=2, type=float, default=(7.2, 6.0))
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def normalize_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def gene_short(value: object) -> str:
    text = normalize_text(value)
    return text.replace("HLA-", "").strip()


def digit_class(label: object) -> str:
    text = normalize_text(label)
    if "*" not in text:
        return "unclassified"
    allele_body = text.split("*", 1)[1]
    field_count = len(allele_body.split(":"))
    if field_count == 2:
        return "2-digit"
    if field_count == 3:
        return "3-digit"
    if field_count == 4:
        return "4-digit"
    return f"{field_count}-field"


def load_and_resolve_labels(args: argparse.Namespace) -> pd.DataFrame:
    input_path = Path(args.input_xlsx).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input xlsx not found: {input_path}")

    all_df = pd.read_excel(input_path, sheet_name=args.all_sheet, dtype=object)
    rec_df = pd.read_excel(input_path, sheet_name=args.rec_sheet, dtype=object)

    required_all = {"Haplotype", "Gene", args.all_label_col}
    required_rec = {"Haplotype", "Gene", args.rec_label_col}
    missing_all = required_all.difference(all_df.columns)
    missing_rec = required_rec.difference(rec_df.columns)
    if missing_all:
        raise ValueError(f"Missing columns in {args.all_sheet}: {sorted(missing_all)}")
    if missing_rec:
        raise ValueError(f"Missing columns in {args.rec_sheet}: {sorted(missing_rec)}")

    all_df = all_df.rename(columns={args.all_label_col: "Immuannot_Candidate_Label"})
    rec_df = rec_df.rename(columns={args.rec_label_col: "Manual_Checked_Label"})

    resolved = all_df.merge(
        rec_df[["Haplotype", "Gene", "Manual_Checked_Label"]],
        on=["Haplotype", "Gene"],
        how="left",
        validate="one_to_one",
    )
    resolved["True_Label"] = resolved["Manual_Checked_Label"].where(
        resolved["Manual_Checked_Label"].notna(),
        resolved["Immuannot_Candidate_Label"],
    )
    resolved["Gene_Short"] = resolved["Gene"].map(gene_short)
    resolved["Digit_Class"] = resolved["True_Label"].map(digit_class)
    resolved["Used_Manual_Checked_Label"] = resolved["Manual_Checked_Label"].notna()
    return resolved


def summarize(resolved: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    unique_labels = (
        resolved.drop_duplicates(["Gene_Short", "True_Label"])
        .sort_values(["Gene_Short", "Digit_Class", "True_Label"])
        .reset_index(drop=True)
    )

    summary = (
        unique_labels.pivot_table(
            index="Gene_Short",
            columns="Digit_Class",
            values="True_Label",
            aggfunc="count",
            fill_value=0,
        )
        .reindex(GENE_ORDER)
        .fillna(0)
    )
    for digit in DIGIT_CLASSES:
        if digit not in summary.columns:
            summary[digit] = 0
    summary = summary[DIGIT_CLASSES].astype(int)
    summary["Total_True_New_Allele_Species"] = summary.sum(axis=1)
    summary = summary.reset_index().rename(columns={"Gene_Short": "Gene"})

    diagnostics_rows = []
    for gene_name, group in resolved.groupby("Gene_Short"):
        row_count = int(len(group))
        unique_count = int(group.drop_duplicates(["Gene_Short", "True_Label"]).shape[0])
        recurrent_row_count = int(group["Used_Manual_Checked_Label"].sum())
        diagnostics_rows.append(
            {
                "Gene": gene_name,
                "Candidate_New_Haplotype_Rows": row_count,
                "True_New_Allele_Species": unique_count,
                "Manual_Label_Rows": recurrent_row_count,
                "Reduction_Implied_By_Manual_Labels": row_count - unique_count,
            }
        )
    diagnostics = pd.DataFrame(diagnostics_rows).sort_values("Gene").reset_index(drop=True)
    return summary, unique_labels, diagnostics


def load_fig4b_comparison(fig4b_summary_path: Path, summary: pd.DataFrame) -> pd.DataFrame:
    if not fig4b_summary_path.exists():
        return pd.DataFrame(
            {
                "Note": [
                    f"Fig.4B summary not found, skipped consistency check: {fig4b_summary_path}"
                ]
            }
        )

    fig4b = pd.read_excel(fig4b_summary_path)
    comparison = summary[["Gene", "Total_True_New_Allele_Species"]].merge(
        fig4b[["Gene", "Final_New_Subtype_Count", "Manual_Duplicate_Reduction"]],
        on="Gene",
        how="outer",
    )
    comparison["Difference_vs_Fig4B"] = (
        comparison["Total_True_New_Allele_Species"] - comparison["Final_New_Subtype_Count"]
    )
    return comparison


def configure_matplotlib(font_family: str) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = font_family
    plt.rcParams["axes.linewidth"] = 1.0


def add_color_grid_legend(axis) -> None:
    from matplotlib.patches import Rectangle

    x_start = 1.03
    y_start = 0.92
    row_height = 0.075
    col_width = 0.105
    label_width = 0.16
    box_width = 0.075
    box_height = 0.045

    for col_index, digit in enumerate(DIGIT_CLASSES):
        axis.text(
            x_start + label_width + col_index * col_width + box_width / 2,
            y_start,
            digit,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontsize=11,
            clip_on=False,
        )

    for row_index, gene in enumerate(GENE_ORDER):
        y_pos = y_start - (row_index + 1) * row_height
        axis.text(
            x_start,
            y_pos + box_height / 2,
            f"HLA-{gene}",
            transform=axis.transAxes,
            ha="left",
            va="center",
            fontsize=11,
            clip_on=False,
        )
        for col_index, digit in enumerate(DIGIT_CLASSES):
            x_pos = x_start + label_width + col_index * col_width
            axis.add_patch(
                Rectangle(
                    (x_pos, y_pos),
                    box_width,
                    box_height,
                    transform=axis.transAxes,
                    facecolor=GENE_DIGIT_COLORS[gene][digit],
                    edgecolor="black",
                    linewidth=0.5,
                    clip_on=False,
                )
            )


def plot_summary(summary: pd.DataFrame, args: argparse.Namespace, output_path: Path, png_path: Path) -> None:
    import matplotlib.pyplot as plt

    configure_matplotlib(args.font_family)

    plot_df = summary.set_index("Gene").reindex(GENE_ORDER)[DIGIT_CLASSES]
    fig, axis = plt.subplots(figsize=tuple(args.figsize))
    fig.subplots_adjust(right=0.73)

    bottom = pd.Series([0] * len(plot_df), index=plot_df.index, dtype=float)
    x_positions = range(len(plot_df.index))
    for digit in DIGIT_CLASSES:
        values = plot_df[digit].astype(int)
        segment_colors = [GENE_DIGIT_COLORS[gene][digit] for gene in plot_df.index]
        axis.bar(
            x_positions,
            values,
            bottom=bottom,
            color=segment_colors,
            width=0.62,
            edgecolor="none",
            linewidth=0,
        )
        for x_pos, value, base in zip(x_positions, values, bottom):
            if value > 0:
                axis.text(
                    x_pos,
                    base + value / 2,
                    str(int(value)),
                    ha="center",
                    va="center",
                    fontsize=13,
                    color="black",
                )
        bottom += values

    for x_pos, total in zip(x_positions, bottom):
        axis.text(x_pos, total + max(bottom) * 0.015, str(int(total)), ha="center", va="bottom", fontsize=13)

    axis.set_xticks(list(x_positions))
    axis.set_xticklabels([f"HLA-{gene}" for gene in plot_df.index], fontsize=18)
    axis.set_ylabel(args.ylabel, fontsize=20)
    axis.set_xlabel(args.xlabel, fontsize=20)
    if args.title:
        axis.set_title(args.title, fontsize=20)
    axis.tick_params(axis="y", labelsize=14)
    add_color_grid_legend(axis)
    axis.spines["top"].set_visible(True)
    axis.spines["right"].set_visible(True)
    axis.set_ylim(0, max(bottom) * 1.12)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


def default_output_path(input_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return script_dir / f"Fig4C_true_new_allele_digit_counts_{input_path.stem}_26-7-16.pdf"


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_xlsx).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve() if args.output else default_output_path(input_path)
    png_path = (
        Path(args.png_output).expanduser().resolve()
        if args.png_output
        else output_path.with_suffix(".png")
    )
    summary_path = (
        Path(args.summary_output).expanduser().resolve()
        if args.summary_output
        else output_path.with_name(output_path.stem + "_summary.xlsx")
    )
    fig4b_summary_path = Path(args.fig4b_summary).expanduser().resolve()

    resolved = load_and_resolve_labels(args)
    summary, unique_labels, diagnostics = summarize(resolved)
    comparison = load_fig4b_comparison(fig4b_summary_path, summary)

    with pd.ExcelWriter(summary_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="summary_counts", index=False)
        diagnostics.to_excel(writer, sheet_name="diagnostics", index=False)
        comparison.to_excel(writer, sheet_name="compare_to_Fig4B", index=False)
        unique_labels.to_excel(writer, sheet_name="unique_true_labels", index=False)
        resolved.to_excel(writer, sheet_name="resolved_haplotype_rows", index=False)

    plot_summary(summary, args, output_path, png_path)

    print("New Fig.4C summary")
    print(summary.to_string(index=False))
    print("\nDiagnostics")
    print(diagnostics.to_string(index=False))
    print("\nComparison to Fig.4B")
    print(comparison.to_string(index=False))
    print(f"\nOutput PDF: {output_path}")
    print(f"Output PNG: {png_path}")
    print(f"Output summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
