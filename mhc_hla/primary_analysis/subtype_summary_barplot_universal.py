# -*- coding: utf-8 -*-
"""Draw the Fig. 4B HLA subtype summary bar plot.

This script uses the same manually QC-filtered A/B/C worksheets as Fig. 4A.
Known subtypes are counted as unique non-empty, non-NULL values that do not
contain ``new``.

Important: the new-subtype count is not fully automatic.  The historical
analysis treated some repeated ``new`` haplotypes as the same new subtype only
after manual sequence review.  Therefore:

new subtype count = new haplotype count - manually reviewed duplicate excess

The duplicate excess must be supplied per gene by
``--new-duplicate-reduction-a/b/c``.  The public-release default worksheet
names are English; historical Chinese worksheet names are accepted as
fallbacks when the English sheets are absent. The exact historical worksheet
mapping is documented in CHINESE_SCHEMA_AND_RENAMING_IMPACT_V3.md.

Example:
python subtype_summary_barplot_universal.py immuannot_scz.xlsx \
    --new-duplicate-reduction-a 1 \
    --new-duplicate-reduction-b 8 \
    --new-duplicate-reduction-c 7
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


PUBLIC_QC_SHEETS = {
    "A": "A_qc_passed",
    "B": "B_qc_passed",
    "C": "C_qc_passed",
}

LEGACY_QC_SHEETS = {
    "A": "A_validated",
    "B": "B_validated",
    "C": "C_validated",
}

DEFAULT_SHEETS = PUBLIC_QC_SHEETS
DEFAULT_COLORS = ["#1f78b4", "#4DBBD5FF"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count Known/New HLA subtypes from manually QC-filtered Immuannot worksheets."
    )
    parser.add_argument("input_xlsx", help="Path to the Immuannot workbook.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output plot path. By default, a PDF is written next to this script.",
    )
    parser.add_argument(
        "--sheet-a",
        default=DEFAULT_SHEETS["A"],
        help=f"Worksheet for HLA-A. Default: {DEFAULT_SHEETS['A']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--sheet-b",
        default=DEFAULT_SHEETS["B"],
        help=f"Worksheet for HLA-B. Default: {DEFAULT_SHEETS['B']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--sheet-c",
        default=DEFAULT_SHEETS["C"],
        help=f"Worksheet for HLA-C. Default: {DEFAULT_SHEETS['C']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--new-duplicate-reduction-a",
        type=int,
        default=0,
        help="Manual duplicate excess to subtract from HLA-A new haplotypes. Default: 0.",
    )
    parser.add_argument(
        "--new-duplicate-reduction-b",
        type=int,
        default=0,
        help="Manual duplicate excess to subtract from HLA-B new haplotypes. Default: 0.",
    )
    parser.add_argument(
        "--new-duplicate-reduction-c",
        type=int,
        default=0,
        help="Manual duplicate excess to subtract from HLA-C new haplotypes. Default: 0.",
    )
    parser.add_argument(
        "--summary-output",
        help="Optional path for the summary workbook.",
    )
    parser.add_argument(
        "--title",
        default="Subtype Number of HLA-A/B/C",
        help="Plot title.",
    )
    parser.add_argument(
        "--xlabel",
        default="HLA Genes",
        help="X-axis title.",
    )
    parser.add_argument(
        "--ylabel",
        default="Number",
        help="Y-axis title.",
    )
    parser.add_argument(
        "--font-family",
        default="Arial",
        help="Matplotlib font family. Default: Arial.",
    )
    parser.add_argument(
        "--figsize",
        nargs=2,
        type=float,
        default=(8, 6),
        metavar=("WIDTH", "HEIGHT"),
        help="Figure size in inches. Default: 8 6.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Output DPI. Default: 600.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot window after saving.",
    )
    return parser.parse_args()


def base_column_name(column_name: object) -> str:
    return re.sub(r"\.\d+$", "", str(column_name)).strip()


def normalize_value(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def is_valid_value(value: str) -> bool:
    return bool(value) and value.upper() not in {"NULL", "UNKNOWN", "(blank)"}


def count_nonempty_values(series: pd.Series) -> int:
    return sum(is_valid_value(normalize_value(value)) for value in series)


def pick_gene_column(dataframe: pd.DataFrame, gene_name: str) -> str:
    candidate_columns = [
        column_name
        for column_name in dataframe.columns
        if base_column_name(column_name).upper() == gene_name.upper()
    ]

    if not candidate_columns:
        raise ValueError(
            f"Could not find gene column {gene_name!r}. Available columns: {list(dataframe.columns)}"
        )

    candidate_columns.sort(
        key=lambda column_name: count_nonempty_values(dataframe[column_name]),
        reverse=True,
    )
    return candidate_columns[0]


def sheet_candidates_for_gene(gene_name: str, preferred_sheet: str) -> list[str]:
    candidates = [preferred_sheet]
    if preferred_sheet == PUBLIC_QC_SHEETS[gene_name]:
        candidates.append(LEGACY_QC_SHEETS[gene_name])
    return list(dict.fromkeys(candidate for candidate in candidates if candidate))


def resolve_sheet_name(data_path: Path, candidates: list[str]) -> str:
    workbook = pd.ExcelFile(data_path)
    for sheet_name in candidates:
        if sheet_name in workbook.sheet_names:
            return sheet_name
    raise ValueError(
        f"None of the requested worksheets were found: {candidates}. "
        f"Available worksheets: {workbook.sheet_names}"
    )


def load_sheet(data_path: Path, sheet_name: str) -> pd.DataFrame:
    dataframe = pd.read_excel(data_path, sheet_name=sheet_name, dtype=object)
    dataframe = dataframe.dropna(how="all")
    return dataframe


def summarize_gene(
    dataframe: pd.DataFrame,
    gene_name: str,
    source_sheet: str,
    duplicate_reduction: int,
) -> dict[str, object]:
    gene_column = pick_gene_column(dataframe, gene_name)
    values = dataframe[gene_column].map(normalize_value)
    values = values[values.map(is_valid_value)].astype(str)

    known_values = values[~values.str.contains("new", case=False, na=False)]
    new_values = values[values.str.contains("new", case=False, na=False)]

    known_subtype_count = int(known_values.nunique())
    new_haplotype_count = int(len(new_values))
    new_unique_string_count = int(new_values.nunique())

    if duplicate_reduction < 0:
        raise ValueError(f"{gene_name} duplicate reduction cannot be negative.")
    if duplicate_reduction > new_haplotype_count:
        raise ValueError(
            f"{gene_name} duplicate reduction {duplicate_reduction} exceeds "
            f"new haplotype count {new_haplotype_count}."
        )

    final_new_subtype_count = new_haplotype_count - duplicate_reduction

    return {
        "Gene": gene_name,
        "SourceSheet": source_sheet,
        "SourceColumn": gene_column,
        "Known_Haplotype_Count": int(len(known_values)),
        "Known_Subtype_Count": known_subtype_count,
        "New_Haplotype_Count": new_haplotype_count,
        "New_Unique_String_Count": new_unique_string_count,
        "Manual_Duplicate_Reduction": duplicate_reduction,
        "Final_New_Subtype_Count": final_new_subtype_count,
    }


def build_summary(data_path: Path, sheet_map: dict[str, str], reduction_map: dict[str, int]) -> pd.DataFrame:
    summary_rows = []
    for gene_name in ("A", "B", "C"):
        source_sheet = resolve_sheet_name(
            data_path,
            sheet_candidates_for_gene(gene_name, sheet_map[gene_name]),
        )
        dataframe = load_sheet(data_path, source_sheet)
        summary_rows.append(
            summarize_gene(
                dataframe=dataframe,
                gene_name=gene_name,
                source_sheet=source_sheet,
                duplicate_reduction=reduction_map[gene_name],
            )
        )
    return pd.DataFrame(summary_rows)


def configure_matplotlib(font_family: str) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = font_family


def annotate_bars(axis, known_values: pd.Series, new_values: pd.Series) -> None:
    for index, (known_count, new_count) in enumerate(zip(known_values, new_values)):
        axis.text(index, known_count / 2, str(int(known_count)), ha="center", va="bottom", fontsize=18, color="black")
        axis.text(index, known_count + new_count / 2, str(int(new_count)), ha="center", va="bottom", fontsize=18, color="black")


def plot_summary(summary_df: pd.DataFrame, args: argparse.Namespace, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    configure_matplotlib(args.font_family)

    data = {
        row["Gene"]: [row["Known_Subtype_Count"], row["Final_New_Subtype_Count"]]
        for _, row in summary_df.iterrows()
    }
    plot_df = pd.DataFrame(data, index=["Known Subtypes", "New Subtypes"])

    fig, axis = plt.subplots(figsize=tuple(args.figsize))

    known_series = plot_df.loc["Known Subtypes"]
    new_series = plot_df.loc["New Subtypes"]

    known_series.plot(kind="bar", stacked=True, color=DEFAULT_COLORS[0], ax=axis, label="Known Subtypes")
    new_series.plot(kind="bar", stacked=True, color=DEFAULT_COLORS[1], ax=axis, bottom=known_series, label="New Subtypes")

    axis.set_xticks(range(len(plot_df.columns)))
    axis.set_xticklabels(plot_df.columns, rotation=0, fontsize=30)
    axis.set_ylabel(args.ylabel, fontsize=30)
    axis.set_xlabel(args.xlabel, fontsize=30)
    axis.set_title(args.title, fontsize=32)
    axis.legend(fontsize=16)

    annotate_bars(axis, known_series, new_series)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=args.dpi, format=output_path.suffix.lstrip("."))

    if args.show:
        plt.show()
    else:
        plt.close()


def default_output_path(input_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return script_dir / f"Subtype_Stat_{input_path.stem}_universal.pdf"


def default_summary_output_path(input_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return script_dir / f"Subtype_Stat_{input_path.stem}_summary_universal.xlsx"


def main() -> int:
    args = parse_args()

    input_path = Path(args.input_xlsx).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input workbook does not exist: {input_path}")

    output_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else default_output_path(input_path)
    )

    summary_output_path = (
        Path(args.summary_output).expanduser().resolve()
        if args.summary_output
        else default_summary_output_path(input_path)
    )

    sheet_map = {
        "A": args.sheet_a,
        "B": args.sheet_b,
        "C": args.sheet_c,
    }
    reduction_map = {
        "A": args.new_duplicate_reduction_a,
        "B": args.new_duplicate_reduction_b,
        "C": args.new_duplicate_reduction_c,
    }

    summary_df = build_summary(
        data_path=input_path,
        sheet_map=sheet_map,
        reduction_map=reduction_map,
    )

    print("Fig.4B summary")
    print(summary_df.to_string(index=False))
    print("\nImportant rule:")
    print("Final_New_Subtype_Count = New_Haplotype_Count - Manual_Duplicate_Reduction")
    print("Manual_Duplicate_Reduction must come from sequence-level manual review.")
    print(f"\nOutput plot: {output_path}")
    print(f"Output summary: {summary_output_path}")

    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_excel(summary_output_path, index=False)

    plot_summary(summary_df, args, output_path)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
