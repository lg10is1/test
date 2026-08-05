# -*- coding: utf-8 -*-
"""Count new HLA typing entries by digit level.

The script reads A/B/C worksheets from a new-subtype workbook and classifies
entries containing ``new`` by the number of colon-delimited fields:

- two fields -> 2-digit
- three fields -> 3-digit
- four fields -> 4-digit

Example:
python novel_typing_digit_distribution_universal.py new_mhc_subtypes.xlsx
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


DEFAULT_SHEETS = {
    "A": "A",
    "B": "B",
    "C": "C",
}

DEFAULT_COLORS = {
    "A": "#47a1a2",
    "B": "#da7271",
    "C": "#1f78b4",
}

CATEGORIES = ["2-digit", "3-digit", "4-digit"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count Fig. 4C digit-level distributions from A/B/C worksheets in a new-HLA-subtype workbook."
    )
    parser.add_argument("input_xlsx", help="Input workbook path.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output plot path. By default, a PDF is written next to this script.",
    )
    parser.add_argument(
        "--summary-output",
        help="Optional path for the summary workbook.",
    )
    parser.add_argument(
        "--sheet-a",
        default=DEFAULT_SHEETS["A"],
        help=f"Worksheet for HLA-A. Default: {DEFAULT_SHEETS['A']}.",
    )
    parser.add_argument(
        "--sheet-b",
        default=DEFAULT_SHEETS["B"],
        help=f"Worksheet for HLA-B. Default: {DEFAULT_SHEETS['B']}.",
    )
    parser.add_argument(
        "--sheet-c",
        default=DEFAULT_SHEETS["C"],
        help=f"Worksheet for HLA-C. Default: {DEFAULT_SHEETS['C']}.",
    )
    parser.add_argument(
        "--title",
        default="Distribution of Novel HLA Alleles by Digit Groups",
        help="Plot title.",
    )
    parser.add_argument(
        "--xlabel",
        default="New Subtype Digits",
        help="X-axis title.",
    )
    parser.add_argument(
        "--ylabel",
        default="Haplotype Count",
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
        default=(10, 6),
        metavar=("WIDTH", "HEIGHT"),
        help="Figure size in inches. Default: 10 6.",
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


def normalize_value(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def base_column_name(column_name: object) -> str:
    return re.sub(r"\.\d+$", "", str(column_name)).strip()


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
    return candidate_columns[0]


def load_sheet(data_path: Path, sheet_name: str) -> pd.DataFrame:
    dataframe = pd.read_excel(data_path, sheet_name=sheet_name, dtype=object)
    dataframe = dataframe.dropna(how="all")
    return dataframe


def classify_digit_group(value: str) -> str | None:
    if "new" not in value.lower():
        return None

    parts = value.split(":")
    if len(parts) == 2:
        return "2-digit"
    if len(parts) == 3:
        return "3-digit"
    if len(parts) == 4:
        return "4-digit"
    return None


def summarize_gene(dataframe: pd.DataFrame, gene_name: str, source_sheet: str) -> dict[str, object]:
    gene_column = pick_gene_column(dataframe, gene_name)
    values = dataframe[gene_column].map(normalize_value)

    counts = {category: 0 for category in CATEGORIES}
    for value in values:
        group = classify_digit_group(value)
        if group:
            counts[group] += 1

    return {
        "Gene": gene_name,
        "SourceSheet": source_sheet,
        "SourceColumn": gene_column,
        **counts,
        "Total_New_Haplotypes": sum(counts.values()),
    }


def build_summary(data_path: Path, sheet_map: dict[str, str]) -> pd.DataFrame:
    summary_rows = []
    for gene_name in ("A", "B", "C"):
        dataframe = load_sheet(data_path, sheet_map[gene_name])
        summary_rows.append(
            summarize_gene(
                dataframe=dataframe,
                gene_name=gene_name,
                source_sheet=sheet_map[gene_name],
            )
        )
    return pd.DataFrame(summary_rows)


def configure_matplotlib(font_family: str) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = font_family


def plot_summary(summary_df: pd.DataFrame, args: argparse.Namespace, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    configure_matplotlib(args.font_family)

    plot_df = summary_df.set_index("Gene")[CATEGORIES]
    hla_types = ["A", "B", "C"]
    x_positions = range(len(CATEGORIES))
    width = 0.25

    plt.figure(figsize=tuple(args.figsize))

    for index, hla_type in enumerate(hla_types):
        heights = [int(plot_df.loc[hla_type, category]) for category in CATEGORIES]
        bars = plt.bar(
            [position + index * width for position in x_positions],
            heights,
            width=width,
            label=f"HLA-{hla_type}",
            color=DEFAULT_COLORS[hla_type],
        )

        for bar in bars:
            y_value = bar.get_height()
            x_value = bar.get_x() + bar.get_width() / 2
            plt.text(x_value, y_value, f"{int(y_value)}", ha="center", va="bottom", fontsize=16, color="black")

    plt.xlabel(args.xlabel, fontsize=30)
    plt.ylabel(args.ylabel, fontsize=30)
    plt.title(args.title, fontsize=27)
    plt.xticks([position + width for position in x_positions], CATEGORIES, fontsize=18)
    plt.legend(fontsize=18)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        output_path,
        format=output_path.suffix.lstrip("."),
        bbox_inches="tight",
        dpi=args.dpi,
        metadata={
            "Creator": "",
            "Producer": "",
            "CreationDate": None,
            "ModDate": None,
            "Title": "Matplotlib PDF",
        },
        facecolor="w",
        edgecolor="none",
    )

    if args.show:
        plt.show()
    else:
        plt.close()


def default_output_path(input_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return script_dir / f"hla_typing_distribution_by_digits_{input_path.stem}_universal.pdf"


def default_summary_output_path(input_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return script_dir / f"hla_typing_distribution_by_digits_{input_path.stem}_summary_universal.xlsx"


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

    summary_df = build_summary(input_path, sheet_map)

    print("Fig.4C summary")
    print(summary_df.to_string(index=False))
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
