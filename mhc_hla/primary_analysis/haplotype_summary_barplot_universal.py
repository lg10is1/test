# -*- coding: utf-8 -*-
"""Draw the Fig. 4A HLA haplotype summary bar plot.

Known and new haplotypes are counted from manually QC-filtered A/B/C
worksheets in the Immuannot workbook.  The public-release default worksheet
names are English; historical Chinese worksheet names are accepted as
fallbacks when the English sheets are absent. The exact historical worksheet
mapping is documented in CHINESE_SCHEMA_AND_RENAMING_IMPACT_V3.md.

Unknown haplotypes are calculated as:
``total_haplotype_count - known_count - new_count``.

Example:
python haplotype_summary_barplot_universal.py immuannot_scz.xlsx

Example with an explicit total:
python haplotype_summary_barplot_universal.py immuannot_scz.xlsx --total-count 420
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
DEFAULT_COLORS = ["#1f78b4", "#4DBBD5FF", "lightgrey"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count Known/New/Unknown HLA haplotypes from manually QC-filtered Immuannot worksheets."
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
        "--total-count",
        type=int,
        help="Total haplotype count. If omitted, the row count is inferred from the total worksheet.",
    )
    parser.add_argument(
        "--total-sheet",
        default="ABC",
        help="Worksheet used to infer total haplotype count. Default: ABC.",
    )
    parser.add_argument(
        "--title",
        default="Haplotype Number of HLA-A/B/C",
        help="Plot title.",
    )
    parser.add_argument(
        "--xlabel",
        default="Haplotype Count",
        help="X-axis title.",
    )
    parser.add_argument(
        "--ylabel",
        default="HLA Genes",
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
        default=(12, 6),
        metavar=("WIDTH", "HEIGHT"),
        help="Figure size in inches. Default: 12 6.",
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


def count_nonempty_values(series: pd.Series) -> int:
    return sum(bool(normalize_value(value)) for value in series)


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


def count_known_new(dataframe: pd.DataFrame, gene_name: str) -> tuple[int, int]:
    gene_column = pick_gene_column(dataframe, gene_name)
    values = dataframe[gene_column].map(normalize_value)

    known_count = 0
    new_count = 0

    for value in values:
        if not value or value.upper() == "NULL":
            continue
        if "new" in value.lower():
            new_count += 1
        else:
            known_count += 1

    return known_count, new_count


def infer_total_count(data_path: Path, explicit_total: int | None, total_sheet: str) -> int:
    if explicit_total is not None:
        return explicit_total

    candidate_sheets = []
    if total_sheet:
        candidate_sheets.append(total_sheet)
    for fallback_sheet in ("ABC", "combined"):
        if fallback_sheet not in candidate_sheets:
            candidate_sheets.append(fallback_sheet)

    last_error: Exception | None = None
    for sheet_name in candidate_sheets:
        try:
            dataframe = load_sheet(data_path, sheet_name)
            return len(dataframe)
        except Exception as error:
            last_error = error

    try:
        workbook = pd.ExcelFile(data_path)
        if workbook.sheet_names:
            first_sheet_name = workbook.sheet_names[0]
            dataframe = load_sheet(data_path, first_sheet_name)
            return len(dataframe)
    except Exception as error:
        last_error = error

    raise ValueError(
        "Could not infer total haplotype count automatically. Please provide --total-count."
    ) from last_error


def build_summary(data_path: Path, sheet_map: dict[str, str], total_count: int) -> pd.DataFrame:
    summary_rows = []

    for gene_name in ("A", "B", "C"):
        sheet_name = resolve_sheet_name(
            data_path,
            sheet_candidates_for_gene(gene_name, sheet_map[gene_name]),
        )
        dataframe = load_sheet(data_path, sheet_name)
        known_count, new_count = count_known_new(dataframe, gene_name)
        unknown_count = total_count - known_count - new_count

        if unknown_count < 0:
            raise ValueError(
                f"Invalid {gene_name} summary: total={total_count}, known={known_count}, "
                f"new={new_count}; unknown would be negative. Check --total-count and input worksheets."
            )

        summary_rows.append(
            {
                "Gene": gene_name,
                "Known": known_count,
                "New": new_count,
                "Unknown": unknown_count,
                "Total": total_count,
                "SourceSheet": sheet_name,
            }
        )

    return pd.DataFrame(summary_rows)


def configure_matplotlib(font_family: str) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = font_family


def annotate_segments(axis) -> None:
    for patch in axis.patches:
        width = patch.get_width()
        if width <= 0:
            continue
        x_position = patch.get_x() + width / 2
        y_position = patch.get_y() + patch.get_height() / 2
        axis.text(
            x_position,
            y_position,
            f"{int(width)}",
            ha="center",
            va="center",
            fontsize=16,
            color="black",
        )


def plot_summary(summary_df: pd.DataFrame, args: argparse.Namespace, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    configure_matplotlib(args.font_family)

    plot_df = summary_df.set_index("Gene")[["Known", "New", "Unknown"]]

    axis = plot_df.plot(
        kind="barh",
        stacked=True,
        color=DEFAULT_COLORS,
        figsize=tuple(args.figsize),
    )

    axis.set_title(args.title, fontsize=32)
    axis.set_xlabel(args.xlabel, fontsize=30)
    axis.set_ylabel(args.ylabel, fontsize=30)
    axis.tick_params(axis="x", labelsize=18)
    axis.tick_params(axis="y", labelsize=24)
    axis.legend(
        labels=["Known", "New", "Unknown"],
        loc="upper left",
        bbox_to_anchor=(1, 0.6),
        ncol=1,
        fontsize=18,
    )

    annotate_segments(axis)

    plt.subplots_adjust(right=0.8)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=args.dpi, format=output_path.suffix.lstrip("."))

    if args.show:
        plt.show()
    else:
        plt.close()


def default_output_path(input_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    file_name = f"Haplotype_Stat_{input_path.stem}_universal.pdf"
    return script_dir / file_name


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

    sheet_map = {
        "A": args.sheet_a,
        "B": args.sheet_b,
        "C": args.sheet_c,
    }

    total_count = infer_total_count(
        data_path=input_path,
        explicit_total=args.total_count,
        total_sheet=args.total_sheet,
    )

    summary_df = build_summary(
        data_path=input_path,
        sheet_map=sheet_map,
        total_count=total_count,
    )

    print("Fig.4A summary")
    print(summary_df.to_string(index=False))
    print(f"\nOutput: {output_path}")

    plot_summary(summary_df, args, output_path)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
