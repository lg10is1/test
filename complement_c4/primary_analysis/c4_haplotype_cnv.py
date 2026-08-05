# -*- coding: utf-8 -*-
"""Count C4 haplotype structures for each input row.

The original Fig. 4F C4 workflow counted occurrences of ``C4AL``, ``C4BL``,
``C4AS``, and ``C4BS`` in a hyphen-delimited ``structure`` column.  This
public-release version keeps the counting rule unchanged and exposes the
input workbook, sheet name, and output workbook as command-line arguments.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_HAPLOTYPES = ("C4AL", "C4BL", "C4AS", "C4BS")
PUBLIC_DEFAULT_SHEET = "mat_pat_format"
LEGACY_DEFAULT_SHEET = "maternal_paternal_format"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count C4 haplotype occurrences in each input row."
    )
    parser.add_argument("input_xlsx", type=Path, help="Input workbook.")
    parser.add_argument(
        "-s",
        "--sheet",
        default=PUBLIC_DEFAULT_SHEET,
        help=(
            "Worksheet containing the structure column. Default: "
            f"{PUBLIC_DEFAULT_SHEET}; legacy fallback documented separately."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output workbook.",
    )
    parser.add_argument(
        "--structure-column",
        default="structure",
        help="Column containing hyphen-delimited C4 structures.",
    )
    parser.add_argument(
        "--delimiter",
        default="-",
        help="Delimiter used in the structure column.",
    )
    return parser.parse_args()


def count_row_haplotypes(
    value: object,
    haplotypes: tuple[str, ...],
    delimiter: str,
) -> dict[str, int]:
    counts = {haplotype: 0 for haplotype in haplotypes}
    if pd.isna(value):
        return counts

    for token in str(value).split(delimiter):
        if token in counts:
            counts[token] += 1
    return counts


def resolve_sheet_name(input_xlsx: Path, requested_sheet: str) -> str:
    workbook = pd.ExcelFile(input_xlsx)
    candidates = [requested_sheet]
    if requested_sheet == PUBLIC_DEFAULT_SHEET:
        candidates.append(LEGACY_DEFAULT_SHEET)
    for sheet_name in dict.fromkeys(candidates):
        if sheet_name in workbook.sheet_names:
            return sheet_name
    raise ValueError(
        f"None of the requested worksheets were found: {candidates}. "
        f"Available worksheets: {workbook.sheet_names}"
    )


def add_haplotype_counts(
    dataframe: pd.DataFrame,
    structure_column: str,
    haplotypes: tuple[str, ...],
    delimiter: str,
) -> pd.DataFrame:
    if structure_column not in dataframe.columns:
        raise ValueError(
            f"Missing required column {structure_column!r}; "
            f"available columns: {list(dataframe.columns)}"
        )

    row_counts = [
        count_row_haplotypes(value, haplotypes, delimiter)
        for value in dataframe[structure_column]
    ]
    result = dataframe.copy()
    for haplotype in haplotypes:
        result[haplotype] = [row[haplotype] for row in row_counts]
    return result


def main() -> None:
    args = parse_args()
    if not args.input_xlsx.is_file():
        raise FileNotFoundError(f"Input workbook not found: {args.input_xlsx}")

    sheet_name = resolve_sheet_name(args.input_xlsx, args.sheet)
    dataframe = pd.read_excel(args.input_xlsx, sheet_name=sheet_name)
    result = add_haplotype_counts(
        dataframe,
        structure_column=args.structure_column,
        haplotypes=DEFAULT_HAPLOTYPES,
        delimiter=args.delimiter,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_excel(args.output, index=False)
    print(f"Rows processed: {len(result)}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
