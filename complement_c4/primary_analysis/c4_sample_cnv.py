# -*- coding: utf-8 -*-
"""Aggregate C4 haplotype structures at the sample level.

Each haplotype row contains a ``haplotype`` identifier and a hyphen-delimited
``structure`` value.  The script extracts the sample identifier from the
haplotype label, groups both haplotypes from the same sample, and reports the
sample-level counts of ``C4AL``, ``C4BL``, ``C4AS``, and ``C4BS``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_HAPLOTYPES = ("C4AL", "C4BL", "C4AS", "C4BS")
DEFAULT_SAMPLE_PATTERNS = {
    "scz": r"^(.*?)-\d+\.\d+\.align\.fasta$",
    "comparison": r"^(.*?)\.\d+\.v0\.9$",
    "public": r"^(.*?)\.\d+$",
}
PUBLIC_DEFAULT_SHEET = "sample_haplotype_map_nonempty"
LEGACY_DEFAULT_SHEET = "renamed_samples_without_empty_haplotypes"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate C4 structure counts from haplotype rows to sample rows."
    )
    parser.add_argument("input_xlsx", type=Path, help="Input workbook.")
    parser.add_argument(
        "-s",
        "--sheet",
        default=PUBLIC_DEFAULT_SHEET,
        help=(
            "Worksheet containing haplotype and structure columns. Default: "
            f"{PUBLIC_DEFAULT_SHEET}; legacy fallback documented separately."
        ),
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output workbook.")
    parser.add_argument(
        "--haplotype-column",
        default="haplotype",
        help="Column containing haplotype identifiers.",
    )
    parser.add_argument(
        "--structure-column",
        default="structure",
        help="Column containing hyphen-delimited C4 structures.",
    )
    parser.add_argument(
        "--cohort-pattern",
        choices=sorted(DEFAULT_SAMPLE_PATTERNS),
        default="public",
        help=(
            "Preset regular expression for extracting sample identifiers. "
            "Use --sample-regex to override this."
        ),
    )
    parser.add_argument(
        "--sample-regex",
        help=(
            "Custom regular expression with one capture group for extracting "
            "sample identifiers from haplotype labels."
        ),
    )
    parser.add_argument(
        "--delimiter",
        default="-",
        help="Delimiter used in the structure column.",
    )
    return parser.parse_args()


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


def count_structures(values: pd.Series, delimiter: str) -> dict[str, int]:
    counts = {haplotype: 0 for haplotype in DEFAULT_HAPLOTYPES}
    for value in values:
        if pd.isna(value):
            continue
        for token in str(value).split(delimiter):
            if token in counts:
                counts[token] += 1
    return counts


def aggregate_sample_counts(
    dataframe: pd.DataFrame,
    haplotype_column: str,
    structure_column: str,
    sample_regex: str,
    delimiter: str,
) -> pd.DataFrame:
    missing_columns = [
        column
        for column in (haplotype_column, structure_column)
        if column not in dataframe.columns
    ]
    if missing_columns:
        raise ValueError(
            f"Missing required columns {missing_columns}; "
            f"available columns: {list(dataframe.columns)}"
        )

    working = dataframe.copy()
    working["sample"] = working[haplotype_column].astype(str).str.extract(sample_regex)[0]
    missing_samples = working.loc[working["sample"].isna(), haplotype_column].head(20).tolist()
    if missing_samples:
        raise ValueError(
            "Failed to extract sample identifiers from some haplotype labels; "
            f"examples: {missing_samples}"
        )

    rows = []
    for sample, group in working.groupby("sample", sort=True):
        counts = count_structures(group[structure_column], delimiter)
        rows.append({"sample": sample, **counts})

    return pd.DataFrame(rows, columns=["sample", *DEFAULT_HAPLOTYPES])


def main() -> None:
    args = parse_args()
    if not args.input_xlsx.is_file():
        raise FileNotFoundError(f"Input workbook not found: {args.input_xlsx}")

    sample_regex = args.sample_regex or DEFAULT_SAMPLE_PATTERNS[args.cohort_pattern]
    sheet_name = resolve_sheet_name(args.input_xlsx, args.sheet)
    dataframe = pd.read_excel(args.input_xlsx, sheet_name=sheet_name)
    result = aggregate_sample_counts(
        dataframe,
        haplotype_column=args.haplotype_column,
        structure_column=args.structure_column,
        sample_regex=sample_regex,
        delimiter=args.delimiter,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_excel(args.output, index=False)
    print(f"Samples processed: {len(result)}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
