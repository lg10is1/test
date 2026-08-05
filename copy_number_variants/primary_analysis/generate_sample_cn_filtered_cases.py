# -*- coding: utf-8 -*-
"""Generate sample-level CNV matrices from haplotype-level copy numbers.

The input matrix is expected to have haplotype identifiers as its index and
genes as columns.  Haplotype identifiers are paired by the historical naming
patterns in the original workflow.  Copy numbers are merged to the sample
level by taking the maximum value per haplotype label and summing the two
haplotypes for each sample.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


INDEX_PATTERNS = (
    (re.compile(r"^(?P<sample>.+?)\.1\.mat_R1_t2t\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>.+?)\.2\.pat_R1_t2t\.2\.scaffold$"), "hap2"),
    (re.compile(r"^(?P<sample>.+?)\.1\.mat\.R1_t2t\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>.+?)\.2\.pat\.R1_t2t\.2\.scaffold$"), "hap2"),
    (re.compile(r"^(?P<sample>.+?)\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>.+?)\.2\.scaffold$"), "hap2"),
)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Convert a haplotype-level CNV matrix to sample-level outputs."
    )
    parser.add_argument(
        "input_xlsx",
        nargs="?",
        type=Path,
        default=script_dir / "haplotype_CN.xlsx",
        help="Input haplotype-level CN matrix.",
    )
    parser.add_argument(
        "--sample-output",
        type=Path,
        default=script_dir / "sample_CN.xlsx",
        help="Output workbook containing the sample-level CN matrix.",
    )
    parser.add_argument(
        "--presence-output",
        type=Path,
        default=script_dir / "sample_CN_presence.xlsx",
        help="Output workbook containing binary sample presence summaries.",
    )
    return parser.parse_args()


def parse_index(index_label: str) -> tuple[str, str]:
    for pattern, haplotype in INDEX_PATTERNS:
        match = pattern.match(str(index_label))
        if match:
            return match.group("sample"), haplotype
    raise ValueError(f"Unrecognized haplotype index format: {index_label}")


def load_haplotype_cn(input_path: Path) -> pd.DataFrame:
    if not input_path.is_file():
        raise FileNotFoundError(f"Input workbook not found: {input_path}")
    matrix = pd.read_excel(input_path, index_col=0)
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def build_metadata(index_labels: pd.Index) -> pd.DataFrame:
    rows = []
    for index_label in index_labels:
        sample, haplotype = parse_index(str(index_label))
        rows.append(
            {
                "raw_index": str(index_label),
                "sample": sample,
                "haplotype": haplotype,
            }
        )
    return pd.DataFrame(rows)


def merge_haplotypes(
    haplotype_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
) -> pd.DataFrame:
    matrix = haplotype_df.copy()
    matrix.index = pd.MultiIndex.from_frame(metadata_df[["sample", "haplotype"]])
    haplotype_level_df = matrix.groupby(level=[0, 1]).max()
    sample_cn_df = haplotype_level_df.groupby(level=0).sum()
    sample_cn_df.index.name = "sample"
    return sample_cn_df.astype(int)


def build_sum_ecn(sample_cn_df: pd.DataFrame) -> pd.DataFrame:
    sum_ecn = sample_cn_df.sum(axis=0)
    return pd.DataFrame({"Gene": sum_ecn.index, "Sum_ECN_true": sum_ecn.values})


def build_binary_sample_outputs(
    sample_cn_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    binary_df = (sample_cn_df > 0).astype(int)
    binary_df.index.name = "sample"

    counts = binary_df.sum(axis=0)
    raw_frequency = (
        counts.rename("Count")
        .reset_index()
        .rename(columns={"index": "Sample", "sample": "Sample"})
        .sort_values(["Count", "Sample"], ascending=[False, True])
        .reset_index(drop=True)
    )
    raw_frequency["Freq"] = raw_frequency["Count"] / binary_df.shape[0]

    count_distribution = (
        raw_frequency["Count"]
        .value_counts()
        .rename_axis("Count")
        .reset_index(name="Sample_N")
        .sort_values("Count")
        .reset_index(drop=True)
    )
    count_distribution["Freq"] = count_distribution["Sample_N"] / binary_df.shape[1]
    return binary_df, raw_frequency, count_distribution


def write_outputs(
    sample_cn_df: pd.DataFrame,
    sample_output: Path,
    presence_output: Path,
) -> None:
    sample_output.parent.mkdir(parents=True, exist_ok=True)
    presence_output.parent.mkdir(parents=True, exist_ok=True)

    sum_ecn_df = build_sum_ecn(sample_cn_df)
    with pd.ExcelWriter(sample_output, engine="openpyxl") as writer:
        sample_cn_df.to_excel(writer, sheet_name="sample_CN")
        sum_ecn_df.to_excel(writer, sheet_name="Sum_ECN", index=False)

    binary_df, raw_frequency_df, count_distribution_df = build_binary_sample_outputs(
        sample_cn_df
    )
    with pd.ExcelWriter(presence_output, engine="openpyxl") as writer:
        binary_df.to_excel(writer, sheet_name="presence_matrix")
        raw_frequency_df.to_excel(
            writer,
            sheet_name="gene_sample_frequency_raw",
            index=False,
        )
        count_distribution_df.to_excel(
            writer,
            sheet_name="gene_sample_frequency_distribution",
            index=False,
        )


def main() -> None:
    args = parse_args()
    haplotype_df = load_haplotype_cn(args.input_xlsx)
    metadata_df = build_metadata(haplotype_df.index)
    sample_cn_df = merge_haplotypes(haplotype_df, metadata_df)
    write_outputs(
        sample_cn_df,
        sample_output=args.sample_output,
        presence_output=args.presence_output,
    )

    print(f"Haplotypes: {len(haplotype_df)}")
    print(f"Samples: {len(sample_cn_df)}")
    print(f"Genes: {sample_cn_df.shape[1]}")
    print(f"Saved: {args.sample_output}")
    print(f"Saved: {args.presence_output}")


if __name__ == "__main__":
    main()
