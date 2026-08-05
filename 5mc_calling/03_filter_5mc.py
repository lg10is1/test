#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def arguments():
    parser = argparse.ArgumentParser(
        description="Filter merged CpG sites by sample missingness."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--max-missing", type=float, default=0.30)
    return parser.parse_args()


def main():
    args = arguments()
    prefix = args.chrom
    args.output_dir.mkdir(parents=True, exist_ok=True)

    col4 = pd.read_csv(
        args.input_dir / f"{prefix}_col4_val.tsv",
        sep="\t",
        dtype=np.float32,
        na_values="NA",
    )
    keep = col4.notna().mean(axis=1) >= (1 - args.max_missing)

    positions = pd.read_csv(args.input_dir / f"{prefix}_positions.tsv", sep="\t")
    positions.loc[keep].to_csv(
        args.output_dir / f"{prefix}_positions.tsv", sep="\t", index=False
    )

    for label in ("col4", "col6", "col9"):
        matrix = pd.read_csv(
            args.input_dir / f"{prefix}_{label}_val.tsv",
            sep="\t",
            dtype=np.float32,
            na_values="NA",
        )
        matrix.loc[keep].to_csv(
            args.output_dir / f"{prefix}_{label}_val.tsv",
            sep="\t",
            index=False,
            na_rep="NA",
            float_format="%.3f",
        )

    print(f"Retained {int(keep.sum())} of {len(keep)} CpG sites")


if __name__ == "__main__":
    main()

