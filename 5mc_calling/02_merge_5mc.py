#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def arguments():
    parser = argparse.ArgumentParser(
        description="Merge per-sample pb-CpG-tools BED values for one chromosome."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chrom", required=True)
    return parser.parse_args()


def sample_name(path):
    name = path.name
    for suffix in (".combined.bed.gz", ".combined.bed"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_columns(path, columns, chrom):
    frame = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        usecols=columns,
        low_memory=False,
    )
    return frame.loc[frame[0].astype(str) == chrom]


def write_matrix(values, samples, path, chunk_size=20000):
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(samples) + "\n")
        for start in range(0, values.shape[0], chunk_size):
            block = values[start : start + chunk_size]
            text = np.char.mod("%.3f", block)
            text[np.isnan(block)] = "NA"
            handle.writelines("\t".join(row) + "\n" for row in text)


def main():
    args = arguments()
    bed_files = sorted(args.input_dir.rglob("*.combined.bed.gz"))
    if not bed_files:
        bed_files = sorted(args.input_dir.rglob("*.combined.bed"))

    samples = [sample_name(path) for path in bed_files]
    coordinates = set()

    for path in bed_files:
        frame = read_columns(path, [0, 1, 2], args.chrom)
        coordinates.update(zip(frame[1], frame[2]))

    master = sorted(coordinates)
    master_index = pd.MultiIndex.from_tuples(master, names=["start", "end"])
    args.output_dir.mkdir(parents=True, exist_ok=True)

    positions = pd.DataFrame(master, columns=["start", "end"])
    positions.insert(0, "chr", args.chrom)
    positions.to_csv(
        args.output_dir / f"{args.chrom}_positions.tsv", sep="\t", index=False
    )

    for column, label in ((3, "col4"), (5, "col6"), (8, "col9")):
        matrix = np.full((len(master), len(bed_files)), np.nan, dtype=np.float32)

        for sample_index, path in enumerate(bed_files):
            frame = read_columns(path, [0, 1, 2, column], args.chrom)
            current = pd.MultiIndex.from_arrays([frame[1], frame[2]])
            row_index = master_index.get_indexer(current)
            valid = row_index >= 0
            matrix[row_index[valid], sample_index] = frame.loc[valid, column].to_numpy(
                dtype=np.float32
            )

        write_matrix(
            matrix,
            samples,
            args.output_dir / f"{args.chrom}_{label}_val.tsv",
        )


if __name__ == "__main__":
    main()

