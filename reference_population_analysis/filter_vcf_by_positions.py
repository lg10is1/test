#!/usr/bin/env python3
"""Filter a VCF by exact CHROM and POS values without loading it into memory."""

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path
from typing import Dict, Set, TextIO


def open_text(path: Path, mode: str) -> TextIO:
    """Open plain or gzip-compressed UTF-8 text based on the filename suffix."""
    if path.name.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def load_positions(position_file: Path) -> Dict[str, Set[int]]:
    """Load a tab-separated CHROM/POS file, allowing blank lines and comments."""
    positions: Dict[str, Set[int]] = {}
    with position_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(f"Expected CHROM and POS at {position_file}:{line_number}")
            chromosome = row[0].strip()
            try:
                position = int(row[1])
            except ValueError as exc:
                raise ValueError(f"Invalid POS at {position_file}:{line_number}: {row[1]!r}") from exc
            if position < 1:
                raise ValueError(f"POS must be positive at {position_file}:{line_number}")
            positions.setdefault(chromosome, set()).add(position)
    if not positions:
        raise ValueError(f"No positions were loaded from {position_file}")
    return positions


def filter_vcf(input_vcf: Path, positions: Dict[str, Set[int]], output_vcf: Path) -> int:
    """Write matching records atomically and return the number retained."""
    if input_vcf.resolve() == output_vcf.resolve():
        raise ValueError("Input and output VCF paths must differ")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".vcf.gz" if output_vcf.name.endswith(".gz") else ".vcf"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_vcf.name}.", suffix=suffix, dir=str(output_vcf.parent)
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    retained = 0
    saw_header = False
    try:
        with open_text(input_vcf, "rt") as source, open_text(temporary_path, "wt") as destination:
            for line_number, line in enumerate(source, start=1):
                if line.startswith("#"):
                    destination.write(line)
                    saw_header = True
                    continue
                if not line.strip():
                    continue
                fields = line.split("\t", 2)
                if len(fields) < 2:
                    raise ValueError(f"Malformed VCF record at {input_vcf}:{line_number}")
                try:
                    position = int(fields[1])
                except ValueError as exc:
                    raise ValueError(f"Invalid VCF POS at {input_vcf}:{line_number}") from exc
                if position in positions.get(fields[0], set()):
                    destination.write(line)
                    retained += 1
        if not saw_header:
            raise ValueError(f"No VCF header was found in {input_vcf}")
        os.replace(str(temporary_path), str(output_vcf))
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    return retained


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter VCF records by exact CHROM/POS pairs.")
    parser.add_argument("-p", "--position-file", required=True, type=Path, help="Tab-separated CHROM and POS file.")
    parser.add_argument("-v", "--vcf-file", required=True, type=Path, help="Input .vcf or .vcf.gz file.")
    parser.add_argument("-o", "--output-vcf", required=True, type=Path, help="Output .vcf or .vcf.gz file.")
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Accepted for source CLI compatibility; filtering is streaming and single-process.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    if not args.position_file.is_file():
        raise FileNotFoundError(f"Position file not found: {args.position_file}")
    if not args.vcf_file.is_file():
        raise FileNotFoundError(f"Input VCF not found: {args.vcf_file}")
    positions = load_positions(args.position_file)
    retained = filter_vcf(args.vcf_file, positions, args.output_vcf)
    print(f"Retained VCF records: {retained}")
    print(f"Filtered VCF: {args.output_vcf}")


if __name__ == "__main__":
    main()
