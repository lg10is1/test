#!/usr/bin/env python3
"""Merge single-sample inversion calls and compare carrier rates.

The clustering and statistical model follow the source script from the EOSCZ
collection. Missing sample/cluster combinations are treated as homozygous
reference, and the association test is a one-sided Fisher exact test for case
enrichment. Confirm those assumptions against the approved analysis plan.
"""

import argparse
import gzip
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd


def calculate_reciprocal_overlap(s1: int, e1: int, s2: int, e2: int) -> float:
    """Return the smaller fraction of each interval covered by the overlap."""
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 <= 0 or len2 <= 0:
        return 0.0
    overlap_len = min(e1, e2) - max(s1, s2)
    if overlap_len <= 0:
        return 0.0
    return min(overlap_len / len1, overlap_len / len2)


def parse_gt(gt_str: str) -> int:
    """Return the number of ALT allele 1 copies in a diploid GT string."""
    if not gt_str or gt_str == "." or "." in gt_str.replace("|", "/").split("/"):
        return 0
    return sum(allele == "1" for allele in gt_str.replace("|", "/").split("/"))


def load_sample_list(file_path: Path) -> List[str]:
    """Load a non-empty, unique, one-sample-per-line manifest."""
    if not file_path.is_file():
        raise FileNotFoundError(f"Sample list not found: {file_path}")
    samples = [line.strip() for line in file_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not samples:
        raise ValueError(f"Sample list is empty: {file_path}")
    duplicates = sorted({sample for sample in samples if samples.count(sample) > 1})
    if duplicates:
        preview = ", ".join(duplicates[:5])
        raise ValueError(f"Duplicate sample IDs in {file_path}: {preview}")
    return samples


def hypergeom_pmf(k: int, population: int, successes: int, draws: int) -> float:
    """Calculate a hypergeometric probability mass."""
    try:
        numerator = math.comb(successes, k) * math.comb(population - successes, draws - k)
        return numerator / math.comb(population, draws)
    except (ValueError, ZeroDivisionError, OverflowError):
        return 0.0


def fisher_exact_onesided_greater(a: int, b: int, c: int, d: int) -> float:
    """Calculate a one-sided Fisher exact p-value for case enrichment."""
    population = a + b + c + d
    case_total = a + b
    carrier_total = a + c
    upper = min(case_total, carrier_total)
    p_value = sum(
        hypergeom_pmf(i, population, carrier_total, case_total)
        for i in range(a, upper + 1)
    )
    return min(p_value, 1.0)


def calculate_fdr(p_values: Sequence[float]) -> List[float]:
    """Apply the Benjamini-Hochberg false-discovery-rate correction."""
    n_values = len(p_values)
    if n_values == 0:
        return []
    ranked = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [0.0] * n_values
    running_min = 1.0
    for rank_index in range(n_values - 1, -1, -1):
        original_index, p_value = ranked[rank_index]
        candidate = p_value * n_values / (rank_index + 1)
        running_min = min(running_min, candidate)
        adjusted[original_index] = min(running_min, 1.0)
    return adjusted


def open_text(path: Path):
    """Open plain-text or gzip-compressed input as UTF-8 text."""
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def sample_name_from_path(path: Path) -> str:
    """Derive the sample name from an *_filtered.vcf[.gz] filename."""
    name = path.name
    for suffix in ("_filtered.vcf.gz", "_filtered.vcf"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Unexpected input filename: {path}")


def parse_end(position: int, info_text: str) -> int:
    """Read END or SVLEN from VCF INFO and return a valid interval end."""
    info: Dict[str, str] = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value.split(",", 1)[0]
    if "END" in info:
        end = int(info["END"])
    elif "SVLEN" in info:
        end = position + abs(int(info["SVLEN"]))
    else:
        raise ValueError("Inversion record has neither END nor SVLEN")
    if end <= position:
        raise ValueError(f"Invalid inversion interval: start={position}, end={end}")
    return end


def load_variants(
    input_dir: Path, allowed_samples: Iterable[str]
) -> Tuple[List[dict], List[str]]:
    """Load inversion records from single-sample *_filtered.vcf[.gz] files."""
    allowed = set(allowed_samples)
    vcf_files = sorted(input_dir.glob("*_filtered.vcf")) + sorted(input_dir.glob("*_filtered.vcf.gz"))
    if not vcf_files:
        raise FileNotFoundError(f"No *_filtered.vcf or *_filtered.vcf.gz files found in {input_dir}")

    all_variants: List[dict] = []
    found_samples: List[str] = []
    seen_samples = set()
    for vcf_path in vcf_files:
        sample_name = sample_name_from_path(vcf_path)
        if sample_name not in allowed:
            continue
        if sample_name in seen_samples:
            raise ValueError(f"More than one input VCF was found for sample {sample_name}")
        seen_samples.add(sample_name)
        found_samples.append(sample_name)

        with open_text(vcf_path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                columns = line.rstrip("\n").split("\t")
                if len(columns) < 10:
                    raise ValueError(f"Expected at least 10 VCF columns in {vcf_path}:{line_number}")
                try:
                    position = int(columns[1])
                    end = parse_end(position, columns[7])
                except ValueError as exc:
                    raise ValueError(f"Invalid record in {vcf_path}:{line_number}: {exc}") from exc
                genotype = columns[9].split(":", 1)[0]
                all_variants.append(
                    {
                        "chrom": columns[0],
                        "start": position,
                        "end": end,
                        "sample": sample_name,
                        "gt": genotype,
                    }
                )
    if not found_samples:
        raise ValueError("No VCF filename matched a sample in the case or control manifests")
    return all_variants, found_samples


def merge_variants(variants: List[dict], max_distance: int, min_overlap: float) -> List[dict]:
    """Cluster variants using breakpoint distance and reciprocal overlap."""
    variants.sort(key=lambda value: (value["chrom"], value["start"], value["end"]))
    clusters: List[dict] = []
    for variant in variants:
        matched = False
        for cluster in reversed(clusters):
            if variant["chrom"] != cluster["chrom"]:
                break
            if variant["start"] - cluster["start"] > max_distance:
                break
            if (
                abs(variant["start"] - cluster["start"]) <= max_distance
                and abs(variant["end"] - cluster["end"]) <= max_distance
                and calculate_reciprocal_overlap(
                    variant["start"], variant["end"], cluster["start"], cluster["end"]
                )
                >= min_overlap
            ):
                cluster["samples"][variant["sample"]] = variant["gt"]
                matched = True
                break
        if not matched:
            clusters.append(
                {
                    "chrom": variant["chrom"],
                    "start": variant["start"],
                    "end": variant["end"],
                    "samples": {variant["sample"]: variant["gt"]},
                }
            )
    return clusters


def build_results(clusters: List[dict], cases: List[str], controls: List[str]) -> pd.DataFrame:
    """Build carrier statistics and one-sided association-test results."""
    rows = []
    n_cases = len(cases)
    n_controls = len(controls)
    for cluster in clusters:
        case_gt = [parse_gt(cluster["samples"].get(sample, "0/0")) for sample in cases]
        control_gt = [parse_gt(cluster["samples"].get(sample, "0/0")) for sample in controls]
        case_carriers = sum(genotype > 0 for genotype in case_gt)
        control_carriers = sum(genotype > 0 for genotype in control_gt)
        case_ac = sum(case_gt)
        control_ac = sum(control_gt)
        p_value = fisher_exact_onesided_greater(
            case_carriers,
            n_cases - case_carriers,
            control_carriers,
            n_controls - control_carriers,
        )
        odds_ratio = (
            (case_carriers + 0.5) * (n_controls - control_carriers + 0.5)
        ) / ((n_cases - case_carriers + 0.5) * (control_carriers + 0.5))
        row = {
            "CHROM": cluster["chrom"],
            "START": cluster["start"],
            "END": cluster["end"],
            "SVLEN": cluster["end"] - cluster["start"],
            "Total_Carriers": case_carriers + control_carriers,
            "CASE_AC": case_ac,
            "CASE_AF": round(case_ac / (n_cases * 2), 6),
            "CTRL_AC": control_ac,
            "CTRL_AF": round(control_ac / (n_controls * 2), 6),
            "OR": round(odds_ratio, 4),
            "P_value": p_value,
        }
        for sample in cases + controls:
            row[sample] = cluster["samples"].get(sample, "0/0")
        rows.append(row)

    base_columns = [
        "CHROM",
        "START",
        "END",
        "SVLEN",
        "Total_Carriers",
        "CASE_AC",
        "CASE_AF",
        "CTRL_AC",
        "CTRL_AF",
        "OR",
        "P_value",
    ]
    frame = pd.DataFrame(rows, columns=base_columns + cases + controls)
    frame.insert(base_columns.index("P_value") + 1, "FDR_P", calculate_fdr(frame["P_value"].tolist()))
    return frame


def no_singleton_path(output_path: Path) -> Path:
    """Return a sibling path for the carrier-count-greater-than-one result."""
    return output_path.with_name(f"{output_path.stem}_no_singletons{output_path.suffix}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge single-sample inversion VCFs and compare carrier rates."
    )
    parser.add_argument("-i", "--input-dir", required=True, type=Path, help="Directory of *_filtered.vcf[.gz] files.")
    parser.add_argument("-case", "--case-file", required=True, type=Path, help="One case sample ID per line.")
    parser.add_argument("-ctrl", "--control-file", required=True, type=Path, help="One control sample ID per line.")
    parser.add_argument("-d", "--dist", type=int, default=1000, help="Maximum breakpoint distance (default: 1000).")
    parser.add_argument("-r", "--overlap", type=float, default=0.5, help="Minimum reciprocal overlap (default: 0.5).")
    parser.add_argument("-o", "--output", type=Path, default=Path("INV_Stats_Full.txt"), help="Full result TSV path.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.dist < 0:
        raise ValueError("--dist must be non-negative")
    if not 0.0 <= args.overlap <= 1.0:
        raise ValueError("--overlap must be between 0 and 1")
    if not args.input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")

    cases = load_sample_list(args.case_file)
    controls = load_sample_list(args.control_file)
    overlap = sorted(set(cases).intersection(controls))
    if overlap:
        raise ValueError(f"Samples occur in both manifests: {', '.join(overlap[:5])}")

    variants, found_samples = load_variants(args.input_dir, cases + controls)
    actual_cases = [sample for sample in cases if sample in found_samples]
    actual_controls = [sample for sample in controls if sample in found_samples]
    if not actual_cases or not actual_controls:
        raise ValueError("At least one case and one control VCF are required")

    clusters = merge_variants(variants, args.dist, args.overlap)
    results = build_results(clusters, actual_cases, actual_controls)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.sort_values(["P_value", "CHROM", "START"]).to_csv(args.output, sep="\t", index=False)

    filtered = results[results["Total_Carriers"] > 1]
    filtered_path = no_singleton_path(args.output)
    filtered.sort_values(["P_value", "CHROM", "START"]).to_csv(filtered_path, sep="\t", index=False)
    print(f"Merged inversion sites: {len(results)}")
    print(f"Non-singleton sites: {len(filtered)}")
    print(f"Full results: {args.output}")
    print(f"Non-singleton results: {filtered_path}")


if __name__ == "__main__":
    main()
