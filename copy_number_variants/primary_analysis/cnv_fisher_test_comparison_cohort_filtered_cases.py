# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


SCRIPT_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class ComparisonConfig:
    output_prefix: str
    scz_dir: Path
    comparison_cohort_dir: Path
    public_reference_dir: Path | None = None


def find_project_dir(anchor_dir: Path) -> Path:
    for parent in (anchor_dir, *anchor_dir.parents):
        cnv_dir = parent / "CNV"
        if not cnv_dir.is_dir():
            continue
        if any(path.is_dir() for path in cnv_dir.glob("*CNV")):
            return parent
    raise FileNotFoundError("Cannot locate project directory containing copy_number_variants/*CNV")


def find_cnv_base_dir(anchor_dir: Path) -> Path:
    project_dir = find_project_dir(anchor_dir)
    cnv_roots = sorted(path for path in (project_dir / "CNV").glob("*CNV") if path.is_dir())
    for cnv_root in cnv_roots:
        matches = sorted(path for path in cnv_root.glob("*protein coding genes*") if path.is_dir())
        if matches:
            return matches[0]
    raise FileNotFoundError("Cannot locate CNV filtered directory matching *protein coding genes*")


def build_comparisons(cnv_base_dir: Path) -> tuple[ComparisonConfig, ...]:
    return (
        ComparisonConfig(
            output_prefix="true_case_cohort_comparison_cohort_haplotype",
            scz_dir=cnv_base_dir / "slurm_scripts_case_cohort",
            comparison_cohort_dir=cnv_base_dir / "slurm_scripts_comparison_site",
        ),
        ComparisonConfig(
            output_prefix="true_case_cohort_comparison_cohort_public_reference_haplotype",
            scz_dir=cnv_base_dir / "slurm_scripts_case_cohort",
            comparison_cohort_dir=cnv_base_dir / "slurm_scripts_comparison_site",
            public_reference_dir=cnv_base_dir / "slurm_scripts_public_reference",
        ),
    )


def load_frequency_table(dataset_dir: Path, count_column: str) -> pd.DataFrame:
    table = pd.read_table(
        dataset_dir / "gene_frequencies_filtered.txt",
        header=None,
        names=["Gene", count_column, f"SourceFrequency_{count_column}"],
    )
    table[count_column] = pd.to_numeric(table[count_column], errors="coerce").fillna(0).astype(int)
    return table[["Gene", count_column]]


def infer_total_from_frequency_file(dataset_dir: Path) -> int:
    table = pd.read_table(
        dataset_dir / "gene_frequencies_filtered.txt",
        header=None,
        names=["Gene", "Count", "Frequency"],
    )
    table["Count"] = pd.to_numeric(table["Count"], errors="coerce")
    table["Frequency"] = pd.to_numeric(table["Frequency"], errors="coerce")
    informative = table[(table["Count"] > 0) & (table["Frequency"] > 0)].copy()
    if informative.empty:
        raise ValueError(f"Cannot infer total from {dataset_dir}")

    approximations = informative["Count"] / informative["Frequency"]
    center = int(round(float(approximations.median())))
    candidates = {center}
    for approx_value in approximations.head(100):
        rounded = int(round(float(approx_value)))
        for candidate in range(max(1, rounded - 3), rounded + 4):
            candidates.add(candidate)

    best_total = None
    best_score = None
    for candidate in sorted(candidates):
        score = (informative["Count"] / candidate - informative["Frequency"]).abs().sum()
        if best_score is None or score < best_score:
            best_total = candidate
            best_score = score

    if best_total is None:
        raise ValueError(f"Failed to infer total from {dataset_dir}")
    return int(best_total)


@lru_cache(maxsize=None)
def infer_total_haplotypes(dataset_dir: Path) -> int:
    haplotype_file = dataset_dir / "haplotype_CN.xlsx"
    if haplotype_file.exists():
        return int(pd.read_excel(haplotype_file, index_col=0).shape[0])
    return infer_total_from_frequency_file(dataset_dir)


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    count = len(values)
    if count == 0:
        return pd.Series(dtype=float, index=p_values.index)

    order = np.argsort(values)
    ranked = values[order]
    ranks = np.arange(1, count + 1, dtype=float)
    adjusted_ranked = np.minimum.accumulate((ranked * count / ranks)[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)

    adjusted = np.empty(count, dtype=float)
    adjusted[order] = adjusted_ranked
    return pd.Series(adjusted, index=p_values.index)


def bonferroni_adjust(p_values: pd.Series) -> pd.Series:
    return pd.Series(np.minimum(p_values.to_numpy(dtype=float) * len(p_values), 1.0), index=p_values.index)


def validate_adjustments(table: pd.DataFrame) -> tuple[float | None, float | None]:
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        return None, None

    p_values = table["P_value"].to_numpy(dtype=float)
    fdr_reference = multipletests(p_values, method="fdr_bh")[1]
    bonferroni_reference = multipletests(p_values, method="bonferroni")[1]
    fdr_max_abs_diff = float(np.max(np.abs(table["FDR"].to_numpy(dtype=float) - fdr_reference)))
    bonferroni_max_abs_diff = float(
        np.max(np.abs(table["Bonferroni"].to_numpy(dtype=float) - bonferroni_reference))
    )
    if not np.allclose(table["FDR"].to_numpy(dtype=float), fdr_reference):
        raise AssertionError(f"FDR differs from statsmodels fdr_bh; max abs diff={fdr_max_abs_diff}")
    if not np.allclose(table["Bonferroni"].to_numpy(dtype=float), bonferroni_reference):
        raise AssertionError(f"Bonferroni differs from statsmodels; max abs diff={bonferroni_max_abs_diff}")
    return fdr_max_abs_diff, bonferroni_max_abs_diff


def fisher_p_value(row: pd.Series, control_count_column: str, control_other_column: str) -> float:
    contingency_table = np.array(
        [
            [int(row["Count_Scz"]), int(row["SCZ_other"])],
            [int(row[control_count_column]), int(row[control_other_column])],
        ]
    )
    _, p_value = fisher_exact(contingency_table, alternative="greater")
    return float(p_value)


def fisher_odds_ratio(row: pd.Series, control_count_column: str, control_other_column: str) -> float:
    contingency_table = np.array(
        [
            [int(row["Count_Scz"]), int(row["SCZ_other"])],
            [int(row[control_count_column]), int(row[control_other_column])],
        ]
    )
    odds_ratio, _ = fisher_exact(contingency_table, alternative="greater")
    return float(odds_ratio)


def build_result_table(config: ComparisonConfig) -> tuple[pd.DataFrame, dict[str, int]]:
    scz_total = infer_total_haplotypes(config.scz_dir)
    comparison_cohort_total = infer_total_haplotypes(config.comparison_cohort_dir)

    scz_data = load_frequency_table(config.scz_dir, "Count_Scz")
    comparison_cohort_data = load_frequency_table(config.comparison_cohort_dir, "Count_comparison_cohort")
    combined = scz_data.merge(comparison_cohort_data, on="Gene", how="outer")

    control_count_column = "Count_comparison_cohort"
    control_frequency_column = "Frequency_comparison_cohort"
    control_other_column = "comparison_cohort_other"
    control_total = comparison_cohort_total
    totals = {"scz_total": scz_total, "comparison_cohort_total": comparison_cohort_total, "control_total": control_total}

    if config.public_reference_dir is None:
        combined[["Count_Scz", "Count_comparison_cohort"]] = combined[["Count_Scz", "Count_comparison_cohort"]].fillna(0).astype(int)
        combined["Frequency_Scz"] = combined["Count_Scz"] / scz_total
        combined["Frequency_comparison_cohort"] = combined["Count_comparison_cohort"] / comparison_cohort_total
        combined["SCZ_other"] = scz_total - combined["Count_Scz"]
        combined["comparison_cohort_other"] = comparison_cohort_total - combined["Count_comparison_cohort"]
        output_columns = [
            "Gene",
            "Count_Scz",
            "Frequency_Scz",
            "Count_comparison_cohort",
            "Frequency_comparison_cohort",
            "SCZ_other",
            "comparison_cohort_other",
        ]
    else:
        public_reference_total = infer_total_haplotypes(config.public_reference_dir)
        totals["public_reference_total"] = public_reference_total
        control_total = comparison_cohort_total + public_reference_total
        totals["control_total"] = control_total

        public_reference_data = load_frequency_table(config.public_reference_dir, "Count_public_reference")
        combined = combined.merge(public_reference_data, on="Gene", how="outer")
        combined[["Count_Scz", "Count_comparison_cohort", "Count_public_reference"]] = combined[
            ["Count_Scz", "Count_comparison_cohort", "Count_public_reference"]
        ].fillna(0).astype(int)

        control_count_column = "Count_comparison_cohort+public_reference"
        control_frequency_column = "Frequency_comparison_cohort+public_reference"
        control_other_column = "comparison_cohort+public_reference_other"

        combined[control_count_column] = combined["Count_comparison_cohort"] + combined["Count_public_reference"]
        combined["Frequency_Scz"] = combined["Count_Scz"] / scz_total
        combined[control_frequency_column] = combined[control_count_column] / control_total
        combined["Frequency_comparison_cohort"] = combined["Count_comparison_cohort"] / comparison_cohort_total
        combined["Frequency_public_reference"] = combined["Count_public_reference"] / public_reference_total
        combined["SCZ_other"] = scz_total - combined["Count_Scz"]
        combined[control_other_column] = control_total - combined[control_count_column]

        output_columns = [
            "Gene",
            "Count_Scz",
            "Frequency_Scz",
            "Count_comparison_cohort+public_reference",
            "Frequency_comparison_cohort+public_reference",
            "Count_comparison_cohort",
            "Frequency_comparison_cohort",
            "Count_public_reference",
            "Frequency_public_reference",
            "SCZ_other",
            "comparison_cohort+public_reference_other",
        ]

    if (combined["Count_Scz"] > scz_total).any():
        raise ValueError("Some SCZ counts exceed the inferred SCZ haplotype total")
    if (combined[control_count_column] > control_total).any():
        raise ValueError("Some control counts exceed the inferred control haplotype total")

    combined["P_value"] = combined.apply(
        fisher_p_value,
        axis=1,
        control_count_column=control_count_column,
        control_other_column=control_other_column,
    )
    combined["FDR"] = benjamini_hochberg(combined["P_value"])
    combined["Bonferroni"] = bonferroni_adjust(combined["P_value"])
    combined["Odds Ratio"] = combined.apply(
        fisher_odds_ratio,
        axis=1,
        control_count_column=control_count_column,
        control_other_column=control_other_column,
    )

    combined = combined.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)
    merge_table = combined[output_columns + ["P_value", "FDR", "Bonferroni", "Odds Ratio"]]
    fdr_diff, bonferroni_diff = validate_adjustments(merge_table)
    totals["tested_genes"] = len(merge_table)
    totals["raw_p_lt_0.05"] = int((merge_table["P_value"] < 0.05).sum())
    totals["fdr_lt_0.05"] = int((merge_table["FDR"] < 0.05).sum())
    totals["bonferroni_lt_0.05"] = int((merge_table["Bonferroni"] < 0.05).sum())
    totals["fdr_reference_max_abs_diff"] = -1 if fdr_diff is None else fdr_diff
    totals["bonferroni_reference_max_abs_diff"] = -1 if bonferroni_diff is None else bonferroni_diff
    return merge_table, totals


def write_outputs(config: ComparisonConfig, output_dir: Path) -> None:
    merge_table, totals = build_result_table(config)
    control_frequency_column = "Frequency_comparison_cohort+public_reference" if config.public_reference_dir is not None else "Frequency_comparison_cohort"
    filtered_columns = merge_table.columns.drop("Odds Ratio")

    significant_all = merge_table.loc[merge_table["P_value"] < 0.05, filtered_columns]
    significant_scz_higher = significant_all[
        significant_all["Frequency_Scz"] > significant_all[control_frequency_column]
    ]

    merge_path = output_dir / f"{config.output_prefix}_merge.xlsx"
    all_path = output_dir / f"{config.output_prefix}_fisher_significant_all.xlsx"
    scz_higher_path = output_dir / f"{config.output_prefix}_fisher_significant_scz_higher.xlsx"

    merge_table.to_excel(merge_path, index=False)
    significant_all.to_excel(all_path, index=False)
    significant_scz_higher.to_excel(scz_higher_path, index=False)

    print(
        f"{config.output_prefix}: "
        f"SCZ={totals['scz_total']}, comparison_cohort={totals['comparison_cohort_total']}"
        + (f", public_reference={totals['public_reference_total']}" if "public_reference_total" in totals else "")
        + f", control={totals['control_total']}, tested_genes={totals['tested_genes']}, "
        f"raw_p<0.05={totals['raw_p_lt_0.05']}, "
        f"FDR<0.05={totals['fdr_lt_0.05']}, "
        f"Bonferroni<0.05={totals['bonferroni_lt_0.05']}, "
        f"FDR_check_max_diff={totals['fdr_reference_max_abs_diff']}, "
        f"Bonf_check_max_diff={totals['bonferroni_reference_max_abs_diff']}"
    )
    print(f"  saved: {merge_path}")
    print(f"  saved: {all_path}")
    print(f"  saved: {scz_higher_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run one-sided Fisher tests for CNV-related gene frequencies in "
            "SCZ versus comparison_cohort and/or comparison_cohort+public_reference."
        )
    )
    parser.add_argument(
        "--cnv-base-dir",
        type=Path,
        help=(
            "Directory containing slurm_scripts_case_cohort, "
            "slurm_scripts_comparison_site, and optionally slurm_scripts_public_reference. If "
            "omitted, the historical project layout is searched."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "outputs",
        help="Directory for result workbooks.",
    )
    parser.add_argument(
        "--comparison",
        choices=("scz_vs_comparison_cohort", "scz_vs_comparison_cohort_public_reference", "both"),
        default="both",
        help="Comparison to run.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cnv_base_dir = args.cnv_base_dir or find_cnv_base_dir(SCRIPT_DIR)
    comparisons = build_comparisons(cnv_base_dir)
    if args.comparison == "scz_vs_comparison_cohort":
        comparisons = comparisons[:1]
    elif args.comparison == "scz_vs_comparison_cohort_public_reference":
        comparisons = comparisons[1:]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"CNV_BASE_DIR={cnv_base_dir}")
    for config in comparisons:
        write_outputs(config, args.output_dir)


if __name__ == "__main__":
    main()
