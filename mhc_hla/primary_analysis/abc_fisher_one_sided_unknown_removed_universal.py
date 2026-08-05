# -*- coding: utf-8 -*-
"""Run Fig. 4F HLA-A/B/C one-sided frequency tests.

The script compares comparison_site/control and SCZ HLA haplotype frequencies.  For each
gene it first tries to read manually QC-filtered public-release sheets
``A_qc_passed``, ``B_qc_passed``, and ``C_qc_passed``; if those are absent, it
tries the historical Chinese worksheet names and then falls back to a raw
summary sheet containing gene columns.

Unknown, NULL, and blank values are excluded before testing.  The default
Fisher direction is ``alternative="less"`` for the comparison_site-vs-SCZ table, matching
the historical Fig. 4F script and testing whether the comparison_site proportion is lower
than the SCZ proportion.

Example:
python abc_fisher_one_sided_unknown_removed_universal.py \
    immuannot_comparison.xlsx immuannot_scz.xlsx
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


PUBLIC_GENE_SHEETS = {
    "A": "A_qc_passed",
    "B": "B_qc_passed",
    "C": "C_qc_passed",
}

LEGACY_GENE_SHEETS = {
    "A": "A_validated",
    "B": "B_validated",
    "C": "C_validated",
}

DEFAULT_GENE_SHEETS = PUBLIC_GENE_SHEETS
DEFAULT_RAW_SHEET_CANDIDATES = (
    "ABC",
    "Sheet1",
    "combined",
)

INVALID_TOKENS = {
    "",
    "NULL",
    "UNKNOWN",
    "(blank)",
}


@dataclass
class GeneData:
    gene: str
    dataset_label: str
    source_sheet: str
    source_column: str
    values: pd.Series

    @property
    def total(self) -> int:
        return int(len(self.values))

    @property
    def counts(self) -> pd.Series:
        return self.values.value_counts()


@dataclass
class ExtraControlSpec:
    path: Path
    raw_sheet: str | None
    label: str


class WorkbookReader:
    def __init__(self, path: Path):
        self.path = path
        self.excel = pd.ExcelFile(path)
        self.sheet_names = set(self.excel.sheet_names)
        self._sheet_cache: dict[str, pd.DataFrame] = {}

    def has_sheet(self, sheet_name: str | None) -> bool:
        return bool(sheet_name) and sheet_name in self.sheet_names

    def load_sheet(self, sheet_name: str) -> pd.DataFrame:
        if sheet_name not in self._sheet_cache:
            dataframe = pd.read_excel(self.path, sheet_name=sheet_name, dtype=object)
            dataframe = dataframe.dropna(how="all")
            self._sheet_cache[sheet_name] = dataframe
        return self._sheet_cache[sheet_name]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run HLA-A/B/C frequency tests after removing unknown/NULL values."
    )
    parser.add_argument("comparison_xlsx", help="Control/comparison_site Immuannot workbook path.")
    parser.add_argument("scz_xlsx", help="SCZ Immuannot workbook path.")
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output directory. By default, a directory is generated next to this script.",
    )
    parser.add_argument(
        "--result-tag",
        help="Result filename tag, for example: new_case_cohort_ragtag_actual_totals_26-4-26.",
    )
    parser.add_argument(
        "--control-extra-xlsx",
        action="append",
        default=[],
        help="Additional workbook to merge into the control group. May be provided multiple times.",
    )
    parser.add_argument(
        "--control-extra-raw-sheet",
        action="append",
        default=[],
        help="Raw worksheet for each --control-extra-xlsx entry. May be provided multiple times.",
    )
    parser.add_argument(
        "--control-extra-label",
        action="append",
        default=[],
        help="Label for each --control-extra-xlsx entry. Defaults to the workbook filename.",
    )
    parser.add_argument(
        "--method",
        choices=("fisher", "chi2"),
        default="fisher",
        help="Statistical method. Default: fisher.",
    )
    parser.add_argument(
        "--alternative",
        choices=("less", "greater", "two-sided"),
        default="less",
        help="Fisher exact-test alternative. Default: less.",
    )
    parser.add_argument(
        "--comparison-raw-sheet",
        default="Sheet1",
        help="Fallback raw worksheet for the control/comparison_site workbook. Default: Sheet1.",
    )
    parser.add_argument(
        "--scz-raw-sheet",
        default="ABC",
        help="Fallback raw worksheet for the SCZ workbook. Default: ABC.",
    )
    parser.add_argument(
        "--comparison-sheet-a",
        default=DEFAULT_GENE_SHEETS["A"],
        help=f"Preferred HLA-A worksheet for control/comparison_site. Default: {DEFAULT_GENE_SHEETS['A']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--comparison-sheet-b",
        default=DEFAULT_GENE_SHEETS["B"],
        help=f"Preferred HLA-B worksheet for control/comparison_site. Default: {DEFAULT_GENE_SHEETS['B']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--comparison-sheet-c",
        default=DEFAULT_GENE_SHEETS["C"],
        help=f"Preferred HLA-C worksheet for control/comparison_site. Default: {DEFAULT_GENE_SHEETS['C']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--scz-sheet-a",
        default=DEFAULT_GENE_SHEETS["A"],
        help=f"Preferred HLA-A worksheet for SCZ. Default: {DEFAULT_GENE_SHEETS['A']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--scz-sheet-b",
        default=DEFAULT_GENE_SHEETS["B"],
        help=f"Preferred HLA-B worksheet for SCZ. Default: {DEFAULT_GENE_SHEETS['B']}; legacy fallback documented separately.",
    )
    parser.add_argument(
        "--scz-sheet-c",
        default=DEFAULT_GENE_SHEETS["C"],
        help=f"Preferred HLA-C worksheet for SCZ. Default: {DEFAULT_GENE_SHEETS['C']}; legacy fallback documented separately.",
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


def is_valid_haplotype(value: str) -> bool:
    if not value:
        return False
    return value.strip().upper() not in INVALID_TOKENS


def count_nonempty_values(series: pd.Series) -> int:
    return sum(is_valid_haplotype(normalize_value(value)) for value in series)


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


def sanitize_label(text: str) -> str:
    sanitized = re.sub(r"[^\w\-]+", "_", text.strip())
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    return sanitized or "dataset"


def simplify_immuannot_stem(path: Path) -> str:
    stem = path.stem
    stem = re.sub(r"^Immuannot_+", "", stem)
    stem = re.sub(r"results$", "", stem)
    stem = re.sub(r"_?results?$", "", stem, flags=re.IGNORECASE)
    return sanitize_label(stem)


def choose_raw_sheet(reader: WorkbookReader, preferred_sheet: str | None) -> str:
    candidates = []
    if preferred_sheet:
        candidates.append(preferred_sheet)
    for fallback_sheet in DEFAULT_RAW_SHEET_CANDIDATES:
        if fallback_sheet not in candidates:
            candidates.append(fallback_sheet)
    for sheet_name in candidates:
        if reader.has_sheet(sheet_name):
            return sheet_name
    return reader.excel.sheet_names[0]


def build_extra_control_specs(args: argparse.Namespace) -> list[ExtraControlSpec]:
    extra_paths = args.control_extra_xlsx or []
    extra_sheets = args.control_extra_raw_sheet or []
    extra_labels = args.control_extra_label or []

    if extra_sheets and len(extra_sheets) != len(extra_paths):
        raise ValueError("--control-extra-raw-sheet must be provided once per --control-extra-xlsx.")
    if extra_labels and len(extra_labels) != len(extra_paths):
        raise ValueError("--control-extra-label must be provided once per --control-extra-xlsx.")

    specs: list[ExtraControlSpec] = []
    for index, raw_path in enumerate(extra_paths):
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Additional control workbook does not exist: {path}")
        raw_sheet = extra_sheets[index] if index < len(extra_sheets) else None
        label = extra_labels[index] if index < len(extra_labels) else simplify_immuannot_stem(path)
        specs.append(
            ExtraControlSpec(
                path=path,
                raw_sheet=raw_sheet,
                label=label,
            )
        )
    return specs


def load_gene_data(
    reader: WorkbookReader,
    dataset_label: str,
    gene_name: str,
    preferred_gene_sheet: str | None,
    preferred_raw_sheet: str | None,
) -> GeneData:
    if reader.has_sheet(preferred_gene_sheet):
        source_sheet = preferred_gene_sheet
    elif (
        preferred_gene_sheet == PUBLIC_GENE_SHEETS[gene_name]
        and reader.has_sheet(LEGACY_GENE_SHEETS[gene_name])
    ):
        source_sheet = LEGACY_GENE_SHEETS[gene_name]
    else:
        source_sheet = choose_raw_sheet(reader, preferred_raw_sheet)

    dataframe = reader.load_sheet(source_sheet)
    gene_column = pick_gene_column(dataframe, gene_name)
    values = dataframe[gene_column].map(normalize_value)
    values = values[values.map(is_valid_haplotype)].astype(str)

    return GeneData(
        gene=gene_name,
        dataset_label=dataset_label,
        source_sheet=source_sheet,
        source_column=gene_column,
        values=values.reset_index(drop=True),
    )


def combine_gene_data(gene_name: str, dataset_label: str, gene_datasets: list[GeneData]) -> GeneData:
    if not gene_datasets:
        raise ValueError(f"No mergeable {gene_name} data for {dataset_label}.")

    combined_values = pd.concat([gene_data.values for gene_data in gene_datasets], ignore_index=True)
    combined_source_sheet = "; ".join(
        f"{gene_data.dataset_label}:{gene_data.source_sheet}" for gene_data in gene_datasets
    )
    combined_source_column = "; ".join(
        f"{gene_data.dataset_label}:{gene_data.source_column}" for gene_data in gene_datasets
    )

    return GeneData(
        gene=gene_name,
        dataset_label=dataset_label,
        source_sheet=combined_source_sheet,
        source_column=combined_source_column,
        values=combined_values,
    )


def perform_statistical_test(
    comparison_count: int,
    scz_count: int,
    comparison_total: int,
    scz_total: int,
    method: str,
    alternative: str,
) -> tuple[float, np.ndarray]:
    contingency_table = np.array(
        [
            [comparison_count, scz_count],
            [comparison_total - comparison_count, scz_total - scz_count],
        ]
    )

    if np.any(contingency_table < 0):
        raise ValueError(f"Negative value in contingency table; check inputs: {contingency_table.tolist()}")

    if method == "fisher":
        _, p_value = fisher_exact(contingency_table, alternative=alternative)
    elif method == "chi2":
        _, p_value, _, _ = chi2_contingency(contingency_table, correction=True)
    else:
        raise ValueError(f"Unsupported statistical method: {method}")

    return float(p_value), contingency_table


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    if p_values.empty:
        return pd.Series(dtype=float)

    numeric_p = pd.to_numeric(p_values, errors="coerce")
    result = pd.Series(np.nan, index=numeric_p.index, dtype=float)

    valid_mask = numeric_p.notna()
    valid_p = numeric_p[valid_mask]
    if valid_p.empty:
        return result

    order = np.argsort(valid_p.to_numpy())
    sorted_p = valid_p.to_numpy()[order]
    total = len(sorted_p)

    adjusted = np.empty(total, dtype=float)
    running_min = 1.0
    for reverse_index in range(total - 1, -1, -1):
        rank = reverse_index + 1
        candidate = sorted_p[reverse_index] * total / rank
        running_min = min(running_min, candidate)
        adjusted[reverse_index] = min(running_min, 1.0)

    valid_indices = valid_p.index.to_numpy()[order]
    result.loc[valid_indices] = adjusted
    return result


def build_gene_results(
    gene_name: str,
    comparison_data: GeneData,
    scz_data: GeneData,
    method: str,
    alternative: str,
) -> pd.DataFrame:
    all_genotypes = sorted(set(comparison_data.counts.index).union(set(scz_data.counts.index)))
    results = []

    for genotype in all_genotypes:
        comparison_count = int(comparison_data.counts.get(genotype, 0))
        scz_count = int(scz_data.counts.get(genotype, 0))
        p_value, contingency_table = perform_statistical_test(
            comparison_count=comparison_count,
            scz_count=scz_count,
            comparison_total=comparison_data.total,
            scz_total=scz_data.total,
            method=method,
            alternative=alternative,
        )

        comparison_prop = comparison_count / comparison_data.total if comparison_data.total else 0.0
        scz_prop = scz_count / scz_data.total if scz_data.total else 0.0

        results.append(
            {
                "Gene": gene_name,
                "Genotype": genotype,
                "HC_Count": comparison_count,
                "HC_Total": comparison_data.total,
                "HC_Proportion": comparison_prop,
                "SCZ_Count": scz_count,
                "SCZ_Total": scz_data.total,
                "SCZ_Proportion": scz_prop,
                "SCZ_minus_HC_Proportion": scz_prop - comparison_prop,
                "P-Value": p_value,
                "Method": method,
                "Alternative": alternative if method == "fisher" else "",
                "Contingency_Table": contingency_table.tolist(),
            }
        )

    result_df = pd.DataFrame(results)
    if result_df.empty:
        return result_df

    result_df["FDR_BH"] = benjamini_hochberg(result_df["P-Value"])
    result_df["Significant_FDR_0.05"] = result_df["FDR_BH"] <= 0.05

    return result_df.sort_values(
        by=["FDR_BH", "P-Value", "SCZ_minus_HC_Proportion", "SCZ_Count", "Genotype"],
        ascending=[True, True, False, False, True],
    ).reset_index(drop=True)


def default_output_dir(comparison_path: Path, scz_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    comparison_label = sanitize_label(comparison_path.stem)
    scz_label = sanitize_label(scz_path.stem)
    return script_dir / f"universal_results_{scz_label}_vs_{comparison_label}"


def infer_result_tag(
    explicit_tag: str | None,
    output_dir: Path,
    scz_path: Path,
) -> str:
    if explicit_tag:
        return explicit_tag

    output_dir_name = output_dir.name
    match = re.match(r"^(?P<label>.+)_(?P<date>\d{2}-\d{1,2}-\d{1,2})$", output_dir_name)
    if match:
        label = sanitize_label(match.group("label"))
        date_text = match.group("date")
        return f"new_{label}_actual_totals_{date_text}"

    scz_label = simplify_immuannot_stem(scz_path)
    return f"new_{scz_label}_actual_totals"


def write_outputs(
    output_dir: Path,
    method: str,
    result_tag: str,
    gene_results: dict[str, pd.DataFrame],
    summary_rows: list[dict[str, object]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / f"gene_totals_summary_{result_tag}.xlsx"
    summary_df.to_excel(summary_path, index=False)

    for gene_name, result_df in gene_results.items():
        output_path = output_dir / f"{gene_name}_{method}_results_{result_tag}.xlsx"
        result_df.to_excel(output_path, index=False)


def main() -> int:
    args = parse_args()

    comparison_path = Path(args.comparison_xlsx).expanduser().resolve()
    scz_path = Path(args.scz_xlsx).expanduser().resolve()

    if not comparison_path.exists():
        raise FileNotFoundError(f"Control/comparison_site input workbook does not exist: {comparison_path}")
    if not scz_path.exists():
        raise FileNotFoundError(f"SCZ input workbook does not exist: {scz_path}")

    extra_control_specs = build_extra_control_specs(args)

    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else default_output_dir(comparison_path, scz_path)
    )
    result_tag = infer_result_tag(args.result_tag, output_dir, scz_path)

    comparison_reader = WorkbookReader(comparison_path)
    scz_reader = WorkbookReader(scz_path)
    extra_control_readers = [
        (spec, WorkbookReader(spec.path))
        for spec in extra_control_specs
    ]

    comparison_sheet_map = {
        "A": args.comparison_sheet_a,
        "B": args.comparison_sheet_b,
        "C": args.comparison_sheet_c,
    }
    scz_sheet_map = {
        "A": args.scz_sheet_a,
        "B": args.scz_sheet_b,
        "C": args.scz_sheet_c,
    }

    summary_rows: list[dict[str, object]] = []
    gene_results: dict[str, pd.DataFrame] = {}

    for gene_name in ("A", "B", "C"):
        primary_control_data = load_gene_data(
            reader=comparison_reader,
            dataset_label="comparison_site",
            gene_name=gene_name,
            preferred_gene_sheet=comparison_sheet_map[gene_name],
            preferred_raw_sheet=args.comparison_raw_sheet,
        )
        control_gene_datasets = [primary_control_data]
        for spec, reader in extra_control_readers:
            extra_control_data = load_gene_data(
                reader=reader,
                dataset_label=spec.label,
                gene_name=gene_name,
                preferred_gene_sheet=DEFAULT_GENE_SHEETS.get(gene_name),
                preferred_raw_sheet=spec.raw_sheet,
            )
            control_gene_datasets.append(extra_control_data)

        comparison_data = combine_gene_data(
            gene_name=gene_name,
            dataset_label="comparison_site+extras" if extra_control_readers else "comparison_site",
            gene_datasets=control_gene_datasets,
        )
        scz_data = load_gene_data(
            reader=scz_reader,
            dataset_label="SCZ",
            gene_name=gene_name,
            preferred_gene_sheet=scz_sheet_map[gene_name],
            preferred_raw_sheet=args.scz_raw_sheet,
        )

        summary_rows.append(
            {
                "Gene": gene_name,
                "HC_SourceSheet": comparison_data.source_sheet,
                "HC_SourceColumn": comparison_data.source_column,
                "HC_Total": comparison_data.total,
                "SCZ_SourceSheet": scz_data.source_sheet,
                "SCZ_SourceColumn": scz_data.source_column,
                "SCZ_Total": scz_data.total,
                "Control_Components": "; ".join(
                    f"{gene_data.dataset_label}:{gene_data.total}" for gene_data in control_gene_datasets
                ),
                "Method": args.method,
                "Alternative": args.alternative if args.method == "fisher" else "",
            }
        )

        gene_results[gene_name] = build_gene_results(
            gene_name=gene_name,
            comparison_data=comparison_data,
            scz_data=scz_data,
            method=args.method,
            alternative=args.alternative,
        )

    write_outputs(
        output_dir=output_dir,
        method=args.method,
        result_tag=result_tag,
        gene_results=gene_results,
        summary_rows=summary_rows,
    )

    print("Fig.4F summary")
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print(f"Result tag: {result_tag}")
    print(f"\nOutput dir: {output_dir}")
    for gene_name, result_df in gene_results.items():
        print(f"{gene_name}: {len(result_df)} genotypes")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
