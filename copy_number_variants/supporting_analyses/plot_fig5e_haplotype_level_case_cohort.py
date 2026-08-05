# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


DATE_TAG = "26-7-20"
CASE_LABEL = "case_cohort"
comparison_cohort_LABEL = "comparison_cohort"
ALLPUB_LABEL = "public_reference"

COLORS = {
    CASE_LABEL: "#da7271",
    comparison_cohort_LABEL: "#1f78b4",
    ALLPUB_LABEL: "#d9a627",
}

DISPLAY_LABELS = {
    CASE_LABEL: "SCZ",
    comparison_cohort_LABEL: "Control",
    ALLPUB_LABEL: "Public",
}


@dataclass(frozen=True)
class Dataset:
    label: str
    directory_name: str
    sheet_name: str = "Sheet1"


DATASETS = {
    CASE_LABEL: Dataset(CASE_LABEL, "slurm_scripts_case_cohort"),
    comparison_cohort_LABEL: Dataset(comparison_cohort_LABEL, "slurm_scripts_comparison_site"),
    ALLPUB_LABEL: Dataset(ALLPUB_LABEL, "slurm_scripts_public_reference"),
}

COMPARISONS = {
    "case_cohort_vs_comparison_cohort": [CASE_LABEL, comparison_cohort_LABEL],
    "case_cohort_vs_comparison_cohort_vs_public_reference": [CASE_LABEL, comparison_cohort_LABEL, ALLPUB_LABEL],
}

COMPARISON_SHORT_NAMES = {
    "case_cohort_vs_comparison_cohort": "SCZ_Control",
    "case_cohort_vs_comparison_cohort_vs_public_reference": "SCZ_Control_Public",
}

SCRIPT_DIR = Path(__file__).resolve().parent


def find_cnv_base_dir() -> Path:
    for parent in SCRIPT_DIR.parents:
        cnv_root = parent / "CNV"
        if not cnv_root.is_dir():
            continue
        for candidate in cnv_root.glob("*CNV"):
            if not candidate.is_dir():
                continue
            matches = sorted(path for path in candidate.glob("*protein coding genes*") if path.is_dir())
            if matches:
                return matches[0]
    raise FileNotFoundError("Cannot locate filtered CNV directory")


CNV_BASE_DIR = find_cnv_base_dir()


def dataset_dir(label: str) -> Path:
    return CNV_BASE_DIR / DATASETS[label].directory_name


def load_gene_haplotype_counts(label: str) -> pd.DataFrame:
    dataset = DATASETS[label]
    file_path = dataset_dir(label) / "haplotype_CN.xlsx"
    matrix = pd.read_excel(file_path, sheet_name=dataset.sheet_name, index_col=0)
    matrix = matrix.dropna(how="all")
    if "[xxx]" in matrix.columns:
        matrix = matrix.drop(columns=["[xxx]"])
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    counts = (matrix > 0).sum(axis=0).astype(int)
    counts = counts[counts > 0].sort_values(ascending=False)
    return pd.DataFrame(
        {
            "Gene": counts.index.astype(str),
            "Haplotype_Count": counts.to_numpy(),
            "Dataset": label,
        }
    )


def group_haplotype_original(values: pd.Series) -> pd.Series:
    labels = [str(i) for i in range(1, 10)] + ["10-49", "50+"]

    def group(value: int) -> str:
        if 1 <= value <= 9:
            return str(value)
        if 10 <= value <= 49:
            return "10-49"
        return "50+"

    return values.map(group).value_counts().reindex(labels, fill_value=0) / len(values)


def group_exact(values: pd.Series, upper: int) -> pd.Series:
    filtered = values[(values >= 1) & (values <= upper)]
    labels = [str(i) for i in range(1, upper + 1)]
    if filtered.empty:
        return pd.Series(0.0, index=labels)
    return filtered.astype(str).value_counts().reindex(labels, fill_value=0) / len(filtered)


def exact_distribution_full_denominator(values: pd.Series) -> pd.Series:
    if values.empty:
        return pd.Series(dtype=float)
    max_value = int(values.max())
    labels = [str(i) for i in range(1, max_value + 1)]
    counts = values.astype(str).value_counts().reindex(labels, fill_value=0)
    return counts / len(values)


def crop_distribution_for_plot(distribution: pd.Series, upper: int) -> pd.Series:
    labels = [str(i) for i in range(1, upper + 1)]
    return distribution.reindex(labels, fill_value=0.0)


def write_distribution_table(distributions: dict[str, pd.Series], output_path: Path) -> None:
    table = pd.DataFrame(distributions)
    if all(str(index_value).isdigit() for index_value in table.index):
        table = table.loc[sorted(table.index, key=lambda index_value: int(str(index_value)))]
    table.index.name = "Bin"
    table.to_csv(output_path, sep="\t")


def plot_grouped_bars(
    proportions: dict[str, pd.Series],
    output_prefix: Path,
    xlabel: str,
    ylabel: str = "Proportion of CNV-related Genes",
    figsize: tuple[float, float] = (10.5, 6),
    xtick_step: int = 1,
) -> None:
    labels = list(next(iter(proportions.values())).index)
    x = np.arange(len(labels))
    group_labels = list(proportions)
    width = min(0.8 / len(group_labels), 0.35)

    fig, ax = plt.subplots(figsize=figsize)
    offset_start = -width * (len(group_labels) - 1) / 2
    for group_index, group_label in enumerate(group_labels):
        ax.bar(
            x + offset_start + group_index * width,
            proportions[group_label].values,
            width=width,
            label=DISPLAY_LABELS.get(group_label, group_label),
            color=COLORS[group_label],
            edgecolor="black",
            linewidth=0.4,
        )

    tick_positions = x[::xtick_step]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([labels[i] for i in tick_positions], fontsize=10, rotation=0 if xtick_step == 1 else 45)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.5)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def summarize_counts(count_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for label, table in count_tables.items():
        values = table["Haplotype_Count"]
        rows.append(
            {
                "Dataset": label,
                "Input_File": str(dataset_dir(label) / "haplotype_CN.xlsx"),
                "Sheet": DATASETS[label].sheet_name,
                "CNV_Related_Genes": int(len(values)),
                "Min_Haplotype_Count": int(values.min()),
                "Q1": float(values.quantile(0.25)),
                "Median": float(values.median()),
                "Mean": float(values.mean()),
                "Q3": float(values.quantile(0.75)),
                "Max_Haplotype_Count": int(values.max()),
                "Total_Gene_Haplotype_Presences": int(values.sum()),
            }
        )
    return pd.DataFrame(rows)


def write_direct_crop_note(output_dir: Path, comparison_name: str) -> None:
    note = "\n".join(
        [
            "SupFig.16 haplotype-count 1-10 direct-crop note",
            f"Comparison: {comparison_name}",
            "The 1-10 figure is a direct crop from the full positive CNV-related gene distribution.",
            "Proportions are NOT recalculated within bins 1-10.",
            "The corresponding *_full_data.tsv keeps the full exact distribution, not only bins 1-10.",
            "The plot is generated by truncating that full distribution to bins 1-10.",
            "Denominator = all genes present in at least one haplotype in each cohort.",
            "Therefore, proportions across bins 1-10 are expected to sum to less than 1 when genes with haplotype counts >10 exist.",
            "Plot legends use display labels SCZ, Control, and Public.",
            "",
        ]
    )
    (output_dir / "SupFig16_HAPLOTYPE_COUNT_1_10_DIRECT_CROP_FULL_DENOMINATOR_NOTE.txt").write_text(
        note,
        encoding="utf-8",
    )


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    count_tables = {label: load_gene_haplotype_counts(label) for label in DATASETS}
    all_counts = pd.concat(count_tables.values(), ignore_index=True)
    summary = summarize_counts(count_tables)

    output_root = SCRIPT_DIR / f"output_SupFig16_hapcount_{DATE_TAG}"
    output_root.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(
        output_root / f"SupFig16_haplotype_count_gene_counts_case_cohort_{DATE_TAG}.xlsx",
        engine="openpyxl",
    ) as writer:
        all_counts.to_excel(writer, index=False, sheet_name="gene_haplotype_counts")
        summary.to_excel(writer, index=False, sheet_name="summary")

    all_counts.to_csv(
        output_root / f"SupFig16_haplotype_count_gene_counts_case_cohort_{DATE_TAG}.tsv",
        sep="\t",
        index=False,
    )
    summary.to_csv(
        output_root / f"SupFig16_haplotype_count_summary_case_cohort_{DATE_TAG}.tsv",
        sep="\t",
        index=False,
    )

    for comparison_name, group_labels in COMPARISONS.items():
        comparison_short_name = COMPARISON_SHORT_NAMES[comparison_name]
        output_dir = output_root / comparison_short_name
        output_dir.mkdir(parents=True, exist_ok=True)

        original = {
            label: group_haplotype_original(count_tables[label]["Haplotype_Count"])
            for label in group_labels
        }
        plot_grouped_bars(
            original,
            output_dir / f"SupFig16_HapCount_{comparison_short_name}_grouped_original",
            xlabel="Haplotype Count",
        )
        write_distribution_table(
            original,
            output_dir / f"SupFig16_HapCount_{comparison_short_name}_grouped_original.tsv",
        )

        exact_1_50 = {
            label: group_exact(count_tables[label]["Haplotype_Count"], 50)
            for label in group_labels
        }
        plot_grouped_bars(
            exact_1_50,
            output_dir / f"SupFig16_HapCount_{comparison_short_name}_1_50",
            xlabel="Haplotype Count",
            figsize=(16, 6),
            xtick_step=5,
        )
        write_distribution_table(
            exact_1_50,
            output_dir / f"SupFig16_HapCount_{comparison_short_name}_1_50.tsv",
        )

        full_exact = {
            label: exact_distribution_full_denominator(count_tables[label]["Haplotype_Count"])
            for label in group_labels
        }
        exact_1_10 = {label: crop_distribution_for_plot(full_exact[label], 10) for label in group_labels}
        plot_grouped_bars(
            exact_1_10,
            output_dir / f"SupFig16_HapCount_{comparison_short_name}_1_10_direct_crop_full_denominator",
            xlabel="Haplotype Count",
        )
        write_distribution_table(
            full_exact,
            output_dir / f"SupFig16_HapCount_{comparison_short_name}_1_10_direct_crop_full_denominator_full_data.tsv",
        )
        write_direct_crop_note(output_dir, comparison_name)
        print(f"Saved SupFig16 haplotype-count outputs: {output_dir}")

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
