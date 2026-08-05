# -*- coding: utf-8 -*-

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


SCRIPT_DIR = Path(__file__).resolve().parent
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


DATASETS = {
    CASE_LABEL: Dataset(CASE_LABEL, "slurm_scripts_case_cohort"),
    comparison_cohort_LABEL: Dataset(comparison_cohort_LABEL, "slurm_scripts_comparison_site"),
    ALLPUB_LABEL: Dataset(ALLPUB_LABEL, "slurm_scripts_public_reference"),
}

COMPARISONS = {
    "case_cohort_vs_comparison_cohort": [CASE_LABEL, comparison_cohort_LABEL],
    "case_cohort_vs_comparison_cohort_vs_public_reference": [CASE_LABEL, comparison_cohort_LABEL, ALLPUB_LABEL],
}


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


def load_extra_copy_numbers(label: str) -> pd.Series:
    table = pd.read_excel(dataset_dir(label) / "sample_CN.xlsx", sheet_name="Sum_ECN")
    counts = pd.to_numeric(table["Sum_ECN_true"], errors="coerce").fillna(0).astype(int)
    return counts[counts > 0].reset_index(drop=True)


def find_presence_file(directory: Path) -> Path:
    files = [path for path in directory.iterdir() if path.name.startswith("sample_CN_") and path.suffix == ".xlsx"]
    if not files:
        raise FileNotFoundError(f"Cannot find sample_CN_*.xlsx in {directory}")
    return sorted(files, key=lambda path: path.name)[0]


def load_gene_sample_counts(label: str) -> pd.Series:
    file_path = find_presence_file(dataset_dir(label))
    matrix = pd.read_excel(file_path, sheet_name=0, index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    counts = (matrix > 0).sum(axis=0).astype(int)
    counts = counts[counts > 0]
    counts.index = counts.index.astype(str)
    return counts.reset_index(drop=True)


def group_extra_original(values: pd.Series) -> pd.Series:
    labels = [str(i) for i in range(1, 10)] + ["10-99", "100+"]

    def group(value: int) -> str:
        if 1 <= value <= 9:
            return str(value)
        if 10 <= value <= 99:
            return "10-99"
        return "100+"

    return values.map(group).value_counts().reindex(labels, fill_value=0) / len(values)


def group_sample_original(values: pd.Series) -> pd.Series:
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


def plot_grouped_bars(
    proportions: dict[str, pd.Series],
    output_prefix: Path,
    xlabel: str,
    ylabel: str = "Proportion of CNV-related Genes",
    figsize: tuple[float, float] = (9, 5.5),
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


def write_distribution_table(distributions: dict[str, pd.Series], output_path: Path) -> None:
    table = pd.DataFrame(distributions)
    if all(str(index_value).isdigit() for index_value in table.index):
        table = table.loc[sorted(table.index, key=lambda index_value: int(str(index_value)))]
    table.index.name = "Bin"
    table.to_csv(output_path, sep="\t")


def write_direct_crop_note(output_dir: Path, comparison_name: str) -> None:
    note = "\n".join(
        [
            "Fig.5D/Fig.5E 1-10 direct-crop note",
            f"Comparison: {comparison_name}",
            "The 1-10 figures are direct crops from the full positive CNV-related gene distribution.",
            "Proportions are NOT recalculated within bins 1-10.",
            "The corresponding *_full_data.tsv keeps the full exact distribution, not only bins 1-10.",
            "The plot is generated by truncating that full distribution to bins 1-10.",
            "Denominator = all genes with positive extra copy number for Fig.5D, and all genes present in at least one sample for Fig.5E.",
            "Therefore, proportions across bins 1-10 are expected to sum to less than 1 when genes with values >10 exist.",
            "",
        ]
    )
    (output_dir / "Fig5D_Fig5E_1_10_DIRECT_CROP_FULL_DENOMINATOR_NOTE.txt").write_text(note, encoding="utf-8")


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    extra_values = {label: load_extra_copy_numbers(label) for label in DATASETS}
    sample_values = {label: load_gene_sample_counts(label) for label in DATASETS}

    for comparison_name, group_labels in COMPARISONS.items():
        output_dir = SCRIPT_DIR / f"output_Fig5D_Fig5E_{comparison_name}_26-6-16"
        output_dir.mkdir(parents=True, exist_ok=True)

        extra_original = {label: group_extra_original(extra_values[label]) for label in group_labels}
        sample_original = {label: group_sample_original(sample_values[label]) for label in group_labels}

        plot_grouped_bars(
            extra_original,
            output_dir / f"Fig5D_Extra_Copy_Number_{comparison_name}_grouped_original",
            xlabel="Extra Copy Number",
            figsize=(10.5, 6),
        )
        write_distribution_table(
            extra_original,
            output_dir / f"Fig5D_Extra_Copy_Number_{comparison_name}_grouped_original.tsv",
        )

        plot_grouped_bars(
            sample_original,
            output_dir / f"Fig5E_Sample_Count_{comparison_name}_grouped_original",
            xlabel="Sample Count",
            figsize=(10.5, 6),
        )
        write_distribution_table(
            sample_original,
            output_dir / f"Fig5E_Sample_Count_{comparison_name}_grouped_original.tsv",
        )

        extra_1_100 = {label: group_exact(extra_values[label], 100) for label in group_labels}
        plot_grouped_bars(
            extra_1_100,
            output_dir / f"Fig5D_Extra_Copy_Number_{comparison_name}_1_100",
            xlabel="Extra Copy Number",
            figsize=(20, 6),
            xtick_step=5,
        )
        write_distribution_table(extra_1_100, output_dir / f"Fig5D_Extra_Copy_Number_{comparison_name}_1_100.tsv")

        extra_full_exact = {label: exact_distribution_full_denominator(extra_values[label]) for label in group_labels}
        extra_1_10 = {label: crop_distribution_for_plot(extra_full_exact[label], 10) for label in group_labels}
        plot_grouped_bars(
            extra_1_10,
            output_dir / f"Fig5D_Extra_Copy_Number_{comparison_name}_1_10_direct_crop_full_denominator",
            xlabel="Extra Copy Number",
            figsize=(10.5, 6),
        )
        write_distribution_table(
            extra_full_exact,
            output_dir / f"Fig5D_Extra_Copy_Number_{comparison_name}_1_10_direct_crop_full_denominator_full_data.tsv",
        )

        sample_1_50 = {label: group_exact(sample_values[label], 50) for label in group_labels}
        plot_grouped_bars(
            sample_1_50,
            output_dir / f"Fig5E_Sample_Count_{comparison_name}_1_50",
            xlabel="Sample Count",
            figsize=(16, 6),
            xtick_step=5,
        )
        write_distribution_table(sample_1_50, output_dir / f"Fig5E_Sample_Count_{comparison_name}_1_50.tsv")

        sample_full_exact = {label: exact_distribution_full_denominator(sample_values[label]) for label in group_labels}
        sample_1_10 = {label: crop_distribution_for_plot(sample_full_exact[label], 10) for label in group_labels}
        plot_grouped_bars(
            sample_1_10,
            output_dir / f"Fig5E_Sample_Count_{comparison_name}_1_10_direct_crop_full_denominator",
            xlabel="Sample Count",
            figsize=(10.5, 6),
        )
        write_distribution_table(
            sample_full_exact,
            output_dir / f"Fig5E_Sample_Count_{comparison_name}_1_10_direct_crop_full_denominator_full_data.tsv",
        )

        write_direct_crop_note(output_dir, comparison_name)

        print(f"Saved Fig5D/Fig5E outputs: {output_dir}")

    print(f"CNV_BASE_DIR={CNV_BASE_DIR}")
    for label in DATASETS:
        print(
            f"{label}: extra_copy_genes={len(extra_values[label])}, "
            f"sample_count_genes={len(sample_values[label])}"
        )


if __name__ == "__main__":
    main()
