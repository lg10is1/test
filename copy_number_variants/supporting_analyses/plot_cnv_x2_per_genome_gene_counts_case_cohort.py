# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATE_TAG = "26-6-18"
CASE_LABEL = "SCZ"
comparison_cohort_LABEL = "comparison_cohort"
PUB_LABEL = "PUB"

COLORS = {
    CASE_LABEL: "#da7271",
    comparison_cohort_LABEL: "#1f78b4",
    PUB_LABEL: "#d9a627",
}


@dataclass(frozen=True)
class Dataset:
    cohort: str
    input_label: str
    relative_file: Path
    sheet_name: str
    legacy_relative_file: Path | None = None
    legacy_sheet_name: str | None = None


CNV_BASE = Path(os.environ.get("EOSCZ_CNV_BASE_DIR", Path(__file__).resolve().parent / "cnv_base"))
OUTPUT_DIR = Path(__file__).resolve().parent

PUBLIC_BINARY_SAMPLE_CN = "sample_CN_binary.xlsx"
LEGACY_BINARY_SAMPLE_CN = "sample_copy_number_presence.xlsx"
PUBLIC_BINARY_SHEET = "binary_matrix"
LEGACY_BINARY_SHEET = "matrix"

DATASETS = (
    Dataset(
        cohort=CASE_LABEL,
        input_label="case_cohort",
        relative_file=Path("slurm_scripts_case_cohort") / PUBLIC_BINARY_SAMPLE_CN,
        sheet_name=PUBLIC_BINARY_SHEET,
        legacy_relative_file=Path("slurm_scripts_case_cohort") / LEGACY_BINARY_SAMPLE_CN,
        legacy_sheet_name=LEGACY_BINARY_SHEET,
    ),
    Dataset(
        cohort=comparison_cohort_LABEL,
        input_label="comparison_cohort",
        relative_file=Path("slurm_scripts_comparison_site") / PUBLIC_BINARY_SAMPLE_CN,
        sheet_name=PUBLIC_BINARY_SHEET,
        legacy_relative_file=Path("slurm_scripts_comparison_site") / LEGACY_BINARY_SAMPLE_CN,
        legacy_sheet_name=LEGACY_BINARY_SHEET,
    ),
    Dataset(
        cohort=PUB_LABEL,
        input_label="public_reference",
        relative_file=Path("slurm_scripts_public_reference") / PUBLIC_BINARY_SAMPLE_CN,
        sheet_name="Sheet1",
        legacy_relative_file=Path("slurm_scripts_public_reference") / LEGACY_BINARY_SAMPLE_CN,
        legacy_sheet_name="Sheet1",
    ),
)


def load_matrix(dataset: Dataset) -> pd.DataFrame:
    path = CNV_BASE / dataset.relative_file
    if not path.exists() and dataset.legacy_relative_file is not None:
        path = CNV_BASE / dataset.legacy_relative_file
    if not path.exists():
        raise FileNotFoundError(f"CNV sample binary matrix not found: {CNV_BASE / dataset.relative_file}")

    workbook = pd.ExcelFile(path)
    if dataset.sheet_name in workbook.sheet_names:
        sheet_name = dataset.sheet_name
    elif dataset.legacy_sheet_name and dataset.legacy_sheet_name in workbook.sheet_names:
        sheet_name = dataset.legacy_sheet_name
    else:
        sheet_name = workbook.sheet_names[0]
    matrix = pd.read_excel(path, sheet_name=sheet_name, index_col=0)
    matrix = matrix.dropna(how="all")
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return (matrix > 0).astype(int)


def summarize(values: pd.Series) -> dict[str, float | int | str]:
    return {
        "N_Assemblies": int(values.size),
        "Mean": float(values.mean()),
        "Median": float(values.median()),
        "SD": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        "Min": int(values.min()),
        "Q1": float(values.quantile(0.25)),
        "Q3": float(values.quantile(0.75)),
        "Max": int(values.max()),
    }


def plot_counts(per_genome: pd.DataFrame, output_prefix: Path) -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    cohorts = [CASE_LABEL, comparison_cohort_LABEL, PUB_LABEL]
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    values_by_group = [
        per_genome.loc[per_genome["Cohort"] == cohort, "CNV_Related_ProteinCoding_Genes"].to_numpy()
        for cohort in cohorts
    ]

    box = ax.boxplot(
        values_by_group,
        positions=np.arange(1, len(cohorts) + 1),
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.1},
        boxprops={"linewidth": 0.8},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
    )
    for patch, cohort in zip(box["boxes"], cohorts):
        patch.set_facecolor(COLORS[cohort])
        patch.set_alpha(0.48)
        patch.set_edgecolor("black")

    rng = np.random.default_rng(20260618)
    for x, cohort, values in zip(range(1, len(cohorts) + 1), cohorts, values_by_group):
        jitter = rng.uniform(-0.15, 0.15, size=len(values))
        ax.scatter(
            np.full(len(values), x) + jitter,
            values,
            s=11,
            alpha=0.65,
            color=COLORS[cohort],
            edgecolor="white",
            linewidth=0.25,
        )

    ax.set_xticks(np.arange(1, len(cohorts) + 1))
    ax.set_xticklabels(cohorts, fontsize=11)
    ax.set_ylabel("CNV-related protein-coding genes per assembly", fontsize=11)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    per_genome_frames = []
    manifest_rows = []

    for dataset in DATASETS:
        matrix = load_matrix(dataset)
        counts = matrix.sum(axis=1).astype(int)
        frame = pd.DataFrame(
            {
                "Cohort": dataset.cohort,
                "Input_Label": dataset.input_label,
                "Sample": counts.index.astype(str),
                "CNV_Related_ProteinCoding_Genes": counts.to_numpy(),
            }
        )
        per_genome_frames.append(frame)
        manifest_rows.append(
            {
                "Cohort": dataset.cohort,
                "Input_Label": dataset.input_label,
                "Input_File": str(CNV_BASE / dataset.relative_file),
                "Sheet": dataset.sheet_name,
                "N_Assemblies": matrix.shape[0],
                "N_Genes": matrix.shape[1],
            }
        )

    per_genome = pd.concat(per_genome_frames, ignore_index=True)
    summary_rows = []
    for (cohort, input_label), group in per_genome.groupby(["Cohort", "Input_Label"], sort=False):
        row = {"Cohort": cohort, "Input_Label": input_label}
        row.update(summarize(group["CNV_Related_ProteinCoding_Genes"]))
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    manifest = pd.DataFrame(manifest_rows)

    output_prefix = OUTPUT_DIR / f"CNV_X2_per_genome_gene_counts_case_cohort_comparison_cohort_PUB_{DATE_TAG}"
    output_xlsx = output_prefix.with_suffix(".xlsx")
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        per_genome.to_excel(writer, index=False, sheet_name="per_genome_counts")
        summary.to_excel(writer, index=False, sheet_name="summary")
        manifest.to_excel(writer, index=False, sheet_name="input_manifest")

    per_genome.to_csv(output_prefix.with_name(output_prefix.name + "_per_genome_counts.tsv"), sep="\t", index=False)
    summary.to_csv(output_prefix.with_name(output_prefix.name + "_summary.tsv"), sep="\t", index=False)
    manifest.to_csv(output_prefix.with_name(output_prefix.name + "_input_manifest.tsv"), sep="\t", index=False)
    plot_counts(per_genome, output_prefix)

    print(summary.to_string(index=False))
    print(f"Saved: {output_xlsx}")
    print(f"Saved: {output_prefix.with_suffix('.pdf')}")
    print(f"Saved: {output_prefix.with_suffix('.png')}")


if __name__ == "__main__":
    main()
