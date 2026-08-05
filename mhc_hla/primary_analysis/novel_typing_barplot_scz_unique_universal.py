from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


GROUP_ORDER = ["C", "B", "A"]
COLORS = {
    "A": "#47a1a2",
    "B": "#da7271",
    "C": "#1f78b4",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot SCZ-unique new HLA subtype counts from manually checked *_ge2 sheets."
    )
    parser.add_argument("--input-xlsx", required=True, help="Input workbook containing A_ge2/B_ge2/C_ge2.")
    parser.add_argument("--output-prefix", required=True, help="Output prefix without extension.")
    parser.add_argument(
        "--title",
        default="",
        help="Optional plot title.",
    )
    return parser.parse_args()


def find_sczu_col(df: pd.DataFrame) -> str:
    matches = [col for col in df.columns if "SCZ unique" in str(col)]
    if not matches:
        raise ValueError("Could not find final SCZ-unique decision column.")
    return matches[0]


def normalize_yes(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().upper() == "Y"


def load_counts(path: Path, gene: str) -> pd.DataFrame:
    sheet = f"{gene}_ge2"
    df = pd.read_excel(path, sheet_name=sheet)
    if gene not in df.columns:
        raise ValueError(f"{sheet} missing {gene} column")

    final_col = find_sczu_col(df)
    sub = df[[gene, final_col]].copy()
    sub = sub[sub[gene].notna()].copy()
    sub[gene] = sub[gene].astype(str).str.strip()
    sub = sub[sub[gene] != ""].copy()
    sub = sub[sub[gene].str.contains(r"^[ABC]\*", regex=True)].copy()
    sub = sub[sub[final_col].map(normalize_yes)].copy()

    counts = sub[gene].value_counts().rename_axis("Genotype").reset_index(name="Count")
    counts["Group"] = gene
    counts["Source_Sheet"] = sheet
    counts["Selection_Rule"] = "Final sequence-level manual decision: SCZ unique == Y"
    return counts


def build_summary(path: Path) -> pd.DataFrame:
    frames = [load_counts(path, gene) for gene in ["A", "B", "C"]]
    combined = pd.concat(frames, ignore_index=True)
    if combined.empty:
        raise ValueError("No SCZ-unique subtypes were selected.")

    combined["Group"] = pd.Categorical(combined["Group"], categories=GROUP_ORDER, ordered=True)
    combined = combined.sort_values(["Group", "Count", "Genotype"], ascending=[True, True, True]).reset_index(drop=True)
    return combined


def plot_counts(summary: pd.DataFrame, output_prefix: Path, title: str) -> None:
    plt.figure(figsize=(12, 10))

    for group in GROUP_ORDER:
        group_data = summary[summary["Group"] == group]
        if group_data.empty:
            continue
        plt.barh(
            group_data["Genotype"],
            group_data["Count"],
            color=COLORS[group],
            label=group,
            height=0.5,
        )

    plt.xlabel("Haplotype Count of SCZ-unique New Subtypes", fontsize=28)
    plt.ylabel("Subtypes", fontsize=28)
    if title:
        plt.title(title, fontsize=22)

    handles, labels = plt.gca().get_legend_handles_labels()
    legend_order = ["A", "B", "C"]
    new_handles = [handles[labels.index(group)] for group in legend_order if group in labels]
    new_labels = [group for group in legend_order if group in labels]
    if new_handles:
        plt.legend(new_handles, new_labels, fontsize=24)

    plt.tick_params(axis="both", labelsize=24)
    plt.tight_layout()
    plt.savefig(output_prefix.with_suffix(".pdf"), dpi=600, format="pdf")
    plt.savefig(output_prefix.with_suffix(".png"), dpi=600, format="png")
    plt.close()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_xlsx)
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary = build_summary(input_path)
    plot_counts(summary, output_prefix=output_prefix, title=args.title)

    summary.to_excel(output_prefix.with_name(output_prefix.name + "_summary.xlsx"), index=False)

    print(f"Plot PDF: {output_prefix.with_suffix('.pdf')}")
    print(f"Plot PNG: {output_prefix.with_suffix('.png')}")
    print(f"Summary XLSX: {output_prefix.with_name(output_prefix.name + '_summary.xlsx')}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
