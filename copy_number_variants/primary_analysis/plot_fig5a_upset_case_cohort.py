# -*- coding: utf-8 -*-

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet


SCRIPT_DIR = Path(__file__).resolve().parent
CASE_LABEL = "case_cohort"
comparison_cohort_LABEL = "comparison_cohort"
ALLPUB_LABEL = "public_reference"

OUTPUT_PDF = SCRIPT_DIR / "Fig5A_upset_case_cohort_comparison_cohort_public_reference_26-6-16.pdf"
OUTPUT_PNG = SCRIPT_DIR / "Fig5A_upset_case_cohort_comparison_cohort_public_reference_26-6-16.png"
OUTPUT_OVERLAP = SCRIPT_DIR / "Fig5A_overlap_genes_case_cohort_comparison_cohort_public_reference_26-6-16.txt"


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
INPUTS = {
    CASE_LABEL: CNV_BASE_DIR / "slurm_scripts_case_cohort" / "gene_frequencies_filtered.txt",
    comparison_cohort_LABEL: CNV_BASE_DIR / "slurm_scripts_comparison_site" / "gene_frequencies_filtered.txt",
    ALLPUB_LABEL: CNV_BASE_DIR / "slurm_scripts_public_reference" / "gene_frequencies_filtered.txt",
}


def read_gene_frequency_table(file_path: Path) -> pd.DataFrame:
    table = pd.read_table(file_path, header=None, names=["Gene", "Count", "Frequency"])
    table["Gene"] = table["Gene"].astype(str)
    table["Count"] = pd.to_numeric(table["Count"], errors="coerce").fillna(0).astype(int)
    table["Frequency"] = pd.to_numeric(table["Frequency"], errors="coerce").fillna(0.0)
    return table


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"

    tables = {label: read_gene_frequency_table(path) for label, path in INPUTS.items()}
    gene_sets = {label: set(table["Gene"]) for label, table in tables.items()}

    memberships = Counter()
    all_genes = sorted(set().union(*gene_sets.values()))
    for gene in all_genes:
        memberships[tuple(gene in gene_sets[label] for label in INPUTS)] += 1

    upset_data = pd.Series(memberships, dtype=int)
    upset_data.index = pd.MultiIndex.from_tuples(upset_data.index, names=list(INPUTS))

    fig = plt.figure(figsize=(11, 6.5))
    upset = UpSet(upset_data, show_counts=True, sort_by="cardinality", element_size=48)
    upset.plot(fig=fig)
    fig.suptitle("CNV-related gene overlap", fontsize=16)
    fig.savefig(OUTPUT_PDF, dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)

    overlap_rows = []
    for gene in sorted(set.intersection(*gene_sets.values())):
        row = {"Gene": gene}
        for label, table in tables.items():
            values = table.loc[table["Gene"] == gene].iloc[0]
            row[f"{label}_Count"] = int(values["Count"])
            row[f"{label}_Frequency"] = float(values["Frequency"])
        overlap_rows.append(row)
    pd.DataFrame(overlap_rows).to_csv(OUTPUT_OVERLAP, sep="\t", index=False)

    print(f"CNV_BASE_DIR={CNV_BASE_DIR}")
    for label, genes in gene_sets.items():
        print(f"{label}: {len(genes)} genes")
    print(f"Union: {len(all_genes)} genes")
    print(f"Shared by all three groups: {len(overlap_rows)} genes")
    print(f"Saved: {OUTPUT_PDF}")
    print(f"Saved: {OUTPUT_PNG}")
    print(f"Saved: {OUTPUT_OVERLAP}")


if __name__ == "__main__":
    main()
