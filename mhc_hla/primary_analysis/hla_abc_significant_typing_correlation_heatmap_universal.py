from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


mpl.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "Arial",
        "axes.unicode_minus": False,
    }
)

GENE_ORDER = ["A", "B", "C"]
REVERSE_GENE_ORDER = ["C", "B", "A"]
PUBLIC_QC_SHEETS = {
    "A": "A_qc_passed",
    "B": "B_qc_passed",
    "C": "C_qc_passed",
}
LEGACY_QC_SHEETS = {
    "A": "A_validated",
    "B": "B_validated",
    "C": "C_validated",
}


@dataclass
class GeneSelection:
    gene: str
    alleles: list[str]
    source: str
    result_table: pd.DataFrame
    manual_alleles: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one combined ABC correlation heatmap using all selected significant alleles."
    )
    parser.add_argument("--scz-xlsx", required=True)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--result-tag", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--sample-label", default="case_cohort_ragtag")
    parser.add_argument("--comparison-label", default="")
    parser.add_argument("--min-scz-count", type=int, default=5)
    parser.add_argument("--allow-nominal-fallback", action="store_true")
    parser.add_argument("--nominal-p-threshold", type=float, default=0.05)
    parser.add_argument("--label-digits", type=int, default=2)
    parser.add_argument("--a-alleles", default="")
    parser.add_argument("--b-alleles", default="")
    parser.add_argument("--c-alleles", default="")
    return parser.parse_args()


def normalize_values(values: Iterable) -> list[str | None]:
    normalized: list[str | None] = []
    for value in values:
        if pd.isna(value):
            normalized.append(None)
            continue
        text = str(value).strip()
        if not text:
            normalized.append(None)
            continue
        if text.upper() == "NULL":
            normalized.append(None)
            continue
        if text.lower() == "unknown":
            normalized.append(None)
            continue
        normalized.append(text)
    return normalized


def split_csv(text: str) -> list[str]:
    if not text:
        return []
    return [part.strip() for part in text.split(",") if part.strip()]


def bh_fdr(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    prev = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        value = ranked[idx] * n / rank
        prev = min(prev, value)
        adjusted[idx] = min(prev, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out.tolist()


def load_gene_calls(xlsx_path: Path, gene: str) -> pd.DataFrame:
    xls = pd.ExcelFile(xlsx_path)
    candidates = [PUBLIC_QC_SHEETS[gene], LEGACY_QC_SHEETS[gene]]
    sheet = next((candidate for candidate in candidates if candidate in xls.sheet_names), xls.sheet_names[0])
    df = pd.read_excel(xlsx_path, sheet_name=sheet)
    if "Sample_name" not in df.columns or gene not in df.columns:
        raise ValueError(f"Sheet {sheet} missing Sample_name/{gene} columns")
    out = df[["Sample_name", gene]].copy()
    out.columns = ["Sample_name", "Genotype"]
    out["Sample_name"] = out["Sample_name"].astype(str)
    out["Genotype"] = normalize_values(out["Genotype"])
    out["Gene"] = gene
    out["Source_Sheet"] = sheet
    return out.drop_duplicates()


def read_result_file(result_dir: Path, gene: str, result_tag: str) -> pd.DataFrame:
    path = result_dir / f"{gene}_fisher_results_new_{result_tag}.xlsx"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_excel(path)
    if "Genotype" not in df.columns:
        raise ValueError(f"Missing Genotype column in {path}")
    df["Genotype"] = normalize_values(df["Genotype"])
    return df


def select_gene_alleles(
    gene: str,
    result_dir: Path,
    result_tag: str,
    min_scz_count: int,
    allow_nominal_fallback: bool,
    nominal_p_threshold: float,
    manual_alleles_text: str,
) -> GeneSelection:
    df = read_result_file(result_dir, gene, result_tag)
    working = df.copy()
    if "Significant_FDR_0.05" in working.columns:
        selected = working[working["Significant_FDR_0.05"] == True].copy()
    elif "FDR_BH" in working.columns:
        selected = working[working["FDR_BH"] <= 0.05].copy()
    else:
        selected = working.iloc[0:0].copy()

    if "SCZ_Count" in selected.columns:
        selected = selected[selected["SCZ_Count"] >= min_scz_count].copy()

    source = "fdr_or_significance_column"
    if selected.empty and allow_nominal_fallback and "P-Value" in working.columns:
        selected = working.copy()
        if "SCZ_Count" in selected.columns:
            selected = selected[selected["SCZ_Count"] >= min_scz_count].copy()
        selected = selected[selected["P-Value"] <= nominal_p_threshold].copy()
        source = f"nominal_p<={nominal_p_threshold}"

    if not selected.empty:
        sort_cols = [col for col in ["SCZ_Count", "FDR_BH", "P-Value", "Genotype"] if col in selected.columns]
        ascending = []
        for col in sort_cols:
            if col == "SCZ_Count":
                ascending.append(False)
            else:
                ascending.append(True)
        selected = selected.sort_values(sort_cols, ascending=ascending).copy()

    auto_alleles = [x for x in selected["Genotype"].tolist() if x]
    manual_alleles = [x for x in split_csv(manual_alleles_text) if x]

    final_alleles: list[str] = []
    for allele in auto_alleles + manual_alleles:
        if allele and allele not in final_alleles:
            final_alleles.append(allele)

    if not final_alleles and manual_alleles:
        source = "manual_only"

    return GeneSelection(
        gene=gene,
        alleles=final_alleles,
        source=source,
        result_table=selected.copy(),
        manual_alleles=manual_alleles,
    )


def compute_pearson_correlation(vec1: pd.Series, vec2: pd.Series) -> float | None:
    array1 = vec1.astype(int).to_numpy()
    array2 = vec2.astype(int).to_numpy()
    if array1.std() == 0 or array2.std() == 0:
        return None
    return float(np.corrcoef(array1, array2)[0, 1])


def compute_cross_gene_stats(
    calls_by_gene: dict[str, pd.DataFrame],
    selection_by_gene: dict[str, GeneSelection],
) -> pd.DataFrame:
    records: list[dict] = []
    p_values: list[float] = []

    for gene1, gene2 in [("A", "B"), ("A", "C"), ("B", "C")]:
        alleles1 = selection_by_gene[gene1].alleles
        alleles2 = selection_by_gene[gene2].alleles
        if not alleles1 or not alleles2:
            continue

        df1 = calls_by_gene[gene1][["Sample_name", "Genotype"]].rename(columns={"Genotype": "Genotype1"})
        df2 = calls_by_gene[gene2][["Sample_name", "Genotype"]].rename(columns={"Genotype": "Genotype2"})
        merged = df1.merge(df2, on="Sample_name", how="inner")
        merged = merged.dropna(subset=["Genotype1", "Genotype2"]).copy()
        pair_total = len(merged)
        if pair_total == 0:
            continue

        for allele1 in alleles1:
            vec1 = merged["Genotype1"] == allele1
            count1 = int(vec1.sum())
            for allele2 in alleles2:
                vec2 = merged["Genotype2"] == allele2
                count2 = int(vec2.sum())
                both = int((vec1 & vec2).sum())
                gene1_only = int((vec1 & ~vec2).sum())
                gene2_only = int((~vec1 & vec2).sum())
                neither = int((~vec1 & ~vec2).sum())
                contingency = [[both, gene1_only], [gene2_only, neither]]
                odds_ratio, p_value = fisher_exact(contingency, alternative="greater")
                pearson_correlation = compute_pearson_correlation(vec1, vec2)
                expected_both = (count1 * count2 / pair_total) if pair_total else np.nan
                enrichment = (both / expected_both) if pd.notna(expected_both) and expected_both > 0 else np.nan

                record = {
                    "Gene1": gene1,
                    "Allele1": allele1,
                    "Gene2": gene2,
                    "Allele2": allele2,
                    "Pair_Label": f"{gene1}:{allele1}__{gene2}:{allele2}",
                    "Pair_Total": pair_total,
                    "Count1": count1,
                    "Count2": count2,
                    "Both_Count": both,
                    "Gene1_Only": gene1_only,
                    "Gene2_Only": gene2_only,
                    "Neither": neither,
                    "Frequency1": count1 / pair_total,
                    "Frequency2": count2 / pair_total,
                    "Both_Frequency": both / pair_total,
                    "Both_in_Gene1": (both / count1) if count1 else np.nan,
                    "Both_in_Gene2": (both / count2) if count2 else np.nan,
                    "Expected_Both": expected_both,
                    "Enrichment_vs_Expected": enrichment,
                    "Pearson_Correlation": pearson_correlation,
                    "Odds_Ratio": odds_ratio,
                    "P_Value": p_value,
                    "Contingency_Table": str(contingency),
                }
                records.append(record)
                p_values.append(p_value)

    stats = pd.DataFrame(records)
    if stats.empty:
        return stats

    stats["FDR_BH"] = bh_fdr(p_values)
    stats["Significant_FDR_0.05"] = stats["FDR_BH"] <= 0.05
    return stats.sort_values(["FDR_BH", "P_Value", "Gene1", "Allele1", "Gene2", "Allele2"]).reset_index(drop=True)


def build_axis_items(selection_by_gene: dict[str, GeneSelection], reverse: bool = False) -> list[tuple[str, str]]:
    gene_order = REVERSE_GENE_ORDER if reverse else GENE_ORDER
    items: list[tuple[str, str]] = []
    for gene in gene_order:
        alleles = selection_by_gene[gene].alleles
        if reverse:
            alleles = list(reversed(alleles))
        for allele in alleles:
            items.append((gene, allele))
    return items


def stats_lookup(stats_df: pd.DataFrame) -> dict[tuple[str, str, str, str], dict]:
    lookup: dict[tuple[str, str, str, str], dict] = {}
    if stats_df.empty:
        return lookup
    for _, row in stats_df.iterrows():
        key = (row["Gene1"], row["Allele1"], row["Gene2"], row["Allele2"])
        lookup[key] = row.to_dict()
        reverse_key = (row["Gene2"], row["Allele2"], row["Gene1"], row["Allele1"])
        lookup[reverse_key] = row.to_dict()
    return lookup


def build_heatmap_tables(
    selection_by_gene: dict[str, GeneSelection],
    stats_df: pd.DataFrame,
    label_digits: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    col_items = build_axis_items(selection_by_gene, reverse=False)
    row_items = build_axis_items(selection_by_gene, reverse=False)
    lookup = stats_lookup(stats_df)

    matrix_rows: list[list[float]] = []
    label_rows: list[list[str]] = []
    long_records: list[dict] = []

    for row_gene, row_allele in row_items:
        value_row: list[float] = []
        label_row: list[str] = []
        for col_gene, col_allele in col_items:
            label = ""
            significant = False

            if row_gene == col_gene:
                if row_allele == col_allele:
                    value = 1.0
                    label = f"{1.0:.{label_digits}f}"
                    cell_type = "self"
                else:
                    value = np.nan
                    label = ""
                    cell_type = "same_gene_offdiag"
            else:
                stat = lookup.get((col_gene, col_allele, row_gene, row_allele))
                if stat is None:
                    value = np.nan
                    cell_type = "missing"
                else:
                    value = stat["Pearson_Correlation"] if pd.notna(stat["Pearson_Correlation"]) else np.nan
                    significant = bool(stat.get("Significant_FDR_0.05", False))
                    if pd.notna(value):
                        label = f"{value:.{label_digits}f}" + ("*" if significant else "")
                    cell_type = "cross_gene"

            value_row.append(value)
            label_row.append(label)
            long_records.append(
                {
                    "Row_Gene": row_gene,
                    "Row_Allele": row_allele,
                    "Col_Gene": col_gene,
                    "Col_Allele": col_allele,
                    "Pearson_Correlation": value,
                    "Label": label,
                    "Cell_Type": cell_type,
                    "Significant_FDR_0.05": significant,
                }
            )

        matrix_rows.append(value_row)
        label_rows.append(label_row)

    matrix_df = pd.DataFrame(
        matrix_rows,
        index=[allele for _, allele in row_items],
        columns=[allele for _, allele in col_items],
    )
    labels_df = pd.DataFrame(
        label_rows,
        index=[allele for _, allele in row_items],
        columns=[allele for _, allele in col_items],
    )
    long_df = pd.DataFrame(long_records)
    return matrix_df, labels_df, long_df


def add_group_lines(ax: plt.Axes, selection_by_gene: dict[str, GeneSelection], n_rows: int, n_cols: int) -> None:
    col_boundaries = []
    running = 0
    for gene in GENE_ORDER:
        running += len(selection_by_gene[gene].alleles)
        if 0 < running < n_cols:
            col_boundaries.append(running - 0.5)

    row_boundaries = []
    running = 0
    for gene in GENE_ORDER:
        running += len(selection_by_gene[gene].alleles)
        if 0 < running < n_rows:
            row_boundaries.append(running - 0.5)

    for x in col_boundaries:
        ax.axvline(x=x, color="black", linewidth=1.2)
    for y in row_boundaries:
        ax.axhline(y=y, color="black", linewidth=1.2)


def fallback_output_path(path: Path) -> Path:
    if not path.exists():
        return path
    for index in range(1, 100):
        suffix = "_updated" if index == 1 else f"_updated{index}"
        candidate = path.with_name(f"{path.stem}{suffix}{path.suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Unable to allocate fallback path for {path}")


def plot_heatmap(
    matrix_df: pd.DataFrame,
    labels_df: pd.DataFrame,
    selection_by_gene: dict[str, GeneSelection],
    output_pdf: Path,
    output_png: Path,
    title: str,
) -> tuple[Path, Path]:
    n_rows, n_cols = matrix_df.shape
    figsize = (max(8, 1.2 * n_cols + 2), max(7, 1.0 * n_rows + 2.5))

    fig, ax = plt.subplots(figsize=figsize)
    cmap = LinearSegmentedColormap.from_list(
        "hla_red_white_blue",
        ["#1f78b4", "white", "#da7271"],
    )
    cmap.set_bad("#d9d9d9")
    norm = Normalize(vmin=-1, vmax=1)

    values = matrix_df.to_numpy(dtype=float)
    for row_idx in range(n_rows):
        for col_idx in range(n_cols):
            value = values[row_idx, col_idx]
            color = "#d9d9d9" if np.isnan(value) else cmap(norm(value))
            ax.add_patch(
                Rectangle(
                    (col_idx - 0.5, row_idx - 0.5),
                    1,
                    1,
                    facecolor=color,
                    edgecolor="none",
                    linewidth=0,
                )
            )

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_aspect("auto")

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(matrix_df.columns.tolist(), rotation=15, ha="right", fontsize=9)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(matrix_df.index.tolist(), fontsize=9)

    for row_idx in range(n_rows):
        for col_idx in range(n_cols):
            label = labels_df.iat[row_idx, col_idx]
            if label:
                ax.text(col_idx, row_idx, label, ha="center", va="center", fontsize=8)

    add_group_lines(ax, selection_by_gene, n_rows=n_rows, n_cols=n_cols)

    ax.set_title(title, fontsize=16, fontweight="bold", pad=18)
    scalar_mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Pearson Correlation", fontsize=8)
    colorbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    actual_pdf = output_pdf
    try:
        fig.savefig(actual_pdf)
    except PermissionError:
        actual_pdf = fallback_output_path(output_pdf)
        fig.savefig(actual_pdf)

    actual_png = output_png
    try:
        fig.savefig(actual_png, dpi=300, bbox_inches="tight")
    except PermissionError:
        actual_png = fallback_output_path(output_png)
        fig.savefig(actual_png, dpi=300, bbox_inches="tight")

    plt.close(fig)
    return actual_pdf, actual_png


def build_selected_alleles_sheet(selection_by_gene: dict[str, GeneSelection]) -> pd.DataFrame:
    rows = []
    for gene in GENE_ORDER:
        selection = selection_by_gene[gene]
        for allele in selection.alleles:
            matched = selection.result_table[selection.result_table["Genotype"] == allele]
            row = {
                "Gene": gene,
                "Allele": allele,
                "Selection_Source": selection.source,
                "Selected_Manually": allele in selection.manual_alleles,
            }
            if not matched.empty:
                first = matched.iloc[0].to_dict()
                row.update({f"Result_{key}": value for key, value in first.items()})
            rows.append(row)
    return pd.DataFrame(rows)


def write_summary_xlsx(
    output_xlsx: Path,
    settings_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    cross_gene_stats_df: pd.DataFrame,
    matrix_long_df: pd.DataFrame,
    matrix_df: pd.DataFrame,
) -> Path:
    actual_xlsx = output_xlsx
    try:
        writer_context = pd.ExcelWriter(actual_xlsx, engine="openpyxl")
    except PermissionError:
        actual_xlsx = fallback_output_path(output_xlsx)
        writer_context = pd.ExcelWriter(actual_xlsx, engine="openpyxl")

    with writer_context as writer:
        settings_df.to_excel(writer, sheet_name="settings", index=False)
        selected_df.to_excel(writer, sheet_name="selected_alleles", index=False)
        cross_gene_stats_df.to_excel(writer, sheet_name="cross_gene_stats", index=False)
        matrix_long_df.to_excel(writer, sheet_name="matrix_long", index=False)
        matrix_df.reset_index().rename(columns={"index": "Row_Allele"}).to_excel(writer, sheet_name="matrix_values", index=False)
    return actual_xlsx


def main() -> None:
    args = parse_args()
    scz_xlsx = Path(args.scz_xlsx)
    result_dir = Path(args.result_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    calls_by_gene = {gene: load_gene_calls(scz_xlsx, gene) for gene in GENE_ORDER}

    manual_map = {
        "A": args.a_alleles,
        "B": args.b_alleles,
        "C": args.c_alleles,
    }
    selection_by_gene = {
        gene: select_gene_alleles(
            gene=gene,
            result_dir=result_dir,
            result_tag=args.result_tag,
            min_scz_count=args.min_scz_count,
            allow_nominal_fallback=args.allow_nominal_fallback,
            nominal_p_threshold=args.nominal_p_threshold,
            manual_alleles_text=manual_map[gene],
        )
        for gene in GENE_ORDER
    }

    for gene in GENE_ORDER:
        observed = set(calls_by_gene[gene]["Genotype"].dropna())
        selection_by_gene[gene].alleles = [allele for allele in selection_by_gene[gene].alleles if allele in observed]

    total_selected = sum(len(selection_by_gene[gene].alleles) for gene in GENE_ORDER)
    if total_selected == 0:
        raise ValueError("No alleles selected for combined heatmap.")

    cross_gene_stats_df = compute_cross_gene_stats(calls_by_gene, selection_by_gene)
    matrix_df, labels_df, matrix_long_df = build_heatmap_tables(
        selection_by_gene=selection_by_gene,
        stats_df=cross_gene_stats_df,
        label_digits=args.label_digits,
    )

    comparison_text = f", {args.comparison_label}" if args.comparison_label else ""
    title = f"HLA-A/B/C combined correlation heatmap ({args.sample_label}{comparison_text})"

    output_pdf = output_dir / f"{args.output_prefix}_heatmap.pdf"
    output_png = output_dir / f"{args.output_prefix}_heatmap.png"
    actual_pdf, actual_png = plot_heatmap(
        matrix_df=matrix_df,
        labels_df=labels_df,
        selection_by_gene=selection_by_gene,
        output_pdf=output_pdf,
        output_png=output_png,
        title=title,
    )

    settings_df = pd.DataFrame(
        {
            "Parameter": [
                "scz_xlsx",
                "result_dir",
                "result_tag",
                "sample_label",
                "comparison_label",
                "min_scz_count",
                "allow_nominal_fallback",
                "nominal_p_threshold",
                "selected_A",
                "selected_B",
                "selected_C",
                "source_A",
                "source_B",
                "source_C",
            ],
            "Value": [
                str(scz_xlsx),
                str(result_dir),
                args.result_tag,
                args.sample_label,
                args.comparison_label,
                args.min_scz_count,
                args.allow_nominal_fallback,
                args.nominal_p_threshold,
                len(selection_by_gene["A"].alleles),
                len(selection_by_gene["B"].alleles),
                len(selection_by_gene["C"].alleles),
                selection_by_gene["A"].source,
                selection_by_gene["B"].source,
                selection_by_gene["C"].source,
            ],
        }
    )
    selected_df = build_selected_alleles_sheet(selection_by_gene)
    output_xlsx = output_dir / f"{args.output_prefix}_summary.xlsx"
    actual_xlsx = write_summary_xlsx(
        output_xlsx=output_xlsx,
        settings_df=settings_df,
        selected_df=selected_df,
        cross_gene_stats_df=cross_gene_stats_df,
        matrix_long_df=matrix_long_df,
        matrix_df=matrix_df,
    )

    print(f"Heatmap PDF: {actual_pdf}")
    print(f"Heatmap PNG: {actual_png}")
    print(f"Summary XLSX: {actual_xlsx}")


if __name__ == "__main__":
    main()
