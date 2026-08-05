# -*- coding: utf-8 -*-

from __future__ import annotations

import re
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

CNV_BASE = Path(os.environ.get("EOSCZ_CNV_BASE_DIR", Path(__file__).resolve().parent / "cnv_base"))
OUTPUT_DIR = Path(os.environ.get("EOSCZ_CNV_SUPPORT_OUTPUT_DIR", Path(__file__).resolve().parent / "outputs"))
QC_XLSX = OUTPUT_DIR / "sample_CN_vs_haplotype_CN_QC_26-7-14.xlsx"
RESULT_PREFIX = "true_case_cohort_comparison_cohort_public_reference_sample_common_CTRLgt1pct_SCZrecurrent2"
PUBLIC_BINARY_SAMPLE_CN = "sample_CN_binary.xlsx"
LEGACY_BINARY_SAMPLE_CN = "sample_copy_number_presence.xlsx"
PUBLIC_BINARY_SHEET = "binary_matrix"
LEGACY_BINARY_SHEET = "matrix"

COHORTS = {
    "SCZ_filtered_cases": CNV_BASE / "slurm_scripts_case_cohort",
    "comparison_cohort": CNV_BASE / "slurm_scripts_comparison_site",
    "public_reference": CNV_BASE / "slurm_scripts_public_reference",
}

SCZ_PATTERNS = (
    (re.compile(r"^(?P<sample>.+?)\.1\.mat_R1_t2t\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>.+?)\.2\.pat_R1_t2t\.2\.scaffold$"), "hap2"),
    (re.compile(r"^(?P<sample>.+?)\.1\.mat\.R1_t2t\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>.+?)\.2\.pat\.R1_t2t\.2\.scaffold$"), "hap2"),
    (re.compile(r"^(?P<sample>.+?)\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>.+?)\.2\.scaffold$"), "hap2"),
)

ALLPUB_PATTERNS = (
    (re.compile(r"^(?P<sample>[^.]+)\.1\.scaffold$"), "hap1"),
    (re.compile(r"^(?P<sample>[^.]+)\.2\.scaffold$"), "hap2"),
    (re.compile(r"^(?P<sample>[^.]+)\.1\.(?:mat_R1_rmN_singleline|mat\.R1_rmN_singleline)$"), "hap1"),
    (re.compile(r"^(?P<sample>[^.]+)\.2\.(?:pat_R1_rmN_singleline|pat\.R1_rmN_singleline)$"), "hap2"),
    (re.compile(r"^\d+_(?P<sample>[^.]+)\.(?P<hap>mat|pat|hap1|hap2)\.fa(?:\..+)?$"), None),
    (re.compile(r"^(?P<sample>[^.]+)\..*?-hap(?P<hap>[12])$"), None),
    (re.compile(r"^(?P<sample>[^.]+)\.(?P<hap>[12])$"), None),
    (re.compile(r"^v\d+_(?P<sample>[^._]+).*\.h(?P<hap>[12])(?:[-.].*)?$"), None),
    (re.compile(r"^\d+_(?P<sample>[^.]+)\.pri\.fa(?:\..+)?$"), "single"),
)

HAP_MAP = {"1": "hap1", "2": "hap2", "mat": "hap1", "pat": "hap2", "hap1": "hap1", "hap2": "hap2"}


def load_hap_matrix(cohort_dir: Path) -> pd.DataFrame:
    df = pd.read_excel(cohort_dir / "haplotype_CN.xlsx", index_col=0)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def parse_scz(label: str) -> tuple[str, str]:
    for pattern, hap in SCZ_PATTERNS:
        match = pattern.match(str(label))
        if match:
            return match.group("sample"), hap
    raise ValueError(f"Unrecognized SCZ haplotype index: {label}")


def parse_comparison_cohort(label: str) -> tuple[str, str]:
    match = re.match(r"^(?P<sample>.+?)_(?P<hap>Mat|Pat)\.v0\.9$", str(label))
    if not match:
        raise ValueError(f"Unrecognized comparison_cohort haplotype index: {label}")
    return match.group("sample"), HAP_MAP[match.group("hap").lower()]


def parse_public_reference(label: str) -> tuple[str, str]:
    text = str(label)
    for pattern, fixed_hap in ALLPUB_PATTERNS:
        match = pattern.match(text)
        if not match:
            continue
        sample = match.group("sample")
        hap = fixed_hap if fixed_hap is not None else HAP_MAP[match.group("hap")]
        return sample, hap
    raise ValueError(f"Unrecognized public_reference haplotype index: {label}")


def build_metadata(hap_df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    parser = {"SCZ_filtered_cases": parse_scz, "comparison_cohort": parse_comparison_cohort, "public_reference": parse_public_reference}[cohort]
    rows = []
    for raw_index in hap_df.index:
        sample, hap = parser(str(raw_index))
        rows.append({"raw_index": str(raw_index), "sample": sample, "haplotype_group": hap})
    return pd.DataFrame(rows)


def rebuild_sample_cn(hap_df: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    matrix = hap_df.copy()
    matrix.index = pd.MultiIndex.from_frame(metadata[["sample", "haplotype_group"]])
    hap_level = matrix.groupby(level=[0, 1]).max()
    sample = hap_level.groupby(level=0).sum()
    sample.index.name = "sample"
    return sample.astype(int)


def read_existing_sample_table(cohort_dir: Path, cohort: str) -> tuple[pd.DataFrame, str]:
    xls = pd.ExcelFile(cohort_dir / "sample_CN.xlsx")
    preferred = {"SCZ_filtered_cases": "sample_CN", "comparison_cohort": "CN_true", "public_reference": "sample_CN"}[cohort]
    sheet = preferred if preferred in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(cohort_dir / "sample_CN.xlsx", sheet_name=sheet, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int), sheet


def read_existing_binary_table(cohort_dir: Path, cohort: str) -> tuple[pd.DataFrame, str]:
    public_path = cohort_dir / PUBLIC_BINARY_SAMPLE_CN
    legacy_path = cohort_dir / LEGACY_BINARY_SAMPLE_CN
    binary_path = public_path if public_path.exists() else legacy_path
    xls = pd.ExcelFile(binary_path)
    if PUBLIC_BINARY_SHEET in xls.sheet_names:
        sheet = PUBLIC_BINARY_SHEET
    elif LEGACY_BINARY_SHEET in xls.sheet_names:
        sheet = LEGACY_BINARY_SHEET
    elif "Sheet1" in xls.sheet_names:
        sheet = "Sheet1"
    else:
        sheet = xls.sheet_names[0]
    df = pd.read_excel(binary_path, sheet_name=sheet, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int), f"{binary_path.name}:{sheet}"


def compare_matrices(label: str, rebuilt: pd.DataFrame, existing: pd.DataFrame, table_name: str, sheet: str) -> tuple[dict, pd.DataFrame]:
    rebuilt = rebuilt.copy()
    existing = existing.copy()
    rebuilt.index = rebuilt.index.astype(str)
    rebuilt.columns = rebuilt.columns.astype(str)
    existing.index = existing.index.astype(str)
    existing.columns = existing.columns.astype(str)

    common_rows = sorted(set(rebuilt.index) & set(existing.index))
    common_cols = sorted(set(rebuilt.columns) & set(existing.columns))
    row_missing_in_existing = sorted(set(rebuilt.index) - set(existing.index))
    row_extra_in_existing = sorted(set(existing.index) - set(rebuilt.index))
    col_missing_in_existing = sorted(set(rebuilt.columns) - set(existing.columns))
    col_extra_in_existing = sorted(set(existing.columns) - set(rebuilt.columns))

    mismatch_count = np.nan
    max_abs_diff = np.nan
    examples = pd.DataFrame()
    if common_rows and common_cols:
        a = rebuilt.loc[common_rows, common_cols]
        b = existing.loc[common_rows, common_cols]
        diff = a - b
        mismatch_mask = diff.ne(0)
        mismatch_count = int(mismatch_mask.to_numpy().sum())
        max_abs_diff = int(diff.abs().to_numpy().max()) if mismatch_count else 0
        if mismatch_count:
            coords = np.argwhere(mismatch_mask.to_numpy())[:200]
            examples = pd.DataFrame(
                [
                    {
                        "Cohort": label,
                        "Table": table_name,
                        "Sheet": sheet,
                        "Sample": common_rows[i],
                        "Gene": common_cols[j],
                        "Rebuilt": int(a.iat[i, j]),
                        "Existing": int(b.iat[i, j]),
                        "Diff": int(diff.iat[i, j]),
                    }
                    for i, j in coords
                ]
            )

    summary = {
        "Cohort": label,
        "Table": table_name,
        "Sheet": sheet,
        "RebuiltRows": rebuilt.shape[0],
        "ExistingRows": existing.shape[0],
        "CommonRows": len(common_rows),
        "RowsMissingInExisting": len(row_missing_in_existing),
        "RowsExtraInExisting": len(row_extra_in_existing),
        "RebuiltCols": rebuilt.shape[1],
        "ExistingCols": existing.shape[1],
        "CommonCols": len(common_cols),
        "ColsMissingInExisting": len(col_missing_in_existing),
        "ColsExtraInExisting": len(col_extra_in_existing),
        "ValueMismatchCount": mismatch_count,
        "MaxAbsDiff": max_abs_diff,
        "MissingRowsExample": ", ".join(row_missing_in_existing[:20]),
        "ExtraRowsExample": ", ".join(row_extra_in_existing[:20]),
        "MissingColsExample": ", ".join(col_missing_in_existing[:20]),
        "ExtraColsExample": ", ".join(col_extra_in_existing[:20]),
    }
    return summary, examples


def cohort_qc() -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rebuilt_tables = {}
    summaries = []
    metadata_rows = []
    mismatch_examples = []
    hap_group_rows = []

    for cohort, cohort_dir in COHORTS.items():
        hap = load_hap_matrix(cohort_dir)
        metadata = build_metadata(hap, cohort)
        rebuilt = rebuild_sample_cn(hap, metadata)
        rebuilt_tables[cohort] = rebuilt

        hap_counts = metadata.groupby("sample").size().reset_index(name="haplotype_rows")
        hap_group_counts = (
            metadata.groupby(["sample", "haplotype_group"]).size().reset_index(name="rows_per_sample_hap_group")
        )
        for _, row in hap_counts.iterrows():
            hap_group_rows.append({"Cohort": cohort, **row.to_dict()})
        metadata_rows.append(
            {
                "Cohort": cohort,
                "HaplotypeRows": hap.shape[0],
                "Genes": hap.shape[1],
                "ParsedRows": metadata.shape[0],
                "RebuiltSamples": rebuilt.shape[0],
                "SamplesWith1HaplotypeRow": int((hap_counts["haplotype_rows"] == 1).sum()),
                "SamplesWith2HaplotypeRows": int((hap_counts["haplotype_rows"] == 2).sum()),
                "SamplesWithMoreThan2HaplotypeRows": int((hap_counts["haplotype_rows"] > 2).sum()),
                "SampleHapGroupsWithDuplicatedRows": int((hap_group_counts["rows_per_sample_hap_group"] > 1).sum()),
            }
        )

        existing_sample, sample_sheet = read_existing_sample_table(cohort_dir, cohort)
        summary, examples = compare_matrices(cohort, rebuilt, existing_sample, "sample_CN.xlsx", sample_sheet)
        summaries.append(summary)
        if not examples.empty:
            mismatch_examples.append(examples)

        existing_binary, binary_sheet = read_existing_binary_table(cohort_dir, cohort)
        rebuilt_binary = (rebuilt > 0).astype(int)
        summary, examples = compare_matrices(cohort, rebuilt_binary, existing_binary, "sample_CN_binary_or_legacy.xlsx", binary_sheet)
        summaries.append(summary)
        if not examples.empty:
            mismatch_examples.append(examples)

    summary_df = pd.DataFrame(summaries)
    metadata_df = pd.DataFrame(metadata_rows)
    mismatch_df = pd.concat(mismatch_examples, ignore_index=True) if mismatch_examples else pd.DataFrame()
    hap_group_df = pd.DataFrame(hap_group_rows)
    return rebuilt_tables, summary_df, metadata_df, mismatch_df, hap_group_df


def sample_level_common_fisher(rebuilt_tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scz = rebuilt_tables["SCZ_filtered_cases"]
    comparison_cohort = rebuilt_tables["comparison_cohort"]
    public_reference = rebuilt_tables["public_reference"]

    genes = sorted(set(scz.columns) | set(comparison_cohort.columns) | set(public_reference.columns))
    scz = scz.reindex(columns=genes, fill_value=0)
    comparison_cohort = comparison_cohort.reindex(columns=genes, fill_value=0)
    public_reference = public_reference.reindex(columns=genes, fill_value=0)
    ctrl = pd.concat([comparison_cohort, public_reference], axis=0)

    scz_binary = (scz > 0).astype(int)
    ctrl_binary = (ctrl > 0).astype(int)
    comparison_cohort_binary = (comparison_cohort > 0).astype(int)
    public_reference_binary = (public_reference > 0).astype(int)

    result = pd.DataFrame({"Gene": genes})
    result["Count_Scz_sample"] = [int(scz_binary[g].sum()) for g in genes]
    result["Frequency_Scz_sample"] = result["Count_Scz_sample"] / scz_binary.shape[0]
    result["Count_comparison_cohort+public_reference_sample"] = [int(ctrl_binary[g].sum()) for g in genes]
    result["Frequency_comparison_cohort+public_reference_sample"] = result["Count_comparison_cohort+public_reference_sample"] / ctrl_binary.shape[0]
    result["Count_comparison_cohort_sample"] = [int(comparison_cohort_binary[g].sum()) for g in genes]
    result["Frequency_comparison_cohort_sample"] = result["Count_comparison_cohort_sample"] / comparison_cohort_binary.shape[0]
    result["Count_public_reference_sample"] = [int(public_reference_binary[g].sum()) for g in genes]
    result["Frequency_public_reference_sample"] = result["Count_public_reference_sample"] / public_reference_binary.shape[0]
    result["SCZ_other_sample"] = scz_binary.shape[0] - result["Count_Scz_sample"]
    result["comparison_cohort+public_reference_other_sample"] = ctrl_binary.shape[0] - result["Count_comparison_cohort+public_reference_sample"]

    common = result[
        (result["Frequency_comparison_cohort+public_reference_sample"] > 0.01)
        & (result["Count_Scz_sample"] >= 2)
    ].copy()

    p_values = []
    odds_ratios = []
    for _, row in common.iterrows():
        contingency = np.array(
            [
                [int(row["Count_Scz_sample"]), int(row["SCZ_other_sample"])],
                [int(row["Count_comparison_cohort+public_reference_sample"]), int(row["comparison_cohort+public_reference_other_sample"])],
            ]
        )
        odds_ratio, p_value = fisher_exact(contingency, alternative="greater")
        odds_ratios.append(float(odds_ratio))
        p_values.append(float(p_value))
    common["P_value"] = p_values
    common["Odds Ratio"] = odds_ratios
    common["FDR"] = multipletests(common["P_value"].to_numpy(dtype=float), method="fdr_bh")[1]
    common["Bonferroni"] = multipletests(common["P_value"].to_numpy(dtype=float), method="bonferroni")[1]
    common["Common_definition"] = "Frequency_comparison_cohort+public_reference_sample > 0.01 and Count_Scz_sample >= 2"
    common = common.sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

    raw_sig = common[common["P_value"] < 0.05].copy()
    fdr_sig = common[common["FDR"] < 0.05].copy()
    summary = pd.DataFrame(
        [
            {"Metric": "SCZ_samples", "Value": scz_binary.shape[0]},
            {"Metric": "comparison_cohort_samples", "Value": comparison_cohort_binary.shape[0]},
            {"Metric": "public_reference_samples", "Value": public_reference_binary.shape[0]},
            {"Metric": "comparison_cohort_public_reference_samples", "Value": ctrl_binary.shape[0]},
            {"Metric": "all_tested_genes_union", "Value": len(result)},
            {"Metric": "common_genes", "Value": len(common)},
            {"Metric": "common_definition", "Value": "Frequency_comparison_cohort+public_reference_sample > 0.01 and Count_Scz_sample >= 2"},
            {"Metric": "raw_p_lt_0.05", "Value": len(raw_sig)},
            {"Metric": "FDR_lt_0.05", "Value": len(fdr_sig)},
            {"Metric": "Bonferroni_lt_0.05", "Value": int((common["Bonferroni"] < 0.05).sum())},
            {"Metric": "min_p_value", "Value": common["P_value"].min() if len(common) else np.nan},
            {"Metric": "min_FDR", "Value": common["FDR"].min() if len(common) else np.nan},
        ]
    )
    return common, raw_sig, summary


def main() -> None:
    rebuilt_tables, qc_summary, hap_parse_summary, mismatch_examples, hap_group_counts = cohort_qc()
    common, raw_sig, fisher_summary = sample_level_common_fisher(rebuilt_tables)
    fdr_sig = common[common["FDR"] < 0.05].copy()
    bonf_sig = common[common["Bonferroni"] < 0.05].copy()

    with pd.ExcelWriter(QC_XLSX, engine="openpyxl") as writer:
        hap_parse_summary.to_excel(writer, sheet_name="haplotype_parse_summary", index=False)
        qc_summary.to_excel(writer, sheet_name="sample_CN_compare_summary", index=False)
        hap_group_counts.to_excel(writer, sheet_name="haplotype_rows_per_sample", index=False)
        if not mismatch_examples.empty:
            mismatch_examples.to_excel(writer, sheet_name="mismatch_examples", index=False)
        fisher_summary.to_excel(writer, sheet_name="sample_fisher_summary", index=False)
        common.to_excel(writer, sheet_name="sample_common_all", index=False)
        raw_sig.to_excel(writer, sheet_name="sample_common_rawP", index=False)
        fdr_sig.to_excel(writer, sheet_name="sample_common_FDR", index=False)

    out_merge = OUTPUT_DIR / f"{RESULT_PREFIX}_tested_genes.xlsx"
    out_raw = OUTPUT_DIR / f"{RESULT_PREFIX}_nominal_p_lt_0.05_all.xlsx"
    out_fdr = OUTPUT_DIR / f"{RESULT_PREFIX}_fdr_lt_0.05.xlsx"
    out_summary = OUTPUT_DIR / f"{RESULT_PREFIX}_summary.xlsx"
    common.to_excel(out_merge, index=False)
    raw_sig.to_excel(out_raw, index=False)
    fdr_sig.to_excel(out_fdr, index=False)
    with pd.ExcelWriter(out_summary, engine="openpyxl") as writer:
        fisher_summary.to_excel(writer, sheet_name="Summary", index=False)
        common.head(50).to_excel(writer, sheet_name="Top50_by_P", index=False)
        raw_sig.to_excel(writer, sheet_name="RawP_lt_0.05", index=False)
        fdr_sig.to_excel(writer, sheet_name="FDR_lt_0.05", index=False)
        bonf_sig.to_excel(writer, sheet_name="Bonf_lt_0.05", index=False)
        qc_summary.to_excel(writer, sheet_name="QC_sample_CN", index=False)

    print("QC summary:")
    print(qc_summary.to_string(index=False))
    print("\nHaplotype parse summary:")
    print(hap_parse_summary.to_string(index=False))
    print("\nSample-level Fisher summary:")
    print(fisher_summary.to_string(index=False))
    print("\nTop sample-level common Fisher:")
    cols = [
        "Gene", "Count_Scz_sample", "Frequency_Scz_sample", "Count_comparison_cohort+public_reference_sample",
        "Frequency_comparison_cohort+public_reference_sample", "P_value", "FDR", "Bonferroni", "Odds Ratio"
    ]
    print(common[cols].head(20).to_string(index=False))
    if len(fdr_sig):
        print("\nFDR significant genes:")
        print(fdr_sig[cols].to_string(index=False))
    else:
        print("\nNo FDR significant genes in sample-level restricted common-gene Fisher test.")
    print(f"\nSaved QC: {QC_XLSX}")
    print(f"Saved: {out_merge}")
    print(f"Saved: {out_raw}")
    print(f"Saved: {out_fdr}")
    print(f"Saved: {out_summary}")


if __name__ == "__main__":
    main()

