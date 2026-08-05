import os
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

CNV_BASE = Path(os.environ.get("EOSCZ_CNV_BASE_DIR", Path(__file__).resolve().parent / "cnv_base"))
OUTDIR = Path(os.environ.get("EOSCZ_CNV_SUPPORT_OUTPUT_DIR", Path(__file__).resolve().parent / "outputs"))
SCZ_DIR = CNV_BASE / "slurm_scripts_case_cohort"
comparison_cohort_DIR = CNV_BASE / "slurm_scripts_comparison_site"
EAS_DIR = CNV_BASE / "slurm_scripts_public_reference_east_asian_subset"
PREFIX = "true_case_cohort_comparison_cohort_public_reference_east_asian_subset_haplotype_corrected"

OUTDIR.mkdir(parents=True, exist_ok=True)

def read_freq_file(path: Path, count_name: str) -> pd.DataFrame:
    df = pd.read_table(path, header=None, names=["Gene", count_name, f"SourceFrequency_{count_name}"])
    df["Gene"] = df["Gene"].astype(str)
    df[count_name] = pd.to_numeric(df[count_name], errors="raise").astype(int)
    return df[["Gene", count_name]]

def infer_total_haplotypes(dataset_dir: Path, preferred: str = "haplotype_CN.xlsx") -> int:
    hap = dataset_dir / preferred
    if hap.exists():
        return int(pd.read_excel(hap, index_col=0).shape[0])
    fallback = dataset_dir / "gene_frequencies_filtered.txt"
    df = pd.read_table(fallback, header=None, names=["Gene", "Count", "Frequency"])
    usable = df[(df["Frequency"] > 0) & df["Count"].notna()]
    totals = (usable["Count"] / usable["Frequency"]).round().astype(int)
    if totals.empty:
        raise RuntimeError(f"Cannot infer total for {dataset_dir}")
    return int(totals.mode().iloc[0])

def build_eas_haplotype_frequency() -> tuple[pd.DataFrame, int, Path]:
    hap_tsv = EAS_DIR / "haplotype_CN.east_asian_subset.tsv"
    hap_xlsx = EAS_DIR / "haplotype_CN.east_asian_subset.xlsx"
    if hap_tsv.exists():
        hap = pd.read_table(hap_tsv, index_col=0)
    elif hap_xlsx.exists():
        hap = pd.read_excel(hap_xlsx, index_col=0)
    else:
        raise FileNotFoundError("Missing haplotype_CN.east_asian_subset.tsv/xlsx")
    hap = hap.apply(pd.to_numeric, errors="coerce").fillna(0)
    total = int(hap.shape[0])
    counts = (hap > 0).sum(axis=0).astype(int)
    freq = counts / total
    out = pd.DataFrame({"Gene": counts.index.astype(str), "Count_public_reference_east_asian_subset": counts.values, "Frequency_public_reference_east_asian_subset": freq.values})
    freq_path = EAS_DIR / "gene_frequencies_filtered.haplotype_east_asian_subset.txt"
    out[["Gene", "Count_public_reference_east_asian_subset", "Frequency_public_reference_east_asian_subset"]].to_csv(freq_path, sep="\t", header=False, index=False, float_format="%.12g")
    return out[["Gene", "Count_public_reference_east_asian_subset"]], total, freq_path

def main():
    eas_freq, eas_total, eas_freq_path = build_eas_haplotype_frequency()
    scz_total = infer_total_haplotypes(SCZ_DIR)
    comparison_cohort_total = infer_total_haplotypes(comparison_cohort_DIR)
    control_total = comparison_cohort_total + eas_total

    scz = read_freq_file(SCZ_DIR / "gene_frequencies_filtered.txt", "Count_Scz")
    comparison_cohort = read_freq_file(comparison_cohort_DIR / "gene_frequencies_filtered.txt", "Count_comparison_cohort")

    # Corrected logic: keep the full SCZ/comparison_cohort gene universe used previously.
    # Missing genes in comparison_cohort or public_reference_east_asian_subset are true zero-count controls, not genes to drop.
    df = scz.merge(comparison_cohort, on="Gene", how="outer").merge(eas_freq, on="Gene", how="left")
    for col in ["Count_Scz", "Count_comparison_cohort", "Count_public_reference_east_asian_subset"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    df["Count_comparison_cohort+public_reference_east_asian_subset"] = df["Count_comparison_cohort"] + df["Count_public_reference_east_asian_subset"]

    df["Frequency_Scz"] = df["Count_Scz"] / scz_total
    df["Frequency_comparison_cohort"] = df["Count_comparison_cohort"] / comparison_cohort_total
    df["Frequency_public_reference_east_asian_subset"] = df["Count_public_reference_east_asian_subset"] / eas_total
    df["Frequency_comparison_cohort+public_reference_east_asian_subset"] = df["Count_comparison_cohort+public_reference_east_asian_subset"] / control_total
    df["SCZ_other"] = scz_total - df["Count_Scz"]
    df["comparison_cohort+public_reference_east_asian_subset_other"] = control_total - df["Count_comparison_cohort+public_reference_east_asian_subset"]
    df["SCZ_frequency_higher"] = df["Frequency_Scz"] > df["Frequency_comparison_cohort+public_reference_east_asian_subset"]
    df["Gene_present_in_public_reference_east_asian_subset_matrix"] = df["Gene"].isin(set(eas_freq["Gene"]))

    pvals = []
    odds = []
    for _, row in df.iterrows():
        a = int(row["Count_Scz"])
        c = int(row["Count_comparison_cohort+public_reference_east_asian_subset"])
        b = scz_total - a
        d = control_total - c
        odds_ratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds.append(float(odds_ratio))
        pvals.append(float(p))
    df["P_value"] = pvals
    _, fdr, _, _ = multipletests(df["P_value"].to_numpy(), method="fdr_bh")
    df["FDR"] = fdr
    df["Bonferroni"] = np.minimum(df["P_value"] * len(df), 1.0)
    df["Odds Ratio"] = odds

    cols = [
        "Gene", "Count_Scz", "Frequency_Scz", "Count_comparison_cohort+public_reference_east_asian_subset", "Frequency_comparison_cohort+public_reference_east_asian_subset",
        "Count_comparison_cohort", "Frequency_comparison_cohort", "Count_public_reference_east_asian_subset", "Frequency_public_reference_east_asian_subset",
        "SCZ_other", "comparison_cohort+public_reference_east_asian_subset_other", "P_value", "FDR", "Bonferroni", "Odds Ratio",
        "SCZ_frequency_higher", "Gene_present_in_public_reference_east_asian_subset_matrix"
    ]
    df = df[cols].sort_values(["P_value", "Gene"], kind="mergesort").reset_index(drop=True)

    merge_path = OUTDIR / f"{PREFIX}_all_tested_genes.xlsx"
    p_all_path = OUTDIR / f"{PREFIX}_one_sided_fisher_nominal_p_lt_0.05_all.xlsx"
    p_scz_path = OUTDIR / f"{PREFIX}_one_sided_fisher_nominal_p_lt_0.05_scz_higher_frequency.xlsx"
    top_path = OUTDIR / f"{PREFIX}_top50_by_pvalue.xlsx"
    summary_path = OUTDIR / "CNV_fisher_SCZ_vs_comparison_cohort_public_reference_east_asian_subset_one_sided_corrected_summary.csv"
    audit_path = OUTDIR / "CNV_fisher_SCZ_vs_comparison_cohort_public_reference_east_asian_subset_one_sided_corrected_audit.json"

    df.to_excel(merge_path, index=False)
    df[df["P_value"] < 0.05].to_excel(p_all_path, index=False)
    df[(df["P_value"] < 0.05) & df["SCZ_frequency_higher"]].to_excel(p_scz_path, index=False)
    df.head(50).to_excel(top_path, index=False)

    summary = pd.DataFrame([
        {"metric": "scz_haplotypes", "value": scz_total},
        {"metric": "comparison_cohort_haplotypes", "value": comparison_cohort_total},
        {"metric": "public_reference_east_asian_subset_haplotypes", "value": eas_total},
        {"metric": "control_haplotypes_comparison_cohort_plus_public_reference_east_asian_subset", "value": control_total},
        {"metric": "tested_genes_outer_SCZ_or_comparison_cohort", "value": len(df)},
        {"metric": "genes_present_in_public_reference_east_asian_subset_matrix", "value": int(df["Gene_present_in_public_reference_east_asian_subset_matrix"].sum())},
        {"metric": "genes_missing_in_public_reference_east_asian_subset_treated_as_zero", "value": int((~df["Gene_present_in_public_reference_east_asian_subset_matrix"]).sum())},
        {"metric": "p_lt_0.05", "value": int((df["P_value"] < 0.05).sum())},
        {"metric": "p_lt_0.05_and_SCZ_frequency_higher", "value": int(((df["P_value"] < 0.05) & df["SCZ_frequency_higher"]).sum())},
        {"metric": "fdr_lt_0.05", "value": int((df["FDR"] < 0.05).sum())},
        {"metric": "bonferroni_lt_0.05", "value": int((df["Bonferroni"] < 0.05).sum())},
    ])
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    audit = {
        "alternative": "greater (SCZ > comparison_cohort+public_reference_east_asian_subset)",
        "corrected_logic": "outer join on SCZ/comparison_cohort gene universe; public_reference_east_asian_subset missing genes are treated as Count_public_reference_east_asian_subset=0",
        "previous_error": "inner join restricted tests to genes present in public_reference_east_asian_subset matrix",
        "frequency_level": "haplotype-level CNV presence: count haplotypes with CN > 0 per gene",
        "generated_public_reference_east_asian_subset_haplotype_frequency": str(eas_freq_path),
        "outputs": [str(merge_path), str(p_all_path), str(p_scz_path), str(top_path), str(summary_path)],
        "totals": {"SCZ": scz_total, "comparison_cohort": comparison_cohort_total, "public_reference_east_asian_subset": eas_total, "comparison_cohort+public_reference_east_asian_subset": control_total},
        "summary": summary.to_dict(orient="records"),
        "top10": df.head(10).replace({np.nan: None}).to_dict(orient="records"),
    }
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    print("DONE_CORRECTED")
    print(f"OUTDIR={OUTDIR}")
    print(f"haplotype_freq={eas_freq_path}")
    print(summary.to_string(index=False))
    print("\nTOP20")
    print(df.head(20).to_string(index=False))

if __name__ == "__main__":
    main()

