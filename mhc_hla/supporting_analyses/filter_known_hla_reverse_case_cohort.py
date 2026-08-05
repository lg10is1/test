# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

import pandas as pd


TAG = "known_reverse_case_cohort_vs_comparison_cohort_public_reference_actual_totals_26-6-18"
GENES = ("A", "B", "C")


def is_known(genotype: object) -> bool:
    return "new" not in str(genotype).lower()


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    all_known_results = []
    all_candidates = []
    summary_rows = []

    for gene in GENES:
        result_path = output_dir / f"{gene}_fisher_results_{TAG}.xlsx"
        df = pd.read_excel(result_path, sheet_name="Sheet1")
        df["Is_Known"] = df["Genotype"].map(is_known)

        known = df[df["Is_Known"]].copy()
        candidates = known[
            (known["HC_Count"] >= 2)
            & (known["SCZ_Count"] == 0)
        ].copy()
        candidates["Selection_Criteria"] = "known; comparison_cohort+public_reference_count>=2; SCZ_count=0"
        candidates = candidates.sort_values(["P-Value", "HC_Count"], ascending=[True, False])

        all_known_results.append(known)
        all_candidates.append(candidates)
        summary_rows.append(
            {
                "Gene": gene,
                "Total_Tested_Genotypes": len(df),
                "Known_Tested_Genotypes": len(known),
                "Known_ControlSpecific_HCge2_SCZ0": len(candidates),
                "Known_ControlSpecific_RawP_lt_0.05": int((candidates["P-Value"] < 0.05).sum()),
                "Known_ControlSpecific_FDR_lt_0.05": int((candidates["FDR_BH"] < 0.05).sum()),
                "Min_Candidate_P": candidates["P-Value"].min() if len(candidates) else pd.NA,
            }
        )

    summary = pd.DataFrame(summary_rows)
    known_all = pd.concat(all_known_results, ignore_index=True)
    candidates_all = pd.concat(all_candidates, ignore_index=True)

    output_xlsx = output_dir / f"known_HLA_reverse_candidates_{TAG}.xlsx"
    output_tsv = output_dir / f"known_HLA_reverse_candidates_{TAG}.tsv"
    summary_tsv = output_dir / f"known_HLA_reverse_summary_{TAG}.tsv"

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, index=False, sheet_name="summary")
        candidates_all.to_excel(writer, index=False, sheet_name="candidates_all")
        for gene in GENES:
            gene_candidates = candidates_all[candidates_all["Gene"] == gene]
            gene_candidates.to_excel(writer, index=False, sheet_name=f"{gene}_candidates")
        known_all.to_excel(writer, index=False, sheet_name="known_all_results")

    candidates_all.to_csv(output_tsv, sep="\t", index=False)
    summary.to_csv(summary_tsv, sep="\t", index=False)

    print(summary.to_string(index=False))
    print(f"Candidates: {len(candidates_all)}")
    print(f"Saved: {output_xlsx}")
    print(f"Saved: {output_tsv}")
    print(f"Saved: {summary_tsv}")


if __name__ == "__main__":
    main()
