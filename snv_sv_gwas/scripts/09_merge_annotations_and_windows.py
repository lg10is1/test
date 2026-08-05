#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_TABLE_DIR = Path(
    "/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public/tables"
)
PANGENIE_SOURCES = {"set00"}


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"File does not exist: {path}")
    return pd.read_csv(path, sep="\t", dtype=str)


def first_existing(columns: list[str], df: pd.DataFrame) -> str | None:
    for col in columns:
        if col in df.columns:
            return col
    return None


def clean_chr(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    text = text.replace("chr", "").replace("CHR", "").strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def collapse_values(values: pd.Series) -> str:
    out = []
    seen = set()
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if not text or text == "NA" or text == ".":
            continue
        for part in text.split(";"):
            part = part.strip()
            if part and part not in seen:
                seen.add(part)
                out.append(part)
    return ";".join(out)


def normalize_variant_group(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().lower()
    if text == "sv":
        return "sv"
    if text in {"snv_indel", "snv", "indel", "snv/indel", "snv_indel_lt50bp"}:
        return "snv_indel"
    if "sv" in text and "snv" not in text:
        return "sv"
    if "snv" in text or "indel" in text:
        return "snv_indel"
    return text or "unknown"


def merge_annotations(lead: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    key_cols = ["lead_id", "source_set", "variant_type"]
    for col in key_cols:
        if col not in lead.columns:
            raise ValueError(f"lead table missing required column: {col}")
        if col not in catalog.columns:
            raise ValueError(f"catalog table missing required column: {col}")

    catalog_unique = catalog.drop_duplicates(subset=key_cols, keep="first").copy()

    annotation_cols = [
        "lead_id",
        "source_set",
        "variant_type",
        "START",
        "END",
        "REF",
        "ALT",
        "original_REF_or_A2",
        "original_ALT_or_A1",
        "Func_refGene",
        "Gene_refGene",
        "GeneDetail_refGene",
        "ExonicFunc_refGene",
        "AAChange_refGene",
        "GWASCatalog_SCZ_related",
        "GWASCatalog_SCZ_related_genes",
        "GWASCatalog_SCZ_related_traits",
        "GWASCatalog_SCZ_related_association_count",
        "GWASCatalog_SCZ_related_study_accessions",
    ]
    annotation_cols = [col for col in annotation_cols if col in catalog_unique.columns]
    catalog_for_merge = catalog_unique[annotation_cols].copy()

    rename_map = {}
    for col in annotation_cols:
        if col in key_cols:
            continue
        if col in lead.columns:
            rename_map[col] = f"anno_{col}"

    if rename_map:
        catalog_for_merge = catalog_for_merge.rename(columns=rename_map)

    merged = lead.merge(
        catalog_for_merge,
        on=key_cols,
        how="left",
        validate="many_to_one",
    )

    chr_col = first_existing(["lead_chr", "CHR", "#CHROM", "CHROM"], merged)
    pos_col = first_existing(["lead_pos", "POS", "BP"], merged)
    p_col = first_existing(["lead_p", "P", "P_noSPA"], merged)

    if chr_col is None or pos_col is None or p_col is None:
        raise ValueError(
            "Cannot resolve lead chr/pos/p columns. "
            f"chr_col={chr_col}, pos_col={pos_col}, p_col={p_col}"
        )

    merged["final_chr"] = merged[chr_col].map(clean_chr)
    merged["final_pos"] = to_numeric(merged[pos_col])
    merged["final_p"] = to_numeric(merged[p_col])

    preferred_cols = [
        "source_set",
        "variant_type",
        "lead_id",
        "final_chr",
        "final_pos",
        "final_p",
        "lead_chr",
        "lead_pos",
        "lead_p",
        "found_in_gwas",
        "lead_source",
        "lead_source_file",
        "gwas_file",
        "SNP",
        "CHR",
        "POS",
        "A1",
        "A2",
        "REF",
        "ALT",
        "anno_REF",
        "anno_ALT",
        "original_REF_or_A2",
        "original_ALT_or_A1",
        "START",
        "END",
        "anno_START",
        "anno_END",
        "Func_refGene",
        "Gene_refGene",
        "GeneDetail_refGene",
        "ExonicFunc_refGene",
        "AAChange_refGene",
        "GWASCatalog_SCZ_related",
        "GWASCatalog_SCZ_related_genes",
        "GWASCatalog_SCZ_related_traits",
        "GWASCatalog_SCZ_related_association_count",
        "GWASCatalog_SCZ_related_study_accessions",
        "BETA",
        "SE",
        "T",
        "OR",
        "Z_STAT",
        "A1_FREQ",
        "AF1",
        "N",
        "OBS_CT",
        "P",
        "P_noSPA",
    ]

    preferred_cols = [col for col in preferred_cols if col in merged.columns]
    other_cols = [
        col
        for col in merged.columns
        if col not in preferred_cols
        and not col.endswith(".x")
        and not col.endswith(".y")
        and not col.endswith("_x")
        and not col.endswith("_y")
        and col not in {"duplicate_index_in_gwas"}
    ]

    return merged[preferred_cols + other_cols]


def make_window_table(final_df: pd.DataFrame, window_kb: int) -> pd.DataFrame:
    window_bp = int(window_kb) * 1000
    df = final_df.copy()
    df = df[df["final_chr"].notna() & df["final_pos"].notna() & df["final_p"].notna()].copy()
    df = df[df["final_chr"].astype(str) != ""].copy()

    if df.empty:
        return pd.DataFrame()

    if "variant_type" not in df.columns:
        raise ValueError("final table missing required column for separated windows: variant_type")

    df["window_variant_group"] = df["variant_type"].map(normalize_variant_group)
    df["window_start"] = ((df["final_pos"].astype(int) - 1) // window_bp) * window_bp + 1
    df["window_end"] = df["window_start"] + window_bp - 1
    df["window_region_id"] = (
        "chr"
        + df["final_chr"].astype(str)
        + ":"
        + df["window_start"].astype(str)
        + "-"
        + df["window_end"].astype(str)
    )
    if "source_set" not in df.columns or "lead_id" not in df.columns:
        raise ValueError("final table missing source_set or lead_id for canonical window rules")
    is_paragraph_sv = (
        (df["window_variant_group"] == "sv")
        & (df["source_set"].astype(str) == "paragraph")
    )
    unexpected_sv_sources = set(
        df.loc[
            (df["window_variant_group"] == "sv")
            & ~df["source_set"].astype(str).isin(PANGENIE_SOURCES | {"paragraph"}),
            "source_set",
        ].astype(str)
    )
    if unexpected_sv_sources:
        raise ValueError(
            "Unexpected SV source(s) for window policy: "
            + ", ".join(sorted(unexpected_sv_sources))
        )
    df["window_merge_scope"] = "cross_source_fixed_window"
    df.loc[is_paragraph_sv, "window_merge_scope"] = "paragraph_clumped_lead_independent"
    df["window_id"] = df["window_region_id"] + "|" + df["window_variant_group"]
    df.loc[is_paragraph_sv, "window_id"] = (
        df.loc[is_paragraph_sv, "window_id"]
        + "|paragraph_clumped|"
        + df.loc[is_paragraph_sv, "lead_id"].astype(str)
    )

    rows = []
    for (chrom, window_start, window_end, region_id, variant_group, merge_scope, window_id), sub in df.groupby(
        [
            "final_chr",
            "window_start",
            "window_end",
            "window_region_id",
            "window_variant_group",
            "window_merge_scope",
            "window_id",
        ],
        dropna=False,
        sort=False,
    ):
        # Representative row uses the smallest P value, i.e. the strongest association.
        idx = sub["final_p"].astype(float).idxmin()
        rep = sub.loc[idx].copy()
        rep["window_chr"] = chrom
        rep["window_start"] = window_start
        rep["window_end"] = window_end
        rep["window_region_id"] = region_id
        rep["window_variant_group"] = variant_group
        rep["window_merge_scope"] = merge_scope
        rep["window_id"] = window_id
        rep["window_n_signals"] = len(sub)
        rep["window_lead_ids"] = collapse_values(sub["lead_id"])
        rep["window_source_sets"] = collapse_values(sub["source_set"])
        rep["window_variant_types"] = collapse_values(sub["variant_type"])
        rep["window_genes"] = collapse_values(sub.get("Gene_refGene", pd.Series(dtype=str)))
        rep["window_catalog_scz_related"] = (
            "yes"
            if any(str(x).lower() == "yes" for x in sub.get("GWASCatalog_SCZ_related", []))
            else "no"
        )
        rep["window_catalog_traits"] = collapse_values(
            sub.get("GWASCatalog_SCZ_related_traits", pd.Series(dtype=str))
        )
        rep["representative_rule"] = "min_final_p_most_significant"
        rows.append(rep)

    out = pd.DataFrame(rows)
    front = [
        "window_id",
        "window_chr",
        "window_start",
        "window_end",
        "window_region_id",
        "window_variant_group",
        "window_merge_scope",
        "window_n_signals",
        "window_lead_ids",
        "window_source_sets",
        "window_variant_types",
        "window_genes",
        "window_catalog_scz_related",
        "window_catalog_traits",
        "representative_rule",
    ]
    front = [col for col in front if col in out.columns]
    rest = [col for col in out.columns if col not in front]
    return out[front + rest]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge GWAS lead, ANNOVAR/refGene, GWAS Catalog annotations, and make 1Mb window summary."
    )
    parser.add_argument("--table-dir", default=str(DEFAULT_TABLE_DIR))
    parser.add_argument("--window-kb", type=int, default=1000)
    parser.add_argument("--lead-table")
    parser.add_argument("--catalog-table")
    parser.add_argument("--out-final")
    parser.add_argument("--out-window")
    parser.add_argument("--out-window-snv-indel")
    parser.add_argument("--out-window-sv")
    parser.add_argument(
        "--expected-sv", type=int, default=-1,
        help="Optional batch-specific check; negative disables the check (default).",
    )
    parser.add_argument(
        "--expected-snv-indel", type=int, default=-1,
        help="Optional batch-specific check; negative disables the check (default).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    table_dir = Path(args.table_dir)

    lead_table = Path(args.lead_table) if args.lead_table else table_dir / "lead_sig_from_gwas.canonical_1000kb.tsv"
    catalog_table = (
        Path(args.catalog_table)
        if args.catalog_table
        else table_dir / "lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.gwas_scz_related.tsv"
    )
    out_final = (
        Path(args.out_final)
        if args.out_final
        else table_dir / "lead_sig_from_gwas.final_merged.cleaned.tsv"
    )
    out_window = (
        Path(args.out_window)
        if args.out_window
        else table_dir / f"lead_sig_from_gwas.final_merged.cleaned.window_{args.window_kb}kb.tsv"
    )
    out_window_snv_indel = (
        Path(args.out_window_snv_indel)
        if args.out_window_snv_indel
        else table_dir / f"lead_sig_from_gwas.final_merged.cleaned.window_{args.window_kb}kb.snv_indel.tsv"
    )
    out_window_sv = (
        Path(args.out_window_sv)
        if args.out_window_sv
        else table_dir / f"lead_sig_from_gwas.final_merged.cleaned.window_{args.window_kb}kb.sv.tsv"
    )

    lead = read_table(lead_table)
    catalog = read_table(catalog_table)

    final_df = merge_annotations(lead, catalog)
    window_df = make_window_table(final_df, args.window_kb)

    out_final.parent.mkdir(parents=True, exist_ok=True)
    out_window.parent.mkdir(parents=True, exist_ok=True)
    out_window_snv_indel.parent.mkdir(parents=True, exist_ok=True)
    out_window_sv.parent.mkdir(parents=True, exist_ok=True)
    if "window_variant_group" in window_df.columns:
        window_snv_indel = window_df[window_df["window_variant_group"] == "snv_indel"]
        window_sv = window_df[window_df["window_variant_group"] == "sv"]
    else:
        window_snv_indel = window_df
        window_sv = window_df
    if args.expected_sv >= 0 and len(window_sv) != args.expected_sv:
        raise ValueError(f"Expected {args.expected_sv} SV windows/leads, observed {len(window_sv)}")
    if args.expected_snv_indel >= 0 and len(window_snv_indel) != args.expected_snv_indel:
        raise ValueError(
            f"Expected {args.expected_snv_indel} SNV/indel windows, observed {len(window_snv_indel)}"
        )
    final_df.to_csv(out_final, sep="\t", index=False)
    window_df.to_csv(out_window, sep="\t", index=False)
    window_snv_indel.to_csv(
        out_window_snv_indel,
        sep="\t",
        index=False,
    )
    window_sv.to_csv(
        out_window_sv,
        sep="\t",
        index=False,
    )

    print(f"lead_table={lead_table}")
    print(f"catalog_table={catalog_table}")
    print(f"final_rows={len(final_df)}")
    print(f"window_rows={len(window_df)}")
    print(f"out_final={out_final}")
    print(f"out_window={out_window}")
    print(f"out_window_snv_indel={out_window_snv_indel}")
    print(f"out_window_sv={out_window_sv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
