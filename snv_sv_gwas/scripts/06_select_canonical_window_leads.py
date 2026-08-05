#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_TABLE_DIR = Path(
    "/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public/tables"
)
PANGENIE_SOURCES = {"set00"}
SCRIPT_DIR = Path(__file__).resolve().parent


def clean_chr(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.endswith(".0"):
        text = text[:-2]
    return text


def variant_group(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().lower()
    if text == "sv":
        return "sv"
    if text in {"snv_indel", "snv", "indel", "snv/indel", "sgv"}:
        return "snv_indel"
    raise ValueError(f"Unsupported variant_type for canonical selection: {value!r}")


def collapse(values: pd.Series) -> str:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = "" if pd.isna(value) else str(value).strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return ";".join(out)


def row_key(source_set: object, variant_type: object, lead_id: object) -> str:
    return "||".join(
        [
            str(source_set).strip().lower(),
            variant_group(variant_type),
            str(lead_id).strip(),
        ]
    )


def ensure_lead_columns(df: pd.DataFrame, label: str) -> pd.DataFrame:
    df = df.copy()
    if "lead_chr" not in df.columns:
        for col in ("CHR", "#CHROM", "CHROM"):
            if col in df.columns:
                df["lead_chr"] = df[col]
                break
    if "lead_pos" not in df.columns:
        for col in ("POS", "BP"):
            if col in df.columns:
                df["lead_pos"] = df[col]
                break
    if "lead_p" not in df.columns and "P" in df.columns:
        df["lead_p"] = df["P"]

    required = {"source_set", "variant_type", "lead_id", "lead_chr", "lead_pos", "lead_p"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{label} table missing required columns: {', '.join(missing)}")
    return df


def add_canonical_columns(df: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    df = df.copy()
    if "found_in_gwas" in df.columns:
        found = df["found_in_gwas"].astype(str).str.lower().isin({"true", "1", "t"})
        df = df[found].copy()

    df["source_set"] = df["source_set"].astype(str).str.strip()
    df["lead_id"] = df["lead_id"].astype(str).str.strip()
    df["canonical_chr"] = df["lead_chr"].map(clean_chr)
    df["canonical_pos"] = pd.to_numeric(df["lead_pos"], errors="coerce")
    df["canonical_p"] = pd.to_numeric(df["lead_p"], errors="coerce")
    df["canonical_variant_group"] = df["variant_type"].map(variant_group)
    df = df[
        df["canonical_chr"].isin([str(x) for x in range(1, 23)])
        & df["canonical_pos"].notna()
        & df["canonical_p"].notna()
        & (df["lead_id"] != "")
    ].copy()

    df["canonical_window_start"] = (
        ((df["canonical_pos"].astype(int) - 1) // window_bp) * window_bp + 1
    )
    df["canonical_window_end"] = df["canonical_window_start"] + window_bp - 1
    df["canonical_region_id"] = (
        "chr"
        + df["canonical_chr"]
        + ":"
        + df["canonical_window_start"].astype(str)
        + "-"
        + df["canonical_window_end"].astype(str)
    )

    def make_group_id(row: pd.Series) -> str:
        group = row["canonical_variant_group"]
        source = str(row["source_set"])
        if group == "sv" and source == "paragraph":
            return f"sv|paragraph_clumped|{row['lead_id']}"
        if group == "sv" and source not in PANGENIE_SOURCES:
            raise ValueError(f"Unexpected non-Paragraph SV source: {source}")
        return f"{group}|{row['canonical_region_id']}"

    df["canonical_group_id"] = df.apply(make_group_id, axis=1)
    df["manual_exact_key"] = [
        row_key(source, group, lead)
        for source, group, lead in zip(
            df["source_set"], df["canonical_variant_group"], df["lead_id"]
        )
    ]
    return df


def read_removed_rows(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    out = pd.read_csv(path, sep="\t", dtype=str)
    if out.empty:
        return out
    if "manual_exclusion_key" in out.columns:
        out = out.drop(columns=["manual_exclusion_key"])
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select canonical 1 Mb GWAS leads. Manual exclusions are exact-row "
            "exclusions by source_set + variant_type + lead_id. If an excluded row "
            "would have been the canonical representative, the group is dropped and "
            "no replacement representative is selected."
        )
    )
    parser.add_argument("--table-dir", default=str(DEFAULT_TABLE_DIR))
    parser.add_argument("--input")
    parser.add_argument("--removed-input")
    parser.add_argument("--output")
    parser.add_argument("--mapping-output")
    parser.add_argument("--summary-output")
    parser.add_argument(
        "--manual-exclusions",
        default=str(SCRIPT_DIR / "manual_excluded_leads.tsv"),
        help="TSV with exact source_set, variant_type, lead_id exclusions.",
    )
    parser.add_argument("--manual-audit-output")
    parser.add_argument("--window-kb", type=int, default=1000)
    parser.add_argument("--expected-sv", type=int, default=-1)
    parser.add_argument("--expected-snv-indel", type=int, default=-1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    table_dir = Path(args.table_dir)
    input_file = Path(args.input) if args.input else table_dir / "lead_sig_from_gwas.all.tsv"
    removed_file = (
        Path(args.removed_input)
        if args.removed_input
        else table_dir / "lead_sig_from_gwas.manual_exclusion.audit.removed.tsv"
    )
    output_file = (
        Path(args.output)
        if args.output
        else table_dir / f"lead_sig_from_gwas.canonical_{args.window_kb}kb.tsv"
    )
    mapping_file = (
        Path(args.mapping_output)
        if args.mapping_output
        else table_dir / f"lead_sig_from_gwas.canonical_{args.window_kb}kb.mapping.tsv"
    )
    summary_file = (
        Path(args.summary_output)
        if args.summary_output
        else table_dir / f"lead_sig_from_gwas.canonical_{args.window_kb}kb.summary.tsv"
    )
    manual_audit_file = (
        Path(args.manual_audit_output)
        if args.manual_audit_output
        else table_dir / "lead_sig_from_gwas.manual_excluded_exact_representatives.tsv"
    )

    df_current_raw = ensure_lead_columns(
        pd.read_csv(input_file, sep="\t", dtype=str),
        label="Lead",
    )

    window_bp = int(args.window_kb) * 1000
    current = add_canonical_columns(df_current_raw, window_bp)
    current["manual_exact_excluded"] = False

    removed_raw = read_removed_rows(removed_file)
    if not removed_raw.empty:
        removed_raw = ensure_lead_columns(removed_raw, label="Removed-row")
        removed = add_canonical_columns(removed_raw, window_bp)
        removed["manual_exact_excluded"] = True
    else:
        removed = current.iloc[0:0].copy()
        removed["manual_exact_excluded"] = pd.Series(dtype=bool)

    original = pd.concat([current, removed], ignore_index=True, sort=False)
    original = original.drop_duplicates(
        subset=["source_set", "canonical_variant_group", "lead_id", "canonical_group_id"],
        keep="first",
    )
    original = original.sort_values(
        ["canonical_p", "source_set", "lead_id"], kind="mergesort"
    ).copy()

    representatives: list[pd.Series] = []
    mapping_parts: list[pd.DataFrame] = []
    for _, sub in original.groupby("canonical_group_id", sort=False, dropna=False):
        sub = sub.copy()
        rep = sub.iloc[0].copy()
        rep["canonical_n_source_signals"] = len(sub)
        rep["canonical_member_lead_ids"] = collapse(sub["lead_id"])
        rep["canonical_member_source_sets"] = collapse(sub["source_set"])
        rep["canonical_representative_rule"] = "minimum_lead_p"

        rep_excluded = bool(rep.get("manual_exact_excluded", False))
        rep["manual_excluded_representative"] = rep_excluded
        if not rep_excluded:
            representatives.append(rep)

        sub["canonical_representative_source_set"] = rep["source_set"]
        sub["canonical_representative_lead_id"] = rep["lead_id"]
        sub["canonical_representative_manual_exact_excluded"] = rep_excluded
        sub["canonical_selected"] = (
            (~sub["manual_exact_excluded"].astype(bool))
            & (sub["source_set"] == rep["source_set"])
            & (sub["lead_id"] == rep["lead_id"])
            & (not rep_excluded)
        )
        mapping_parts.append(sub)

    canonical = pd.DataFrame(representatives)
    mapping = (
        pd.concat(mapping_parts, ignore_index=True, sort=False)
        if mapping_parts
        else original.iloc[0:0].copy()
    )

    if canonical.empty:
        raise ValueError("No canonical leads remain after exact manual exclusions.")

    canonical = canonical.sort_values(
        ["canonical_variant_group", "canonical_chr", "canonical_pos", "source_set", "lead_id"],
        kind="mergesort",
    )

    manual_audit = mapping[mapping["manual_exact_excluded"].astype(bool)].copy()
    if manual_audit.empty:
        manual_audit = pd.DataFrame(
            columns=[
                "source_set",
                "variant_type",
                "lead_id",
                "canonical_group_id",
                "canonical_representative_source_set",
                "canonical_representative_lead_id",
                "canonical_representative_manual_exact_excluded",
            ]
        )

    counts = canonical.groupby("canonical_variant_group").size().to_dict()
    n_sv = int(counts.get("sv", 0))
    n_snv = int(counts.get("snv_indel", 0))
    if args.expected_sv >= 0 and n_sv != args.expected_sv:
        raise ValueError(f"Expected {args.expected_sv} canonical SV leads, observed {n_sv}")
    if args.expected_snv_indel >= 0 and n_snv != args.expected_snv_indel:
        raise ValueError(
            f"Expected {args.expected_snv_indel} canonical SNV/indel leads, observed {n_snv}"
        )

    n_manual_rows = int(mapping["manual_exact_excluded"].astype(bool).sum())
    n_manual_reps = int(mapping["canonical_representative_manual_exact_excluded"].astype(bool).sum())
    n_dropped_groups = int(
        mapping.loc[
            mapping["canonical_representative_manual_exact_excluded"].astype(bool),
            "canonical_group_id",
        ].nunique()
    )

    summary = pd.DataFrame(
        [
            {"metric": "input_source_level_leads_after_exact_manual_exclusion", "value": len(current)},
            {"metric": "manual_exact_excluded_input_rows", "value": n_manual_rows},
            {"metric": "manual_exact_excluded_canonical_groups_no_replacement", "value": n_dropped_groups},
            {"metric": "canonical_total_leads", "value": len(canonical)},
            {"metric": "canonical_sv_leads", "value": n_sv},
            {"metric": "canonical_snv_indel_leads", "value": n_snv},
            {
                "metric": "sv_rule",
                "value": "Pangenie set00-set02 merged by 1Mb window; Paragraph clumped leads retained independently",
            },
            {
                "metric": "snv_indel_rule",
                "value": "established cross-source 1Mb window merge; minimum P representative",
            },
            {
                "metric": "manual_exclusion_rule",
                "value": "exact source_set + variant_type + lead_id rows are removed; if the removed row was the original canonical representative, no replacement representative is selected",
            },
        ]
    )

    for path in (output_file, mapping_file, summary_file, manual_audit_file):
        path.parent.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(output_file, sep="\t", index=False)
    mapping.to_csv(mapping_file, sep="\t", index=False)
    summary.to_csv(summary_file, sep="\t", index=False)
    manual_audit.to_csv(manual_audit_file, sep="\t", index=False)

    print(f"input={input_file}")
    print(f"removed_rows={removed_file}")
    print(f"canonical={output_file}")
    print(f"mapping={mapping_file}")
    print(f"summary={summary_file}")
    print(f"manual_exact_audit={manual_audit_file}")
    print(f"canonical_total={len(canonical)} canonical_sv={n_sv} canonical_snv_indel={n_snv}")
    print(f"manual_exact_rows={n_manual_rows}; dropped_canonical_groups={n_dropped_groups}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
