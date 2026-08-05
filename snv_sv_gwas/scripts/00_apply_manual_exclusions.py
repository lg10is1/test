#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".clumped":
        return pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=False, engine="python")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_table(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def candidate_id_columns(df: pd.DataFrame) -> list[str]:
    cols: list[str] = []
    exact = {
        "lead_id",
        "SNP",
        "ID",
        "MarkerName",
        "MARKERNAME",
        "variant_id",
        "VARIANT_ID",
        "unified_SNP",
    }
    for col in df.columns:
        if col in exact or col.startswith("SNP_") or col.startswith("ID_"):
            cols.append(col)
    return cols


def file_context(path: Path) -> tuple[str | None, str | None]:
    text = path.as_posix().lower()
    source = next((s for s in ("set00", "set01", "set02", "paragraph", "deepvariant") if s in text), None)
    variant_type = None
    if ".sv" in text or "/sv" in text or "sv_" in text or "sv." in text:
        variant_type = "sv"
    if "snv_indel" in text or "snv" in text:
        variant_type = "snv_indel"
    return source, variant_type


def row_mask_for_exclusion(df: pd.DataFrame, path: Path, exclusion: pd.Series) -> pd.Series:
    source = str(exclusion["source_set"]).strip()
    variant_type = str(exclusion["variant_type"]).strip().lower()
    lead_id = str(exclusion["lead_id"]).strip()
    if not lead_id:
        return pd.Series(False, index=df.index)

    id_cols = candidate_id_columns(df)
    if not id_cols:
        return pd.Series(False, index=df.index)

    mask = pd.Series(False, index=df.index)
    for col in id_cols:
        col_mask = df[col].astype(str) == lead_id
        if col.endswith(f"_{source}") or col in {"SNP", "ID", "lead_id", "MarkerName", "MARKERNAME"}:
            mask |= col_mask
        elif col_mask.any():
            mask |= col_mask

    if not mask.any():
        return mask

    file_source, file_type = file_context(path)
    if "source_set" in df.columns:
        mask &= df["source_set"].astype(str).str.lower() == source.lower()
    elif "set" in df.columns:
        mask &= df["set"].astype(str).str.lower() == source.lower()
    elif file_source is not None:
        mask &= file_source == source.lower()

    if "variant_type" in df.columns:
        mask &= df["variant_type"].astype(str).str.lower().isin({variant_type, variant_type.upper()})
    elif file_type is not None:
        mask &= file_type == variant_type

    return mask


def filter_file(path: Path, exclusions: pd.DataFrame) -> dict[str, object] | None:
    try:
        df = read_table(path)
    except Exception as exc:
        return {"file": str(path), "status": "READ_ERROR", "rows_before": "NA", "rows_after": "NA", "removed": 0, "note": str(exc)}

    if df.empty:
        return {"file": str(path), "status": "EMPTY", "rows_before": 0, "rows_after": 0, "removed": 0, "note": ""}

    remove = pd.Series(False, index=df.index)
    for _, exclusion in exclusions.iterrows():
        remove |= row_mask_for_exclusion(df, path, exclusion)

    n_remove = int(remove.sum())
    if n_remove:
        out = df.loc[~remove].copy()
        write_table(out, path)
        return {
            "file": str(path),
            "status": "FILTERED",
            "rows_before": len(df),
            "rows_after": len(out),
            "removed": n_remove,
            "note": "",
        }
    return {"file": str(path), "status": "UNCHANGED", "rows_before": len(df), "rows_after": len(df), "removed": 0, "note": ""}


def refresh_clumping_summary(gwas_dir: Path) -> None:
    clump_dir = gwas_dir / "clumping_by_set_subtype"
    summary_file = clump_dir / "clumping_summary.by_set_subtype.tsv"
    if not summary_file.exists():
        return
    summary = read_table(summary_file)
    if not {"set", "subtype"}.issubset(summary.columns):
        return
    for idx, row in summary.iterrows():
        set_name = str(row["set"])
        subtype = str(row["subtype"])
        prefix = f"{set_name}.{subtype}.clump_p1_5e-06.r2_0.01.kb_1000"
        independent = clump_dir / f"{prefix}.independent_lead_signals.tsv"
        clumped = clump_dir / f"{prefix}.clumped"
        lead_n = None
        if independent.exists():
            lead_n = max(len(read_table(independent)), 0)
        elif clumped.exists():
            lead_n = max(len(read_table(clumped)), 0)
        if lead_n is not None and "lead_n" in summary.columns:
            summary.loc[idx, "lead_n"] = str(lead_n)
    write_table(summary, summary_file)


def main() -> int:
    parser = argparse.ArgumentParser(description="Remove manually excluded leads from reused public GWAS/clumping outputs.")
    parser.add_argument("--root", required=True, help="GWAS output directory to filter in place.")
    parser.add_argument("--exclusions", required=True, help="TSV with source_set, variant_type, lead_id.")
    parser.add_argument("--audit", help="Audit output TSV.")
    args = parser.parse_args()

    root = Path(args.root)
    exclusions = read_table(Path(args.exclusions))
    required = {"source_set", "variant_type", "lead_id"}
    missing = required.difference(exclusions.columns)
    if missing:
        raise ValueError(f"Exclusion table missing columns: {', '.join(sorted(missing))}")

    records = []
    paths = sorted(set(root.rglob("*.tsv")) | set(root.rglob("*.clumped")))
    for path in paths:
        records.append(filter_file(path, exclusions))
    refresh_clumping_summary(root)

    audit = pd.DataFrame([r for r in records if r is not None])
    audit_file = Path(args.audit) if args.audit else root / "manual_exclusion_filter.audit.tsv"
    audit_file.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(audit_file, sep="\t", index=False)
    print(f"[DONE] Manual exclusion audit: {audit_file}")
    print(audit.groupby("status")["file"].count().to_string() if not audit.empty else "No TSV files found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
