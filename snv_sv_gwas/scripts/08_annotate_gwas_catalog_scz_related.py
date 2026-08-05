#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path


DEFAULT_INPUT_BASENAME = "lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.tsv"
DEFAULT_TRAIT_KEYWORDS = [
    "schizophrenia",
    "schizoaffective",
    "psychosis",
    "psychotic",
    "hallucination",
    "hallucinations",
    "delusion",
    "delusions",
    "antipsychotic",
]


def split_input_genes(value: str) -> list[str]:
    genes = []
    for part in (value or "").replace(",", ";").split(";"):
        gene = part.strip()
        if gene and gene != ".":
            genes.append(gene)
    return genes


def split_gwas_genes(value: str) -> list[str]:
    if not value or value in {"NR", "NA"}:
        return []
    parts = re.split(r"\s+-\s+|[,;]", value)
    return [part.strip() for part in parts if part.strip() and part.strip() not in {"NR", "NA"}]


def trait_is_related(row: dict[str, str], keywords: list[str]) -> bool:
    text = " | ".join([row.get("MAPPED_TRAIT", ""), row.get("DISEASE/TRAIT", "")]).lower()
    return any(keyword in text for keyword in keywords)


def default_input_candidates() -> list[Path]:
    return [
        Path("/path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public/tables")
        / DEFAULT_INPUT_BASENAME,
        Path(r"/path/to/local/figure") / DEFAULT_INPUT_BASENAME,
        Path.cwd() / DEFAULT_INPUT_BASENAME,
    ]


def resolve_input_path(value: str | None) -> Path:
    if value:
        return Path(value)
    for candidate in default_input_candidates():
        if candidate.exists():
            return candidate
    candidates = "\n  ".join(str(path) for path in default_input_candidates())
    raise FileNotFoundError(
        "Input TSV was not provided and no default candidate exists.\n"
        "Use --input PATH.\n"
        f"Default candidates:\n  {candidates}"
    )


def gwas_zip_candidates(script_dir: Path, input_path: Path) -> list[Path]:
    candidates = []
    env_path = os.environ.get("GWAS_CATALOG_ASSOC_ZIP")
    if env_path:
        candidates.append(Path(env_path))
    candidates.extend(
        [
            script_dir / "gwas_catalog_v1.0.2-associations_full.zip",
            script_dir / "gwas_catalog_v1.0.2-associations_e115_r2026-06-01_full.zip",
            input_path.parent / "gwas_catalog_v1.0.2-associations_full.zip",
            input_path.parent / "gwas_catalog_v1.0.2-associations_e115_r2026-06-01_full.zip",
        ]
    )
    return candidates


def validate_gwas_zip(path: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        tsv_names = [name for name in archive.namelist() if name.lower().endswith(".tsv")]
        if not tsv_names:
            raise ValueError(f"No TSV file found inside {path}")


def resolve_gwas_zip(value: str | None, script_dir: Path, input_path: Path) -> Path:
    requested = Path(value) if value else None
    if requested is not None:
        if requested.exists():
            validate_gwas_zip(requested)
            return requested
        raise FileNotFoundError(f"GWAS Catalog zip not found: {requested}")

    for candidate in gwas_zip_candidates(script_dir, input_path):
        if candidate.exists():
            validate_gwas_zip(candidate)
            return candidate

    raise FileNotFoundError(
        "A local GWAS Catalog associations zip is required. "
        "Use --gwas-zip PATH or GWAS_CATALOG_ASSOC_ZIP. "
        "This public-release script does not download data."
    )


def output_paths(input_path: Path, out_dir: Path | None) -> tuple[Path, Path, Path]:
    output_dir = out_dir or input_path.parent
    stem = input_path.stem
    return (
        output_dir / f"{stem}.gwas_scz_related.tsv",
        output_dir / f"{stem}.gwas_scz_related_gene_summary.tsv",
        output_dir / f"{stem}.gwas_scz_related_matching_associations.tsv",
    )


def read_input_rows(input_path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if "Gene_refGene" not in fieldnames:
        raise ValueError(f"Input TSV must contain Gene_refGene. Found columns: {', '.join(fieldnames)}")
    return rows, fieldnames


def annotate(args: argparse.Namespace) -> None:
    script_dir = Path(__file__).resolve().parent
    input_path = resolve_input_path(args.input)
    gwas_zip = resolve_gwas_zip(args.gwas_zip, script_dir, input_path)
    out_dir = Path(args.out_dir) if args.out_dir else None
    output, summary, matching_assoc = output_paths(input_path, out_dir)
    output.parent.mkdir(parents=True, exist_ok=True)

    keywords = [keyword.strip().lower() for keyword in args.keywords.split(",") if keyword.strip()]
    rows, input_fieldnames = read_input_rows(input_path)
    input_genes = {
        gene
        for row in rows
        for gene in split_input_genes(row.get("Gene_refGene", ""))
    }

    traits_by_gene: dict[str, set[str]] = defaultdict(set)
    accessions_by_gene: dict[str, set[str]] = defaultdict(set)
    association_count_by_gene: dict[str, int] = defaultdict(int)
    matching_rows: list[dict[str, str]] = []
    total_gwas_rows = 0
    related_gwas_rows = 0

    with zipfile.ZipFile(gwas_zip) as archive:
        tsv_names = [name for name in archive.namelist() if name.lower().endswith(".tsv")]
        if not tsv_names:
            raise ValueError(f"No TSV file found inside {gwas_zip}")
        with archive.open(tsv_names[0]) as raw_handle:
            text_handle = (line.decode("utf-8", "replace") for line in raw_handle)
            reader = csv.DictReader(text_handle, delimiter="\t")
            required_cols = {"MAPPED_GENE", "MAPPED_TRAIT", "DISEASE/TRAIT"}
            missing_cols = required_cols - set(reader.fieldnames or [])
            if missing_cols:
                raise ValueError(f"GWAS TSV missing required columns: {', '.join(sorted(missing_cols))}")
            for gwas_row in reader:
                total_gwas_rows += 1
                if not trait_is_related(gwas_row, keywords):
                    continue
                related_gwas_rows += 1

                mapped_genes = split_gwas_genes(gwas_row.get("MAPPED_GENE", ""))
                matched_genes = sorted(set(mapped_genes) & input_genes)
                if not matched_genes:
                    continue

                mapped_trait = gwas_row.get("MAPPED_TRAIT", "").strip()
                disease_trait = gwas_row.get("DISEASE/TRAIT", "").strip()
                accession = gwas_row.get("STUDY ACCESSION", "").strip()
                trait_label = mapped_trait or disease_trait
                for gene in matched_genes:
                    association_count_by_gene[gene] += 1
                    if trait_label:
                        traits_by_gene[gene].add(trait_label)
                    if accession:
                        accessions_by_gene[gene].add(accession)

                matching_rows.append(
                    {
                        "matched_input_genes": ";".join(matched_genes),
                        "MAPPED_GENE": gwas_row.get("MAPPED_GENE", ""),
                        "MAPPED_TRAIT": mapped_trait,
                        "DISEASE/TRAIT": disease_trait,
                        "MAPPED_TRAIT_URI": gwas_row.get("MAPPED_TRAIT_URI", ""),
                        "STUDY ACCESSION": accession,
                        "PUBMEDID": gwas_row.get("PUBMEDID", ""),
                        "P-VALUE": gwas_row.get("P-VALUE", ""),
                        "SNPS": gwas_row.get("SNPS", ""),
                    }
                )

    added_fields = [
        "GWASCatalog_SCZ_related",
        "GWASCatalog_SCZ_related_genes",
        "GWASCatalog_SCZ_related_traits",
        "GWASCatalog_SCZ_related_association_count",
        "GWASCatalog_SCZ_related_study_accessions",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=input_fieldnames + added_fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            row_genes = split_input_genes(row.get("Gene_refGene", ""))
            matched_genes = sorted(gene for gene in row_genes if gene in traits_by_gene)
            matched_traits = sorted({trait for gene in matched_genes for trait in traits_by_gene[gene]})
            matched_accessions = sorted({acc for gene in matched_genes for acc in accessions_by_gene[gene]})
            row["GWASCatalog_SCZ_related"] = "yes" if matched_genes else "no"
            row["GWASCatalog_SCZ_related_genes"] = ";".join(matched_genes)
            row["GWASCatalog_SCZ_related_traits"] = ";".join(matched_traits)
            row["GWASCatalog_SCZ_related_association_count"] = (
                str(sum(association_count_by_gene[gene] for gene in matched_genes))
                if matched_genes
                else "0"
            )
            row["GWASCatalog_SCZ_related_study_accessions"] = ";".join(matched_accessions)
            writer.writerow(row)

    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "traits", "association_count", "study_accessions"])
        for gene in sorted(traits_by_gene):
            writer.writerow(
                [
                    gene,
                    ";".join(sorted(traits_by_gene[gene])),
                    association_count_by_gene[gene],
                    ";".join(sorted(accessions_by_gene[gene])),
                ]
            )

    with matching_assoc.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "matched_input_genes",
            "MAPPED_GENE",
            "MAPPED_TRAIT",
            "DISEASE/TRAIT",
            "MAPPED_TRAIT_URI",
            "STUDY ACCESSION",
            "PUBMEDID",
            "P-VALUE",
            "SNPS",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(matching_rows)

    matched_input_rows = sum(
        1
        for row in rows
        if any(gene in traits_by_gene for gene in split_input_genes(row.get("Gene_refGene", "")))
    )
    print(f"input={input_path}")
    print(f"gwas_zip={gwas_zip}")
    print(f"input_rows={len(rows)}")
    print(f"input_unique_genes={len(input_genes)}")
    print(f"gwas_rows_scanned={total_gwas_rows}")
    print(f"gwas_related_rows={related_gwas_rows}")
    print(f"matching_association_rows={len(matching_rows)}")
    print(f"matched_genes={len(traits_by_gene)}")
    print(f"matched_input_rows={matched_input_rows}")
    print(f"output={output}")
    print(f"summary={summary}")
    print(f"matching_associations={matching_assoc}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate ANNOVAR Gene_refGene rows with GWAS Catalog schizophrenia-related associations."
    )
    parser.add_argument("--input", help="Input ANNOVAR TSV. Must contain Gene_refGene.")
    parser.add_argument("--gwas-zip", help="GWAS Catalog all-associations v1.0.2 zip.")
    parser.add_argument("--out-dir", help="Output directory. Defaults to the input TSV directory.")
    parser.add_argument(
        "--keywords",
        default=",".join(DEFAULT_TRAIT_KEYWORDS),
        help="Comma-separated lower-case trait keywords used against MAPPED_TRAIT and DISEASE/TRAIT.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    try:
        annotate(parse_args(sys.argv[1:] if argv is None else argv))
    except Exception as exc:  # noqa: BLE001
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
