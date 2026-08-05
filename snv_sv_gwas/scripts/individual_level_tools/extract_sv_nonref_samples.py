#!/usr/bin/env python3
# PRIVACY WARNING: This optional tool writes individual sample identifiers.
# Use it only on authorized controlled data and never commit its outputs.
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path("/path/to/local/project")
DEFAULT_INPUT = (
    ROOT
    / "Figure3_manual/summary_maintext_public/01_GWAS_leads"
    / "lead_sig_from_gwas.final_merged.cleaned.noloc.compact.tsv"
)
DEFAULT_OUTDIR = ROOT / "Figure3_manual/summary_maintext_public/01_GWAS_leads"
LRS_PANGENIE_VCF = (
    ROOT
    / "TGS_SV_merge_SCZ/truvari_single_sample"
    / "truvari_merged_sort_pP0.5.sv_len_gt50.sorted.vcf.gz"
)
LRS_PARAGRAPH_VCF = (
    ROOT
    / "TGS_SV_merge_SCZ/truvari_single_sample/merge_2caller_pav"
    / "merge_high_conf.2share.vcf.gz"
)
PANGENIE_GWAS = ROOT / "TGS_callset/Pangenie_v3/06.gwas"
NGS_BFILES = {
    # Public example: set00 only. Add set01/set02 if available.
    "set00": PANGENIE_GWAS / "set00/NGS.QCsite.QCind",
    "paragraph": (
        ROOT
        / "GWAS/Deepvariant_paragraph"
        / "chr_all2.strict_step2_genimi.common_samples.merged"
    ),
}
PANGENIE_SOURCES = {"set00"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Add NGS and LRS heterozygous/homozygous-alternate carrier sample IDs "
            "to the canonical compact SV table."
        )
    )
    p.add_argument("--input", default=str(DEFAULT_INPUT))
    p.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    p.add_argument("--plink", default="plink")
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--lrs-pangenie-vcf", default=str(LRS_PANGENIE_VCF))
    p.add_argument("--lrs-paragraph-vcf", default=str(LRS_PARAGRAPH_VCF))
    p.add_argument("--window-bp", type=int, default=1000)
    return p.parse_args()


def run(cmd: list[str], *, capture: bool = True) -> str:
    result = subprocess.run(
        cmd,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )
    return result.stdout if capture else ""


def check_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def check_bfile(prefix: Path) -> None:
    missing = [str(prefix) + ext for ext in (".bed", ".bim", ".fam") if not Path(str(prefix) + ext).is_file()]
    if missing:
        raise FileNotFoundError("Missing bfile component(s): " + ", ".join(missing))


def clean_chr(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    return text[:-2] if text.endswith(".0") else text


def split_alt(value: str) -> list[str]:
    return [x.upper() for x in str(value).split(",") if x]


def sample_list(values: list[str]) -> str:
    return ";".join(sorted(set(str(x) for x in values if str(x))))


def get_vcf_samples(vcf: Path, bcftools: str) -> list[str]:
    samples = [x for x in run([bcftools, "query", "-l", str(vcf)]).splitlines() if x]
    if not samples:
        raise RuntimeError(f"No samples in VCF: {vcf}")
    return samples


def get_vcf_contigs(vcf: Path, bcftools: str) -> set[str]:
    header = run([bcftools, "view", "-h", str(vcf)])
    return set(re.findall(r"^##contig=<ID=([^,>]+)", header, flags=re.MULTILINE))


def vcf_chr(chrom: str, contigs: set[str]) -> str:
    plain = clean_chr(chrom)
    for candidate in (f"chr{plain}", plain):
        if candidate in contigs:
            return candidate
    raise ValueError(f"Cannot map chromosome {chrom!r} to VCF contigs")


def query_region(vcf: Path, chrom: str, start: int, end: int, bcftools: str) -> list[list[str]]:
    region = f"{chrom}:{max(1, start)}-{max(1, end)}"
    selector = "-r" if Path(str(vcf) + ".tbi").is_file() or Path(str(vcf) + ".csi").is_file() else "-t"
    text = run([bcftools, "view", "-H", selector, region, str(vcf)])
    return [line.split("\t") for line in text.splitlines() if line]


def parse_record(fields: list[str], samples: list[str]) -> dict[str, object]:
    if len(fields) < 9:
        raise ValueError("VCF record has fewer than 9 columns")
    fmt = fields[8].split(":")
    if "GT" not in fmt:
        raise ValueError(f"VCF record has no GT FORMAT: {fields[2]}")
    gt_index = fmt.index("GT")
    sample_fields = fields[9:]
    if len(sample_fields) != len(samples):
        raise ValueError(
            f"VCF sample count mismatch for {fields[2]}: fields={len(sample_fields)}, header={len(samples)}"
        )
    genotypes: dict[str, str] = {}
    for sample, value in zip(samples, sample_fields):
        parts = value.split(":")
        genotypes[sample] = parts[gt_index] if gt_index < len(parts) else "."
    info = fields[7]
    return {
        "chrom": fields[0],
        "pos": int(fields[1]),
        "id": fields[2],
        "ref": fields[3],
        "alt": fields[4],
        "info": info,
        "genotypes": genotypes,
    }


def classify_gt(gt: str) -> tuple[str, int]:
    """Return carrier class and number of missing alleles imputed as reference."""
    gt = str(gt).split(":", 1)[0]
    if gt in {"", "."}:
        return "missing", 0
    alleles = re.split(r"[/|]", gt)
    if not alleles:
        return "missing", 0
    n_imputed = sum(a == "." for a in alleles)
    alleles = ["0" if a == "." else a for a in alleles]
    if any(not re.fullmatch(r"\d+", a) for a in alleles):
        return "missing", n_imputed
    nonref = sum(int(a) > 0 for a in alleles)
    if nonref == 0:
        return "0/0", n_imputed
    if nonref == 1:
        return "0/1", n_imputed
    return "1/1", n_imputed


def carriers_from_record(record: dict[str, object]) -> tuple[str, str, dict[str, int]]:
    het: list[str] = []
    hom_alt: list[str] = []
    counts = {"0/0": 0, "0/1": 0, "1/1": 0, "missing": 0, "imputed_alleles": 0}
    for sample, gt in record["genotypes"].items():  # type: ignore[union-attr]
        label, n_imputed = classify_gt(str(gt))
        counts[label] += 1
        counts["imputed_alleles"] += n_imputed
        if label == "0/1":
            het.append(str(sample))
        elif label == "1/1":
            hom_alt.append(str(sample))
    return sample_list(het), sample_list(hom_alt), counts


def choose_pangenie_lrs_match(
    lead: pd.Series,
    records: list[dict[str, object]],
    window_bp: int,
) -> tuple[dict[str, object] | None, str, int | None]:
    a1 = str(lead["A1"]).upper()
    a2 = str(lead["A2"]).upper()
    lead_id = str(lead["lead_id"])
    lead_pos = int(float(lead["POS"]))
    ranked: list[tuple[tuple[int, int, int, int], dict[str, object], str, int]] = []
    for record in records:
        distance = abs(int(record["pos"]) - lead_pos)
        ref = str(record["ref"]).upper()
        alts = split_alt(str(record["alt"]))
        exact = distance == 0 and ((ref == a1 and a2 in alts) or (ref == a2 and a1 in alts))
        fuzzy = distance <= window_bp and ref in {a1, a2}
        id_same = str(record["id"]) == lead_id
        match_type = "exact_ref_alt" if exact else "fuzzy_ref_same" if fuzzy else "none"
        rank = (int(exact), int(fuzzy), int(id_same), -distance)
        ranked.append((rank, record, match_type, distance))
    if not ranked:
        return None, "none", None
    ranked.sort(key=lambda x: x[0], reverse=True)
    _, record, match_type, distance = ranked[0]
    if match_type == "none":
        return None, "none", distance
    return record, match_type, distance


def extract_lrs(
    sv: pd.DataFrame,
    pangenie_vcf: Path,
    paragraph_vcf: Path,
    bcftools: str,
    window_bp: int,
) -> tuple[dict[tuple[str, str], dict[str, object]], list[dict[str, object]]]:
    resources: dict[str, dict[str, object]] = {}
    for name, vcf in {"pangenie_lrs": pangenie_vcf, "paragraph_lrs": paragraph_vcf}.items():
        check_file(vcf, name + " VCF")
        resources[name] = {
            "vcf": vcf,
            "samples": get_vcf_samples(vcf, bcftools),
            "contigs": get_vcf_contigs(vcf, bcftools),
        }

    result: dict[tuple[str, str], dict[str, object]] = {}
    audit: list[dict[str, object]] = []
    for _, lead in sv.iterrows():
        source = str(lead["source_set"])
        lead_id = str(lead["lead_id"])
        key = (source, lead_id)
        resource_name = "paragraph_lrs" if source == "paragraph" else "pangenie_lrs"
        resource = resources[resource_name]
        vcf = resource["vcf"]  # type: ignore[assignment]
        samples = resource["samples"]  # type: ignore[assignment]
        contigs = resource["contigs"]  # type: ignore[assignment]
        chrom = vcf_chr(str(lead["CHR"]), contigs)  # type: ignore[arg-type]
        pos = int(float(lead["POS"]))
        raw_records = query_region(vcf, chrom, pos - window_bp, pos + window_bp, bcftools)  # type: ignore[arg-type]
        records = [parse_record(x, samples) for x in raw_records]  # type: ignore[arg-type]

        if source == "paragraph":
            hits = [r for r in records if str(r["id"]) == lead_id]
            if len(hits) > 1:
                raise RuntimeError(f"Multiple direct Paragraph LRS ID matches: {lead_id}")
            record = hits[0] if hits else None
            match_type = "direct_id" if record else "none"
            distance = abs(int(record["pos"]) - pos) if record else None
        elif source in PANGENIE_SOURCES:
            record, match_type, distance = choose_pangenie_lrs_match(lead, records, window_bp)
        else:
            raise ValueError(f"Unexpected SV source_set: {source}")

        if record is None:
            result[key] = {"het": "", "hom_alt": ""}
            audit.append(
                {
                    "platform": "LRS",
                    "source_set": source,
                    "lead_id": lead_id,
                    "status": "NO_MATCH",
                    "resource": str(vcf),
                    "match_type": match_type,
                    "n_region_candidates": len(records),
                }
            )
            continue

        het, hom_alt, counts = carriers_from_record(record)
        result[key] = {"het": het, "hom_alt": hom_alt}
        audit.append(
            {
                "platform": "LRS",
                "source_set": source,
                "lead_id": lead_id,
                "status": "PASS",
                "resource": str(vcf),
                "match_type": match_type,
                "matched_variant_id": record["id"],
                "matched_chr": record["chrom"],
                "matched_pos": record["pos"],
                "distance_bp": distance,
                "matched_ref": record["ref"],
                "matched_alt": record["alt"],
                "n_region_candidates": len(records),
                "n_samples": len(samples),  # type: ignore[arg-type]
                "n_gt_0_0": counts["0/0"],
                "n_gt_0_1": counts["0/1"],
                "n_gt_1_1": counts["1/1"],
                "n_gt_missing": counts["missing"],
                "n_missing_alleles_imputed_as_ref": counts["imputed_alleles"],
            }
        )
    return result, audit


def read_bim(prefix: Path) -> pd.DataFrame:
    return pd.read_csv(
        str(prefix) + ".bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        dtype=str,
    )


def extract_ngs(
    sv: pd.DataFrame,
    plink: str,
    workdir: Path,
) -> tuple[dict[tuple[str, str], dict[str, object]], list[dict[str, object]]]:
    result: dict[tuple[str, str], dict[str, object]] = {}
    audit: list[dict[str, object]] = []
    workdir.mkdir(parents=True, exist_ok=True)
    for source, source_leads in sv.groupby("source_set", sort=False):
        source = str(source)
        if source not in NGS_BFILES:
            raise ValueError(f"No NGS bfile configured for source_set={source}")
        prefix = NGS_BFILES[source]
        check_bfile(prefix)
        bim = read_bim(prefix)
        bim_by_id = bim.drop_duplicates("SNP").set_index("SNP", drop=False)
        ids = [str(x) for x in source_leads["lead_id"]]
        found = [x for x in ids if x in bim_by_id.index]
        missing = [x for x in ids if x not in bim_by_id.index]
        for lead_id in missing:
            result[(source, lead_id)] = {"het": "", "hom_alt": ""}
            audit.append(
                {
                    "platform": "NGS",
                    "source_set": source,
                    "lead_id": lead_id,
                    "status": "ABSENT_FROM_BIM",
                    "resource": str(prefix),
                }
            )
        if not found:
            continue

        ids_file = workdir / f"{source}.sv_leads.ids"
        ids_file.write_text("\n".join(found) + "\n", encoding="utf-8")
        out_prefix = workdir / f"{source}.sv_leads"
        run(
            [
                plink,
                "--bfile", str(prefix),
                "--extract", str(ids_file),
                "--recode", "A",
                "--keep-allele-order",
                "--allow-no-sex",
                "--threads", "1",
                "--out", str(out_prefix),
            ],
            capture=True,
        )
        raw_file = Path(str(out_prefix) + ".raw")
        check_file(raw_file, source + " PLINK raw")
        raw = pd.read_csv(raw_file, sep=r"\s+", dtype=str)
        if "IID" not in raw.columns:
            raise ValueError(f"PLINK raw lacks IID: {raw_file}")

        for lead_id in found:
            b = bim_by_id.loc[lead_id]
            expected = f"{lead_id}_{b['A1']}"
            candidates = [c for c in raw.columns if c == expected or c.startswith(lead_id + "_")]
            if expected in raw.columns:
                column = expected
            elif len(candidates) == 1:
                column = candidates[0]
            else:
                raise ValueError(f"Cannot identify unique PLINK dosage column for {source}/{lead_id}: {candidates}")
            dosage = pd.to_numeric(raw[column], errors="coerce")
            het = raw.loc[dosage.round() == 1, "IID"].astype(str).tolist()
            hom_alt = raw.loc[dosage.round() == 2, "IID"].astype(str).tolist()
            result[(source, lead_id)] = {"het": sample_list(het), "hom_alt": sample_list(hom_alt)}
            audit.append(
                {
                    "platform": "NGS",
                    "source_set": source,
                    "lead_id": lead_id,
                    "status": "PASS",
                    "resource": str(prefix),
                    "plink_raw_column": column,
                    "plink_counted_allele_A1": b["A1"],
                    "plink_other_allele_A2": b["A2"],
                    "n_samples": len(raw),
                    "n_gt_0_0": int((dosage.round() == 0).sum()),
                    "n_gt_0_1": len(het),
                    "n_gt_1_1": len(hom_alt),
                    "n_gt_missing": int(dosage.isna().sum()),
                }
            )
    return result, audit


def main() -> int:
    args = parse_args()
    input_file = Path(args.input)
    outdir = Path(args.outdir)
    workdir = outdir / "work"
    check_file(input_file, "compact lead table")
    if args.window_bp < 0:
        raise ValueError("--window-bp must be non-negative")

    compact = pd.read_csv(input_file, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_set", "variant_type", "lead_id", "CHR", "POS", "A1", "A2"}
    missing_columns = sorted(required.difference(compact.columns))
    if missing_columns:
        raise ValueError("Compact table missing: " + ", ".join(missing_columns))
    sv = compact[compact["variant_type"].str.lower() == "sv"].copy()
    sv = sv[sv["CHR"].map(clean_chr).isin([str(x) for x in range(1, 23)])].copy()
    if sv.empty:
        raise ValueError("No autosomal SV rows in compact table")
    if sv.duplicated(["source_set", "lead_id"]).any():
        raise ValueError("Duplicate source_set/lead_id rows in compact SV table")
    unexpected = sorted(set(sv["source_set"]) - (PANGENIE_SOURCES | {"paragraph"}))
    if unexpected:
        raise ValueError("Unexpected SV source_set values: " + ", ".join(unexpected))

    outdir.mkdir(parents=True, exist_ok=True)
    lrs, lrs_audit = extract_lrs(
        sv,
        Path(args.lrs_pangenie_vcf),
        Path(args.lrs_paragraph_vcf),
        args.bcftools,
        args.window_bp,
    )
    ngs, ngs_audit = extract_ngs(sv, args.plink, workdir)

    sv["NGS_GT_0_1_sample_ids"] = [ngs.get((r.source_set, r.lead_id), {}).get("het", "") for r in sv.itertuples()]
    sv["NGS_GT_1_1_sample_ids"] = [ngs.get((r.source_set, r.lead_id), {}).get("hom_alt", "") for r in sv.itertuples()]
    sv["LRS_GT_0_1_sample_ids"] = [lrs.get((r.source_set, r.lead_id), {}).get("het", "") for r in sv.itertuples()]
    sv["LRS_GT_1_1_sample_ids"] = [lrs.get((r.source_set, r.lead_id), {}).get("hom_alt", "") for r in sv.itertuples()]

    stem = "lead_sig_from_gwas.final_merged.cleaned.noloc.compact.SV_sample"
    out_tsv = outdir / f"{stem}.tsv"
    out_csv = outdir / f"{stem}.csv"
    audit_file = outdir / "SV_sample_extraction.audit.tsv"
    config_file = outdir / "SV_sample_extraction.config.tsv"
    sv.to_csv(out_tsv, sep="\t", index=False)
    sv.to_csv(out_csv, index=False)
    pd.DataFrame(lrs_audit + ngs_audit).to_csv(audit_file, sep="\t", index=False)
    pd.DataFrame(
        [
            {"parameter": "input", "value": str(input_file)},
            {"parameter": "n_autosomal_sv", "value": len(sv)},
            {"parameter": "lrs_pangenie_vcf", "value": args.lrs_pangenie_vcf},
            {"parameter": "lrs_paragraph_vcf", "value": args.lrs_paragraph_vcf},
            {"parameter": "lrs_pangenie_match", "value": f"exact_ref_alt then fuzzy_ref_same within +/-{args.window_bp} bp"},
            {"parameter": "lrs_paragraph_match", "value": "direct VCF ID match within coordinate query window"},
            {"parameter": "gt_missing_allele_policy", "value": ". allele -> reference 0 before carrier classification"},
            {"parameter": "carrier_definition", "value": "one non-reference allele=0/1; two non-reference alleles=1/1"},
            {"parameter": "ngs_plink_rule", "value": "PLINK --recode A; dosage of BIM A1: 1=0/1, 2=1/1"},
        ]
    ).to_csv(config_file, sep="\t", index=False)

    print(f"n_autosomal_sv={len(sv)}")
    print(f"output_tsv={out_tsv}")
    print(f"output_csv={out_csv}")
    print(f"audit={audit_file}")
    print(f"config={config_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
