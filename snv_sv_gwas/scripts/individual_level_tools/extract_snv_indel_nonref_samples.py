# PRIVACY WARNING: This optional tool writes individual sample identifiers.
# Use it only on authorized controlled data and never commit its outputs.
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path("/path/to/local/project")
GWAS_LEAD_DIR = ROOT / "Figure3_manual/summary_maintext_public/01_GWAS_leads"
DEFAULT_INPUT = GWAS_LEAD_DIR / "lead_sig_from_gwas.final_merged.cleaned.noloc.compact.tsv"
DEFAULT_OUTDIR = GWAS_LEAD_DIR
PANGENIE_GWAS = ROOT / "TGS_callset/Pangenie_v3/06.gwas"
PANGENIE_SOURCES = ("set00",)
DEFAULT_TGS_SAMPLES = PANGENIE_GWAS / "tgs_sample.txt"
DEFAULT_IGV_OUTDIR = ROOT / "Figure3_manual/summary_public/10_igv_check"
DEFAULT_LRS_BAM_DIR = Path(
    "/path/to/shared_resources/project/authorized_long_read_data/align_t2t"
)
DEFAULT_IGV_GENOME = Path(
    "/path/to/shared_resources/references/GCF_009914755.1_T2T-CHM13v2.0_genomic_chr.fasta"
)
DEFAULT_IGV_GTF = Path(
    "/path/to/shared_resources/references/T2T-CHM13V2.0_UCSC.sorted.gtf"
)
NGS_TGS_BFILES = {
    source: PANGENIE_GWAS / source / "NGS_TGS.QCsite"
    for source in PANGENIE_SOURCES
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Add merged NGS+TGS heterozygous/homozygous-alternate carrier "
            "sample IDs to the public compact SNV/indel lead table."
        )
    )
    p.add_argument("--input", default=str(DEFAULT_INPUT))
    p.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    p.add_argument("--plink", default="plink")
    p.add_argument("--tgs-samples", default=str(DEFAULT_TGS_SAMPLES))
    p.add_argument("--igv-outdir", default=str(DEFAULT_IGV_OUTDIR))
    p.add_argument("--lrs-bam-dir", default=str(DEFAULT_LRS_BAM_DIR))
    p.add_argument("--igv-genome", default=str(DEFAULT_IGV_GENOME))
    p.add_argument("--igv-gtf", default=str(DEFAULT_IGV_GTF))
    p.add_argument("--igv-window-bp", type=int, default=1000)
    p.add_argument("--igv-max-samples-per-genotype", type=int, default=5)
    p.add_argument("--skip-igv-scripts", action="store_true")
    return p.parse_args()


def run(cmd: list[str], *, capture: bool = True) -> str:
    result = subprocess.run(
        cmd,
        check=False,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )
    if result.returncode != 0:
        msg = "Command failed: " + " ".join(cmd)
        if capture:
            msg += "\nSTDOUT:\n" + (result.stdout or "")
            msg += "\nSTDERR:\n" + (result.stderr or "")
        raise RuntimeError(msg)
    return result.stdout or ""


def check_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def check_bfile(prefix: Path) -> None:
    for suffix in (".bed", ".bim", ".fam"):
        check_file(Path(str(prefix) + suffix), f"PLINK bfile {suffix}")


def clean_chr(value: object) -> str:
    return str(value).replace("chr", "").replace("CHR", "")


def sample_list(values: list[str]) -> str:
    return ";".join(sorted(set(str(x) for x in values if str(x))))


def split_sample_list(value: object) -> list[str]:
    text = str(value or "")
    return [x for x in text.split(";") if x]


def safe_token(value: object) -> str:
    text = str(value)
    out = []
    for ch in text:
        if ch.isalnum() or ch in ("-", "_", "."):
            out.append(ch)
        else:
            out.append("_")
    token = "".join(out).strip("_")
    return token or "NA"


def chr_label(value: object) -> str:
    chrom = clean_chr(value)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def read_tgs_samples(path: Path) -> set[str]:
    check_file(path, "TGS sample list")
    dt = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, comment="#")
    if dt.empty:
        raise ValueError(f"TGS sample list is empty: {path}")
    sample_col = 1 if dt.shape[1] >= 2 else 0
    samples = set(dt.iloc[:, sample_col].dropna().astype(str))
    if not samples:
        raise ValueError(f"No sample IDs parsed from TGS sample list: {path}")
    return samples


def read_bim(prefix: Path) -> pd.DataFrame:
    return pd.read_csv(
        str(prefix) + ".bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        dtype=str,
    )


def write_igv_batch_scripts(
    snv_indel: pd.DataFrame,
    igv_outdir: Path,
    bam_dir: Path,
    genome: Path,
    gtf: Path,
    window_bp: int,
    max_samples: int,
) -> tuple[Path, pd.DataFrame]:
    if window_bp < 0:
        raise ValueError("--igv-window-bp must be non-negative")
    if max_samples < 1:
        raise ValueError("--igv-max-samples-per-genotype must be positive")

    igv_outdir.mkdir(parents=True, exist_ok=True)
    script_file = igv_outdir / "snv_indel_lrs_all.igv"

    genotype_specs = [
        ("0_0", "LRS_GT_0_0_sample_ids", "GT0"),
        ("0_1", "LRS_GT_0_1_sample_ids", "GT1"),
        ("1_1", "LRS_GT_1_1_sample_ids", "GT2"),
    ]
    manifest_rows: list[dict[str, object]] = []
    lines = [
        "setSleepInterval 50",
        "maxPanelHeight 3000",
    ]

    for idx, row in enumerate(snv_indel.itertuples(index=False), start=1):
        lead_id = getattr(row, "lead_id")
        chrom = chr_label(getattr(row, "CHR"))
        pos = int(float(getattr(row, "POS")))
        start = max(1, pos - window_bp)
        end = max(1, pos + window_bp)
        lead_token = safe_token(lead_id)
        variant_stem = f"{idx:04d}_{lead_token}_{chrom}_{pos}"

        for genotype_label, sample_col, suffix in genotype_specs:
            samples = split_sample_list(getattr(row, sample_col, ""))[:max_samples]
            bam_paths = [
                bam_dir
                / (
                    f"{sample}.align_chm13.bam"
                    if re.fullmatch(r"SPECIAL_SAMPLE_[0-9]+", sample)
                    else f"{sample}.align.bam"
                )
                for sample in samples
            ]
            snapshot_name = f"{variant_stem}_{suffix}.png"
            lines.append(f"# {variant_stem} {genotype_label} samples={','.join(samples) if samples else 'NONE'}")
            lines.append("new")
            lines.append(f"genome {genome}")
            for bam_path in bam_paths:
                lines.append(f"load {bam_path}")
            lines.append(f"load {gtf}")
            lines.append("collapse")
            lines.append("expand")
            lines.append(f"snapshotDirectory {igv_outdir}")
            lines.append(f"goto {chrom}:{start}-{end}")
            lines.append("expand")
            lines.append("sort base")
            lines.append(f"snapshot {snapshot_name}")
            manifest_rows.append(
                {
                    "script": str(script_file),
                    "snapshot": str(igv_outdir / snapshot_name),
                    "variant_index": f"{idx:04d}",
                    "source_set": getattr(row, "source_set"),
                    "lead_id": lead_id,
                    "chr": chrom,
                    "pos": pos,
                    "genotype": genotype_label.replace("_", "/"),
                    "n_samples_loaded": len(samples),
                    "samples_loaded": ";".join(samples),
                    "bam_paths_loaded": ";".join(str(path) for path in bam_paths),
                }
            )

    script_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return script_file, pd.DataFrame(manifest_rows)


def extract_ngs_tgs(
    snv_indel: pd.DataFrame,
    plink: str,
    workdir: Path,
    tgs_samples: set[str],
) -> tuple[dict[tuple[str, str], dict[str, object]], list[dict[str, object]]]:
    result: dict[tuple[str, str], dict[str, object]] = {}
    audit: list[dict[str, object]] = []
    workdir.mkdir(parents=True, exist_ok=True)

    for source, source_leads in snv_indel.groupby("source_set", sort=False):
        source = str(source)
        if source not in NGS_TGS_BFILES:
            raise ValueError(f"No NGS_TGS bfile configured for source_set={source}")
        prefix = NGS_TGS_BFILES[source]
        check_bfile(prefix)
        bim = read_bim(prefix)
        bim_by_id = bim.drop_duplicates("SNP").set_index("SNP", drop=False)

        ids = [str(x) for x in source_leads["lead_id"]]
        found = [x for x in ids if x in bim_by_id.index]
        missing = [x for x in ids if x not in bim_by_id.index]
        for lead_id in missing:
            result[(source, lead_id)] = {
                "lrs_ref": "",
                "ngs_het": "",
                "ngs_hom_alt": "",
                "lrs_het": "",
                "lrs_hom_alt": "",
            }
            audit.append(
                {
                    "platform": "NGS_TGS_split",
                    "source_set": source,
                    "lead_id": lead_id,
                    "status": "ABSENT_FROM_BIM",
                    "resource": str(prefix),
                }
            )
        if not found:
            continue

        ids_file = workdir / f"{source}.snv_indel_leads.ids"
        ids_file.write_text("\n".join(found) + "\n", encoding="utf-8")
        out_prefix = workdir / f"{source}.snv_indel_leads"
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
                raise ValueError(
                    f"Cannot identify unique PLINK dosage column for {source}/{lead_id}: {candidates}"
                )

            dosage = pd.to_numeric(raw[column], errors="coerce")
            rounded = dosage.round()
            iid = raw["IID"].astype(str)
            is_lrs = iid.isin(tgs_samples)
            is_ngs = ~is_lrs
            lrs_ref = raw.loc[(rounded == 0) & is_lrs, "IID"].astype(str).tolist()
            ngs_het = raw.loc[(rounded == 1) & is_ngs, "IID"].astype(str).tolist()
            ngs_hom_alt = raw.loc[(rounded == 2) & is_ngs, "IID"].astype(str).tolist()
            lrs_het = raw.loc[(rounded == 1) & is_lrs, "IID"].astype(str).tolist()
            lrs_hom_alt = raw.loc[(rounded == 2) & is_lrs, "IID"].astype(str).tolist()
            result[(source, lead_id)] = {
                "lrs_ref": sample_list(lrs_ref),
                "ngs_het": sample_list(ngs_het),
                "ngs_hom_alt": sample_list(ngs_hom_alt),
                "lrs_het": sample_list(lrs_het),
                "lrs_hom_alt": sample_list(lrs_hom_alt),
            }
            audit.append(
                {
                    "platform": "NGS_TGS_split",
                    "source_set": source,
                    "lead_id": lead_id,
                    "status": "PASS",
                    "resource": str(prefix),
                    "plink_raw_column": column,
                    "plink_counted_allele_A1": b["A1"],
                    "plink_other_allele_A2": b["A2"],
                    "n_samples": len(raw),
                    "n_ngs_samples": int(is_ngs.sum()),
                    "n_lrs_samples": int(is_lrs.sum()),
                    "n_gt_0_0": int((rounded == 0).sum()),
                    "n_gt_0_1": int((rounded == 1).sum()),
                    "n_gt_1_1": int((rounded == 2).sum()),
                    "n_gt_missing": int(dosage.isna().sum()),
                    "n_ngs_gt_0_0": int(((rounded == 0) & is_ngs).sum()),
                    "n_ngs_gt_0_1": len(ngs_het),
                    "n_ngs_gt_1_1": len(ngs_hom_alt),
                    "n_ngs_gt_missing": int((dosage.isna() & is_ngs).sum()),
                    "n_lrs_gt_0_0": int(((rounded == 0) & is_lrs).sum()),
                    "n_lrs_gt_0_1": len(lrs_het),
                    "n_lrs_gt_1_1": len(lrs_hom_alt),
                    "n_lrs_gt_missing": int((dosage.isna() & is_lrs).sum()),
                }
            )

    return result, audit


def main() -> int:
    args = parse_args()
    input_file = Path(args.input)
    outdir = Path(args.outdir)
    workdir = outdir / "work_snv_indel_sample"
    check_file(input_file, "compact lead table")

    compact = pd.read_csv(input_file, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_set", "variant_type", "lead_id", "CHR", "POS", "A1", "A2"}
    missing_columns = sorted(required.difference(compact.columns))
    if missing_columns:
        raise ValueError("Compact table missing: " + ", ".join(missing_columns))

    snv_indel = compact[
        (compact["variant_type"].str.lower() == "snv_indel")
        & compact["source_set"].isin(PANGENIE_SOURCES)
    ].copy()
    snv_indel = snv_indel[
        snv_indel["CHR"].map(clean_chr).isin([str(x) for x in range(1, 23)])
    ].copy()
    if snv_indel.empty:
        raise ValueError(f"No autosomal {','.join(PANGENIE_SOURCES)} SNV/indel rows in compact table")
    if snv_indel.duplicated(["source_set", "lead_id"]).any():
        raise ValueError("Duplicate source_set/lead_id rows in compact SNV/indel table")

    outdir.mkdir(parents=True, exist_ok=True)
    tgs_samples = read_tgs_samples(Path(args.tgs_samples))
    carriers, audit = extract_ngs_tgs(snv_indel, args.plink, workdir, tgs_samples)

    snv_indel["LRS_GT_0_0_sample_ids"] = [
        carriers.get((r.source_set, r.lead_id), {}).get("lrs_ref", "")
        for r in snv_indel.itertuples()
    ]
    snv_indel["NGS_GT_0_1_sample_ids"] = [
        carriers.get((r.source_set, r.lead_id), {}).get("ngs_het", "")
        for r in snv_indel.itertuples()
    ]
    snv_indel["NGS_GT_1_1_sample_ids"] = [
        carriers.get((r.source_set, r.lead_id), {}).get("ngs_hom_alt", "")
        for r in snv_indel.itertuples()
    ]
    snv_indel["LRS_GT_0_1_sample_ids"] = [
        carriers.get((r.source_set, r.lead_id), {}).get("lrs_het", "")
        for r in snv_indel.itertuples()
    ]
    snv_indel["LRS_GT_1_1_sample_ids"] = [
        carriers.get((r.source_set, r.lead_id), {}).get("lrs_hom_alt", "")
        for r in snv_indel.itertuples()
    ]

    stem = "lead_sig_from_gwas.final_merged.cleaned.noloc.compact.SNV_INDEL_sample"
    out_tsv = outdir / f"{stem}.tsv"
    out_csv = outdir / f"{stem}.csv"
    audit_file = outdir / "SNV_INDEL_sample_extraction.audit.tsv"
    config_file = outdir / "SNV_INDEL_sample_extraction.config.tsv"
    igv_manifest_file = Path(args.igv_outdir) / "snv_indel_lrs_igv_manifest.tsv"
    snv_indel.to_csv(out_tsv, sep="\t", index=False)
    snv_indel.to_csv(out_csv, index=False)
    pd.DataFrame(audit).to_csv(audit_file, sep="\t", index=False)

    igv_batch = ""
    if not args.skip_igv_scripts:
        igv_batch_path, igv_manifest = write_igv_batch_scripts(
            snv_indel,
            Path(args.igv_outdir),
            Path(args.lrs_bam_dir),
            Path(args.igv_genome),
            Path(args.igv_gtf),
            args.igv_window_bp,
            args.igv_max_samples_per_genotype,
        )
        igv_manifest.to_csv(igv_manifest_file, sep="\t", index=False)
        igv_batch = str(igv_batch_path)

    pd.DataFrame(
        [
            {"parameter": "input", "value": str(input_file)},
            {"parameter": "n_autosomal_snv_indel", "value": len(snv_indel)},
            {"parameter": "source_sets", "value": ",".join(PANGENIE_SOURCES)},
            *[
                {"parameter": f"{source}_bfile", "value": str(NGS_TGS_BFILES[source])}
                for source in PANGENIE_SOURCES
            ],
            {"parameter": "tgs_sample_file", "value": str(args.tgs_samples)},
            {"parameter": "n_tgs_samples_in_sample_file", "value": len(tgs_samples)},
            {"parameter": "excluded_source_set", "value": "deepvariant"},
            {"parameter": "carrier_definition", "value": "PLINK A1 dosage 1=0/1; dosage 2=1/1"},
            {"parameter": "sample_split_rule", "value": "IID in tgs_sample.txt -> LRS/TGS; all other IID -> NGS"},
            {"parameter": "lrs_reference_column", "value": "LRS_GT_0_0_sample_ids records only LRS/TGS 0/0 samples; NGS 0/0 samples are intentionally omitted"},
            {"parameter": "plink_rule", "value": "PLINK --recode A with --keep-allele-order"},
            {"parameter": "igv_scripts_enabled", "value": not args.skip_igv_scripts},
            {"parameter": "igv_outdir", "value": str(args.igv_outdir)},
            {"parameter": "igv_manifest", "value": str(igv_manifest_file) if not args.skip_igv_scripts else ""},
            {"parameter": "igv_batch", "value": igv_batch},
            {"parameter": "igv_lrs_bam_dir", "value": str(args.lrs_bam_dir)},
            {"parameter": "igv_genome", "value": str(args.igv_genome)},
            {"parameter": "igv_gtf", "value": str(args.igv_gtf)},
            {"parameter": "igv_window_bp", "value": args.igv_window_bp},
            {"parameter": "igv_max_samples_per_genotype", "value": args.igv_max_samples_per_genotype},
        ]
    ).to_csv(config_file, sep="\t", index=False)

    print(f"n_autosomal_snv_indel={len(snv_indel)}")
    print(f"output_tsv={out_tsv}")
    print(f"output_csv={out_csv}")
    print(f"audit={audit_file}")
    print(f"config={config_file}")
    if not args.skip_igv_scripts:
        print(f"igv_manifest={igv_manifest_file}")
        print(f"igv_batch={igv_batch}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
