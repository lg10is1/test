from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


DEFAULT_PANELS = [("A", "B"), ("A", "C"), ("B", "C")]


@dataclass
class ComparisonSpec:
    output_tag: str
    result_dir: Path
    result_tag: str
    output_dir: Path
    label: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-generate HLA A-B, A-C, and B-C correlation heatmaps."
    )
    parser.add_argument("--scz-xlsx", required=True, help="SCZ Immuannot xlsx path.")
    parser.add_argument(
        "--comparison",
        action="append",
        required=True,
        help=(
            "Comparison spec with 5 fields separated by '|': "
            "output_tag|result_dir|result_tag|output_dir|label"
        ),
    )
    parser.add_argument(
        "--panel-script",
        default=str(Path(__file__).with_name("hla_abc_correlation_heatmap_universal.R")),
        help="Single-panel R script path.",
    )
    parser.add_argument(
        "--rscript-path",
        default="",
        help="Optional explicit Rscript.exe path.",
    )
    parser.add_argument(
        "--sample-label",
        default="case_cohort_ragtag",
        help="Sample label used in plot titles.",
    )
    parser.add_argument(
        "--min-scz-count",
        type=int,
        default=5,
        help="Minimum SCZ count threshold when selecting alleles from Fig.4F outputs.",
    )
    parser.add_argument(
        "--label-digits",
        type=int,
        default=2,
        help="Decimal digits shown in heatmap labels.",
    )
    parser.add_argument(
        "--allow-nominal-fallback",
        action="store_true",
        help="If a gene has no FDR-significant alleles, fall back to nominal P threshold.",
    )
    parser.add_argument(
        "--nominal-p-threshold",
        type=float,
        default=0.05,
        help="Nominal P-value threshold used when fallback is enabled.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Reuse existing panel outputs instead of overwriting them.",
    )
    parser.add_argument("--a-alleles", default="", help="Optional extra A alleles, comma-separated.")
    parser.add_argument("--b-alleles", default="", help="Optional extra B alleles, comma-separated.")
    parser.add_argument("--c-alleles", default="", help="Optional extra C alleles, comma-separated.")
    return parser.parse_args()


def parse_comparison(raw: str) -> ComparisonSpec:
    parts = [part.strip() for part in raw.split("|")]
    if len(parts) != 5:
        raise ValueError(
            "Each --comparison must be: output_tag|result_dir|result_tag|output_dir|label"
        )
    output_tag, result_dir, result_tag, output_dir, label = parts
    return ComparisonSpec(
        output_tag=output_tag,
        result_dir=Path(result_dir),
        result_tag=result_tag,
        output_dir=Path(output_dir),
        label=label,
    )


def resolve_rscript(explicit: str) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit))

    script_dir = Path(sys.executable).resolve().parent
    candidates.extend(
        [
            script_dir / "Rscript.exe",
            script_dir / "Scripts" / "Rscript.exe",
        ]
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate

    found = shutil.which("Rscript")
    if found:
        return Path(found)

    raise FileNotFoundError("Unable to locate Rscript.exe")


def result_file(result_dir: Path, gene: str, result_tag: str) -> Path:
    path = result_dir / f"{gene}_fisher_results_new_{result_tag}.xlsx"
    if not path.exists():
        raise FileNotFoundError(f"Missing result file: {path}")
    return path


def manual_alleles_for_gene(gene: str, args: argparse.Namespace) -> str:
    mapping = {
        "A": args.a_alleles,
        "B": args.b_alleles,
        "C": args.c_alleles,
    }
    return mapping.get(gene, "")


def run_panel(
    rscript_path: Path,
    panel_script: Path,
    scz_xlsx: Path,
    comparison: ComparisonSpec,
    sample_label: str,
    min_scz_count: int,
    label_digits: int,
    allow_nominal_fallback: bool,
    nominal_p_threshold: float,
    skip_existing: bool,
    col_gene: str,
    row_gene: str,
    col_alleles: str,
    row_alleles: str,
) -> dict:
    comparison.output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = f"{col_gene}{row_gene}_correlation_{comparison.output_tag}"
    title = f"HLA-{col_gene} × HLA-{row_gene} correlation heatmap ({sample_label}, {comparison.label})"
    pdf_path = comparison.output_dir / f"{output_prefix}_heatmap.pdf"
    png_path = comparison.output_dir / f"{output_prefix}_heatmap.png"
    summary_csv_dir = comparison.output_dir / f"{output_prefix}_summary_csv"
    summary_xlsx = comparison.output_dir / f"{output_prefix}_summary.xlsx"

    if skip_existing and pdf_path.exists() and png_path.exists() and (summary_xlsx.exists() or summary_csv_dir.exists()):
        if summary_csv_dir.exists() and not summary_xlsx.exists():
            convert_summary_csv_dir_to_xlsx(summary_csv_dir, summary_xlsx)
        return {
            "panel": f"{col_gene}x{row_gene}",
            "output_prefix": output_prefix,
            "pdf": pdf_path,
            "png": png_path,
            "summary_xlsx": summary_xlsx if summary_xlsx.exists() else "",
            "summary_csv_dir": summary_csv_dir if summary_csv_dir.exists() else "",
        }

    command = [
        str(rscript_path),
        str(panel_script),
        f"--scz-xlsx={scz_xlsx}",
        f"--col-gene={col_gene}",
        f"--row-gene={row_gene}",
        f"--col-result={result_file(comparison.result_dir, col_gene, comparison.result_tag)}",
        f"--row-result={result_file(comparison.result_dir, row_gene, comparison.result_tag)}",
        f"--col-alleles={col_alleles}",
        f"--row-alleles={row_alleles}",
        f"--min-scz-count={min_scz_count}",
        f"--label-digits={label_digits}",
        f"--allow-nominal-fallback={'TRUE' if allow_nominal_fallback else 'FALSE'}",
        f"--nominal-p-threshold={nominal_p_threshold}",
        f"--output-dir={comparison.output_dir}",
        f"--output-prefix={output_prefix}",
        f"--plot-title={title}",
    ]

    subprocess.run(command, check=True)

    if summary_csv_dir.exists() and not summary_xlsx.exists():
        convert_summary_csv_dir_to_xlsx(summary_csv_dir, summary_xlsx)

    return {
        "panel": f"{col_gene}x{row_gene}",
        "output_prefix": output_prefix,
        "pdf": comparison.output_dir / f"{output_prefix}_heatmap.pdf",
        "png": comparison.output_dir / f"{output_prefix}_heatmap.png",
        "summary_xlsx": summary_xlsx if summary_xlsx.exists() else "",
        "summary_csv_dir": summary_csv_dir if summary_csv_dir.exists() else "",
    }


def convert_summary_csv_dir_to_xlsx(summary_dir: Path, output_xlsx: Path) -> None:
    csv_files = sorted(summary_dir.glob("*.csv"))
    if not csv_files:
        return

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        for csv_file in csv_files:
            df = pd.read_csv(csv_file)
            df.to_excel(writer, sheet_name=csv_file.stem[:31], index=False)


def write_manifest(output_dir: Path, output_tag: str, rows: list[dict]) -> None:
    manifest = pd.DataFrame(rows)
    csv_path = output_dir / f"panel_manifest_{output_tag}.csv"
    xlsx_path = output_dir / f"panel_manifest_{output_tag}.xlsx"
    manifest.to_csv(csv_path, index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        manifest.to_excel(writer, sheet_name="manifest", index=False)


def main() -> None:
    args = parse_args()
    scz_xlsx = Path(args.scz_xlsx)
    panel_script = Path(args.panel_script)
    rscript_path = resolve_rscript(args.rscript_path)

    if not scz_xlsx.exists():
        raise FileNotFoundError(f"Missing SCZ xlsx: {scz_xlsx}")
    if not panel_script.exists():
        raise FileNotFoundError(f"Missing panel script: {panel_script}")

    comparisons = [parse_comparison(raw) for raw in args.comparison]

    for comparison in comparisons:
        manifest_rows = []
        for col_gene, row_gene in DEFAULT_PANELS:
            result = run_panel(
                rscript_path=rscript_path,
                panel_script=panel_script,
                scz_xlsx=scz_xlsx,
                comparison=comparison,
                sample_label=args.sample_label,
                min_scz_count=args.min_scz_count,
                label_digits=args.label_digits,
                allow_nominal_fallback=args.allow_nominal_fallback,
                nominal_p_threshold=args.nominal_p_threshold,
                skip_existing=args.skip_existing,
                col_gene=col_gene,
                row_gene=row_gene,
                col_alleles=manual_alleles_for_gene(col_gene, args),
                row_alleles=manual_alleles_for_gene(row_gene, args),
            )
            manifest_rows.append(result)

        write_manifest(comparison.output_dir, comparison.output_tag, manifest_rows)
        print(f"Completed comparison: {comparison.output_tag}")


if __name__ == "__main__":
    main()
