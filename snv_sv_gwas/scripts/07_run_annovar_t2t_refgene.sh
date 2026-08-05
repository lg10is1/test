#!/usr/bin/env bash
set -euo pipefail

cd /path/to/EOSCZ_PROJECT/figure_analysis/SV_SNV_LD/LD_decay_public/tables

INPUT="lead_sig_from_gwas.all.tsv"
OUT_PREFIX="lead_sig_from_gwas.all.annovar"
OUT_PREFIX_NOLOC="lead_sig_from_gwas.all.annovar_noloc"
ANNOVAR_DB_T2T="${ANNOVAR_DB_T2T:-$HOME/software/annovar/test}"
ANNOVAR_DB_T2T_NOLOC="${ANNOVAR_DB_T2T_NOLOC:-$HOME/software/annovar/t2t_noloc}"

check_dir() {
  local dir="$1"
  local label="$2"
  if [[ ! -d "${dir}" ]]; then
    echo "[ERROR] ${label} ANNOVAR database directory not found: ${dir}" >&2
    exit 1
  fi
}

check_dir "${ANNOVAR_DB_T2T}" "standard t2t"
check_dir "${ANNOVAR_DB_T2T_NOLOC}" "t2t_noloc"

python3 - <<'PY'
import pandas as pd
import re

input_file = "lead_sig_from_gwas.all.tsv"
out_avinput = "lead_sig_from_gwas.all.annovar.avinput"
out_map = "lead_sig_from_gwas.all.annovar.input_map.tsv"

df = pd.read_csv(input_file, sep="\t", dtype=str)

def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

chr_col = pick_col(df, ["lead_chr", "CHR", "#CHROM", "CHROM"])
pos_col = pick_col(df, ["lead_pos", "POS", "BP", "Start", "START"])

ref_col = pick_col(df, ["REF", "Ref"])
alt_col = pick_col(df, ["ALT", "Alt"])

# fastGWA commonly uses A1/A2. Fallback: A2 as REF, A1 as ALT.
# refGene annotation mainly depends on coordinates; alleles keep avinput valid.
if ref_col is None or alt_col is None:
    ref_col = pick_col(df, ["A2", "OTHER_ALLELE"])
    alt_col = pick_col(df, ["A1", "EFFECT_ALLELE"])

if chr_col is None or pos_col is None or ref_col is None or alt_col is None:
    raise RuntimeError(
        "Cannot find required columns.\n"
        f"chr_col={chr_col}, pos_col={pos_col}, ref_col={ref_col}, alt_col={alt_col}\n"
        f"Available columns: {list(df.columns)}"
    )

records = []
maps = []

for _, row in df.iterrows():
    chrom = str(row[chr_col])
    pos_raw = str(row[pos_col])
    ref = str(row[ref_col])
    alt = str(row[alt_col])

    if chrom in ["nan", "NA", ".", ""]:
        continue
    if pos_raw in ["nan", "NA", ".", ""]:
        continue
    if ref in ["nan", "NA", ".", ""]:
        continue
    if alt in ["nan", "NA", ".", ""]:
        continue

    chrom = re.sub(r"^chr", "", chrom)
    pos = int(float(pos_raw))

    ref = ref.upper()
    alt = alt.upper()

    if "," in ref or "," in alt:
        continue
    if ";" in ref or ";" in alt:
        continue

    if len(ref) == len(alt):
        start = pos
        end = pos + len(ref) - 1
        av_ref = ref
        av_alt = alt
    elif len(ref) < len(alt) and alt.startswith(ref):
        ins_seq = alt[len(ref):]
        start = pos + len(ref) - 1
        end = start
        av_ref = "-"
        av_alt = ins_seq
    elif len(ref) > len(alt) and ref.startswith(alt):
        del_seq = ref[len(alt):]
        start = pos + len(alt)
        end = pos + len(ref) - 1
        av_ref = del_seq
        av_alt = "-"
    else:
        start = pos
        end = pos + len(ref) - 1
        av_ref = ref
        av_alt = alt

    lead_id = row["lead_id"] if "lead_id" in df.columns else f"{chrom}:{pos}:{ref}:{alt}"
    source_set = row["source_set"] if "source_set" in df.columns else "."
    variant_type = row["variant_type"] if "variant_type" in df.columns else "."

    records.append([chrom, start, end, av_ref, av_alt])
    maps.append([
        chrom, start, end, av_ref, av_alt,
        lead_id, source_set, variant_type,
        ref, alt
    ])

av = pd.DataFrame(records, columns=["CHR", "START", "END", "REF", "ALT"])
mp = pd.DataFrame(
    maps,
    columns=[
        "CHR", "START", "END", "REF", "ALT",
        "lead_id", "source_set", "variant_type",
        "original_REF_or_A2", "original_ALT_or_A1"
    ]
)

av.to_csv(out_avinput, sep="\t", header=False, index=False)
mp.to_csv(out_map, sep="\t", index=False)

print(f"[INFO] avinput written: {out_avinput}")
print(f"[INFO] map written: {out_map}")
print(f"[INFO] variants for ANNOVAR: {av.shape[0]}")
PY

run_annovar_refgene() {
  local out_prefix="$1"
  local annovar_db="$2"
  local label="$3"
  local buildver="$4"

  echo
  echo "============================================================"
  echo "[RUN] ANNOVAR ${label}"
  echo "============================================================"
  echo "[INFO] database: ${annovar_db}"

  "${ANNOVAR_TABLE:?Set ANNOVAR_TABLE to table_annovar.pl}" \
    "${OUT_PREFIX}.avinput" \
    "${annovar_db}" \
    -buildver "${buildver}" \
    -out "${out_prefix}" \
    -remove \
    -protocol refGene \
    -operation g \
    -nastring . \
    -csvout \
    -polish \
    -thread 40 \
    -maxgenethread 40
}

run_annovar_refgene "${OUT_PREFIX}" "${ANNOVAR_DB_T2T}" "t2t/refGene" "t2t"
run_annovar_refgene "${OUT_PREFIX_NOLOC}" "${ANNOVAR_DB_T2T_NOLOC}" "t2t_noloc/refGene" "t2t_noloc"

python3 - <<'PY'
import pandas as pd

prefixes = [
    ("lead_sig_from_gwas.all.annovar", "t2t"),
    ("lead_sig_from_gwas.all.annovar_noloc", "t2t_noloc"),
]
map_file = "lead_sig_from_gwas.all.annovar.input_map.tsv"

for out_prefix, buildver in prefixes:
    anno_file = out_prefix + f".{buildver}_multianno.csv"

    anno = pd.read_csv(anno_file, dtype=str)
    mp = pd.read_csv(map_file, sep="\t", dtype=str)

    anno = anno.rename(columns={
        "Chr": "CHR",
        "Start": "START",
        "End": "END",
        "Ref": "REF",
        "Alt": "ALT",
        "Func.refGene": "Func_refGene",
        "Gene.refGene": "Gene_refGene",
        "GeneDetail.refGene": "GeneDetail_refGene",
        "ExonicFunc.refGene": "ExonicFunc_refGene",
        "AAChange.refGene": "AAChange_refGene"
    })

    key_cols = ["CHR", "START", "END", "REF", "ALT"]

    for c in key_cols:
        anno[c] = anno[c].astype(str)
        mp[c] = mp[c].astype(str)

    anno_unique = anno.drop_duplicates(subset=key_cols, keep="first").copy()

    merged = mp.merge(
        anno_unique,
        on=key_cols,
        how="left",
        validate="many_to_one"
    )

    front_cols = [
        "lead_id",
        "source_set",
        "variant_type",
        "CHR",
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
        "AAChange_refGene"
    ]

    front_cols = [c for c in front_cols if c in merged.columns]
    other_cols = [c for c in merged.columns if c not in front_cols]
    merged = merged[front_cols + other_cols]

    merged.to_csv(out_prefix + ".refGene.annotated.fixed.tsv", sep="\t", index=False)
    merged.to_csv(out_prefix + ".refGene.annotated.fixed.csv", index=False)

    print("[INFO] prefix:", out_prefix)
    print("[INFO] input_map rows:", mp.shape[0])
    print("[INFO] annovar rows:", anno.shape[0])
    print("[INFO] annovar unique coordinate rows:", anno_unique.shape[0])
    print("[INFO] fixed annotated rows:", merged.shape[0])
    print("[INFO] output TSV:", out_prefix + ".refGene.annotated.fixed.tsv")
    print("[INFO] output CSV:", out_prefix + ".refGene.annotated.fixed.csv")
PY

echo
echo "Final outputs:"
echo "$(pwd)/lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.tsv"
echo "$(pwd)/lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.csv"
echo "$(pwd)/lead_sig_from_gwas.all.annovar_noloc.refGene.annotated.fixed.tsv"
echo "$(pwd)/lead_sig_from_gwas.all.annovar_noloc.refGene.annotated.fixed.csv"
echo
echo "Line counts:"
wc -l lead_sig_from_gwas.all.tsv \
      lead_sig_from_gwas.all.annovar.avinput \
      lead_sig_from_gwas.all.annovar.input_map.tsv \
      lead_sig_from_gwas.all.annovar.t2t_multianno.csv \
      lead_sig_from_gwas.all.annovar.refGene.annotated.fixed.tsv \
      lead_sig_from_gwas.all.annovar_noloc.t2t_noloc_multianno.csv \
      lead_sig_from_gwas.all.annovar_noloc.refGene.annotated.fixed.tsv
