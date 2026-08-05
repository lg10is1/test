#!/usr/bin/env bash

# Convert the merged VCF to PLINK files, calculate genotype consistency, filter
# variants and samples, run PCA/projection PCA, and build GCTA GRMs.
#
# Usage:
#   bash 04_run_bfile_qc.sh PROJECT_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PROJECT_DIR" >&2
    exit 2
fi

PROJECT_DIR="$(cd "$1" && pwd)"
BASE_WORK_DIR="${PROJECT_DIR}/06.gwas"
VCF_BASE="${PROJECT_DIR}/05.merge"
DATA_DIR="set00"
SUPPORT_DIR="${SUPPORT_DIR:-${BASE_WORK_DIR}}"
TGS_SAMPLE_FILE="${TGS_SAMPLE_FILE:-${SUPPORT_DIR}/tgs_sample.txt}"
BAD_SAMPLE_FILE="${BAD_SAMPLE_FILE:-${SUPPORT_DIR}/bad.sample}"
EAS_SAMPLE_FILE="${EAS_SAMPLE_FILE:-${SUPPORT_DIR}/eas.sample}"
EXCLUDE_COMPLEX_FILE="${EXCLUDE_COMPLEX_FILE:-${SUPPORT_DIR}/exclude_complex}"
MAP_FILE="${MAP_FILE:-${PROJECT_DIR}/01.split_vcf/sample_map.txt}"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-${BASE_WORK_DIR}/genotype_consistency.py}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
R2_THRESHOLD="${R2_THRESHOLD:-0.4}"
CONCORDANCE_THRESHOLD="${CONCORDANCE_THRESHOLD:-0.7}"
MAX_CHROMOSOME_JOBS="${MAX_CHROMOSOME_JOBS:-12}"
THREADS_PER_JOB="${THREADS_PER_JOB:-2}"
CHUNK_SIZE="${CHUNK_SIZE:-50000}"
CHROMOSOMES="${CHROMOSOMES:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24}"
KEEP_CONSISTENCY_TEMP="${KEEP_CONSISTENCY_TEMP:-0}"
GCTA_THREADS="${GCTA_THREADS:-20}"
SKIP_PLOTS="${SKIP_PLOTS:-0}"

for value in "${MAX_CHROMOSOME_JOBS}" "${THREADS_PER_JOB}" "${CHUNK_SIZE}" "${GCTA_THREADS}"; do
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: thread, job, and chunk parameters must be positive integers." >&2
        exit 2
    fi
done

for tool in plink plink2 gcta64 awk "${PYTHON_BIN}"; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: required command is not available in PATH: ${tool}" >&2
        exit 1
    }
done

mkdir -p "${BASE_WORK_DIR}"

require_file() {
    local file="$1"
    local description="$2"
    [[ -f "${file}" ]] || {
        echo "ERROR: ${description} not found: ${file}" >&2
        return 1
    }
}

have_bfile() {
    local prefix="$1"
    [[ -f "${prefix}.bed" && -f "${prefix}.bim" && -f "${prefix}.fam" ]]
}

# Run chromosome-level consistency jobs and merge their reports. This logic is
# kept here because consistency is part of this QC pipeline, not a separate
# user-facing step.
calculate_consistency() (
    set -euo pipefail

    local bfile_prefix="$1"
    local final_out="$2"
    local work_dir
    local tmp_dir
    local tmp_out
    local completed=0
    local failed=0
    local found=0
    local wrote_header=0
    local chr
    local result
    local active_pids=()
    local active_chrs=()

    work_dir="$(dirname "${bfile_prefix}")"

    for suffix in bed bim fam; do
        require_file "${bfile_prefix}.${suffix}" "PLINK consistency input"
    done
    require_file "${MAP_FILE}" "sample map"
    require_file "${PYTHON_SCRIPT}" "genotype-consistency program"

    tmp_dir="$(mktemp -d "${work_dir}/.consistency.XXXXXX")"
    tmp_out="${final_out}.tmp.$$"

    cleanup_consistency() {
        rm -f -- "${tmp_out}"
        if [[ "${completed}" -eq 1 && "${KEEP_CONSISTENCY_TEMP}" -eq 0 ]]; then
            rm -rf -- "${tmp_dir}"
        else
            echo "Consistency temporary files retained in: ${tmp_dir}" >&2
        fi
    }
    trap cleanup_consistency EXIT

    run_chromosome() {
        local chromosome="$1"
        local prefix="${tmp_dir}/chr${chromosome}"
        local chromosome_result="${tmp_dir}/res_chr${chromosome}.txt"
        local chromosome_log="${tmp_dir}/chr${chromosome}.log"

        {
            plink \
                --bfile "${bfile_prefix}" \
                --chr "${chromosome}" \
                --make-bed \
                --out "${prefix}" \
                --keep-allele-order \
                --noweb

            "${PYTHON_BIN}" "${PYTHON_SCRIPT}" \
                --bfile "${prefix}" \
                --map "${MAP_FILE}" \
                --out "${chromosome_result}" \
                --threads "${THREADS_PER_JOB}" \
                --chunk "${CHUNK_SIZE}"
        } >"${chromosome_log}" 2>&1

        [[ -s "${chromosome_result}" ]]
    }

    wait_for_first_chromosome() {
        local pid="${active_pids[0]}"
        local chromosome="${active_chrs[0]}"

        if ! wait "${pid}"; then
            echo "ERROR: chromosome ${chromosome} consistency failed; log follows:" >&2
            sed -n '1,160p' "${tmp_dir}/chr${chromosome}.log" >&2 || true
            failed=1
        else
            echo "[DONE] Consistency chromosome ${chromosome}"
        fi

        active_pids=("${active_pids[@]:1}")
        active_chrs=("${active_chrs[@]:1}")
    }

    echo "[START] Genotype consistency"

    for chr in ${CHROMOSOMES}; do
        if ! awk -v chr="${chr}" \
            '$1 == chr { found=1; exit } END { exit !found }' \
            "${bfile_prefix}.bim"; then
            continue
        fi

        found=$((found + 1))
        echo "[RUN] Consistency chromosome ${chr}"
        run_chromosome "${chr}" &
        active_pids+=("$!")
        active_chrs+=("${chr}")

        if (( ${#active_pids[@]} >= MAX_CHROMOSOME_JOBS )); then
            wait_for_first_chromosome
        fi
    done

    while (( ${#active_pids[@]} > 0 )); do
        wait_for_first_chromosome
    done

    if [[ "${found}" -eq 0 ]]; then
        echo "ERROR: no requested chromosomes were found in ${bfile_prefix}.bim." >&2
        exit 1
    fi

    if [[ "${failed}" -ne 0 ]]; then
        echo "ERROR: one or more chromosome consistency jobs failed." >&2
        exit 1
    fi

    for chr in ${CHROMOSOMES}; do
        result="${tmp_dir}/res_chr${chr}.txt"
        [[ -s "${result}" ]] || continue

        if [[ "${wrote_header}" -eq 0 ]]; then
            head -n 1 "${result}" > "${tmp_out}"
            wrote_header=1
        fi
        tail -n +2 "${result}" >> "${tmp_out}"
    done

    if [[ "${wrote_header}" -eq 0 ]]; then
        echo "ERROR: no chromosome consistency result was available for merging." >&2
        exit 1
    fi

    mv -f -- "${tmp_out}" "${final_out}"
    completed=1
    echo "[DONE] Combined consistency report: ${final_out}"
)

plot_pca_outputs() {
    if [[ "${SKIP_PLOTS}" -eq 1 ]]; then
        echo "[PLOT] Skipped because SKIP_PLOTS=1."
        return
    fi

    if ! "${PYTHON_BIN}" -c 'import pandas, plotly' >/dev/null 2>&1; then
        echo "[PLOT] pandas or plotly is unavailable; HTML plots were skipped." >&2
        return
    fi

    "${PYTHON_BIN}" <<'PY'
import os

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def smart_plot(file_path, out_name, title):
    if not os.path.exists(file_path):
        return

    preview = pd.read_csv(file_path, sep=r"\s+", nrows=5)
    column_text = " ".join(str(column).upper() for column in preview.columns)
    has_header = any(token in column_text for token in ("PC", "FID", "SCORE"))

    if has_header:
        frame = pd.read_csv(file_path, sep=r"\s+")
    else:
        names = ["FID", "IID"] + [f"PC{i}" for i in range(1, preview.shape[1] - 1)]
        frame = pd.read_csv(file_path, sep=r"\s+", names=names)

    columns = frame.columns.tolist()
    hover_id = columns[1] if len(columns) > 1 else columns[0]

    def find_pc(number):
        target = f"PC{number}"
        return next((column for column in columns if target in str(column).upper()), None)

    pairs = [(1, 2), (3, 4), (5, 6)]
    figure = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=[f"PC{first} vs PC{second}" for first, second in pairs],
        horizontal_spacing=0.05,
    )

    for index, (first, second) in enumerate(pairs, start=1):
        x_column, y_column = find_pc(first), find_pc(second)
        if x_column is None or y_column is None:
            continue

        figure.add_trace(
            go.Scatter(
                x=frame[x_column],
                y=frame[y_column],
                mode="markers",
                marker={"size": 5, "opacity": 0.6},
                text=frame[hover_id],
                hovertemplate=(
                    f"ID: %{{text}}<br>PC{first}: %{{x}}"
                    f"<br>PC{second}: %{{y}}<extra></extra>"
                ),
            ),
            row=1,
            col=index,
        )
        figure.update_xaxes(title_text=f"PC{first}", row=1, col=index)
        figure.update_yaxes(title_text=f"PC{second}", row=1, col=index)

    figure.update_layout(
        title_text=title,
        template="plotly_white",
        height=500,
        width=1500,
        showlegend=False,
    )
    figure.write_html(f"{out_name}.html")


smart_plot("pca.ngs_tgs.eigenvec", "plot_pca_ngs_tgs", "NGS+TGS PCA")
smart_plot("pca.ngs.eigenvec", "plot_pca_ngs", "NGS PCA")
smart_plot("pca.ngs.eas.eigenvec", "plot_pca_eas", "EAS PCA")
smart_plot(
    "pca.eas.proj.all.sample.clean.sscore",
    "plot_pca_projection",
    "EAS projection",
)
PY
}

run_qc_pipeline() {
    local target_dir="${BASE_WORK_DIR}/${DATA_DIR}"
    local input_vcf="${VCF_BASE}/${DATA_DIR}/step1.split.vcf.gz"
    local final_consistency="${target_dir}/consistency_report_final.txt"

    mkdir -p "${target_dir}"
    cd "${target_dir}"
    echo "[$(date +'%H:%M:%S')] [START] BFILE QC pipeline"

    if ! have_bfile "NGS_TGS"; then
        echo "[1/8] Converting VCF to PLINK BED files."
        require_file "${input_vcf}" "merged input VCF"
        plink \
            --vcf "${input_vcf}" \
            --make-bed \
            --out NGS_TGS \
            --const-fid \
            --vcf-half-call m \
            --keep-allele-order \
            --noweb
        awk 'BEGIN { OFS="\t" } { $2=$1 ":" $4 ":" NR; print }' \
            NGS_TGS.bim > NGS_TGS.bim.tmp
        mv -f NGS_TGS.bim.tmp NGS_TGS.bim
    fi

    if [[ ! -s "${final_consistency}" ]]; then
        echo "[2/8] Calculating genotype consistency."
        calculate_consistency "${target_dir}/NGS_TGS" "${final_consistency}"
    fi

    if ! have_bfile "NGS_TGS.QCsite"; then
        echo "[3/8] Filtering variants using consistency and excluded regions."
        awk \
            -v r2="${R2_THRESHOLD}" \
            -v concordance="${CONCORDANCE_THRESHOLD}" \
            'NR > 1 && $18 > r2 && $13 > concordance { print $1 }' \
            "${final_consistency}" > SNP_r2conc.txt

        [[ -s SNP_r2conc.txt ]] || {
            echo "ERROR: no variant passed the consistency thresholds." >&2
            return 1
        }

        plink \
            --bfile NGS_TGS \
            --extract SNP_r2conc.txt \
            --geno 0.5 \
            --make-bed \
            --out NGS_TGS.QCsite0 \
            --keep-allele-order \
            --noweb

        if [[ -f "${EXCLUDE_COMPLEX_FILE}" ]]; then
            plink \
                --bfile NGS_TGS.QCsite0 \
                --exclude range "${EXCLUDE_COMPLEX_FILE}" \
                --make-bed \
                --out NGS_TGS.QCsite \
                --keep-allele-order \
                --noweb
        else
            echo "WARNING: exclude_complex was not found; range exclusion skipped." >&2
            for suffix in bed bim fam; do
                cp -f "NGS_TGS.QCsite0.${suffix}" "NGS_TGS.QCsite.${suffix}"
            done
        fi
    fi

    if [[ ! -s pca.ngs_tgs.eigenvec ]]; then
        echo "[4/8] Calculating NGS+TGS PCA."
        plink --bfile NGS_TGS.QCsite --mind 0.2 --make-bed \
            --out NGS_TGS.QCsite.QCind --keep-allele-order --noweb
        plink --bfile NGS_TGS.QCsite.QCind --geno 0.2 --maf 0.05 \
            --indep-pairwise 50 5 0.1 --out prune.ngs_tgs \
            --keep-allele-order --noweb
        plink --bfile NGS_TGS.QCsite.QCind --extract prune.ngs_tgs.prune.in \
            --pca 20 --chr 1-22 --out pca.ngs_tgs --keep-allele-order --noweb
    fi

    if [[ ! -s pca.ngs.eigenvec ]]; then
        echo "[5/8] Creating the NGS-only QC dataset and PCA."
        require_file "${TGS_SAMPLE_FILE}" "TGS sample list"
        require_file "${BAD_SAMPLE_FILE}" "bad-sample list"
        plink --bfile NGS_TGS.QCsite --remove "${TGS_SAMPLE_FILE}" --make-bed \
            --out NGS.QCsite --keep-allele-order --noweb
        plink --bfile NGS.QCsite --mind 0.2 --make-bed \
            --out NGS.QCsite.QCind --keep-allele-order --noweb
        plink --bfile NGS.QCsite.QCind --geno 0.2 --maf 0.05 \
            --indep-pairwise 50 5 0.1 --out prune.ngs \
            --keep-allele-order --noweb
        plink --bfile NGS.QCsite.QCind --remove "${BAD_SAMPLE_FILE}" \
            --extract prune.ngs.prune.in --pca 20 --chr 1-22 --out pca.ngs \
            --keep-allele-order --noweb
    fi

    if [[ ! -s pca.ngs.eas.eigenvec ]]; then
        echo "[6/8] Calculating EAS PCA."
        require_file "${BAD_SAMPLE_FILE}" "bad-sample list"
        require_file "${EAS_SAMPLE_FILE}" "EAS sample list"
        plink --bfile NGS.QCsite.QCind --remove "${BAD_SAMPLE_FILE}" \
            --keep "${EAS_SAMPLE_FILE}" --geno 0.5 --maf 0.05 \
            --indep-pairwise 50 5 0.1 --out prune.ngs.eas \
            --keep-allele-order --noweb
        plink --bfile NGS.QCsite.QCind --remove "${BAD_SAMPLE_FILE}" \
            --keep "${EAS_SAMPLE_FILE}" --geno 0.2 \
            --extract prune.ngs.eas.prune.in --pca 10 --chr 1-22 \
            --out pca.ngs.eas --keep-allele-order --noweb
    fi

    if [[ ! -s pca.eas.proj.all.sample.clean.sscore ]]; then
        echo "[7/8] Calculating EAS projection PCA."
        require_file "${BAD_SAMPLE_FILE}" "bad-sample list"
        require_file "${EAS_SAMPLE_FILE}" "EAS sample list"
        plink2 --bfile NGS.QCsite.QCind --remove "${BAD_SAMPLE_FILE}" \
            --keep "${EAS_SAMPLE_FILE}" --geno 0.2 \
            --extract prune.ngs.eas.prune.in --pca allele-wts --freq counts \
            --out pca.ngs.eas.proj --keep-allele-order
        plink2 --bfile NGS.QCsite.QCind \
            --read-freq pca.ngs.eas.proj.acount \
            --score pca.ngs.eas.proj.eigenvec.allele 2 6 header-read \
            no-mean-imputation variance-standardize --score-col-nums 7-16 \
            --out pca.eas.proj.all.sample
        awk '{ printf "%s\t%s", $1, $2; for (i=4; i<=NF; i++) printf "\t%s", $i; print "" }' \
            pca.eas.proj.all.sample.sscore > pca.eas.proj.all.sample.clean.sscore
    fi

    if [[ ! -s sp_grm.grm.sp || ! -s grm.grm.bin ]]; then
        echo "[8/8] Calculating GCTA GRMs."
        plink --bfile NGS.QCsite.QCind --maf 0.05 --geno 0.1 \
            --indep-pairwise 50 10 0.6 --out NGS_for_grm \
            --keep-allele-order --noweb
        plink --bfile NGS.QCsite.QCind --extract NGS_for_grm.prune.in \
            --make-bed --out NGS_for_grm --keep-allele-order --noweb
        gcta64 --bfile NGS_for_grm --make-grm --sparse-cutoff 0.05 \
            --thread-num "${GCTA_THREADS}" --out sp_grm
        gcta64 --bfile NGS_for_grm --autosome --make-grm \
            --thread-num "${GCTA_THREADS}" --out grm
    fi

    plot_pca_outputs
    echo "[$(date +'%H:%M:%S')] [DONE] BFILE QC pipeline"
}

run_qc_pipeline
