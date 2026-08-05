set -euo pipefail

THREADS=32
ALIGN_THREADS=20
SORT_THREADS=8
SORT_MEM="8G"

OFFSET=${OFFSET:-0}
REAL_ID=$((SLURM_ARRAY_TASK_ID + OFFSET))

hg38_ref="GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF="CHM13v2.0_genomic_chr.fasta"

OUT_DIR="OUT_DIR"
REMOTE_DIR="REMOTE_DIR"

TMP_ROOT="tmp"

LIST="t2t_final_path.list"

mkdir -p "$OUT_DIR" task_file "$TMP_ROOT" "${OUT_DIR}/doneflag" "${OUT_DIR}/log" "${OUT_DIR}/checkpoint"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

die() {
    log "[ERROR] $*"
    exit 1
}

check_cram() {
    local cram="$1"
    [[ -s "$cram" ]] && [[ -s "${cram}.crai" ]] && samtools quickcheck "$cram" 2>/dev/null
}

check_gzip() {
    local file="$1"
    [[ -s "$file" ]] && gzip -t "$file" 2>/dev/null
}

mark_stage() {
    local stage_name="$1"
    printf '%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" > "${state_dir}/${stage_name}.ok"
}

clear_stage() {
    local stage_name="$1"
    rm -f "${state_dir}/${stage_name}.ok"
}

has_stage() {
    local stage_name="$1"
    [[ -f "${state_dir}/${stage_name}.ok" ]]
}

cleanup() {
    local exit_code=$?

    if [[ $exit_code -eq 0 ]]; then
        rm -rf "$tmp_dir"
    else
        log "[WARN] Job failed, keeping tmp dir for resume/debug: $tmp_dir"
    fi

    exit "$exit_code"
}

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "SLURM_ARRAY_TASK_ID is not set"
[[ -f "$LIST" ]] || die "Sample list not found: $LIST"

TOTAL=$(wc -l < "$LIST")

if [[ "$REAL_ID" -gt "$TOTAL" ]]; then
    log "[SKIP] REAL_ID $REAL_ID exceeds list"
    exit 0
fi

cram_file=$(sed -n "${REAL_ID}p" "$LIST")


[[ -n "${cram_file:-}" ]] || die "No sample found at line ${SLURM_ARRAY_TASK_ID} in $LIST"
[[ -f "$cram_file" ]] || die "Input CRAM not found: $cram_file"

sample=$(basename "$cram_file" .final.cram)
[[ -n "$sample" ]] || die "Failed to parse sample name from: $cram_file"

out_cram="${OUT_DIR}/${sample}.t2t.cram"
done_flag="${OUT_DIR}/doneflag/${sample}.done"
log_file="${OUT_DIR}/log/${sample}_${SLURM_ARRAY_TASK_ID}.resume.log"
remote_cram="${REMOTE_DIR}/${sample}.t2t.cram"
tmp_dir="${TMP_ROOT}/${sample}_${SLURM_ARRAY_TASK_ID}"
state_dir="${OUT_DIR}/checkpoint/${sample}"
collated_bam="${tmp_dir}/collated.bam"
r1_fq="${tmp_dir}/R1.fq.gz"
r2_fq="${tmp_dir}/R2.fq.gz"

mkdir -p "$tmp_dir" "$state_dir"
trap cleanup EXIT

exec >"$log_file" 2>&1

log "========================================"
log "[START] Sample: $sample, ArrayID: $SLURM_ARRAY_TASK_ID"
log "========================================"

module load miniconda3
source activate bwa
module load samtools

for cmd in samtools bwa-mem2 gzip; do
    command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
done

if [[ -f "$remote_cram" ]]; then
    log "[SKIP] Remote CRAM exists: $remote_cram"
    touch "$done_flag"
    mark_stage "remote_exists"
    exit 0
fi

if check_cram "$out_cram"; then
    log "[SKIP] Local CRAM is valid: $out_cram"
    touch "$done_flag"
    mark_stage "fastq"
    mark_stage "align"
    mark_stage "index"
    exit 0
fi

if [[ -f "$done_flag" ]]; then
    log "[WARN] Found done flag without a valid local CRAM, removing stale flag"
    rm -f "$done_flag"
fi

if has_stage "fastq"; then
    if check_gzip "$r1_fq" && check_gzip "$r2_fq"; then
        log "[RESUME] FASTQ checkpoint found, reusing existing FASTQ files"
    else
        log "[WARN] FASTQ checkpoint exists but FASTQ files are invalid, rebuilding"
        clear_stage "fastq"
    fi
fi

if has_stage "align"; then
    if [[ -s "$out_cram" ]]; then
        log "[RESUME] Alignment checkpoint found, output CRAM exists"
    else
        log "[WARN] Alignment checkpoint exists but CRAM is missing, rerunning alignment"
        clear_stage "align"
        clear_stage "index"
    fi
fi

if has_stage "index"; then
    if check_cram "$out_cram"; then
        log "[RESUME] Indexed CRAM is already valid"
        touch "$done_flag"
        exit 0
    else
        log "[WARN] Index checkpoint exists but CRAM validation failed, rebuilding index/alignment"
        clear_stage "index"
        rm -f "${out_cram}.crai"
    fi
fi

start_time=$(date +%s)

if has_stage "fastq" && check_gzip "$r1_fq" && check_gzip "$r2_fq"; then
    log "[STEP 1] FASTQ stage already complete"
else
    log "[STEP 1] Processing FASTQ extraction"
    clear_stage "fastq"
    clear_stage "align"
    clear_stage "index"

    if check_gzip "$r1_fq" && check_gzip "$r2_fq"; then
        log "[STEP 1] Reusing existing FASTQ files"
    else
        rm -f "$r1_fq" "$r2_fq" "$collated_bam"

        log "[STEP 1] Running samtools collate"
        samtools collate \
            -@ "$THREADS" \
            -u \
            -o "$collated_bam" \
            --reference "$hg38_ref" \
            "$cram_file" \
            "${tmp_dir}/collate_tmp"

        [[ -s "$collated_bam" ]] || die "collated.bam was not created"

        log "[STEP 1] Running samtools fastq"
        samtools fastq \
            -@ "$THREADS" \
            --reference "$hg38_ref" \
            -F 0x900 \
            -n \
            -1 "$r1_fq" \
            -2 "$r2_fq" \
            -0 /dev/null \
            -s /dev/null \
            "$collated_bam"

        check_gzip "$r1_fq" || die "R1 FASTQ is invalid: $r1_fq"
        check_gzip "$r2_fq" || die "R2 FASTQ is invalid: $r2_fq"
        rm -f "$collated_bam"
    fi

    mark_stage "fastq"
    log "[STEP 1] FASTQ ready: R1=$(du -h "$r1_fq" | cut -f1), R2=$(du -h "$r2_fq" | cut -f1)"
fi

if has_stage "align" && [[ -s "$out_cram" ]]; then
    log "[STEP 2] Alignment stage already complete"
else
    log "[STEP 2] Starting alignment"
    clear_stage "align"
    clear_stage "index"
    rm -f "$out_cram" "${out_cram}.crai"

    rg_header="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:WGS"

    bwa-mem2 mem \
        -t "$ALIGN_THREADS" \
        -K 50000000 \
        -R "$rg_header" \
        "$REF" \
        "$r1_fq" \
        "$r2_fq" | \
        samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEM" \
            --reference "$REF" \
            -T "${tmp_dir}/st_sort" \
            -O CRAM \
            -o "$out_cram"

    [[ -s "$out_cram" ]] || die "Alignment output CRAM was not created"
    mark_stage "align"
    log "[STEP 2] Alignment complete: $(du -h "$out_cram" | cut -f1)"
fi

if has_stage "index" && check_cram "$out_cram"; then
    log "[STEP 3] Index stage already complete"
else
    log "[STEP 3] Creating index"
    clear_stage "index"
    rm -f "${out_cram}.crai"
    samtools index -@ "$THREADS" "$out_cram"

    if ! check_cram "$out_cram"; then
        rm -f "$out_cram" "${out_cram}.crai"
        clear_stage "align"
        die "CRAM validation failed after indexing"
    fi

    mark_stage "index"
fi

touch "$done_flag"

end_time=$(date +%s)
duration=$((end_time - start_time))
log "[SUCCESS] Completed in $((duration / 60)) minutes"
log "[SUCCESS] Checkpoints saved in: $state_dir"
log "========================================"
