#!/usr/bin/env bash
set -euo pipefail

# Required environment variables:
#   SOURCE_REFERENCE   GRCh38 reference used by the input CRAM files.
#   TARGET_REFERENCE   CHM13 reference for realignment.
#   INPUT_LIST         One local input CRAM path per line.
#   OUT_DIR            Output directory.
# Optional variables:
#   TMP_ROOT, REMOTE_DIR, THREADS, ALIGN_THREADS, SORT_THREADS, SORT_MEM, OFFSET.
#
# The Slurm array index is treated as one-based. For a zero-based array, set
# OFFSET=1. Activate the required software environment before submitting.

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
: "${SOURCE_REFERENCE:?SOURCE_REFERENCE is required}"
: "${TARGET_REFERENCE:?TARGET_REFERENCE is required}"
: "${INPUT_LIST:?INPUT_LIST is required}"
: "${OUT_DIR:?OUT_DIR is required}"

THREADS=${THREADS:-32}
ALIGN_THREADS=${ALIGN_THREADS:-20}
SORT_THREADS=${SORT_THREADS:-8}
SORT_MEM=${SORT_MEM:-8G}
OFFSET=${OFFSET:-0}
TMP_ROOT=${TMP_ROOT:-tmp}
REMOTE_DIR=${REMOTE_DIR:-}

for value_name in SLURM_ARRAY_TASK_ID THREADS ALIGN_THREADS SORT_THREADS; do
    value=${!value_name}
    [[ "$value" =~ ^[0-9]+$ ]] || die "$value_name must be a non-negative integer"
done
[[ "$OFFSET" =~ ^-?[0-9]+$ ]] || die "OFFSET must be an integer"
[[ -r "$SOURCE_REFERENCE" ]] || die "Source reference is not readable: $SOURCE_REFERENCE"
[[ -r "$TARGET_REFERENCE" ]] || die "Target reference is not readable: $TARGET_REFERENCE"
[[ -r "$INPUT_LIST" ]] || die "Input list is not readable: $INPUT_LIST"
[[ -n "$TMP_ROOT" && "$TMP_ROOT" != "/" ]] || die "TMP_ROOT must not be empty or the filesystem root"

for command_name in samtools bwa-mem2 gzip; do
    command -v "$command_name" >/dev/null 2>&1 || die "Required command not found: $command_name"
done

line_number=$((SLURM_ARRAY_TASK_ID + OFFSET))
[[ "$line_number" -ge 1 ]] || die "Resolved input-list line must be at least 1; got $line_number"
total=$(awk 'END {print NR}' "$INPUT_LIST")
if [[ "$line_number" -gt "$total" ]]; then
    log "SKIP: line $line_number exceeds the $total-line input list"
    exit 0
fi

cram_file=$(sed -n "${line_number}p" "$INPUT_LIST")
[[ -n "$cram_file" ]] || die "No CRAM path at line $line_number"
[[ -r "$cram_file" ]] || die "Input CRAM is not readable: $cram_file"

sample=$(basename -- "$cram_file")
sample=${sample%.cram}
sample=${sample%.final}
[[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || die "Unsafe or empty sample label derived from: $cram_file"

mkdir -p -- "$OUT_DIR" "$TMP_ROOT" "$OUT_DIR/doneflag" "$OUT_DIR/log" "$OUT_DIR/checkpoint"
tmp_root_abs=$(cd -- "$TMP_ROOT" && pwd -P)
[[ -n "$tmp_root_abs" && "$tmp_root_abs" != "/" ]] || die "Resolved TMP_ROOT must not be the filesystem root"
out_cram="${OUT_DIR}/${sample}.t2t.cram"
done_flag="${OUT_DIR}/doneflag/${sample}.done"
log_file="${OUT_DIR}/log/${sample}_${SLURM_ARRAY_TASK_ID}.resume.log"
tmp_dir="${tmp_root_abs}/${sample}_${SLURM_ARRAY_TASK_ID}"
state_dir="${OUT_DIR}/checkpoint/${sample}"
collated_bam="${tmp_dir}/collated.bam"
r1_fq="${tmp_dir}/R1.fq.gz"
r2_fq="${tmp_dir}/R2.fq.gz"
remote_cram=""
if [[ -n "$REMOTE_DIR" ]]; then
    remote_cram="${REMOTE_DIR%/}/${sample}.t2t.cram"
fi
mkdir -p -- "$tmp_dir" "$state_dir"

check_cram() {
    local cram=$1
    [[ -s "$cram" && -s "${cram}.crai" ]] && samtools quickcheck "$cram" 2>/dev/null
}

check_gzip() {
    local file=$1
    [[ -s "$file" ]] && gzip -t "$file" 2>/dev/null
}

mark_stage() {
    local stage_name=$1
    printf '%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" > "${state_dir}/${stage_name}.ok"
}

clear_stage() {
    rm -f -- "${state_dir}/$1.ok"
}

has_stage() {
    [[ -f "${state_dir}/$1.ok" ]]
}

cleanup() {
    local exit_code=$?
    if [[ "$exit_code" -eq 0 ]]; then
        case "$tmp_dir" in
            "${tmp_root_abs}/"*) rm -rf -- "$tmp_dir" ;;
            *) log "WARN: refusing to remove unexpected temporary path: $tmp_dir" ;;
        esac
    else
        log "WARN: job failed; retaining temporary files for review: $tmp_dir"
    fi
    exit "$exit_code"
}
trap cleanup EXIT

exec > "$log_file" 2>&1
log "START: sample=$sample array_id=$SLURM_ARRAY_TASK_ID list_line=$line_number"

if [[ -n "$remote_cram" && -s "$remote_cram" ]]; then
    log "SKIP: remote CRAM exists: $remote_cram"
    touch "$done_flag"
    mark_stage remote_exists
    exit 0
fi

if check_cram "$out_cram"; then
    log "SKIP: local CRAM is already valid: $out_cram"
    touch "$done_flag"
    mark_stage fastq
    mark_stage align
    mark_stage index
    exit 0
fi

if [[ -f "$done_flag" ]]; then
    log "WARN: removing stale done flag because the local CRAM is not valid"
    rm -f -- "$done_flag"
fi

if has_stage fastq && ! { check_gzip "$r1_fq" && check_gzip "$r2_fq"; }; then
    log "WARN: invalid FASTQ checkpoint; rebuilding"
    clear_stage fastq
fi
if has_stage align && [[ ! -s "$out_cram" ]]; then
    log "WARN: alignment checkpoint has no output; rebuilding"
    clear_stage align
    clear_stage index
fi
if has_stage index && ! check_cram "$out_cram"; then
    log "WARN: index checkpoint is invalid; rebuilding"
    clear_stage index
    rm -f -- "${out_cram}.crai"
fi

start_time=$(date +%s)
if has_stage fastq && check_gzip "$r1_fq" && check_gzip "$r2_fq"; then
    log "STEP 1: reusing checkpointed FASTQ files"
else
    clear_stage fastq
    clear_stage align
    clear_stage index
    rm -f -- "$r1_fq" "$r2_fq" "$collated_bam"

    log "STEP 1: collating input CRAM"
    samtools collate \
        -@ "$THREADS" \
        -u \
        -o "$collated_bam" \
        --reference "$SOURCE_REFERENCE" \
        "$cram_file" \
        "${tmp_dir}/collate_tmp"
    [[ -s "$collated_bam" ]] || die "Collated BAM was not created"

    log "STEP 1: extracting paired FASTQ"
    samtools fastq \
        -@ "$THREADS" \
        --reference "$SOURCE_REFERENCE" \
        -F 0x900 \
        -n \
        -1 "$r1_fq" \
        -2 "$r2_fq" \
        -0 /dev/null \
        -s /dev/null \
        "$collated_bam"
    check_gzip "$r1_fq" || die "Invalid R1 FASTQ: $r1_fq"
    check_gzip "$r2_fq" || die "Invalid R2 FASTQ: $r2_fq"
    rm -f -- "$collated_bam"
    mark_stage fastq
fi

if has_stage align && [[ -s "$out_cram" ]]; then
    log "STEP 2: reusing checkpointed alignment"
else
    clear_stage align
    clear_stage index
    rm -f -- "$out_cram" "${out_cram}.crai"
    read_group="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:WGS"
    log "STEP 2: aligning reads to the target reference"
    bwa-mem2 mem \
        -t "$ALIGN_THREADS" \
        -K 50000000 \
        -R "$read_group" \
        "$TARGET_REFERENCE" \
        "$r1_fq" \
        "$r2_fq" | \
        samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEM" \
            --reference "$TARGET_REFERENCE" \
            -T "${tmp_dir}/st_sort" \
            -O CRAM \
            -o "$out_cram"
    [[ -s "$out_cram" ]] || die "Alignment CRAM was not created"
    mark_stage align
fi

if has_stage index && check_cram "$out_cram"; then
    log "STEP 3: reusing valid index"
else
    clear_stage index
    rm -f -- "${out_cram}.crai"
    samtools index -@ "$THREADS" "$out_cram"
    if ! check_cram "$out_cram"; then
        rm -f -- "$out_cram" "${out_cram}.crai"
        clear_stage align
        die "CRAM validation failed after indexing"
    fi
    mark_stage index
fi

touch "$done_flag"
end_time=$(date +%s)
duration=$((end_time - start_time))
log "SUCCESS: completed in $((duration / 60)) minutes"
log "Checkpoints: $state_dir"
