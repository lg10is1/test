#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  download_igsr_ont.sh --workdir DIR [options]

Options:
  --workdir DIR       Destination for downloaded files and logs (required).
  --url-list FILE     HTTPS URL manifest (default: file_list.txt beside script).
  --md5-file FILE     md5sum manifest (default: md5_checksums.txt beside script).
  --jobs N            Concurrent files (default: 12).
  --connections N     aria2 connections per file (default: 8).
  --max-rounds N      Maximum checksum/download rounds (default: 3).
  --max-tries N       aria2 attempts per invocation (default: 5).
  --execute           Perform downloads. Without this flag, only validate and plan.
  -h, --help          Show this help.

The retained manifest contains about one thousand public CRAM/CRAI pairs and can
require substantial bandwidth and storage. Review the plan before using --execute.
EOF
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
workdir=""
url_list="${script_dir}/file_list.txt"
md5_file="${script_dir}/md5_checksums.txt"
jobs=12
connections=8
max_rounds=3
max_tries=5
execute=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --workdir)
            [[ $# -ge 2 ]] || { echo "Missing value for --workdir" >&2; exit 2; }
            workdir=$2
            shift 2
            ;;
        --url-list)
            [[ $# -ge 2 ]] || { echo "Missing value for --url-list" >&2; exit 2; }
            url_list=$2
            shift 2
            ;;
        --md5-file)
            [[ $# -ge 2 ]] || { echo "Missing value for --md5-file" >&2; exit 2; }
            md5_file=$2
            shift 2
            ;;
        --jobs)
            [[ $# -ge 2 ]] || { echo "Missing value for --jobs" >&2; exit 2; }
            jobs=$2
            shift 2
            ;;
        --connections)
            [[ $# -ge 2 ]] || { echo "Missing value for --connections" >&2; exit 2; }
            connections=$2
            shift 2
            ;;
        --max-rounds)
            [[ $# -ge 2 ]] || { echo "Missing value for --max-rounds" >&2; exit 2; }
            max_rounds=$2
            shift 2
            ;;
        --max-tries)
            [[ $# -ge 2 ]] || { echo "Missing value for --max-tries" >&2; exit 2; }
            max_tries=$2
            shift 2
            ;;
        --execute)
            execute=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n "$workdir" ]] || { echo "--workdir is required" >&2; exit 2; }
[[ -r "$url_list" ]] || { echo "URL list is not readable: $url_list" >&2; exit 1; }
[[ -r "$md5_file" ]] || { echo "MD5 manifest is not readable: $md5_file" >&2; exit 1; }
for value_name in jobs connections max_rounds max_tries; do
    value=${!value_name}
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
        echo "--${value_name//_/-} must be a positive integer" >&2
        exit 2
    }
done

url_count=$(awk 'NF {count++} END {print count+0}' "$url_list")
md5_count=$(awk 'NF {count++} END {print count+0}' "$md5_file")
[[ "$url_count" -gt 0 ]] || { echo "URL list is empty" >&2; exit 1; }
[[ "$url_count" -eq "$md5_count" ]] || {
    echo "Manifest counts differ: URLs=$url_count MD5=$md5_count" >&2
    exit 1
}

bad_url=$(awk 'NF && $0 !~ /^https:\/\/ftp\.1000genomes\.ebi\.ac\.uk\// {print; exit}' "$url_list")
[[ -z "$bad_url" ]] || { echo "Unexpected URL outside the approved IGSR HTTPS host: $bad_url" >&2; exit 1; }
bad_md5=$(awk 'NF && $0 !~ /^[[:xdigit:]]{32}[[:space:]]+[^[:space:]]+$/ {print; exit}' "$md5_file")
[[ -z "$bad_md5" ]] || { echo "Malformed MD5 manifest row: $bad_md5" >&2; exit 1; }

printf 'IGSR ONT download plan\n'
printf '  Destination: %s\n' "$workdir"
printf '  Files: %s\n' "$url_count"
printf '  Parallel files: %s\n' "$jobs"
printf '  Connections per file: %s\n' "$connections"
printf '  Maximum rounds: %s\n' "$max_rounds"

if [[ "$execute" != true ]]; then
    echo "Validation complete; no network request was made. Add --execute to download."
    exit 0
fi

for command_name in aria2c md5sum xargs awk grep; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "Required command not found: $command_name" >&2
        exit 127
    }
done

mkdir -p -- "$workdir" "$workdir/logs"
workdir=$(cd -- "$workdir" && pwd)
result_file="${workdir}/md5_check_results.txt"
ok_file="${workdir}/ok_files.txt"
failed_file="${workdir}/failed_download_list.txt"
retry_file="${workdir}/download_failed_again.txt"
touch "$result_file"
cd -- "$workdir"

extract_ok_files() {
    awk -F': ' '$2 == "OK" {print $1}' "$result_file" | sort -u > "$ok_file"
}

generate_failed_list() {
    : > "$failed_file"
    while IFS= read -r url; do
        [[ -n "$url" ]] || continue
        filename=$(basename -- "$url")
        if ! grep -Fxq -- "$filename" "$ok_file"; then
            printf '%s\n' "$url" >> "$failed_file"
        fi
    done < "$url_list"
}

download_one() {
    local url=$1
    local filename md5_line
    filename=$(basename -- "$url")
    if [[ ! "$filename" =~ ^[A-Za-z0-9._-]+$ ]]; then
        printf 'Rejected unsafe filename: %s\n' "$filename" >&2
        printf '%s\n' "$url" >> "$retry_file"
        return 1
    fi

    printf 'Downloading %s\n' "$filename"
    if ! aria2c \
        -x "$connections" \
        -s "$connections" \
        -k 4M \
        -c \
        --file-allocation=none \
        --console-log-level=warn \
        --retry-wait=10 \
        --max-tries="$max_tries" \
        "$url" > "logs/${filename}.log" 2>&1; then
        printf '%s\n' "$url" >> "$retry_file"
        return 1
    fi

    md5_line=$(awk -v name="$filename" '$2 == name {print; exit}' "$md5_file")
    if [[ -z "$md5_line" ]]; then
        printf 'No MD5 entry: %s\n' "$filename" >&2
        printf '%s\n' "$url" >> "$retry_file"
        return 1
    fi
    if printf '%s\n' "$md5_line" | md5sum -c - > "logs/${filename}.md5" 2>&1; then
        printf '%s: OK\n' "$filename" >> "$result_file"
        printf 'Verified %s\n' "$filename"
    else
        printf 'MD5 failed: %s\n' "$filename" >&2
        rm -f -- "$filename"
        printf '%s\n' "$url" >> "$retry_file"
        return 1
    fi
}

export -f download_one
export connections max_tries md5_file result_file retry_file

for ((round=1; round<=max_rounds; round++)); do
    printf '\nRound %s of %s\n' "$round" "$max_rounds"
    extract_ok_files
    generate_failed_list
    need=$(awk 'NF {count++} END {print count+0}' "$failed_file")
    if [[ "$need" -eq 0 ]]; then
        echo "All files are present and match the MD5 manifest."
        exit 0
    fi
    printf 'Files requiring download or verification: %s\n' "$need"
    : > "$retry_file"
    xargs -r -P "$jobs" -I '{}' bash -c 'download_one "$1"' _ '{}' < "$failed_file" || true
done

extract_ok_files
generate_failed_list
remaining=$(awk 'NF {count++} END {print count+0}' "$failed_file")
if [[ "$remaining" -gt 0 ]]; then
    printf 'Download incomplete after %s rounds; remaining files: %s\n' "$max_rounds" "$remaining" >&2
    exit 1
fi
echo "All files are present and match the MD5 manifest."
