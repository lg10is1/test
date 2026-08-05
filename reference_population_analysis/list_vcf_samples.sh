#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: list_vcf_samples.sh INPUT.vcf[.gz] OUTPUT.txt

Write the VCF sample names, one per line, using bcftools query -l.
EOF
}

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 2
fi

input_vcf=$1
output_file=$2

command -v bcftools >/dev/null 2>&1 || {
    echo "Required command not found: bcftools" >&2
    exit 127
}
[[ -r "$input_vcf" ]] || {
    echo "Input VCF is not readable: $input_vcf" >&2
    exit 1
}

output_parent=$(dirname -- "$output_file")
[[ -d "$output_parent" ]] || {
    echo "Output directory does not exist: $output_parent" >&2
    exit 1
}

temporary_file=$(mktemp "${output_file}.tmp.XXXXXX")
cleanup() {
    rm -f -- "$temporary_file"
}
trap cleanup EXIT

bcftools query -l "$input_vcf" > "$temporary_file"
mv -f -- "$temporary_file" "$output_file"
trap - EXIT
printf 'Sample list: %s\n' "$output_file"
