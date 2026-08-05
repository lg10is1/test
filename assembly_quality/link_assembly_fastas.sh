#!/usr/bin/env bash
set -euo pipefail
: "${ASSEMBLY_SOURCE_DIR:?Set ASSEMBLY_SOURCE_DIR}"
: "${ASSEMBLY_LINK_DIR:?Set ASSEMBLY_LINK_DIR}"
mkdir -p "$ASSEMBLY_LINK_DIR"

shopt -s nullglob
for file in "$ASSEMBLY_SOURCE_DIR"/*.fa; do
  filename=$(basename "$file")
  new_filename=$(printf '%s' "$filename" | sed -E 's/-(1|2|3|4)//g')
  ln -sfn "$file" "$ASSEMBLY_LINK_DIR/$new_filename"
done
