# PRIVACY WARNING: outputs contain individual sample identifiers.
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

python3 "$(dirname "$0")/extract_snv_indel_nonref_samples.py" \
  --plink "${PLINK_BIN:-plink}" \
  "$@"
