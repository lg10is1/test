#!/usr/bin/env bash
set -euo pipefail

## Public full workflow runner.
## Edit placeholder paths in scripts/*.R, scripts/*.sh, and scripts/*.py before running.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}/scripts"

bash run_pipeline.sh "$@"
