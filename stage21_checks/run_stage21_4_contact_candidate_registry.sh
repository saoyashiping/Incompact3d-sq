#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
# DECOMP2D_ROOT is compatibility-only; this wrapper never cd's into it.
export DECOMP2D_ROOT
mkdir -p "${REPO_ROOT}/stage21_outputs"
python3 "${SCRIPT_DIR}/assert_stage21_4_contact_candidate_registry.py" "$@"
