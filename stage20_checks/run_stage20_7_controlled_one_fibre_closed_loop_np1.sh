#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
# DECOMP2D_ROOT is accepted only as a compatibility variable; this wrapper does not cd into it.
export DECOMP2D_ROOT
mkdir -p "${REPO_ROOT}/stage20_outputs"
python3 "${SCRIPT_DIR}/assert_stage20_7_controlled_one_fibre_closed_loop_np1.py" "$@"
