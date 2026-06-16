#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)
# DECOMP2D_ROOT is intentionally not used to change directories; it is accepted
# only as a compatibility variable for environments that export it.
: "${DECOMP2D_ROOT:=}"

mkdir -p "${REPO_ROOT}/stage20_outputs"
python3 "${SCRIPT_DIR}/assert_stage20_0_preflight_boundary.py"
