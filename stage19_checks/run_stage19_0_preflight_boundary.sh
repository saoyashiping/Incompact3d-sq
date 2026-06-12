#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Interface compatibility only.  Stage 19.0 must not cd into DECOMP2D_ROOT.
: "${DECOMP2D_ROOT:=${REPO_ROOT}}"

: "${STAGE19_0_ENABLE:=1}"
: "${STAGE19_0_PREFLIGHT_ENABLE:=1}"
: "${STAGE19_0_REQUIRE_STAGE18_CLOSED:=1}"
: "${STAGE19_0_REQUIRE_STAGE18_12_PASS:=1}"
: "${STAGE19_0_REQUIRE_STAGE18_OUTPUTS:=1}"
: "${STAGE19_0_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE19_0_DIAGNOSTIC_ONLY:=1}"
: "${STAGE19_0_RERUN_PRIOR_STAGES:=0}"
: "${STAGE19_0_ZERO_TOL:=1.0e-14}"
: "${STAGE19_0_AUDIT_TOL:=1.0e-12}"
: "${STAGE19_0_TEST_CASE:=stage18_closure_stage19_preflight_boundary}"

export STAGE19_0_ENABLE STAGE19_0_PREFLIGHT_ENABLE
export STAGE19_0_REQUIRE_STAGE18_CLOSED STAGE19_0_REQUIRE_STAGE18_12_PASS STAGE19_0_REQUIRE_STAGE18_OUTPUTS
export STAGE19_0_SINGLE_FIBRE_ONLY STAGE19_0_DIAGNOSTIC_ONLY STAGE19_0_RERUN_PRIOR_STAGES
export STAGE19_0_ZERO_TOL STAGE19_0_AUDIT_TOL STAGE19_0_TEST_CASE

mkdir -p "${REPO_ROOT}/stage19_outputs"
OUTPUT="${REPO_ROOT}/stage19_outputs/fibre_stage19_0_preflight_boundary.dat"

python3 "${SCRIPT_DIR}/assert_stage19_0_preflight_boundary.py" \
  --repo-root "${REPO_ROOT}" \
  --output "${OUTPUT}"
