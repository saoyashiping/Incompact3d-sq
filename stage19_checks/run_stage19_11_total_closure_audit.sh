#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
export DECOMP2D_ROOT

: "${STAGE19_11_ENABLE:=1}"
: "${STAGE19_11_TOTAL_CLOSURE_AUDIT_ENABLE:=1}"
export STAGE19_11_ENABLE STAGE19_11_TOTAL_CLOSURE_AUDIT_ENABLE

mkdir -p "${REPO_ROOT}/stage19_outputs"
OUTPUT="${REPO_ROOT}/stage19_outputs/fibre_stage19_11_total_closure_audit.dat"
CLOSED_FILE="${REPO_ROOT}/stage19_checks/STAGE19_CLOSED.md"
python3 "${REPO_ROOT}/stage19_checks/assert_stage19_11_total_closure_audit.py" \
  --repo-root "${REPO_ROOT}" \
  --output "${OUTPUT}" \
  --closed-file "${CLOSED_FILE}"
