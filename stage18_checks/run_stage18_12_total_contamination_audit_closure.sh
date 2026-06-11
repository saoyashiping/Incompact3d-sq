#!/usr/bin/env bash
# Stage 18.12 total contamination audit and closure wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_12_ENABLE:=1}"
: "${STAGE18_12_TOTAL_AUDIT_ENABLE:=1}"
: "${STAGE18_12_CLOSURE_ENABLE:=1}"
: "${STAGE18_12_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_12_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_12_REQUIRE_PRIOR_OUTPUTS:=1}"
: "${STAGE18_12_RERUN_PRIOR_STAGES:=0}"
: "${STAGE18_12_REQUIRE_STAGE18_11:=1}"
: "${STAGE18_12_WRITE_CLOSURE_FILE:=1}"
: "${STAGE18_12_ZERO_TOL:=1.0e-14}"
: "${STAGE18_12_AUDIT_TOL:=1.0e-12}"
: "${STAGE18_12_TEST_CASE:=stage18_total_audit_closure}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# Compatibility variables only: Stage 18.12 does not cd into DECOMP2D_ROOT,
# build, invoke MPI, launch production DNS, modify production Fortran, or write
# production restart/statistics/visualisation/RHS/IBM/DNS-core outputs.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_12_total_contamination_audit_closure.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-12-enable "${STAGE18_12_ENABLE}" \
  --total-audit-enable "${STAGE18_12_TOTAL_AUDIT_ENABLE}" \
  --closure-enable "${STAGE18_12_CLOSURE_ENABLE}" \
  --single-fibre-only "${STAGE18_12_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_12_DIAGNOSTIC_ONLY}" \
  --require-prior-outputs "${STAGE18_12_REQUIRE_PRIOR_OUTPUTS}" \
  --rerun-prior-stages "${STAGE18_12_RERUN_PRIOR_STAGES}" \
  --require-stage18-11 "${STAGE18_12_REQUIRE_STAGE18_11}" \
  --write-closure-file "${STAGE18_12_WRITE_CLOSURE_FILE}" \
  --zero-tol "${STAGE18_12_ZERO_TOL}" \
  --audit-tol "${STAGE18_12_AUDIT_TOL}" \
  --test-case "${STAGE18_12_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.12 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: PASS"
  echo "STAGE 18.12 FINAL VERDICT: PASS"
else
  echo "STAGE 18.12 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: FAIL"
  echo "STAGE 18.12 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
