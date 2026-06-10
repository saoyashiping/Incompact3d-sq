#!/usr/bin/env bash
# Stage 18.2 structure-state and geometry-operator wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_2_ENABLE:=1}"
: "${STAGE18_2_STRUCTURE_STATE_GEOMETRY_ENABLE:=1}"
: "${STAGE18_2_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_2_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_2_NPTS:=16}"
: "${STAGE18_2_FIBRE_LENGTH:=1.0}"
: "${STAGE18_2_COMPONENT_DIM:=3}"
: "${STAGE18_2_SINE_EPS:=1.0e-3}"
: "${STAGE18_2_SINE_MODE:=1}"
: "${STAGE18_2_ZERO_TOL:=1.0e-14}"
: "${STAGE18_2_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_2_DERIVATIVE_TOL:=5.0e-2}"
: "${STAGE18_2_ARC_ERROR_TOL:=1.0e-2}"
: "${STAGE18_2_TEST_CASE:=straight_and_sine_geometry}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are interface-compatibility
# variables only. Stage 18.2 intentionally does not cd into DECOMP2D_ROOT, build,
# invoke MPI, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_2_structure_state_geometry_operators.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-2-enable "${STAGE18_2_ENABLE}" \
  --structure-state-geometry-enable "${STAGE18_2_STRUCTURE_STATE_GEOMETRY_ENABLE}" \
  --single-fibre-only "${STAGE18_2_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_2_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_2_NPTS}" \
  --fibre-length "${STAGE18_2_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_2_COMPONENT_DIM}" \
  --sine-eps "${STAGE18_2_SINE_EPS}" \
  --sine-mode "${STAGE18_2_SINE_MODE}" \
  --zero-tol "${STAGE18_2_ZERO_TOL}" \
  --formula-tol "${STAGE18_2_FORMULA_TOL}" \
  --derivative-tol "${STAGE18_2_DERIVATIVE_TOL}" \
  --arc-error-tol "${STAGE18_2_ARC_ERROR_TOL}" \
  --test-case "${STAGE18_2_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.2 STRUCTURE STATE GEOMETRY OPERATORS VERDICT: PASS"
  echo "STAGE 18.2 FINAL VERDICT: PASS"
else
  echo "STAGE 18.2 STRUCTURE STATE GEOMETRY OPERATORS VERDICT: FAIL"
  echo "STAGE 18.2 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
