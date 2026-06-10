#!/usr/bin/env bash
# Stage 18.3 physical bending-force operator wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_3_ENABLE:=1}"
: "${STAGE18_3_BENDING_OPERATOR_ENABLE:=1}"
: "${STAGE18_3_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_3_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_3_NPTS:=32}"
: "${STAGE18_3_FIBRE_LENGTH:=1.0}"
: "${STAGE18_3_COMPONENT_DIM:=3}"
: "${STAGE18_3_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_3_GAMMA:=1.0e-3}"
: "${STAGE18_3_USE_DIMENSIONAL_BENDING:=1}"
: "${STAGE18_3_USE_NONDIMENSIONAL_BENDING:=1}"
: "${STAGE18_3_SINE_EPS:=1.0e-3}"
: "${STAGE18_3_SINE_MODE:=1}"
: "${STAGE18_3_QUADRATIC_EPS:=1.0e-3}"
: "${STAGE18_3_ZERO_TOL:=1.0e-14}"
: "${STAGE18_3_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_3_ENERGY_TOL:=1.0e-12}"
: "${STAGE18_3_TEST_CASE:=straight_sine_quadratic_bending}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are interface-compatibility
# variables only. Stage 18.3 intentionally does not cd into DECOMP2D_ROOT, build,
# invoke MPI, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_3_physical_bending_force_operator.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-3-enable "${STAGE18_3_ENABLE}" \
  --bending-operator-enable "${STAGE18_3_BENDING_OPERATOR_ENABLE}" \
  --single-fibre-only "${STAGE18_3_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_3_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_3_NPTS}" \
  --fibre-length "${STAGE18_3_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_3_COMPONENT_DIM}" \
  --bending-stiffness "${STAGE18_3_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_3_GAMMA}" \
  --use-dimensional-bending "${STAGE18_3_USE_DIMENSIONAL_BENDING}" \
  --use-nondimensional-bending "${STAGE18_3_USE_NONDIMENSIONAL_BENDING}" \
  --sine-eps "${STAGE18_3_SINE_EPS}" \
  --sine-mode "${STAGE18_3_SINE_MODE}" \
  --quadratic-eps "${STAGE18_3_QUADRATIC_EPS}" \
  --zero-tol "${STAGE18_3_ZERO_TOL}" \
  --formula-tol "${STAGE18_3_FORMULA_TOL}" \
  --energy-tol "${STAGE18_3_ENERGY_TOL}" \
  --test-case "${STAGE18_3_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.3 PHYSICAL BENDING FORCE OPERATOR VERDICT: PASS"
  echo "STAGE 18.3 FINAL VERDICT: PASS"
else
  echo "STAGE 18.3 PHYSICAL BENDING FORCE OPERATOR VERDICT: FAIL"
  echo "STAGE 18.3 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
