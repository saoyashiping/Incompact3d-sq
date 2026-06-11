#!/usr/bin/env bash
# Stage 18.4 standalone tension/inextensibility diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_4_ENABLE:=1}"
: "${STAGE18_4_TENSION_CONSTRAINT_ENABLE:=1}"
: "${STAGE18_4_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_4_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_4_NPTS:=64}"
: "${STAGE18_4_FIBRE_LENGTH:=1.0}"
: "${STAGE18_4_COMPONENT_DIM:=3}"
: "${STAGE18_4_RHO_L:=1.0}"
: "${STAGE18_4_RHO_TILDE:=1.0}"
: "${STAGE18_4_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_4_GAMMA:=1.0e-3}"
: "${STAGE18_4_VELOCITY_EPS:=1.0e-3}"
: "${STAGE18_4_ARCLENGTH_STRETCH_EPS:=1.0e-3}"
: "${STAGE18_4_ZERO_TOL:=1.0e-14}"
: "${STAGE18_4_FORMULA_TOL:=1.0e-10}"
: "${STAGE18_4_SOLVE_TOL:=5.0e-3}"
: "${STAGE18_4_ARC_ERROR_TOL:=1.0e-6}"
: "${STAGE18_4_TEST_CASE:=straight_manufactured_velocity_arclength}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are interface-compatibility
# variables only. Stage 18.4 intentionally does not cd into DECOMP2D_ROOT, build,
# invoke MPI, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-4-enable "${STAGE18_4_ENABLE}" \
  --tension-constraint-enable "${STAGE18_4_TENSION_CONSTRAINT_ENABLE}" \
  --single-fibre-only "${STAGE18_4_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_4_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_4_NPTS}" \
  --fibre-length "${STAGE18_4_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_4_COMPONENT_DIM}" \
  --rho-l "${STAGE18_4_RHO_L}" \
  --rho-tilde "${STAGE18_4_RHO_TILDE}" \
  --bending-stiffness "${STAGE18_4_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_4_GAMMA}" \
  --velocity-eps "${STAGE18_4_VELOCITY_EPS}" \
  --arclength-stretch-eps "${STAGE18_4_ARCLENGTH_STRETCH_EPS}" \
  --zero-tol "${STAGE18_4_ZERO_TOL}" \
  --formula-tol "${STAGE18_4_FORMULA_TOL}" \
  --solve-tol "${STAGE18_4_SOLVE_TOL}" \
  --arc-error-tol "${STAGE18_4_ARC_ERROR_TOL}" \
  --test-case "${STAGE18_4_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.4 TENSION INEXTENSIBILITY CONSTRAINT VERDICT: PASS"
  echo "STAGE 18.4 FINAL VERDICT: PASS"
else
  echo "STAGE 18.4 TENSION INEXTENSIBILITY CONSTRAINT VERDICT: FAIL"
  echo "STAGE 18.4 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
