#!/usr/bin/env bash
# Stage 18.6 fluid-on-fibre force input diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_6_ENABLE:=1}"
: "${STAGE18_6_FLUID_FORCE_INPUT_ENABLE:=1}"
: "${STAGE18_6_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_6_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_6_NPTS:=16}"
: "${STAGE18_6_FIBRE_LENGTH:=1.0}"
: "${STAGE18_6_COMPONENT_DIM:=3}"
: "${STAGE18_6_RHO_L:=1.0}"
: "${STAGE18_6_RHO_TILDE:=1.0}"
: "${STAGE18_6_USE_DIMENSIONAL_MASS:=1}"
: "${STAGE18_6_USE_NONDIMENSIONAL_MASS:=1}"
: "${STAGE18_6_FLUID_FORCE_MAG:=1.0e-3}"
: "${STAGE18_6_VELOCITY_MAG:=1.0e-3}"
: "${STAGE18_6_ZERO_TOL:=1.0e-14}"
: "${STAGE18_6_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_6_POWER_TOL:=1.0e-12}"
: "${STAGE18_6_TEST_CASE:=zero_uniform_sign_split_power}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are compatibility
# variables only. Stage 18.6 does not cd into DECOMP2D_ROOT, build, invoke MPI,
# spread force to RHS, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_6_fluid_force_input_physical_structure.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-6-enable "${STAGE18_6_ENABLE}" \
  --fluid-force-input-enable "${STAGE18_6_FLUID_FORCE_INPUT_ENABLE}" \
  --single-fibre-only "${STAGE18_6_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_6_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_6_NPTS}" \
  --fibre-length "${STAGE18_6_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_6_COMPONENT_DIM}" \
  --rho-l "${STAGE18_6_RHO_L}" \
  --rho-tilde "${STAGE18_6_RHO_TILDE}" \
  --use-dimensional-mass "${STAGE18_6_USE_DIMENSIONAL_MASS}" \
  --use-nondimensional-mass "${STAGE18_6_USE_NONDIMENSIONAL_MASS}" \
  --fluid-force-mag "${STAGE18_6_FLUID_FORCE_MAG}" \
  --velocity-mag "${STAGE18_6_VELOCITY_MAG}" \
  --zero-tol "${STAGE18_6_ZERO_TOL}" \
  --formula-tol "${STAGE18_6_FORMULA_TOL}" \
  --power-tol "${STAGE18_6_POWER_TOL}" \
  --test-case "${STAGE18_6_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.6 FLUID FORCE INPUT PHYSICAL STRUCTURE VERDICT: PASS"
  echo "STAGE 18.6 FINAL VERDICT: PASS"
else
  echo "STAGE 18.6 FLUID FORCE INPUT PHYSICAL STRUCTURE VERDICT: FAIL"
  echo "STAGE 18.6 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
