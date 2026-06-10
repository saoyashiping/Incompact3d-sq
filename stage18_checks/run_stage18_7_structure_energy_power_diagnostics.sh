#!/usr/bin/env bash
# Stage 18.7 structure energy/power diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_7_ENABLE:=1}"
: "${STAGE18_7_ENERGY_POWER_DIAGNOSTIC_ENABLE:=1}"
: "${STAGE18_7_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_7_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_7_NPTS:=64}"
: "${STAGE18_7_FIBRE_LENGTH:=1.0}"
: "${STAGE18_7_COMPONENT_DIM:=3}"
: "${STAGE18_7_RHO_L:=1.0}"
: "${STAGE18_7_RHO_TILDE:=1.0}"
: "${STAGE18_7_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_7_GAMMA:=1.0e-3}"
: "${STAGE18_7_USE_DIMENSIONAL_ENERGY:=1}"
: "${STAGE18_7_USE_NONDIMENSIONAL_ENERGY:=1}"
: "${STAGE18_7_VELOCITY_MAG:=1.0e-3}"
: "${STAGE18_7_FLUID_FORCE_MAG:=1.0e-3}"
: "${STAGE18_7_SINE_EPS:=1.0e-3}"
: "${STAGE18_7_SINE_MODE:=1}"
: "${STAGE18_7_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_7_ZERO_TOL:=1.0e-14}"
: "${STAGE18_7_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_7_POWER_TOL:=1.0e-12}"
: "${STAGE18_7_ENERGY_TOL:=1.0e-12}"
: "${STAGE18_7_TEST_CASE:=straight_rest_uniform_sine_power_residual}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are compatibility
# variables only. Stage 18.7 does not cd into DECOMP2D_ROOT, build, invoke MPI,
# write production energy output, modify RHS/IBM/DNS-core, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_7_structure_energy_power_diagnostics.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-7-enable "${STAGE18_7_ENABLE}" \
  --energy-power-diagnostic-enable "${STAGE18_7_ENERGY_POWER_DIAGNOSTIC_ENABLE}" \
  --single-fibre-only "${STAGE18_7_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_7_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_7_NPTS}" \
  --fibre-length "${STAGE18_7_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_7_COMPONENT_DIM}" \
  --rho-l "${STAGE18_7_RHO_L}" \
  --rho-tilde "${STAGE18_7_RHO_TILDE}" \
  --bending-stiffness "${STAGE18_7_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_7_GAMMA}" \
  --use-dimensional-energy "${STAGE18_7_USE_DIMENSIONAL_ENERGY}" \
  --use-nondimensional-energy "${STAGE18_7_USE_NONDIMENSIONAL_ENERGY}" \
  --velocity-mag "${STAGE18_7_VELOCITY_MAG}" \
  --fluid-force-mag "${STAGE18_7_FLUID_FORCE_MAG}" \
  --sine-eps "${STAGE18_7_SINE_EPS}" \
  --sine-mode "${STAGE18_7_SINE_MODE}" \
  --dt-structure "${STAGE18_7_DT_STRUCTURE}" \
  --zero-tol "${STAGE18_7_ZERO_TOL}" \
  --formula-tol "${STAGE18_7_FORMULA_TOL}" \
  --power-tol "${STAGE18_7_POWER_TOL}" \
  --energy-tol "${STAGE18_7_ENERGY_TOL}" \
  --test-case "${STAGE18_7_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.7 STRUCTURE ENERGY POWER DIAGNOSTICS VERDICT: PASS"
  echo "STAGE 18.7 FINAL VERDICT: PASS"
else
  echo "STAGE 18.7 STRUCTURE ENERGY POWER DIAGNOSTICS VERDICT: FAIL"
  echo "STAGE 18.7 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
