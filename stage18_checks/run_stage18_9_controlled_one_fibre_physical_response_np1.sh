#!/usr/bin/env bash
# Stage 18.9 controlled one-fibre physical response np=1 diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_9_ENABLE:=1}"
: "${STAGE18_9_CONTROLLED_RESPONSE_ENABLE:=1}"
: "${STAGE18_9_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_9_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_9_NP:=1}"
: "${STAGE18_9_NPTS:=64}"
: "${STAGE18_9_NSTEPS:=5}"
: "${STAGE18_9_FIBRE_LENGTH:=1.0}"
: "${STAGE18_9_COMPONENT_DIM:=3}"
: "${STAGE18_9_RHO_L:=1.0}"
: "${STAGE18_9_RHO_TILDE:=1.0}"
: "${STAGE18_9_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_9_GAMMA:=1.0e-3}"
: "${STAGE18_9_USE_DIMENSIONAL_RESPONSE:=1}"
: "${STAGE18_9_USE_NONDIMENSIONAL_RESPONSE:=1}"
: "${STAGE18_9_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_9_SINE_EPS:=1.0e-3}"
: "${STAGE18_9_SINE_MODE:=1}"
: "${STAGE18_9_FLUID_FORCE_MAG:=1.0e-3}"
: "${STAGE18_9_INITIAL_VELOCITY_MAG:=0.0}"
: "${STAGE18_9_RESPONSE_BOUND:=1.0e-4}"
: "${STAGE18_9_VELOCITY_BOUND:=1.0e-3}"
: "${STAGE18_9_ACCELERATION_BOUND:=1.0e-2}"
: "${STAGE18_9_ZERO_TOL:=1.0e-14}"
: "${STAGE18_9_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_9_ENERGY_TOL:=1.0e-12}"
: "${STAGE18_9_BOUNDED_TOL:=1.0e-8}"
: "${STAGE18_9_TEST_CASE:=np1_dry_forced_bending_energy_bounded}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# Compatibility variables only: Stage 18.9 does not cd into DECOMP2D_ROOT, build,
# invoke MPI, launch production DNS, spread forces to RHS/IBM, write production
# response output, or update production structure state.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_9_controlled_one_fibre_physical_response_np1.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-9-enable "${STAGE18_9_ENABLE}" \
  --controlled-response-enable "${STAGE18_9_CONTROLLED_RESPONSE_ENABLE}" \
  --single-fibre-only "${STAGE18_9_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_9_DIAGNOSTIC_ONLY}" \
  --np "${STAGE18_9_NP}" \
  --npts "${STAGE18_9_NPTS}" \
  --nsteps "${STAGE18_9_NSTEPS}" \
  --fibre-length "${STAGE18_9_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_9_COMPONENT_DIM}" \
  --rho-l "${STAGE18_9_RHO_L}" \
  --rho-tilde "${STAGE18_9_RHO_TILDE}" \
  --bending-stiffness "${STAGE18_9_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_9_GAMMA}" \
  --use-dimensional-response "${STAGE18_9_USE_DIMENSIONAL_RESPONSE}" \
  --use-nondimensional-response "${STAGE18_9_USE_NONDIMENSIONAL_RESPONSE}" \
  --dt-structure "${STAGE18_9_DT_STRUCTURE}" \
  --sine-eps "${STAGE18_9_SINE_EPS}" \
  --sine-mode "${STAGE18_9_SINE_MODE}" \
  --fluid-force-mag "${STAGE18_9_FLUID_FORCE_MAG}" \
  --initial-velocity-mag "${STAGE18_9_INITIAL_VELOCITY_MAG}" \
  --response-bound "${STAGE18_9_RESPONSE_BOUND}" \
  --velocity-bound "${STAGE18_9_VELOCITY_BOUND}" \
  --acceleration-bound "${STAGE18_9_ACCELERATION_BOUND}" \
  --zero-tol "${STAGE18_9_ZERO_TOL}" \
  --formula-tol "${STAGE18_9_FORMULA_TOL}" \
  --energy-tol "${STAGE18_9_ENERGY_TOL}" \
  --bounded-tol "${STAGE18_9_BOUNDED_TOL}" \
  --test-case "${STAGE18_9_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.9 CONTROLLED ONE FIBRE PHYSICAL RESPONSE NP1 VERDICT: PASS"
  echo "STAGE 18.9 FINAL VERDICT: PASS"
else
  echo "STAGE 18.9 CONTROLLED ONE FIBRE PHYSICAL RESPONSE NP1 VERDICT: FAIL"
  echo "STAGE 18.9 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
