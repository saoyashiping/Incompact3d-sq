#!/usr/bin/env bash
# Stage 18.10 parallel-consistency physical-structure diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_10_ENABLE:=1}"
: "${STAGE18_10_PARALLEL_CONSISTENCY_ENABLE:=1}"
: "${STAGE18_10_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_10_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_10_NP_LIST:=1,2,4}"
: "${STAGE18_10_NPTS:=64}"
: "${STAGE18_10_NSTEPS:=5}"
: "${STAGE18_10_FIBRE_LENGTH:=1.0}"
: "${STAGE18_10_COMPONENT_DIM:=3}"
: "${STAGE18_10_RHO_L:=1.0}"
: "${STAGE18_10_RHO_TILDE:=1.0}"
: "${STAGE18_10_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_10_GAMMA:=1.0e-3}"
: "${STAGE18_10_USE_DIMENSIONAL_RESPONSE:=1}"
: "${STAGE18_10_USE_NONDIMENSIONAL_RESPONSE:=1}"
: "${STAGE18_10_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_10_SINE_EPS:=1.0e-3}"
: "${STAGE18_10_SINE_MODE:=1}"
: "${STAGE18_10_FLUID_FORCE_MAG:=1.0e-3}"
: "${STAGE18_10_INITIAL_VELOCITY_MAG:=0.0}"
: "${STAGE18_10_RESPONSE_BOUND:=1.0e-4}"
: "${STAGE18_10_VELOCITY_BOUND:=1.0e-3}"
: "${STAGE18_10_ACCELERATION_BOUND:=1.0e-2}"
: "${STAGE18_10_ZERO_TOL:=1.0e-14}"
: "${STAGE18_10_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_10_REDUCTION_TOL:=1.0e-12}"
: "${STAGE18_10_ENERGY_TOL:=1.0e-12}"
: "${STAGE18_10_BOUNDED_TOL:=1.0e-8}"
: "${STAGE18_10_TEST_CASE:=np1_np2_np4_partition_reduction_consistency}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# Compatibility variables only: Stage 18.10 does not cd into DECOMP2D_ROOT,
# build, invoke MPI, launch production DNS, spread forces to RHS/IBM, write
# production parallel output, or update production structure state.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_10_parallel_consistency_physical_structure.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-10-enable "${STAGE18_10_ENABLE}" \
  --parallel-consistency-enable "${STAGE18_10_PARALLEL_CONSISTENCY_ENABLE}" \
  --single-fibre-only "${STAGE18_10_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_10_DIAGNOSTIC_ONLY}" \
  --np-list "${STAGE18_10_NP_LIST}" \
  --npts "${STAGE18_10_NPTS}" \
  --nsteps "${STAGE18_10_NSTEPS}" \
  --fibre-length "${STAGE18_10_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_10_COMPONENT_DIM}" \
  --rho-l "${STAGE18_10_RHO_L}" \
  --rho-tilde "${STAGE18_10_RHO_TILDE}" \
  --bending-stiffness "${STAGE18_10_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_10_GAMMA}" \
  --use-dimensional-response "${STAGE18_10_USE_DIMENSIONAL_RESPONSE}" \
  --use-nondimensional-response "${STAGE18_10_USE_NONDIMENSIONAL_RESPONSE}" \
  --dt-structure "${STAGE18_10_DT_STRUCTURE}" \
  --sine-eps "${STAGE18_10_SINE_EPS}" \
  --sine-mode "${STAGE18_10_SINE_MODE}" \
  --fluid-force-mag "${STAGE18_10_FLUID_FORCE_MAG}" \
  --initial-velocity-mag "${STAGE18_10_INITIAL_VELOCITY_MAG}" \
  --response-bound "${STAGE18_10_RESPONSE_BOUND}" \
  --velocity-bound "${STAGE18_10_VELOCITY_BOUND}" \
  --acceleration-bound "${STAGE18_10_ACCELERATION_BOUND}" \
  --zero-tol "${STAGE18_10_ZERO_TOL}" \
  --formula-tol "${STAGE18_10_FORMULA_TOL}" \
  --reduction-tol "${STAGE18_10_REDUCTION_TOL}" \
  --energy-tol "${STAGE18_10_ENERGY_TOL}" \
  --bounded-tol "${STAGE18_10_BOUNDED_TOL}" \
  --test-case "${STAGE18_10_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.10 PARALLEL CONSISTENCY PHYSICAL STRUCTURE VERDICT: PASS"
  echo "STAGE 18.10 FINAL VERDICT: PASS"
else
  echo "STAGE 18.10 PARALLEL CONSISTENCY PHYSICAL STRUCTURE VERDICT: FAIL"
  echo "STAGE 18.10 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
