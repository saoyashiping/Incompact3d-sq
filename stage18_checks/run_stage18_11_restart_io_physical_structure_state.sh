#!/usr/bin/env bash
# Stage 18.11 restart/I/O compatibility diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_11_ENABLE:=1}"
: "${STAGE18_11_RESTART_IO_COMPATIBILITY_ENABLE:=1}"
: "${STAGE18_11_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_11_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_11_NP_LIST:=1,2,4}"
: "${STAGE18_11_NPTS:=64}"
: "${STAGE18_11_NSTEPS:=5}"
: "${STAGE18_11_FIBRE_LENGTH:=1.0}"
: "${STAGE18_11_COMPONENT_DIM:=3}"
: "${STAGE18_11_RHO_L:=1.0}"
: "${STAGE18_11_RHO_TILDE:=1.0}"
: "${STAGE18_11_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_11_GAMMA:=1.0e-3}"
: "${STAGE18_11_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_11_SINE_EPS:=1.0e-3}"
: "${STAGE18_11_SINE_MODE:=1}"
: "${STAGE18_11_FLUID_FORCE_MAG:=1.0e-3}"
: "${STAGE18_11_INITIAL_VELOCITY_MAG:=0.0}"
: "${STAGE18_11_RESPONSE_BOUND:=1.0e-4}"
: "${STAGE18_11_VELOCITY_BOUND:=1.0e-3}"
: "${STAGE18_11_ACCELERATION_BOUND:=1.0e-2}"
: "${STAGE18_11_ZERO_TOL:=1.0e-14}"
: "${STAGE18_11_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_11_RESTART_TOL:=1.0e-12}"
: "${STAGE18_11_REDUCTION_TOL:=1.0e-12}"
: "${STAGE18_11_ENERGY_TOL:=1.0e-12}"
: "${STAGE18_11_TEST_CASE:=snapshot_roundtrip_restart_equivalence_partition_io}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# Compatibility variables only: Stage 18.11 does not cd into DECOMP2D_ROOT,
# build, invoke MPI, launch production DNS, or write production restart,
# statistics, visualisation, RHS, IBM, DNS-core, or structure-state files.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_11_restart_io_physical_structure_state.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-11-enable "${STAGE18_11_ENABLE}" \
  --restart-io-compatibility-enable "${STAGE18_11_RESTART_IO_COMPATIBILITY_ENABLE}" \
  --single-fibre-only "${STAGE18_11_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_11_DIAGNOSTIC_ONLY}" \
  --np-list "${STAGE18_11_NP_LIST}" \
  --npts "${STAGE18_11_NPTS}" \
  --nsteps "${STAGE18_11_NSTEPS}" \
  --fibre-length "${STAGE18_11_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_11_COMPONENT_DIM}" \
  --rho-l "${STAGE18_11_RHO_L}" \
  --rho-tilde "${STAGE18_11_RHO_TILDE}" \
  --bending-stiffness "${STAGE18_11_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_11_GAMMA}" \
  --dt-structure "${STAGE18_11_DT_STRUCTURE}" \
  --sine-eps "${STAGE18_11_SINE_EPS}" \
  --sine-mode "${STAGE18_11_SINE_MODE}" \
  --fluid-force-mag "${STAGE18_11_FLUID_FORCE_MAG}" \
  --initial-velocity-mag "${STAGE18_11_INITIAL_VELOCITY_MAG}" \
  --response-bound "${STAGE18_11_RESPONSE_BOUND}" \
  --velocity-bound "${STAGE18_11_VELOCITY_BOUND}" \
  --acceleration-bound "${STAGE18_11_ACCELERATION_BOUND}" \
  --zero-tol "${STAGE18_11_ZERO_TOL}" \
  --formula-tol "${STAGE18_11_FORMULA_TOL}" \
  --restart-tol "${STAGE18_11_RESTART_TOL}" \
  --reduction-tol "${STAGE18_11_REDUCTION_TOL}" \
  --energy-tol "${STAGE18_11_ENERGY_TOL}" \
  --test-case "${STAGE18_11_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.11 RESTART IO PHYSICAL STRUCTURE STATE VERDICT: PASS"
  echo "STAGE 18.11 FINAL VERDICT: PASS"
else
  echo "STAGE 18.11 RESTART IO PHYSICAL STRUCTURE STATE VERDICT: FAIL"
  echo "STAGE 18.11 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
