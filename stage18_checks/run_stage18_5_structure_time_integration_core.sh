#!/usr/bin/env bash
# Stage 18.5 standalone structure time-integration core diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_5_ENABLE:=1}"
: "${STAGE18_5_TIME_INTEGRATION_CORE_ENABLE:=1}"
: "${STAGE18_5_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_5_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_5_NPTS:=16}"
: "${STAGE18_5_FIBRE_LENGTH:=1.0}"
: "${STAGE18_5_COMPONENT_DIM:=3}"
: "${STAGE18_5_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_5_RHO_L:=1.0}"
: "${STAGE18_5_RHO_TILDE:=1.0}"
: "${STAGE18_5_USE_DIMENSIONAL_MASS:=1}"
: "${STAGE18_5_USE_NONDIMENSIONAL_MASS:=1}"
: "${STAGE18_5_UNIFORM_VELOCITY:=1.0e-3}"
: "${STAGE18_5_CONSTANT_FORCE:=1.0e-3}"
: "${STAGE18_5_ZERO_TOL:=1.0e-14}"
: "${STAGE18_5_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_5_DT_REFINEMENT_TOL:=1.0e-12}"
: "${STAGE18_5_TEST_CASE:=zero_uniform_constant_split_dt_refinement}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are interface-compatibility
# variables only. Stage 18.5 intentionally does not cd into DECOMP2D_ROOT, build,
# invoke MPI, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_5_structure_time_integration_core.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-5-enable "${STAGE18_5_ENABLE}" \
  --time-integration-core-enable "${STAGE18_5_TIME_INTEGRATION_CORE_ENABLE}" \
  --single-fibre-only "${STAGE18_5_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_5_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_5_NPTS}" \
  --fibre-length "${STAGE18_5_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_5_COMPONENT_DIM}" \
  --dt-structure "${STAGE18_5_DT_STRUCTURE}" \
  --rho-l "${STAGE18_5_RHO_L}" \
  --rho-tilde "${STAGE18_5_RHO_TILDE}" \
  --use-dimensional-mass "${STAGE18_5_USE_DIMENSIONAL_MASS}" \
  --use-nondimensional-mass "${STAGE18_5_USE_NONDIMENSIONAL_MASS}" \
  --uniform-velocity "${STAGE18_5_UNIFORM_VELOCITY}" \
  --constant-force "${STAGE18_5_CONSTANT_FORCE}" \
  --zero-tol "${STAGE18_5_ZERO_TOL}" \
  --formula-tol "${STAGE18_5_FORMULA_TOL}" \
  --dt-refinement-tol "${STAGE18_5_DT_REFINEMENT_TOL}" \
  --test-case "${STAGE18_5_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.5 STRUCTURE TIME INTEGRATION CORE VERDICT: PASS"
  echo "STAGE 18.5 FINAL VERDICT: PASS"
else
  echo "STAGE 18.5 STRUCTURE TIME INTEGRATION CORE VERDICT: FAIL"
  echo "STAGE 18.5 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
