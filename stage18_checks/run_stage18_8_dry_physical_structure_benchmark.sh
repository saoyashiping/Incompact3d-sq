#!/usr/bin/env bash
# Stage 18.8 dry physical structure benchmark diagnostic wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_8_ENABLE:=1}"
: "${STAGE18_8_DRY_STRUCTURE_BENCHMARK_ENABLE:=1}"
: "${STAGE18_8_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_8_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_8_NPTS:=64}"
: "${STAGE18_8_FIBRE_LENGTH:=1.0}"
: "${STAGE18_8_COMPONENT_DIM:=3}"
: "${STAGE18_8_RHO_L:=1.0}"
: "${STAGE18_8_RHO_TILDE:=1.0}"
: "${STAGE18_8_BENDING_STIFFNESS:=1.0e-3}"
: "${STAGE18_8_GAMMA:=1.0e-3}"
: "${STAGE18_8_USE_DIMENSIONAL_DRY:=1}"
: "${STAGE18_8_USE_NONDIMENSIONAL_DRY:=1}"
: "${STAGE18_8_VELOCITY_MAG:=1.0e-3}"
: "${STAGE18_8_SINE_EPS:=1.0e-3}"
: "${STAGE18_8_SINE_MODE:=1}"
: "${STAGE18_8_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_8_MAX_DISPLACEMENT_INCREMENT:=1.0e-4}"
: "${STAGE18_8_MAX_VELOCITY_INCREMENT:=1.0e-3}"
: "${STAGE18_8_ZERO_TOL:=1.0e-14}"
: "${STAGE18_8_FORMULA_TOL:=1.0e-12}"
: "${STAGE18_8_ENERGY_TOL:=1.0e-12}"
: "${STAGE18_8_BOUNDED_TOL:=1.0e-8}"
: "${STAGE18_8_TEST_CASE:=dry_straight_translation_sine_bounded_energy}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# Compatibility variables only: Stage 18.8 does not cd into DECOMP2D_ROOT, build,
# invoke MPI, write production dry-benchmark output, modify RHS/IBM/DNS-core, or
# run production structure integration.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_8_dry_physical_structure_benchmark.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-8-enable "${STAGE18_8_ENABLE}" \
  --dry-structure-benchmark-enable "${STAGE18_8_DRY_STRUCTURE_BENCHMARK_ENABLE}" \
  --single-fibre-only "${STAGE18_8_SINGLE_FIBRE_ONLY}" \
  --diagnostic-only "${STAGE18_8_DIAGNOSTIC_ONLY}" \
  --npts "${STAGE18_8_NPTS}" \
  --fibre-length "${STAGE18_8_FIBRE_LENGTH}" \
  --component-dim "${STAGE18_8_COMPONENT_DIM}" \
  --rho-l "${STAGE18_8_RHO_L}" \
  --rho-tilde "${STAGE18_8_RHO_TILDE}" \
  --bending-stiffness "${STAGE18_8_BENDING_STIFFNESS}" \
  --gamma "${STAGE18_8_GAMMA}" \
  --use-dimensional-dry "${STAGE18_8_USE_DIMENSIONAL_DRY}" \
  --use-nondimensional-dry "${STAGE18_8_USE_NONDIMENSIONAL_DRY}" \
  --velocity-mag "${STAGE18_8_VELOCITY_MAG}" \
  --sine-eps "${STAGE18_8_SINE_EPS}" \
  --sine-mode "${STAGE18_8_SINE_MODE}" \
  --dt-structure "${STAGE18_8_DT_STRUCTURE}" \
  --max-displacement-increment "${STAGE18_8_MAX_DISPLACEMENT_INCREMENT}" \
  --max-velocity-increment "${STAGE18_8_MAX_VELOCITY_INCREMENT}" \
  --zero-tol "${STAGE18_8_ZERO_TOL}" \
  --formula-tol "${STAGE18_8_FORMULA_TOL}" \
  --energy-tol "${STAGE18_8_ENERGY_TOL}" \
  --bounded-tol "${STAGE18_8_BOUNDED_TOL}" \
  --test-case "${STAGE18_8_TEST_CASE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.8 DRY PHYSICAL STRUCTURE BENCHMARK VERDICT: PASS"
  echo "STAGE 18.8 FINAL VERDICT: PASS"
else
  echo "STAGE 18.8 DRY PHYSICAL STRUCTURE BENCHMARK VERDICT: FAIL"
  echo "STAGE 18.8 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
