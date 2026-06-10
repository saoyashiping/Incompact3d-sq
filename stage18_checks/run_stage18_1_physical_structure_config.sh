#!/usr/bin/env bash
# Stage 18.1 physical structure dynamics configuration wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_1_ENABLE:=1}"
: "${STAGE18_1_STRUCTURE_DYNAMICS_CONFIG_ENABLE:=1}"
: "${STAGE18_1_SINGLE_FIBRE_ONLY:=1}"
: "${STAGE18_1_PHYSICAL_STRUCTURE_BOUNDARY:=1}"
: "${STAGE18_1_BENDING_CONFIG_ENABLE:=1}"
: "${STAGE18_1_TENSION_CONFIG_ENABLE:=1}"
: "${STAGE18_1_INEXTENSIBILITY_CONFIG_ENABLE:=1}"
: "${STAGE18_1_TIME_INTEGRATION_CONFIG_ENABLE:=1}"
: "${STAGE18_1_ENERGY_DIAGNOSTIC_CONFIG_ENABLE:=1}"
: "${STAGE18_1_DIAGNOSTIC_ONLY:=1}"
: "${STAGE18_1_RHO_S:=1.0}"
: "${STAGE18_1_FIBRE_RADIUS:=1.0e-3}"
: "${STAGE18_1_YOUNG_MODULUS:=1.0}"
: "${STAGE18_1_FIBRE_LENGTH:=1.0}"
: "${STAGE18_1_NPTS:=8}"
: "${STAGE18_1_DT_STRUCTURE:=1.0e-4}"
: "${STAGE18_1_RHO_TILDE:=1.0}"
: "${STAGE18_1_GAMMA:=1.0e-3}"
: "${STAGE18_1_ZERO_TOL:=1.0e-14}"
: "${STAGE18_1_FORMULA_TOL:=1.0e-12}"

mkdir -p "${REPO_ROOT}/stage18_outputs"

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are interface-compatibility
# variables only.  Stage 18.1 intentionally does not cd into DECOMP2D_ROOT, build,
# invoke MPI, or run production physics.
set +e
python3 "${REPO_ROOT}/stage18_checks/assert_stage18_1_physical_structure_config.py" \
  --repo-root "${REPO_ROOT}" \
  --stage18-1-enable "${STAGE18_1_ENABLE}" \
  --structure-dynamics-config-enable "${STAGE18_1_STRUCTURE_DYNAMICS_CONFIG_ENABLE}" \
  --single-fibre-only "${STAGE18_1_SINGLE_FIBRE_ONLY}" \
  --physical-structure-boundary "${STAGE18_1_PHYSICAL_STRUCTURE_BOUNDARY}" \
  --bending-config-enable "${STAGE18_1_BENDING_CONFIG_ENABLE}" \
  --tension-config-enable "${STAGE18_1_TENSION_CONFIG_ENABLE}" \
  --inextensibility-config-enable "${STAGE18_1_INEXTENSIBILITY_CONFIG_ENABLE}" \
  --time-integration-config-enable "${STAGE18_1_TIME_INTEGRATION_CONFIG_ENABLE}" \
  --energy-diagnostic-config-enable "${STAGE18_1_ENERGY_DIAGNOSTIC_CONFIG_ENABLE}" \
  --diagnostic-only "${STAGE18_1_DIAGNOSTIC_ONLY}" \
  --rho-s "${STAGE18_1_RHO_S}" \
  --fibre-radius "${STAGE18_1_FIBRE_RADIUS}" \
  --young-modulus "${STAGE18_1_YOUNG_MODULUS}" \
  --fibre-length "${STAGE18_1_FIBRE_LENGTH}" \
  --npts "${STAGE18_1_NPTS}" \
  --dt-structure "${STAGE18_1_DT_STRUCTURE}" \
  --rho-tilde "${STAGE18_1_RHO_TILDE}" \
  --gamma "${STAGE18_1_GAMMA}" \
  --zero-tol "${STAGE18_1_ZERO_TOL}" \
  --formula-tol "${STAGE18_1_FORMULA_TOL}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.1 PHYSICAL STRUCTURE CONFIG VERDICT: PASS"
  echo "STAGE 18.1 FINAL VERDICT: PASS"
else
  echo "STAGE 18.1 PHYSICAL STRUCTURE CONFIG VERDICT: FAIL"
  echo "STAGE 18.1 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
