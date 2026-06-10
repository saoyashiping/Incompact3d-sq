#!/usr/bin/env bash
# Stage 18.0 preflight boundary wrapper.
# Diagnostic-only: builds nothing, runs no MPI, and runs no production physics.
set -euo pipefail

: "${DECOMP2D_ROOT:=$(pwd)}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE18_0_REQUIRE_STAGE17_CLOSED:=1}"
: "${STAGE18_0_ACCEPT_STAGE17_CLOSED_EVIDENCE:=1}"
: "${STAGE18_0_ENABLE:=1}"
: "${STAGE18_0_SINGLE_FIBRE_STRUCTURE_DYNAMICS_BOUNDARY:=1}"
: "${STAGE18_0_DIAGNOSTIC_ONLY:=1}"

cd "${DECOMP2D_ROOT}"
mkdir -p stage18_outputs

# Intentionally unused runtime knobs are retained for interface compatibility.
# BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are accepted but never invoked here.
set +e
python3 stage18_checks/assert_stage18_0_preflight_boundary.py \
  --require-stage17-closed "${STAGE18_0_REQUIRE_STAGE17_CLOSED}" \
  --accept-stage17-closed-evidence "${STAGE18_0_ACCEPT_STAGE17_CLOSED_EVIDENCE}" \
  --stage18-0-enable "${STAGE18_0_ENABLE}" \
  --single-fibre-structure-dynamics-boundary "${STAGE18_0_SINGLE_FIBRE_STRUCTURE_DYNAMICS_BOUNDARY}" \
  --diagnostic-only "${STAGE18_0_DIAGNOSTIC_ONLY}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 18.0 PREFLIGHT BOUNDARY VERDICT: PASS"
  echo "STAGE 18.0 FINAL VERDICT: PASS"
else
  echo "STAGE 18.0 PREFLIGHT BOUNDARY VERDICT: FAIL"
  echo "STAGE 18.0 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
