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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"
mkdir -p stage18_outputs

# DECOMP2D_ROOT, BUILD_DIR, MPIEXEC, and MPIEXEC_FLAGS are accepted for
# interface compatibility with earlier stages, but Stage 18.0 is a repository
# preflight and must run from the repository root, not from DECOMP2D_ROOT.
set +e
python3 "${SCRIPT_DIR}/assert_stage18_0_preflight_boundary.py" \
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
