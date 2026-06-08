#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE17_9_RUN_STAGE17_8:=0}"
: "${STAGE17_9_REQUIRE_STAGE17_8:=1}"
: "${STAGE17_9_ACCEPT_STAGE17_8_CLOSED_EVIDENCE:=1}"
: "${STAGE17_9_RUN_STAGE16_CLOSED_LOOP:=0}"
: "${STAGE17_9_REQUIRE_STAGE16_CLOSED_LOOP:=1}"
: "${STAGE17_9_ACCEPT_STAGE16_CLOSED_LOOP_EVIDENCE:=1}"
: "${STAGE17_9_ENABLE:=1}"
: "${STAGE17_9_WALL_SAFETY_ENABLE:=1}"
: "${STAGE17_9_CONTACT_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_9_FIBRE_COLLISION_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_9_DIAGNOSTIC_ONLY:=1}"
: "${STAGE17_9_SMALL_LAMBDA:=1.0e-8}"
: "${STAGE17_9_MAX_RHS_INCREMENT:=1.0e-8}"
: "${STAGE17_9_MAX_FLUID_DELTA:=1.0e-8}"
: "${STAGE17_9_MAX_CONTACT_FORCE_NORM:=0.0}"
: "${STAGE17_9_MAX_CONTACT_RHS_NORM:=0.0}"
: "${STAGE17_9_MAX_CONTACT_STRUCTURE_UPDATE_NORM:=0.0}"

mkdir -p stage17_outputs

if [ "${STAGE17_9_RUN_STAGE17_8}" = "1" ]; then
  if [ ! -x stage17_checks/run_stage17_8_fibre_fibre_placeholder_geometry.sh ]; then
    echo "STAGE 17.9 CLOSED-LOOP WALL CONTACT COMPATIBILITY VERDICT: FAIL"
    echo "STAGE 17.9 FINAL VERDICT: FAIL"
    echo "reason stage17_8_wrapper_missing_or_not_executable"
    exit 1
  fi
  STAGE17_8_REQUIRE_STAGE17_7="${STAGE17_9_REQUIRE_STAGE17_8}" \
    stage17_checks/run_stage17_8_fibre_fibre_placeholder_geometry.sh
fi

if [ "${STAGE17_9_RUN_STAGE16_CLOSED_LOOP}" = "1" ]; then
  if [ ! -x stage16_checks/run_stage16_12_total_smoke_closure.sh ]; then
    echo "STAGE 17.9 CLOSED-LOOP WALL CONTACT COMPATIBILITY VERDICT: FAIL"
    echo "STAGE 17.9 FINAL VERDICT: FAIL"
    echo "reason stage16_closed_loop_wrapper_missing_or_not_executable"
    exit 1
  fi
  stage16_checks/run_stage16_12_total_smoke_closure.sh
fi

set +e
python3 stage17_checks/assert_stage17_9_closed_loop_wall_contact_compatibility.py \
  --stage17-9-enable "${STAGE17_9_ENABLE}" \
  --wall-safety-enable "${STAGE17_9_WALL_SAFETY_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_9_CONTACT_PLACEHOLDER_ENABLE}" \
  --fibre-collision-placeholder-enable "${STAGE17_9_FIBRE_COLLISION_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_9_DIAGNOSTIC_ONLY}" \
  --small-lambda "${STAGE17_9_SMALL_LAMBDA}" \
  --max-rhs-increment "${STAGE17_9_MAX_RHS_INCREMENT}" \
  --max-fluid-delta "${STAGE17_9_MAX_FLUID_DELTA}" \
  --max-contact-force-norm "${STAGE17_9_MAX_CONTACT_FORCE_NORM}" \
  --max-contact-rhs-norm "${STAGE17_9_MAX_CONTACT_RHS_NORM}" \
  --max-contact-structure-update-norm "${STAGE17_9_MAX_CONTACT_STRUCTURE_UPDATE_NORM}" \
  --accept-stage17-8-closed-evidence "${STAGE17_9_ACCEPT_STAGE17_8_CLOSED_EVIDENCE}" \
  --accept-stage16-closed-loop-evidence "${STAGE17_9_ACCEPT_STAGE16_CLOSED_LOOP_EVIDENCE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 17.9 CLOSED-LOOP WALL CONTACT COMPATIBILITY VERDICT: PASS"
  echo "STAGE 17.9 FINAL VERDICT: PASS"
else
  echo "STAGE 17.9 CLOSED-LOOP WALL CONTACT COMPATIBILITY VERDICT: FAIL"
  echo "STAGE 17.9 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
