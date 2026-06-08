#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE17_6_RUN_STAGE17_5:=0}"
: "${STAGE17_6_REQUIRE_STAGE17_5:=1}"
: "${STAGE17_6_ACCEPT_STAGE17_5_CLOSED_EVIDENCE:=1}"
: "${STAGE17_6_ENABLE:=1}"
: "${STAGE17_6_WALL_SAFETY_ENABLE:=1}"
: "${STAGE17_6_BOUNDARY_CHECK_ENABLE:=1}"
: "${STAGE17_6_FAIL_CLOSED_ENABLE:=1}"
: "${STAGE17_6_CONTACT_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_6_DIAGNOSTIC_ONLY:=1}"
: "${STAGE17_6_Y_MIN:=-1.0}"
: "${STAGE17_6_Y_MAX:=1.0}"
: "${STAGE17_6_EFFECTIVE_FIBRE_RADIUS:=1.0e-3}"
: "${STAGE17_6_MIN_WALL_CLEARANCE:=1.0e-4}"
: "${STAGE17_6_WARNING_WALL_CLEARANCE:=1.0e-3}"
: "${STAGE17_6_PENETRATION_TOLERANCE:=1.0e-12}"
: "${STAGE17_6_NPTS:=8}"
: "${STAGE17_6_TEST_CASE:=all_segments_clear}"
: "${STAGE17_6_ZERO_TOL:=1.0e-14}"

mkdir -p stage17_outputs

if [ "${STAGE17_6_RUN_STAGE17_5}" = "1" ]; then
  if [ ! -x stage17_checks/run_stage17_5_near_wall_contact_state.sh ]; then
    echo "STAGE 17.6 SEGMENT WALL CLEARANCE SAFETY VERDICT: FAIL"
    echo "STAGE 17.6 FINAL VERDICT: FAIL"
    echo "reason stage17_5_wrapper_missing_or_not_executable"
    exit 1
  fi
  STAGE17_5_REQUIRE_STAGE17_4="${STAGE17_6_REQUIRE_STAGE17_5}" \
    stage17_checks/run_stage17_5_near_wall_contact_state.sh
fi

set +e
python3 stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py \
  --stage17-6-enable "${STAGE17_6_ENABLE}" \
  --wall-safety-enable "${STAGE17_6_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_6_BOUNDARY_CHECK_ENABLE}" \
  --fail-closed-enable "${STAGE17_6_FAIL_CLOSED_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_6_CONTACT_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_6_DIAGNOSTIC_ONLY}" \
  --y-min "${STAGE17_6_Y_MIN}" \
  --y-max "${STAGE17_6_Y_MAX}" \
  --effective-fibre-radius "${STAGE17_6_EFFECTIVE_FIBRE_RADIUS}" \
  --min-wall-clearance "${STAGE17_6_MIN_WALL_CLEARANCE}" \
  --warning-wall-clearance "${STAGE17_6_WARNING_WALL_CLEARANCE}" \
  --penetration-tolerance "${STAGE17_6_PENETRATION_TOLERANCE}" \
  --npts "${STAGE17_6_NPTS}" \
  --test-case "${STAGE17_6_TEST_CASE}" \
  --zero-tol "${STAGE17_6_ZERO_TOL}" \
  --accept-stage17-5-closed-evidence "${STAGE17_6_ACCEPT_STAGE17_5_CLOSED_EVIDENCE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 17.6 SEGMENT WALL CLEARANCE SAFETY VERDICT: PASS"
  echo "STAGE 17.6 FINAL VERDICT: PASS"
else
  echo "STAGE 17.6 SEGMENT WALL CLEARANCE SAFETY VERDICT: FAIL"
  echo "STAGE 17.6 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
