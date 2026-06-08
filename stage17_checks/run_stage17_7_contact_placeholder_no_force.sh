#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE17_7_RUN_STAGE17_6:=0}"
: "${STAGE17_7_REQUIRE_STAGE17_6:=1}"
: "${STAGE17_7_ACCEPT_STAGE17_6_CLOSED_EVIDENCE:=1}"
: "${STAGE17_7_ENABLE:=1}"
: "${STAGE17_7_WALL_SAFETY_ENABLE:=1}"
: "${STAGE17_7_BOUNDARY_CHECK_ENABLE:=1}"
: "${STAGE17_7_FAIL_CLOSED_ENABLE:=1}"
: "${STAGE17_7_CONTACT_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_7_FIBRE_COLLISION_PLACEHOLDER_ENABLE:=0}"
: "${STAGE17_7_DIAGNOSTIC_ONLY:=1}"
: "${STAGE17_7_Y_MIN:=-1.0}"
: "${STAGE17_7_Y_MAX:=1.0}"
: "${STAGE17_7_EFFECTIVE_FIBRE_RADIUS:=1.0e-3}"
: "${STAGE17_7_MIN_WALL_CLEARANCE:=1.0e-4}"
: "${STAGE17_7_WARNING_WALL_CLEARANCE:=1.0e-3}"
: "${STAGE17_7_PENETRATION_TOLERANCE:=1.0e-12}"
: "${STAGE17_7_NPTS:=8}"
: "${STAGE17_7_TEST_CASE:=clear_no_contact_record}"
: "${STAGE17_7_ZERO_TOL:=1.0e-14}"

mkdir -p stage17_outputs

if [ "${STAGE17_7_RUN_STAGE17_6}" = "1" ]; then
  if [ ! -x stage17_checks/run_stage17_6_segment_wall_clearance_safety.sh ]; then
    echo "STAGE 17.7 CONTACT PLACEHOLDER NO-FORCE VERDICT: FAIL"
    echo "STAGE 17.7 FINAL VERDICT: FAIL"
    echo "reason stage17_6_wrapper_missing_or_not_executable"
    exit 1
  fi
  STAGE17_6_REQUIRE_STAGE17_5="${STAGE17_7_REQUIRE_STAGE17_6}" \
    stage17_checks/run_stage17_6_segment_wall_clearance_safety.sh
fi

set +e
python3 stage17_checks/assert_stage17_7_contact_placeholder_no_force.py \
  --stage17-7-enable "${STAGE17_7_ENABLE}" \
  --wall-safety-enable "${STAGE17_7_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_7_BOUNDARY_CHECK_ENABLE}" \
  --fail-closed-enable "${STAGE17_7_FAIL_CLOSED_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_7_CONTACT_PLACEHOLDER_ENABLE}" \
  --fibre-collision-placeholder-enable "${STAGE17_7_FIBRE_COLLISION_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_7_DIAGNOSTIC_ONLY}" \
  --y-min "${STAGE17_7_Y_MIN}" \
  --y-max "${STAGE17_7_Y_MAX}" \
  --effective-fibre-radius "${STAGE17_7_EFFECTIVE_FIBRE_RADIUS}" \
  --min-wall-clearance "${STAGE17_7_MIN_WALL_CLEARANCE}" \
  --warning-wall-clearance "${STAGE17_7_WARNING_WALL_CLEARANCE}" \
  --penetration-tolerance "${STAGE17_7_PENETRATION_TOLERANCE}" \
  --npts "${STAGE17_7_NPTS}" \
  --test-case "${STAGE17_7_TEST_CASE}" \
  --zero-tol "${STAGE17_7_ZERO_TOL}" \
  --accept-stage17-6-closed-evidence "${STAGE17_7_ACCEPT_STAGE17_6_CLOSED_EVIDENCE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 17.7 CONTACT PLACEHOLDER NO-FORCE VERDICT: PASS"
  echo "STAGE 17.7 FINAL VERDICT: PASS"
else
  echo "STAGE 17.7 CONTACT PLACEHOLDER NO-FORCE VERDICT: FAIL"
  echo "STAGE 17.7 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
