#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE17_8_RUN_STAGE17_7:=0}"
: "${STAGE17_8_REQUIRE_STAGE17_7:=1}"
: "${STAGE17_8_ACCEPT_STAGE17_7_CLOSED_EVIDENCE:=1}"
: "${STAGE17_8_ENABLE:=1}"
: "${STAGE17_8_FIBRE_COLLISION_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_8_CONTACT_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_8_DIAGNOSTIC_ONLY:=1}"
: "${STAGE17_8_EFFECTIVE_FIBRE_RADIUS:=1.0e-3}"
: "${STAGE17_8_MIN_FIBRE_FIBRE_CLEARANCE:=1.0e-4}"
: "${STAGE17_8_WARNING_FIBRE_FIBRE_CLEARANCE:=1.0e-3}"
: "${STAGE17_8_OVERLAP_TOLERANCE:=1.0e-12}"
: "${STAGE17_8_NPTS:=8}"
: "${STAGE17_8_TEST_CASE:=mock_fibres_clear}"
: "${STAGE17_8_ZERO_TOL:=1.0e-14}"

mkdir -p stage17_outputs

if [ "${STAGE17_8_RUN_STAGE17_7}" = "1" ]; then
  if [ ! -x stage17_checks/run_stage17_7_contact_placeholder_no_force.sh ]; then
    echo "STAGE 17.8 FIBRE-FIBRE PLACEHOLDER GEOMETRY VERDICT: FAIL"
    echo "STAGE 17.8 FINAL VERDICT: FAIL"
    echo "reason stage17_7_wrapper_missing_or_not_executable"
    exit 1
  fi
  STAGE17_7_REQUIRE_STAGE17_6="${STAGE17_8_REQUIRE_STAGE17_7}" \
    stage17_checks/run_stage17_7_contact_placeholder_no_force.sh
fi

set +e
python3 stage17_checks/assert_stage17_8_fibre_fibre_placeholder_geometry.py \
  --stage17-8-enable "${STAGE17_8_ENABLE}" \
  --fibre-collision-placeholder-enable "${STAGE17_8_FIBRE_COLLISION_PLACEHOLDER_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_8_CONTACT_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_8_DIAGNOSTIC_ONLY}" \
  --effective-fibre-radius "${STAGE17_8_EFFECTIVE_FIBRE_RADIUS}" \
  --min-fibre-fibre-clearance "${STAGE17_8_MIN_FIBRE_FIBRE_CLEARANCE}" \
  --warning-fibre-fibre-clearance "${STAGE17_8_WARNING_FIBRE_FIBRE_CLEARANCE}" \
  --overlap-tolerance "${STAGE17_8_OVERLAP_TOLERANCE}" \
  --npts "${STAGE17_8_NPTS}" \
  --test-case "${STAGE17_8_TEST_CASE}" \
  --zero-tol "${STAGE17_8_ZERO_TOL}" \
  --accept-stage17-7-closed-evidence "${STAGE17_8_ACCEPT_STAGE17_7_CLOSED_EVIDENCE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 17.8 FIBRE-FIBRE PLACEHOLDER GEOMETRY VERDICT: PASS"
  echo "STAGE 17.8 FINAL VERDICT: PASS"
else
  echo "STAGE 17.8 FIBRE-FIBRE PLACEHOLDER GEOMETRY VERDICT: FAIL"
  echo "STAGE 17.8 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
