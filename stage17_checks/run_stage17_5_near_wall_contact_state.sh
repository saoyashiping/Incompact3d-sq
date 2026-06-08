#!/usr/bin/env bash
set -u

# Stage 17.5 near-wall/contact-state classification wrapper.
# This wrapper creates diagnostic output only. It builds no target, runs no MPI,
# applies no force, and runs no production physics.

DECOMP2D_ROOT="${DECOMP2D_ROOT:-}"
BUILD_DIR="${BUILD_DIR:-build_stage9}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPIEXEC_FLAGS="${MPIEXEC_FLAGS:-}"
STAGE17_5_RUN_STAGE17_4="${STAGE17_5_RUN_STAGE17_4:-0}"
STAGE17_5_REQUIRE_STAGE17_4="${STAGE17_5_REQUIRE_STAGE17_4:-1}"
STAGE17_5_ACCEPT_STAGE17_4_CLOSED_EVIDENCE="${STAGE17_5_ACCEPT_STAGE17_4_CLOSED_EVIDENCE:-1}"
STAGE17_5_ENABLE="${STAGE17_5_ENABLE:-1}"
STAGE17_5_WALL_SAFETY_ENABLE="${STAGE17_5_WALL_SAFETY_ENABLE:-1}"
STAGE17_5_BOUNDARY_CHECK_ENABLE="${STAGE17_5_BOUNDARY_CHECK_ENABLE:-1}"
STAGE17_5_FAIL_CLOSED_ENABLE="${STAGE17_5_FAIL_CLOSED_ENABLE:-1}"
STAGE17_5_CONTACT_PLACEHOLDER_ENABLE="${STAGE17_5_CONTACT_PLACEHOLDER_ENABLE:-1}"
STAGE17_5_DIAGNOSTIC_ONLY="${STAGE17_5_DIAGNOSTIC_ONLY:-1}"
STAGE17_5_Y_MIN="${STAGE17_5_Y_MIN:--1.0}"
STAGE17_5_Y_MAX="${STAGE17_5_Y_MAX:-1.0}"
STAGE17_5_EFFECTIVE_FIBRE_RADIUS="${STAGE17_5_EFFECTIVE_FIBRE_RADIUS:-1.0e-3}"
STAGE17_5_MIN_WALL_CLEARANCE="${STAGE17_5_MIN_WALL_CLEARANCE:-1.0e-4}"
STAGE17_5_WARNING_WALL_CLEARANCE="${STAGE17_5_WARNING_WALL_CLEARANCE:-1.0e-3}"
STAGE17_5_PENETRATION_TOLERANCE="${STAGE17_5_PENETRATION_TOLERANCE:-1.0e-12}"
STAGE17_5_NPTS="${STAGE17_5_NPTS:-8}"
STAGE17_5_TEST_CASE="${STAGE17_5_TEST_CASE:-all_clear}"
STAGE17_5_ZERO_TOL="${STAGE17_5_ZERO_TOL:-1.0e-14}"

OUT_DIR="stage17_outputs"
OUT_FILE="${OUT_DIR}/fibre_stage17_5_near_wall_contact_state.dat"
REASONS_FILE="${OUT_DIR}/stage17_5_near_wall_contact_state_reasons.tmp"

mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

# Accept standard project environment variables without using them to build or run
# physics during Stage 17.5.
printf 'stage17_5_build_dir %s\n' "${BUILD_DIR}" > /dev/null
printf 'stage17_5_mpiexec %s\n' "${MPIEXEC}" > /dev/null
printf 'stage17_5_mpiexec_flags %s\n' "${MPIEXEC_FLAGS}" > /dev/null
printf 'stage17_5_decomp2d_root %s\n' "${DECOMP2D_ROOT}" > /dev/null

if [ "${STAGE17_5_RUN_STAGE17_4}" = "1" ]; then
  bash stage17_checks/run_stage17_4_boundary_containment_fail_closed.sh >/dev/null 2>>"${REASONS_FILE}"
  stage17_4_status=$?
  if [ "${stage17_4_status}" -ne 0 ]; then
    printf '%s\n' 'stage17_4_optional_prerun_not_pass' >> "${REASONS_FILE}"
  fi
fi

python3 stage17_checks/assert_stage17_5_near_wall_contact_state.py \
  --output "${OUT_FILE}" \
  --require-stage17-4 "${STAGE17_5_REQUIRE_STAGE17_4}" \
  --accept-stage17-4-closed-evidence "${STAGE17_5_ACCEPT_STAGE17_4_CLOSED_EVIDENCE}" \
  --enable "${STAGE17_5_ENABLE}" \
  --wall-safety-enable "${STAGE17_5_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_5_BOUNDARY_CHECK_ENABLE}" \
  --fail-closed-enable "${STAGE17_5_FAIL_CLOSED_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_5_CONTACT_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_5_DIAGNOSTIC_ONLY}" \
  --y-min "${STAGE17_5_Y_MIN}" \
  --y-max "${STAGE17_5_Y_MAX}" \
  --effective-fibre-radius "${STAGE17_5_EFFECTIVE_FIBRE_RADIUS}" \
  --min-wall-clearance "${STAGE17_5_MIN_WALL_CLEARANCE}" \
  --warning-wall-clearance "${STAGE17_5_WARNING_WALL_CLEARANCE}" \
  --penetration-tolerance "${STAGE17_5_PENETRATION_TOLERANCE}" \
  --npts "${STAGE17_5_NPTS}" \
  --test-case "${STAGE17_5_TEST_CASE}" \
  --zero-tol "${STAGE17_5_ZERO_TOL}" \
  --wrapper-reasons-file "${REASONS_FILE}"
helper_status=$?

rm -f "${REASONS_FILE}"

if [ "${helper_status}" -eq 0 ]; then
  printf '%s\n' 'STAGE 17.5 NEAR-WALL CONTACT STATE VERDICT: PASS'
  printf '%s\n' 'STAGE 17.5 FINAL VERDICT: PASS'
else
  printf '%s\n' 'STAGE 17.5 NEAR-WALL CONTACT STATE VERDICT: FAIL'
  printf '%s\n' 'STAGE 17.5 FINAL VERDICT: FAIL'
fi
exit "${helper_status}"
