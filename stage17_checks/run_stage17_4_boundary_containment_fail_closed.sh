#!/usr/bin/env bash
set -u

# Stage 17.4 boundary-containment and wall-penetration fail-closed wrapper.
# This wrapper creates diagnostic output only. It builds no target, runs no MPI,
# applies no force, classifies no contact state, and runs no production physics.

DECOMP2D_ROOT="${DECOMP2D_ROOT:-}"
BUILD_DIR="${BUILD_DIR:-build_stage9}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPIEXEC_FLAGS="${MPIEXEC_FLAGS:-}"
STAGE17_4_RUN_STAGE17_3="${STAGE17_4_RUN_STAGE17_3:-0}"
STAGE17_4_REQUIRE_STAGE17_3="${STAGE17_4_REQUIRE_STAGE17_3:-1}"
STAGE17_4_ACCEPT_STAGE17_3_CLOSED_EVIDENCE="${STAGE17_4_ACCEPT_STAGE17_3_CLOSED_EVIDENCE:-1}"
STAGE17_4_ENABLE="${STAGE17_4_ENABLE:-1}"
STAGE17_4_WALL_SAFETY_ENABLE="${STAGE17_4_WALL_SAFETY_ENABLE:-1}"
STAGE17_4_BOUNDARY_CHECK_ENABLE="${STAGE17_4_BOUNDARY_CHECK_ENABLE:-1}"
STAGE17_4_FAIL_CLOSED_ENABLE="${STAGE17_4_FAIL_CLOSED_ENABLE:-1}"
STAGE17_4_DIAGNOSTIC_ONLY="${STAGE17_4_DIAGNOSTIC_ONLY:-1}"
STAGE17_4_Y_MIN="${STAGE17_4_Y_MIN:--1.0}"
STAGE17_4_Y_MAX="${STAGE17_4_Y_MAX:-1.0}"
STAGE17_4_EFFECTIVE_FIBRE_RADIUS="${STAGE17_4_EFFECTIVE_FIBRE_RADIUS:-1.0e-3}"
STAGE17_4_PENETRATION_TOLERANCE="${STAGE17_4_PENETRATION_TOLERANCE:-1.0e-12}"
STAGE17_4_NPTS="${STAGE17_4_NPTS:-8}"
STAGE17_4_TEST_CASE="${STAGE17_4_TEST_CASE:-contained_clear}"
STAGE17_4_ZERO_TOL="${STAGE17_4_ZERO_TOL:-1.0e-14}"

OUT_DIR="stage17_outputs"
OUT_FILE="${OUT_DIR}/fibre_stage17_4_boundary_containment_fail_closed.dat"
REASONS_FILE="${OUT_DIR}/stage17_4_boundary_containment_fail_closed_reasons.tmp"

mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

# Accept standard project environment variables without using them to build or run
# physics during Stage 17.4.
printf 'stage17_4_build_dir %s\n' "${BUILD_DIR}" > /dev/null
printf 'stage17_4_mpiexec %s\n' "${MPIEXEC}" > /dev/null
printf 'stage17_4_mpiexec_flags %s\n' "${MPIEXEC_FLAGS}" > /dev/null
printf 'stage17_4_decomp2d_root %s\n' "${DECOMP2D_ROOT}" > /dev/null

if [ "${STAGE17_4_RUN_STAGE17_3}" = "1" ]; then
  bash stage17_checks/run_stage17_3_effective_radius_wall_clearance.sh >/dev/null 2>>"${REASONS_FILE}"
  stage17_3_status=$?
  if [ "${stage17_3_status}" -ne 0 ]; then
    printf '%s\n' 'stage17_3_optional_prerun_not_pass' >> "${REASONS_FILE}"
  fi
fi

python3 stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py \
  --output "${OUT_FILE}" \
  --require-stage17-3 "${STAGE17_4_REQUIRE_STAGE17_3}" \
  --accept-stage17-3-closed-evidence "${STAGE17_4_ACCEPT_STAGE17_3_CLOSED_EVIDENCE}" \
  --enable "${STAGE17_4_ENABLE}" \
  --wall-safety-enable "${STAGE17_4_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_4_BOUNDARY_CHECK_ENABLE}" \
  --fail-closed-enable "${STAGE17_4_FAIL_CLOSED_ENABLE}" \
  --diagnostic-only "${STAGE17_4_DIAGNOSTIC_ONLY}" \
  --y-min "${STAGE17_4_Y_MIN}" \
  --y-max "${STAGE17_4_Y_MAX}" \
  --effective-fibre-radius "${STAGE17_4_EFFECTIVE_FIBRE_RADIUS}" \
  --penetration-tolerance "${STAGE17_4_PENETRATION_TOLERANCE}" \
  --npts "${STAGE17_4_NPTS}" \
  --test-case "${STAGE17_4_TEST_CASE}" \
  --zero-tol "${STAGE17_4_ZERO_TOL}" \
  --wrapper-reasons-file "${REASONS_FILE}"
helper_status=$?

rm -f "${REASONS_FILE}"

if [ "${helper_status}" -eq 0 ]; then
  printf '%s\n' 'STAGE 17.4 BOUNDARY CONTAINMENT FAIL-CLOSED VERDICT: PASS'
  printf '%s\n' 'STAGE 17.4 FINAL VERDICT: PASS'
else
  printf '%s\n' 'STAGE 17.4 BOUNDARY CONTAINMENT FAIL-CLOSED VERDICT: FAIL'
  printf '%s\n' 'STAGE 17.4 FINAL VERDICT: FAIL'
fi
exit "${helper_status}"
