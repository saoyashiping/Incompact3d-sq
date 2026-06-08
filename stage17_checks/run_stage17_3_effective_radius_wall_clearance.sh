#!/usr/bin/env bash
set -u

# Stage 17.3 effective-radius wall-clearance diagnostic wrapper.
# This wrapper creates diagnostic output only. It builds no target, runs no MPI,
# applies no force, classifies no contact state, and runs no production physics.

DECOMP2D_ROOT="${DECOMP2D_ROOT:-}"
BUILD_DIR="${BUILD_DIR:-build_stage9}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPIEXEC_FLAGS="${MPIEXEC_FLAGS:-}"
STAGE17_3_RUN_STAGE17_2="${STAGE17_3_RUN_STAGE17_2:-0}"
STAGE17_3_REQUIRE_STAGE17_2="${STAGE17_3_REQUIRE_STAGE17_2:-1}"
STAGE17_3_ACCEPT_STAGE17_2_CLOSED_EVIDENCE="${STAGE17_3_ACCEPT_STAGE17_2_CLOSED_EVIDENCE:-1}"
STAGE17_3_ENABLE="${STAGE17_3_ENABLE:-1}"
STAGE17_3_WALL_SAFETY_ENABLE="${STAGE17_3_WALL_SAFETY_ENABLE:-1}"
STAGE17_3_BOUNDARY_CHECK_ENABLE="${STAGE17_3_BOUNDARY_CHECK_ENABLE:-1}"
STAGE17_3_DIAGNOSTIC_ONLY="${STAGE17_3_DIAGNOSTIC_ONLY:-1}"
STAGE17_3_Y_MIN="${STAGE17_3_Y_MIN:--1.0}"
STAGE17_3_Y_MAX="${STAGE17_3_Y_MAX:-1.0}"
STAGE17_3_EFFECTIVE_FIBRE_RADIUS="${STAGE17_3_EFFECTIVE_FIBRE_RADIUS:-1.0e-3}"
STAGE17_3_NPTS="${STAGE17_3_NPTS:-8}"
STAGE17_3_TEST_CASE="${STAGE17_3_TEST_CASE:-centered_clear}"
STAGE17_3_ZERO_TOL="${STAGE17_3_ZERO_TOL:-1.0e-14}"

OUT_DIR="stage17_outputs"
OUT_FILE="${OUT_DIR}/fibre_stage17_3_effective_radius_wall_clearance.dat"
REASONS_FILE="${OUT_DIR}/stage17_3_effective_radius_wall_clearance_reasons.tmp"

mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

# Accept standard project environment variables without using them to build or run
# physics during Stage 17.3.
printf 'stage17_3_build_dir %s\n' "${BUILD_DIR}" > /dev/null
printf 'stage17_3_mpiexec %s\n' "${MPIEXEC}" > /dev/null
printf 'stage17_3_mpiexec_flags %s\n' "${MPIEXEC_FLAGS}" > /dev/null
printf 'stage17_3_decomp2d_root %s\n' "${DECOMP2D_ROOT}" > /dev/null

if [ "${STAGE17_3_RUN_STAGE17_2}" = "1" ]; then
  bash stage17_checks/run_stage17_2_channel_wall_domain_boundary.sh >/dev/null 2>>"${REASONS_FILE}"
  stage17_2_status=$?
  if [ "${stage17_2_status}" -ne 0 ]; then
    printf '%s\n' 'stage17_2_optional_prerun_not_pass' >> "${REASONS_FILE}"
  fi
fi

python3 stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py \
  --output "${OUT_FILE}" \
  --require-stage17-2 "${STAGE17_3_REQUIRE_STAGE17_2}" \
  --accept-stage17-2-closed-evidence "${STAGE17_3_ACCEPT_STAGE17_2_CLOSED_EVIDENCE}" \
  --enable "${STAGE17_3_ENABLE}" \
  --wall-safety-enable "${STAGE17_3_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_3_BOUNDARY_CHECK_ENABLE}" \
  --diagnostic-only "${STAGE17_3_DIAGNOSTIC_ONLY}" \
  --y-min "${STAGE17_3_Y_MIN}" \
  --y-max "${STAGE17_3_Y_MAX}" \
  --effective-fibre-radius "${STAGE17_3_EFFECTIVE_FIBRE_RADIUS}" \
  --npts "${STAGE17_3_NPTS}" \
  --test-case "${STAGE17_3_TEST_CASE}" \
  --zero-tol "${STAGE17_3_ZERO_TOL}" \
  --wrapper-reasons-file "${REASONS_FILE}"
helper_status=$?

rm -f "${REASONS_FILE}"

if [ "${helper_status}" -eq 0 ]; then
  printf '%s\n' 'STAGE 17.3 EFFECTIVE-RADIUS WALL CLEARANCE VERDICT: PASS'
  printf '%s\n' 'STAGE 17.3 FINAL VERDICT: PASS'
else
  printf '%s\n' 'STAGE 17.3 EFFECTIVE-RADIUS WALL CLEARANCE VERDICT: FAIL'
  printf '%s\n' 'STAGE 17.3 FINAL VERDICT: FAIL'
fi
exit "${helper_status}"
