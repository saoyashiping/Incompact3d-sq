#!/usr/bin/env bash
set -u

# Stage 17.2 channel wall/domain-boundary metadata wrapper.
# This wrapper creates diagnostic output only. It builds no target, runs no MPI,
# computes no wall clearance, and runs no production physics.

DECOMP2D_ROOT="${DECOMP2D_ROOT:-}"
BUILD_DIR="${BUILD_DIR:-build_stage9}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPIEXEC_FLAGS="${MPIEXEC_FLAGS:-}"
STAGE17_2_RUN_STAGE17_1="${STAGE17_2_RUN_STAGE17_1:-0}"
STAGE17_2_REQUIRE_STAGE17_1="${STAGE17_2_REQUIRE_STAGE17_1:-1}"
STAGE17_2_ACCEPT_STAGE17_1_CLOSED_EVIDENCE="${STAGE17_2_ACCEPT_STAGE17_1_CLOSED_EVIDENCE:-1}"
STAGE17_2_ENABLE="${STAGE17_2_ENABLE:-1}"
STAGE17_2_WALL_SAFETY_ENABLE="${STAGE17_2_WALL_SAFETY_ENABLE:-1}"
STAGE17_2_BOUNDARY_CHECK_ENABLE="${STAGE17_2_BOUNDARY_CHECK_ENABLE:-1}"
STAGE17_2_DIAGNOSTIC_ONLY="${STAGE17_2_DIAGNOSTIC_ONLY:-1}"
STAGE17_2_Y_MIN="${STAGE17_2_Y_MIN:--1.0}"
STAGE17_2_Y_MAX="${STAGE17_2_Y_MAX:-1.0}"
STAGE17_2_WALL_NORMAL_DIRECTION="${STAGE17_2_WALL_NORMAL_DIRECTION:-y}"
STAGE17_2_X_BOUNDARY_POLICY="${STAGE17_2_X_BOUNDARY_POLICY:-periodic}"
STAGE17_2_Y_BOUNDARY_POLICY="${STAGE17_2_Y_BOUNDARY_POLICY:-wall_bounded}"
STAGE17_2_Z_BOUNDARY_POLICY="${STAGE17_2_Z_BOUNDARY_POLICY:-periodic}"

OUT_DIR="stage17_outputs"
OUT_FILE="${OUT_DIR}/fibre_stage17_2_channel_wall_domain_boundary.dat"
REASONS_FILE="${OUT_DIR}/stage17_2_channel_wall_domain_boundary_reasons.tmp"

mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

# Accept standard project environment variables without using them to build or run
# physics during Stage 17.2.
printf 'stage17_2_build_dir %s\n' "${BUILD_DIR}" > /dev/null
printf 'stage17_2_mpiexec %s\n' "${MPIEXEC}" > /dev/null
printf 'stage17_2_mpiexec_flags %s\n' "${MPIEXEC_FLAGS}" > /dev/null
printf 'stage17_2_decomp2d_root %s\n' "${DECOMP2D_ROOT}" > /dev/null

if [ "${STAGE17_2_RUN_STAGE17_1}" = "1" ]; then
  bash stage17_checks/run_stage17_1_wall_contact_safety_config.sh >/dev/null 2>>"${REASONS_FILE}"
  stage17_1_status=$?
  if [ "${stage17_1_status}" -ne 0 ]; then
    printf '%s\n' 'stage17_1_optional_prerun_not_pass' >> "${REASONS_FILE}"
  fi
fi

python3 stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py \
  --output "${OUT_FILE}" \
  --require-stage17-1 "${STAGE17_2_REQUIRE_STAGE17_1}" \
  --accept-stage17-1-closed-evidence "${STAGE17_2_ACCEPT_STAGE17_1_CLOSED_EVIDENCE}" \
  --enable "${STAGE17_2_ENABLE}" \
  --wall-safety-enable "${STAGE17_2_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_2_BOUNDARY_CHECK_ENABLE}" \
  --diagnostic-only "${STAGE17_2_DIAGNOSTIC_ONLY}" \
  --y-min "${STAGE17_2_Y_MIN}" \
  --y-max "${STAGE17_2_Y_MAX}" \
  --wall-normal-direction "${STAGE17_2_WALL_NORMAL_DIRECTION}" \
  --x-boundary-policy "${STAGE17_2_X_BOUNDARY_POLICY}" \
  --y-boundary-policy "${STAGE17_2_Y_BOUNDARY_POLICY}" \
  --z-boundary-policy "${STAGE17_2_Z_BOUNDARY_POLICY}" \
  --wrapper-reasons-file "${REASONS_FILE}"
helper_status=$?

rm -f "${REASONS_FILE}"

if [ "${helper_status}" -eq 0 ]; then
  printf '%s\n' 'STAGE 17.2 CHANNEL WALL DOMAIN BOUNDARY VERDICT: PASS'
  printf '%s\n' 'STAGE 17.2 FINAL VERDICT: PASS'
else
  printf '%s\n' 'STAGE 17.2 CHANNEL WALL DOMAIN BOUNDARY VERDICT: FAIL'
  printf '%s\n' 'STAGE 17.2 FINAL VERDICT: FAIL'
fi
exit "${helper_status}"
