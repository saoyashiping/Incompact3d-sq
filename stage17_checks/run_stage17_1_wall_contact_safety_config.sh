#!/usr/bin/env bash
set -u

# Stage 17.1 wall/contact safety configuration wrapper.
# This wrapper creates diagnostic output only. It builds no target, runs no MPI,
# and runs no production physics.

DECOMP2D_ROOT="${DECOMP2D_ROOT:-}"
BUILD_DIR="${BUILD_DIR:-build_stage9}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPIEXEC_FLAGS="${MPIEXEC_FLAGS:-}"
STAGE17_1_RUN_STAGE17_0="${STAGE17_1_RUN_STAGE17_0:-0}"
STAGE17_1_REQUIRE_STAGE17_0="${STAGE17_1_REQUIRE_STAGE17_0:-1}"
STAGE17_1_ACCEPT_STAGE17_0_CLOSED_EVIDENCE="${STAGE17_1_ACCEPT_STAGE17_0_CLOSED_EVIDENCE:-1}"
STAGE17_1_ENABLE="${STAGE17_1_ENABLE:-1}"
STAGE17_1_WALL_SAFETY_ENABLE="${STAGE17_1_WALL_SAFETY_ENABLE:-1}"
STAGE17_1_BOUNDARY_CHECK_ENABLE="${STAGE17_1_BOUNDARY_CHECK_ENABLE:-1}"
STAGE17_1_FAIL_CLOSED_ENABLE="${STAGE17_1_FAIL_CLOSED_ENABLE:-1}"
STAGE17_1_CONTACT_PLACEHOLDER_ENABLE="${STAGE17_1_CONTACT_PLACEHOLDER_ENABLE:-1}"
STAGE17_1_FIBRE_COLLISION_PLACEHOLDER_ENABLE="${STAGE17_1_FIBRE_COLLISION_PLACEHOLDER_ENABLE:-0}"
STAGE17_1_EFFECTIVE_FIBRE_RADIUS="${STAGE17_1_EFFECTIVE_FIBRE_RADIUS:-1.0e-3}"
STAGE17_1_MIN_WALL_CLEARANCE="${STAGE17_1_MIN_WALL_CLEARANCE:-1.0e-4}"
STAGE17_1_WARNING_WALL_CLEARANCE="${STAGE17_1_WARNING_WALL_CLEARANCE:-1.0e-3}"
STAGE17_1_PENETRATION_TOLERANCE="${STAGE17_1_PENETRATION_TOLERANCE:-1.0e-12}"
STAGE17_1_DIAGNOSTIC_ONLY="${STAGE17_1_DIAGNOSTIC_ONLY:-1}"

OUT_DIR="stage17_outputs"
OUT_FILE="${OUT_DIR}/fibre_stage17_1_wall_contact_safety_config.dat"
REASONS_FILE="${OUT_DIR}/stage17_1_wall_contact_safety_config_reasons.tmp"

mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

# Accept standard project environment variables without using them to build or run
# physics during Stage 17.1.
printf 'stage17_1_build_dir %s\n' "${BUILD_DIR}" > /dev/null
printf 'stage17_1_mpiexec %s\n' "${MPIEXEC}" > /dev/null
printf 'stage17_1_mpiexec_flags %s\n' "${MPIEXEC_FLAGS}" > /dev/null
printf 'stage17_1_decomp2d_root %s\n' "${DECOMP2D_ROOT}" > /dev/null

if [ "${STAGE17_1_RUN_STAGE17_0}" = "1" ]; then
  bash stage17_checks/run_stage17_0_preflight_safety_boundary.sh >/dev/null 2>>"${REASONS_FILE}"
  stage17_0_status=$?
  if [ "${stage17_0_status}" -ne 0 ]; then
    printf '%s\n' 'stage17_0_optional_prerun_not_pass' >> "${REASONS_FILE}"
  fi
fi


python3 stage17_checks/assert_stage17_1_wall_contact_safety_config.py \
  --output "${OUT_FILE}" \
  --require-stage17-0 "${STAGE17_1_REQUIRE_STAGE17_0}" \
  --accept-stage17-0-closed-evidence "${STAGE17_1_ACCEPT_STAGE17_0_CLOSED_EVIDENCE}" \
  --enable "${STAGE17_1_ENABLE}" \
  --wall-safety-enable "${STAGE17_1_WALL_SAFETY_ENABLE}" \
  --boundary-check-enable "${STAGE17_1_BOUNDARY_CHECK_ENABLE}" \
  --fail-closed-enable "${STAGE17_1_FAIL_CLOSED_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_1_CONTACT_PLACEHOLDER_ENABLE}" \
  --fibre-collision-placeholder-enable "${STAGE17_1_FIBRE_COLLISION_PLACEHOLDER_ENABLE}" \
  --effective-fibre-radius "${STAGE17_1_EFFECTIVE_FIBRE_RADIUS}" \
  --min-wall-clearance "${STAGE17_1_MIN_WALL_CLEARANCE}" \
  --warning-wall-clearance "${STAGE17_1_WARNING_WALL_CLEARANCE}" \
  --penetration-tolerance "${STAGE17_1_PENETRATION_TOLERANCE}" \
  --diagnostic-only "${STAGE17_1_DIAGNOSTIC_ONLY}" \
  --wrapper-reasons-file "${REASONS_FILE}"
helper_status=$?

rm -f "${REASONS_FILE}"

if [ "${helper_status}" -eq 0 ]; then
  printf '%s\n' 'STAGE 17.1 WALL CONTACT SAFETY CONFIG VERDICT: PASS'
  printf '%s\n' 'STAGE 17.1 FINAL VERDICT: PASS'
else
  printf '%s\n' 'STAGE 17.1 WALL CONTACT SAFETY CONFIG VERDICT: FAIL'
  printf '%s\n' 'STAGE 17.1 FINAL VERDICT: FAIL'
fi
exit "${helper_status}"
