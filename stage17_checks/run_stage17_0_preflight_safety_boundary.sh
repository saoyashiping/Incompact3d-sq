#!/usr/bin/env bash
set -u

# Stage 17.0 preflight wrapper: evidence aggregation only.
# It intentionally builds nothing, runs no MPI command, and executes no physics.

DECOMP2D_ROOT="${DECOMP2D_ROOT:-}"
BUILD_DIR="${BUILD_DIR:-build_stage9}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPIEXEC_FLAGS="${MPIEXEC_FLAGS:-}"
STAGE17_0_REQUIRE_STAGE16_CLOSED="${STAGE17_0_REQUIRE_STAGE16_CLOSED:-1}"
STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE="${STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE:-1}"
STAGE17_0_ENABLE="${STAGE17_0_ENABLE:-1}"
STAGE17_0_DIAGNOSTIC_ONLY="${STAGE17_0_DIAGNOSTIC_ONLY:-1}"

OUT_DIR="stage17_outputs"
OUT_FILE="${OUT_DIR}/fibre_stage17_0_preflight_safety_boundary.dat"
REASONS_FILE="${OUT_DIR}/stage17_0_preflight_wrapper_reasons.tmp"

mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

# Keep the configured environment visible without using it for build/run work.
# This confirms the wrapper accepts the standard project variables while preserving
# Stage 17.0's no-physics/no-build boundary.
printf 'stage17_0_wrapper_build_dir %s\n' "${BUILD_DIR}" > /dev/null
printf 'stage17_0_wrapper_mpiexec %s\n' "${MPIEXEC}" > /dev/null
printf 'stage17_0_wrapper_mpiexec_flags %s\n' "${MPIEXEC_FLAGS}" > /dev/null
printf 'stage17_0_wrapper_decomp2d_root %s\n' "${DECOMP2D_ROOT}" > /dev/null

if [ "${STAGE17_0_ENABLE}" != "1" ]; then
  printf '%s\n' 'stage17_0_enable_not_set_to_1' >> "${REASONS_FILE}"
fi
if [ "${STAGE17_0_DIAGNOSTIC_ONLY}" != "1" ]; then
  printf '%s\n' 'stage17_0_diagnostic_only_not_set_to_1' >> "${REASONS_FILE}"
fi
# Do not fail here if STAGE16_CLOSED.md is absent in a fresh source archive.
# The helper can accept read-only Stage 16.12 closure machinery when
# STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE=1, without modifying closed files.
if [ "${STAGE17_0_REQUIRE_STAGE16_CLOSED}" = "1" ] && \
   [ "${STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE}" != "1" ] && \
   [ ! -s stage16_checks/STAGE16_CLOSED.md ]; then
  printf '%s\n' 'stage16_closed_file_missing_or_empty_before_helper' >> "${REASONS_FILE}"
fi

python3 stage17_checks/assert_stage17_0_preflight_safety_boundary.py \
  --output "${OUT_FILE}" \
  --require-stage16-closed "${STAGE17_0_REQUIRE_STAGE16_CLOSED}" \
  --accept-stage16-closed-evidence "${STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE}" \
  --enable "${STAGE17_0_ENABLE}" \
  --diagnostic-only "${STAGE17_0_DIAGNOSTIC_ONLY}" \
  --wrapper-reasons-file "${REASONS_FILE}"
helper_status=$?

rm -f "${REASONS_FILE}"

if [ "${helper_status}" -eq 0 ]; then
  printf '%s\n' 'STAGE 17.0 PREFLIGHT SAFETY BOUNDARY VERDICT: PASS'
  printf '%s\n' 'STAGE 17.0 FINAL VERDICT: PASS'
else
  printf '%s\n' 'STAGE 17.0 PREFLIGHT SAFETY BOUNDARY VERDICT: FAIL'
  printf '%s\n' 'STAGE 17.0 FINAL VERDICT: FAIL'
fi
exit "${helper_status}"
