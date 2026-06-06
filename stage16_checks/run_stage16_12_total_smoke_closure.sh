#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

: "${DECOMP2D_ROOT:=}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE16_12_RUN_STAGE16_11:=0}"
: "${STAGE16_12_REQUIRE_STAGE14_CLOSED:=1}"
: "${STAGE16_12_REQUIRE_STAGE15_CLOSED:=1}"
: "${STAGE16_12_REQUIRE_STAGE16_11:=1}"
: "${STAGE16_12_ACCEPT_STAGE16_11_CLOSED_EVIDENCE:=1}"
: "${STAGE16_12_ENABLE:=1}"
: "${STAGE16_12_DIAGNOSTIC_ONLY:=1}"
: "${STAGE16_12_GENERATE_CLOSURE_FILE:=1}"

OUT_DIR="stage16_outputs"
REASONS_FILE="${OUT_DIR}/stage16_12_wrapper_reasons.tmp"
mkdir -p "${OUT_DIR}"
: > "${REASONS_FILE}"

add_reason() {
  printf '%s\n' "$1" >> "${REASONS_FILE}"
}

search_silent() {
  pattern="$1"
  shift
  if command -v rg >/dev/null 2>&1; then
    rg --quiet --fixed-strings -- "$pattern" "$@"
  else
    grep -R -q -F -- "$pattern" "$@"
  fi
}

ensure_build_dir() {
  if [ -d "${BUILD_DIR}" ]; then
    return 0
  fi
  if [ -z "${DECOMP2D_ROOT}" ]; then
    add_reason "build_dir_missing_and_decomp2d_root_unset"
    return 1
  fi
  cmake -S . -B "${BUILD_DIR}" -DDECOMP2D_ROOT="${DECOMP2D_ROOT}"
}

if [ "${STAGE16_12_ENABLE}" != "1" ]; then
  add_reason "stage16_12_enable_not_set"
fi

if [ "${STAGE16_12_DIAGNOSTIC_ONLY}" != "1" ]; then
  add_reason "stage16_12_diagnostic_only_not_set"
fi

# Stage 16.12 is an evidence aggregation/closure stage. It builds nothing and
# runs no physics by default. Users may opt in to refreshing the closed Stage
# 16.11 smoke evidence, but the default path only parses/accepts closed evidence.
if [ "${STAGE16_12_RUN_STAGE16_11}" = "1" ]; then
  if ! bash stage16_checks/run_stage16_11_short_time_stability_smoke.sh; then
    add_reason "stage16_11_refresh_run_failed"
  fi
fi

python3 stage16_checks/assert_stage16_12_total_smoke_closure.py \
  --require-stage14-closed "${STAGE16_12_REQUIRE_STAGE14_CLOSED}" \
  --require-stage15-closed "${STAGE16_12_REQUIRE_STAGE15_CLOSED}" \
  --require-stage16-11 "${STAGE16_12_REQUIRE_STAGE16_11}" \
  --accept-stage16-11-closed-evidence "${STAGE16_12_ACCEPT_STAGE16_11_CLOSED_EVIDENCE}" \
  --enable "${STAGE16_12_ENABLE}" \
  --diagnostic-only "${STAGE16_12_DIAGNOSTIC_ONLY}" \
  --generate-closure-file "${STAGE16_12_GENERATE_CLOSURE_FILE}" \
  --wrapper-reasons-file "${REASONS_FILE}"
