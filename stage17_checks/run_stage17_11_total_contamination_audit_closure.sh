#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE17_11_RUN_STAGE17_10:=0}"
: "${STAGE17_11_REQUIRE_STAGE17_10:=1}"
: "${STAGE17_11_ACCEPT_STAGE17_10_CLOSED_EVIDENCE:=1}"
: "${STAGE17_11_REQUIRE_STAGE16_CLOSED_LOOP:=1}"
: "${STAGE17_11_ACCEPT_STAGE16_CLOSED_LOOP_EVIDENCE:=1}"
: "${STAGE17_11_ENABLE:=1}"
: "${STAGE17_11_WALL_SAFETY_ENABLE:=1}"
: "${STAGE17_11_CONTACT_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_11_FIBRE_COLLISION_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_11_DIAGNOSTIC_ONLY:=1}"
: "${STAGE17_11_REQUIRE_PARALLEL_COMPAT:=1}"
: "${STAGE17_11_REQUIRE_RESTART_IO_COMPAT:=1}"
: "${STAGE17_11_REQUIRE_STATS_VISU_COARSE_IO_COMPAT:=1}"
: "${STAGE17_11_MAX_CONTACT_FORCE_NORM:=0.0}"
: "${STAGE17_11_MAX_CONTACT_RHS_NORM:=0.0}"
: "${STAGE17_11_MAX_CONTACT_STRUCTURE_UPDATE_NORM:=0.0}"

mkdir -p stage17_outputs

if [ "${STAGE17_11_RUN_STAGE17_10}" = "1" ]; then
  if [ ! -x stage17_checks/run_stage17_10_parallel_restart_io_wall_safety.sh ]; then
    echo "STAGE 17.11 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: FAIL"
    echo "STAGE 17.11 FINAL VERDICT: FAIL"
    echo "reason stage17_10_wrapper_missing_or_not_executable"
    exit 1
  fi
  STAGE17_10_REQUIRE_STAGE17_9="${STAGE17_11_REQUIRE_STAGE17_10}" \
    stage17_checks/run_stage17_10_parallel_restart_io_wall_safety.sh
fi

set +e
python3 stage17_checks/assert_stage17_11_total_contamination_audit_closure.py \
  --stage17-11-enable "${STAGE17_11_ENABLE}" \
  --wall-safety-enable "${STAGE17_11_WALL_SAFETY_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_11_CONTACT_PLACEHOLDER_ENABLE}" \
  --fibre-collision-placeholder-enable "${STAGE17_11_FIBRE_COLLISION_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_11_DIAGNOSTIC_ONLY}" \
  --require-parallel-compat "${STAGE17_11_REQUIRE_PARALLEL_COMPAT}" \
  --require-restart-io-compat "${STAGE17_11_REQUIRE_RESTART_IO_COMPAT}" \
  --require-stats-visu-coarse-io-compat "${STAGE17_11_REQUIRE_STATS_VISU_COARSE_IO_COMPAT}" \
  --max-contact-force-norm "${STAGE17_11_MAX_CONTACT_FORCE_NORM}" \
  --max-contact-rhs-norm "${STAGE17_11_MAX_CONTACT_RHS_NORM}" \
  --max-contact-structure-update-norm "${STAGE17_11_MAX_CONTACT_STRUCTURE_UPDATE_NORM}" \
  --accept-stage17-10-closed-evidence "${STAGE17_11_ACCEPT_STAGE17_10_CLOSED_EVIDENCE}" \
  --accept-stage16-closed-loop-evidence "${STAGE17_11_ACCEPT_STAGE16_CLOSED_LOOP_EVIDENCE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 17.11 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: PASS"
  echo "STAGE 17.11 FINAL VERDICT: PASS"
else
  echo "STAGE 17.11 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: FAIL"
  echo "STAGE 17.11 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
