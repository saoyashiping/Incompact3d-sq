#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

: "${DECOMP2D_ROOT:=${REPO_ROOT}}"
: "${BUILD_DIR:=build_stage9}"
: "${MPIEXEC:=mpirun}"
: "${MPIEXEC_FLAGS:=}"
: "${STAGE17_10_RUN_STAGE17_9:=0}"
: "${STAGE17_10_REQUIRE_STAGE17_9:=1}"
: "${STAGE17_10_ACCEPT_STAGE17_9_CLOSED_EVIDENCE:=1}"
: "${STAGE17_10_ENABLE:=1}"
: "${STAGE17_10_WALL_SAFETY_ENABLE:=1}"
: "${STAGE17_10_CONTACT_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_10_FIBRE_COLLISION_PLACEHOLDER_ENABLE:=1}"
: "${STAGE17_10_DIAGNOSTIC_ONLY:=1}"
: "${STAGE17_10_CHECK_NP1:=1}"
: "${STAGE17_10_CHECK_NP2:=1}"
: "${STAGE17_10_CHECK_NP4:=1}"
: "${STAGE17_10_REQUIRE_RESTART_IO_COMPAT:=1}"
: "${STAGE17_10_REQUIRE_STATS_VISU_COARSE_IO_COMPAT:=1}"
: "${STAGE17_10_MAX_PARALLEL_SIGNATURE_DELTA:=1.0e-12}"
: "${STAGE17_10_MAX_RESTART_SIGNATURE_DELTA:=1.0e-12}"
: "${STAGE17_10_MAX_CONTACT_FORCE_NORM:=0.0}"
: "${STAGE17_10_MAX_CONTACT_RHS_NORM:=0.0}"
: "${STAGE17_10_MAX_CONTACT_STRUCTURE_UPDATE_NORM:=0.0}"

mkdir -p stage17_outputs

if [ "${STAGE17_10_RUN_STAGE17_9}" = "1" ]; then
  if [ ! -x stage17_checks/run_stage17_9_closed_loop_wall_contact_compatibility.sh ]; then
    echo "STAGE 17.10 PARALLEL RESTART IO WALL SAFETY VERDICT: FAIL"
    echo "STAGE 17.10 FINAL VERDICT: FAIL"
    echo "reason stage17_9_wrapper_missing_or_not_executable"
    exit 1
  fi
  STAGE17_9_REQUIRE_STAGE17_8="${STAGE17_10_REQUIRE_STAGE17_9}" \
    stage17_checks/run_stage17_9_closed_loop_wall_contact_compatibility.sh
fi

set +e
python3 stage17_checks/assert_stage17_10_parallel_restart_io_wall_safety.py \
  --stage17-10-enable "${STAGE17_10_ENABLE}" \
  --wall-safety-enable "${STAGE17_10_WALL_SAFETY_ENABLE}" \
  --contact-placeholder-enable "${STAGE17_10_CONTACT_PLACEHOLDER_ENABLE}" \
  --fibre-collision-placeholder-enable "${STAGE17_10_FIBRE_COLLISION_PLACEHOLDER_ENABLE}" \
  --diagnostic-only "${STAGE17_10_DIAGNOSTIC_ONLY}" \
  --check-np1 "${STAGE17_10_CHECK_NP1}" \
  --check-np2 "${STAGE17_10_CHECK_NP2}" \
  --check-np4 "${STAGE17_10_CHECK_NP4}" \
  --require-restart-io-compat "${STAGE17_10_REQUIRE_RESTART_IO_COMPAT}" \
  --require-stats-visu-coarse-io-compat "${STAGE17_10_REQUIRE_STATS_VISU_COARSE_IO_COMPAT}" \
  --max-parallel-signature-delta "${STAGE17_10_MAX_PARALLEL_SIGNATURE_DELTA}" \
  --max-restart-signature-delta "${STAGE17_10_MAX_RESTART_SIGNATURE_DELTA}" \
  --max-contact-force-norm "${STAGE17_10_MAX_CONTACT_FORCE_NORM}" \
  --max-contact-rhs-norm "${STAGE17_10_MAX_CONTACT_RHS_NORM}" \
  --max-contact-structure-update-norm "${STAGE17_10_MAX_CONTACT_STRUCTURE_UPDATE_NORM}" \
  --accept-stage17-9-closed-evidence "${STAGE17_10_ACCEPT_STAGE17_9_CLOSED_EVIDENCE}"
helper_status=$?
set -e

if [ "${helper_status}" -eq 0 ]; then
  echo "STAGE 17.10 PARALLEL RESTART IO WALL SAFETY VERDICT: PASS"
  echo "STAGE 17.10 FINAL VERDICT: PASS"
else
  echo "STAGE 17.10 PARALLEL RESTART IO WALL SAFETY VERDICT: FAIL"
  echo "STAGE 17.10 FINAL VERDICT: FAIL"
fi

exit "${helper_status}"
