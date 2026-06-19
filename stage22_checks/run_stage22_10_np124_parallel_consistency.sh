#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${STAGE22_10_NP_VALUES:=1,2,4}"
: "${STAGE22_10_ALLOWED_NP_VALUES:=1,2,4}"
: "${STAGE22_10_G1_ENABLE:=1}"
: "${STAGE22_10_G2_ENABLE:=0}"
: "${STAGE22_10_G3_ENABLE:=0}"
: "${STAGE22_10_PRODUCTION_RESTART_TEST_ALLOWED:=0}"
: "${STAGE22_10_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED:=0}"
: "${STAGE22_10_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED:=0}"
: "${STAGE22_10_PHYSICAL_CASE_SETTINGS_IDENTICAL:=1}"

[[ "${STAGE22_10_NP_VALUES}" == "1,2,4" ]] || { echo "Stage 22.10 fail-closed: only np=1,2,4 are supported" >&2; exit 2; }
[[ "${STAGE22_10_ALLOWED_NP_VALUES}" == "1,2,4" ]] || { echo "Stage 22.10 fail-closed: allowed np values must be 1,2,4" >&2; exit 2; }
[[ "${STAGE22_10_G1_ENABLE}" == "1" ]] || { echo "Stage 22.10 fail-closed: G1 must be enabled" >&2; exit 2; }
[[ "${STAGE22_10_G2_ENABLE}" == "0" && "${STAGE22_10_G3_ENABLE}" == "0" ]] || { echo "Stage 22.10 fail-closed: G2/G3 are forbidden" >&2; exit 2; }
[[ "${STAGE22_10_PRODUCTION_RESTART_TEST_ALLOWED}" == "0" ]] || { echo "Stage 22.10 fail-closed: restart testing is forbidden" >&2; exit 2; }
[[ "${STAGE22_10_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED}" == "0" ]] || { echo "Stage 22.10 fail-closed: uncontrolled multifibre is forbidden" >&2; exit 2; }
[[ "${STAGE22_10_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED}" == "0" ]] || { echo "Stage 22.10 fail-closed: source/schema modification is forbidden" >&2; exit 2; }
[[ "${STAGE22_10_PHYSICAL_CASE_SETTINGS_IDENTICAL}" == "1" ]] || { echo "Stage 22.10 fail-closed: physical settings differ across np values" >&2; exit 2; }

mkdir -p "${REPO_ROOT}/stage22_outputs" \
         "${REPO_ROOT}/stage22_cases/stage22_10_np124_parallel_consistency" \
         "${REPO_ROOT}/stage22_outputs/stage22_10_np124_parallel_consistency"

python3 "${SCRIPT_DIR}/assert_stage22_10_np124_parallel_consistency.py"
