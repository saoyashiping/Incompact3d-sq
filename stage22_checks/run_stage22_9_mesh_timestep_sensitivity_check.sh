#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${STAGE22_9_NP:=1}"
: "${STAGE22_9_NP2_ALLOWED:=0}"
: "${STAGE22_9_NP4_ALLOWED:=0}"
: "${STAGE22_9_G1_ENABLE:=1}"
: "${STAGE22_9_G2_ENABLE:=1}"
: "${STAGE22_9_G3_OPTIONAL_ENABLE:=0}"
: "${STAGE22_9_PRODUCTION_RESTART_TEST_ALLOWED:=0}"
: "${STAGE22_9_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED:=0}"
: "${STAGE22_9_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED:=0}"

[[ "${STAGE22_9_NP}" == "1" ]] || { echo "Stage 22.9 fail-closed: np must be exactly 1" >&2; exit 2; }
[[ "${STAGE22_9_NP2_ALLOWED}" == "0" ]] || { echo "Stage 22.9 fail-closed: np=2 is forbidden" >&2; exit 2; }
[[ "${STAGE22_9_NP4_ALLOWED}" == "0" ]] || { echo "Stage 22.9 fail-closed: np=4 is forbidden" >&2; exit 2; }
[[ "${STAGE22_9_G1_ENABLE}" == "1" && "${STAGE22_9_G2_ENABLE}" == "1" ]] || { echo "Stage 22.9 fail-closed: G1 and G2 must be enabled" >&2; exit 2; }
[[ "${STAGE22_9_G3_OPTIONAL_ENABLE}" == "0" ]] || { echo "Stage 22.9 fail-closed: G3 requires a future explicit instruction" >&2; exit 2; }
[[ "${STAGE22_9_PRODUCTION_RESTART_TEST_ALLOWED}" == "0" ]] || { echo "Stage 22.9 fail-closed: restart testing is forbidden" >&2; exit 2; }
[[ "${STAGE22_9_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED}" == "0" ]] || { echo "Stage 22.9 fail-closed: uncontrolled multifibre is forbidden" >&2; exit 2; }
[[ "${STAGE22_9_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED}" == "0" ]] || { echo "Stage 22.9 fail-closed: source/schema modification is forbidden" >&2; exit 2; }

mkdir -p "${REPO_ROOT}/stage22_outputs" \
         "${REPO_ROOT}/stage22_cases/stage22_9_mesh_timestep_sensitivity_check" \
         "${REPO_ROOT}/stage22_outputs/stage22_9_mesh_timestep_sensitivity_check"

python3 "${SCRIPT_DIR}/assert_stage22_9_mesh_timestep_sensitivity_check.py"
