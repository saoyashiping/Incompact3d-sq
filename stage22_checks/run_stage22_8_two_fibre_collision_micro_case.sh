#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${STAGE22_8_NP:=1}"
: "${STAGE22_8_N_FIBRE:=2}"
: "${STAGE22_8_LAMBDA_CONTACT:=1.0}"
: "${STAGE22_8_COLLISION_FORCE_ENABLE:=1}"
: "${STAGE22_8_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE:=1}"
: "${STAGE22_8_COLLISION_FORCE_APPLY_ENABLE:=1}"
: "${STAGE22_8_WALL_CONTACT_SAFETY_ENABLE:=1}"
: "${STAGE22_8_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED:=0}"

[[ "${STAGE22_8_NP}" == "1" ]] || { echo "Stage 22.8 fail-closed: np must be exactly 1" >&2; exit 2; }
[[ "${STAGE22_8_N_FIBRE}" == "2" ]] || { echo "Stage 22.8 fail-closed: n_fibre must be exactly 2" >&2; exit 2; }
[[ "${STAGE22_8_LAMBDA_CONTACT}" == "1.0" ]] || { echo "Stage 22.8 fail-closed: lambda_contact must be 1.0" >&2; exit 2; }
[[ "${STAGE22_8_COLLISION_FORCE_ENABLE}" == "1" ]] || { echo "Stage 22.8 fail-closed: collision force gate must be enabled" >&2; exit 2; }
[[ "${STAGE22_8_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE}" == "1" ]] || { echo "Stage 22.8 fail-closed: collision candidate gate must be enabled" >&2; exit 2; }
[[ "${STAGE22_8_COLLISION_FORCE_APPLY_ENABLE}" == "1" ]] || { echo "Stage 22.8 fail-closed: collision apply gate must be enabled" >&2; exit 2; }
[[ "${STAGE22_8_WALL_CONTACT_SAFETY_ENABLE}" == "1" ]] || { echo "Stage 22.8 fail-closed: wall contact safety gate must be enabled" >&2; exit 2; }
[[ "${STAGE22_8_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED}" == "0" ]] || { echo "Stage 22.8 fail-closed: uncontrolled multifibre is forbidden" >&2; exit 2; }

mkdir -p "${REPO_ROOT}/stage22_outputs" \
         "${REPO_ROOT}/stage22_cases/stage22_8_two_fibre_collision_micro_case" \
         "${REPO_ROOT}/stage22_outputs/stage22_8_two_fibre_collision_micro_case"

python3 "${SCRIPT_DIR}/assert_stage22_8_two_fibre_collision_micro_case.py"
