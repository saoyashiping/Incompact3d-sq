#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
: "${STAGE22_7_NP:=1}"
: "${STAGE22_7_LAMBDA_CONTACT:=1.0}"
: "${STAGE22_7_N_FIBRE:=1}"
: "${STAGE22_7_WALL_CONTACT_FORCE_ENABLE:=1}"
: "${STAGE22_7_COLLISION_FORCE_ENABLE:=0}"
: "${STAGE22_7_PRODUCTION_MULTIFIBRE_ENABLE:=0}"
if [ "${STAGE22_7_NP}" != "1" ]; then echo "Stage 22.7 fail-closed: np must be exactly 1" >&2; exit 1; fi
if [ "${STAGE22_7_LAMBDA_CONTACT}" != "1.0" ]; then echo "Stage 22.7 fail-closed: lambda_contact must be 1.0" >&2; exit 1; fi
if [ "${STAGE22_7_N_FIBRE}" != "1" ]; then echo "Stage 22.7 fail-closed: n_fibre must be exactly 1" >&2; exit 1; fi
if [ "${STAGE22_7_WALL_CONTACT_FORCE_ENABLE}" != "1" ]; then echo "Stage 22.7 fail-closed: wall contact gate must be enabled" >&2; exit 1; fi
if [ "${STAGE22_7_COLLISION_FORCE_ENABLE}" != "0" ] || [ "${STAGE22_7_PRODUCTION_MULTIFIBRE_ENABLE}" != "0" ]; then echo "Stage 22.7 fail-closed: collision/multifibre gates must be disabled" >&2; exit 1; fi
mkdir -p "${repo_root}/stage22_cases/stage22_7_single_fibre_near_wall_contact_micro_case"
mkdir -p "${repo_root}/stage22_outputs/stage22_7_single_fibre_near_wall_contact_micro_case"
mkdir -p "${repo_root}/stage22_outputs"
python3 "${repo_root}/stage22_checks/assert_stage22_7_single_fibre_near_wall_contact_micro_case.py"
