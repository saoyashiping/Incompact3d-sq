#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${STAGE22_11_NP_PRIMARY:=1}"
: "${STAGE22_11_NP4_ENABLE:=0}"
: "${STAGE22_11_G1_ENABLE:=1}"
: "${STAGE22_11_G2_ENABLE:=0}"
: "${STAGE22_11_G3_ENABLE:=0}"
: "${STAGE22_11_RESTART_TEST_ENABLE:=1}"
: "${STAGE22_11_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED:=0}"
: "${STAGE22_11_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED:=0}"
: "${STAGE22_11_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED:=0}"
: "${STAGE22_11_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED:=0}"
: "${STAGE22_11_FRESH_RESTART_PHYSICAL_SETTINGS_IDENTICAL:=1}"
: "${STAGE22_11_RESTART_FILE_READABLE:=1}"
: "${STAGE22_11_NAN_INF_DETECTED:=0}"
: "${STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATION_ALLOWED:=0}"
: "${STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATED:=0}"

[[ "${STAGE22_11_NP_PRIMARY}" == "1" ]] || { echo "Stage 22.11 fail-closed: np=1 primary is required" >&2; exit 2; }
[[ "${STAGE22_11_NP4_ENABLE}" == "0" ]] || { echo "Stage 22.11 fail-closed: np=4 is disabled by default" >&2; exit 2; }
[[ "${STAGE22_11_G1_ENABLE}" == "1" && "${STAGE22_11_G2_ENABLE}" == "0" && "${STAGE22_11_G3_ENABLE}" == "0" ]] || { echo "Stage 22.11 fail-closed: only G1 is allowed" >&2; exit 2; }
[[ "${STAGE22_11_RESTART_TEST_ENABLE}" == "1" ]] || { echo "Stage 22.11 fail-closed: restart audit must be enabled" >&2; exit 2; }
[[ "${STAGE22_11_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED}" == "0" && "${STAGE22_11_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED}" == "0" && "${STAGE22_11_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED}" == "0" ]] || { echo "Stage 22.11 fail-closed: schema modification is forbidden" >&2; exit 2; }
[[ "${STAGE22_11_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED}" == "0" ]] || { echo "Stage 22.11 fail-closed: source/schema modification attempted" >&2; exit 2; }
[[ "${STAGE22_11_FRESH_RESTART_PHYSICAL_SETTINGS_IDENTICAL}" == "1" ]] || { echo "Stage 22.11 fail-closed: fresh/restart physical settings differ" >&2; exit 2; }
[[ "${STAGE22_11_RESTART_FILE_READABLE}" == "1" ]] || { echo "Stage 22.11 fail-closed: restart file is not readable" >&2; exit 2; }
[[ "${STAGE22_11_NAN_INF_DETECTED}" == "0" ]] || { echo "Stage 22.11 fail-closed: NaN/Inf detected" >&2; exit 2; }
[[ "${STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATION_ALLOWED}" == "0" && "${STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATED}" == "0" ]] || { echo "Stage 22.11 fail-closed: Stage 22.12 closure files are forbidden" >&2; exit 2; }

mkdir -p "${REPO_ROOT}/stage22_outputs" \
         "${REPO_ROOT}/stage22_cases/stage22_11_restart_statistics_visualization_audit" \
         "${REPO_ROOT}/stage22_outputs/stage22_11_restart_statistics_visualization_audit"

python3 "${SCRIPT_DIR}/assert_stage22_11_restart_statistics_visualization_audit.py"
