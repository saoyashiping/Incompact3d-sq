#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${STAGE22_12_ENABLE:=1}"
: "${STAGE22_12_FINAL_TOTAL_CLOSURE_ENABLE:=1}"
: "${STAGE22_12_BUILD_ALLOWED:=0}"
: "${STAGE22_12_DNS_ALLOWED:=0}"
: "${STAGE22_12_MPI_ALLOWED:=0}"
: "${STAGE22_12_RESTART_RERUN_ALLOWED:=0}"
: "${STAGE22_12_STATISTICS_RERUN_ALLOWED:=0}"
: "${STAGE22_12_VISUALIZATION_RERUN_ALLOWED:=0}"
: "${STAGE22_12_MESH_REFINEMENT_ALLOWED:=0}"
: "${STAGE22_12_NEW_PHYSICS_ALLOWED:=0}"
: "${STAGE22_12_SOURCE_MODIFICATION_ALLOWED:=0}"
: "${STAGE22_12_SCHEMA_MODIFICATION_ALLOWED:=0}"
: "${STAGE22_12_CLOSED_STAGE_MODIFICATION_ALLOWED:=0}"
: "${STAGE22_12_CREATE_STAGE22_CLOSED:=1}"
: "${STAGE22_12_CREATE_PROJECT_FINAL_CLOSED:=1}"
: "${STAGE22_12_DECLARE_NO_STAGE23:=1}"

[[ "${STAGE22_12_ENABLE}" == "1" && "${STAGE22_12_FINAL_TOTAL_CLOSURE_ENABLE}" == "1" ]] || { echo "Stage 22.12 fail-closed: closure must be enabled" >&2; exit 2; }
[[ "${STAGE22_12_BUILD_ALLOWED}" == "0" ]] || { echo "Stage 22.12 fail-closed: build is forbidden" >&2; exit 2; }
[[ "${STAGE22_12_DNS_ALLOWED}" == "0" ]] || { echo "Stage 22.12 fail-closed: DNS is forbidden" >&2; exit 2; }
[[ "${STAGE22_12_MPI_ALLOWED}" == "0" ]] || { echo "Stage 22.12 fail-closed: MPI is forbidden" >&2; exit 2; }
[[ "${STAGE22_12_RESTART_RERUN_ALLOWED}" == "0" && "${STAGE22_12_STATISTICS_RERUN_ALLOWED}" == "0" && "${STAGE22_12_VISUALIZATION_RERUN_ALLOWED}" == "0" ]] || { echo "Stage 22.12 fail-closed: reruns are forbidden" >&2; exit 2; }
[[ "${STAGE22_12_MESH_REFINEMENT_ALLOWED}" == "0" && "${STAGE22_12_NEW_PHYSICS_ALLOWED}" == "0" ]] || { echo "Stage 22.12 fail-closed: refinement/new physics forbidden" >&2; exit 2; }
[[ "${STAGE22_12_SOURCE_MODIFICATION_ALLOWED}" == "0" && "${STAGE22_12_SCHEMA_MODIFICATION_ALLOWED}" == "0" && "${STAGE22_12_CLOSED_STAGE_MODIFICATION_ALLOWED}" == "0" ]] || { echo "Stage 22.12 fail-closed: source/schema/closed-stage modification forbidden" >&2; exit 2; }
[[ "${STAGE22_12_CREATE_STAGE22_CLOSED}" == "1" && "${STAGE22_12_CREATE_PROJECT_FINAL_CLOSED}" == "1" && "${STAGE22_12_DECLARE_NO_STAGE23}" == "1" ]] || { echo "Stage 22.12 fail-closed: closure docs and no Stage 23 declaration are required" >&2; exit 2; }

mkdir -p "${REPO_ROOT}/stage22_outputs" \
         "${REPO_ROOT}/stage22_outputs/stage22_12_final_total_closure"

python3 "${SCRIPT_DIR}/assert_stage22_12_final_total_closure.py"
