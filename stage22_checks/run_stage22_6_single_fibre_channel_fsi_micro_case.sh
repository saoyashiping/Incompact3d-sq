#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

: "${STAGE22_6_NP:=1}"
: "${STAGE22_6_LAMBDA_CONTACT:=0.0}"
: "${STAGE22_6_N_FIBRE:=1}"

if [ "${STAGE22_6_NP}" != "1" ]; then echo "Stage 22.6 fail-closed: np must be exactly 1" >&2; exit 1; fi
if [ "${STAGE22_6_LAMBDA_CONTACT}" != "0.0" ]; then echo "Stage 22.6 fail-closed: lambda_contact must be 0.0" >&2; exit 1; fi
if [ "${STAGE22_6_N_FIBRE}" != "1" ]; then echo "Stage 22.6 fail-closed: n_fibre must be exactly 1" >&2; exit 1; fi

mkdir -p "${repo_root}/stage22_cases/stage22_6_single_fibre_channel_fsi_micro_case"
mkdir -p "${repo_root}/stage22_outputs/stage22_6_single_fibre_channel_fsi_micro_case"
mkdir -p "${repo_root}/stage22_outputs"
python3 "${repo_root}/stage22_checks/assert_stage22_6_single_fibre_channel_fsi_micro_case.py"
