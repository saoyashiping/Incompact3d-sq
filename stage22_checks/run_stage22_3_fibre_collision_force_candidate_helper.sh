#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

mkdir -p "${repo_root}/stage22_outputs"
python3 "${repo_root}/stage22_checks/assert_stage22_3_fibre_collision_force_candidate_helper.py"
