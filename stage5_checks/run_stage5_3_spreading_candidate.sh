#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p stage5_outputs build_stage5
cmake -S . -B build_stage5
cmake --build build_stage5 --target fibre_stage5_spreading_candidate_check -j
exe="build_stage5/bin/fibre_stage5_spreading_candidate_check"
if [[ ! -x "$exe" ]]; then
  exe="$(find build_stage5 -type f -name 'fibre_stage5_spreading_candidate_check' | head -n 1)"
fi
"$exe"
cat stage5_outputs/fibre_stage5_spreading_candidate_check.dat
