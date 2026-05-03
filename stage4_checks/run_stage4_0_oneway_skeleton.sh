#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p stage4_outputs build_stage4

cmake -S . -B build_stage4
cmake --build build_stage4 --target fibre_stage4_oneway_skeleton_check -j

exe="build_stage4/bin/fibre_stage4_oneway_skeleton_check"
if [[ -x "$exe" ]]; then
  "$exe"
else
  exe_found="$(find build_stage4 -type f -name 'fibre_stage4_oneway_skeleton_check' | head -n 1)"
  if [[ -z "$exe_found" ]]; then
    echo "Could not locate fibre_stage4_oneway_skeleton_check under build_stage4" >&2
    exit 1
  fi
  "$exe_found"
fi

cat stage4_outputs/fibre_stage4_oneway_skeleton_check.dat
