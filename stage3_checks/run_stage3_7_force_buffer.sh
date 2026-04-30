#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

build_dir="build_stage3"
out_dir="stage3_outputs"

mkdir -p "$out_dir"

cmake -S . -B "$build_dir"
cmake --build "$build_dir" --target fibre_ibm_force_buffer_check

exe="$build_dir/bin/fibre_ibm_force_buffer_check"
if [ ! -x "$exe" ]; then
  exe=$(find "$build_dir" -type f -name fibre_ibm_force_buffer_check -perm -111 | head -n 1 || true)
fi
if [ -z "${exe:-}" ] || [ ! -x "$exe" ]; then
  echo "ERROR: fibre_ibm_force_buffer_check executable not found under $build_dir"
  exit 1
fi

"$exe"

cat "$out_dir"/fibre_ibm_force_buffer_check.dat
