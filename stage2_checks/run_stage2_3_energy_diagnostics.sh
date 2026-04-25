#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

build_dir="build_stage2"
out_dir="stage2_outputs"

mkdir -p "$out_dir"

cmake -S . -B "$build_dir"
cmake --build "$build_dir" --target fibre_energy_diagnostics_check

exe="$build_dir/bin/fibre_energy_diagnostics_check"
if [ ! -x "$exe" ]; then
  exe=$(find "$build_dir" -type f -name fibre_energy_diagnostics_check -perm -111 | head -n 1 || true)
fi
if [ -z "${exe:-}" ] || [ ! -x "$exe" ]; then
  echo "ERROR: fibre_energy_diagnostics_check executable not found under $build_dir"
  exit 1
fi

"$exe"

cat "$out_dir"/fibre_energy_diagnostics_check.dat
