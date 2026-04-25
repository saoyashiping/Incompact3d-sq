#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

build_dir="build_stage2"
out_dir="stage2_outputs"

mkdir -p "$out_dir"

cmake -S . -B "$build_dir"
cmake --build "$build_dir" --target fibre_structure_check
"$build_dir"/src/fibre_structure_check

cat "$out_dir"/fibre_geometry_check.dat
