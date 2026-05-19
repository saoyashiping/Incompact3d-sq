#!/usr/bin/env bash
set -euo pipefail
: "${DECOMP2D_ROOT:=/home/sq/opt/2decomp-fft}"
if [ ! -d "$DECOMP2D_ROOT" ]; then echo "ERROR: DECOMP2D_ROOT=$DECOMP2D_ROOT does not exist."; exit 1; fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"; cd "$ROOT_DIR"
mkdir -p stage7_outputs
cmake -S . -B build_stage7 -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft}"
cmake --build build_stage7 --target fibre_stage7_boundary_safety_check
EXE="build_stage7/bin/fibre_stage7_boundary_safety_check"
if [[ ! -x "$EXE" ]]; then EXE="$(find build_stage7 -type f -name fibre_stage7_boundary_safety_check | head -n 1)"; fi
if [[ -z "${EXE:-}" || ! -x "$EXE" ]]; then echo "ERROR: executable fibre_stage7_boundary_safety_check not found"; exit 1; fi
"$EXE"
if [[ ! -s stage7_outputs/fibre_stage7_boundary_safety_check.dat ]]; then echo "ERROR: stage7_outputs/fibre_stage7_boundary_safety_check.dat was not generated"; exit 1; fi
cat stage7_outputs/fibre_stage7_boundary_safety_check.dat
