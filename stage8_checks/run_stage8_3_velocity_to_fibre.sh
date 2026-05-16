#!/usr/bin/env bash
set -euo pipefail
: "${DECOMP2D_ROOT:=/home/sq/opt/2decomp-fft}"
if [ ! -d "$DECOMP2D_ROOT" ]; then
  echo "ERROR: DECOMP2D_ROOT=$DECOMP2D_ROOT does not exist."
  exit 1
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"
mkdir -p stage8_outputs
cmake -S . -B build_stage8 -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft}"
cmake --build build_stage8 --target fibre_stage8_velocity_to_fibre_check
EXE="build_stage8/bin/fibre_stage8_velocity_to_fibre_check"
if [[ ! -x "$EXE" ]]; then EXE="$(find build_stage8 -type f -name fibre_stage8_velocity_to_fibre_check | head -n 1)"; fi
if [[ -z "${EXE:-}" || ! -x "$EXE" ]]; then echo "ERROR: executable fibre_stage8_velocity_to_fibre_check not found"; exit 1; fi
"$EXE"
if [[ ! -s stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat ]]; then echo "ERROR: stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat was not generated"; exit 1; fi
cat stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat
