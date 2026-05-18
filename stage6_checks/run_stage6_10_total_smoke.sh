#!/usr/bin/env bash
set -euo pipefail

: "${DECOMP2D_ROOT:=/home/sq/opt/2decomp-fft}"
if [ ! -d "$DECOMP2D_ROOT" ]; then
  echo "ERROR: DECOMP2D_ROOT=$DECOMP2D_ROOT does not exist."
  echo "Install 2decomp-fft once, or export DECOMP2D_ROOT=/path/to/2decomp-fft/install."
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p stage6_outputs
rm -f stage6_outputs/STAGE6_CLOSED.md
cmake -S . -B build_stage6 \
  -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft}"
cmake --build build_stage6 --target fibre_stage6_total_smoke_check
EXE="build_stage6/bin/fibre_stage6_total_smoke_check"
if [[ ! -x "$EXE" ]]; then
  EXE=$(find build_stage6 -type f -name fibre_stage6_total_smoke_check | head -n 1)
fi
"$EXE"
cat stage6_outputs/fibre_stage6_total_smoke_check.dat
