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

mkdir -p stage5_outputs build_stage5
cmake -S . -B build_stage5 \
  -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft}"
cmake --build build_stage5 --target fibre_stage5_coupled_order_check -j
exe="build_stage5/bin/fibre_stage5_coupled_order_check"
if [[ ! -x "$exe" ]]; then
  exe="$(find build_stage5 -type f -name 'fibre_stage5_coupled_order_check' | head -n 1)"
fi
"$exe"
cat stage5_outputs/fibre_stage5_coupled_order_check.dat
