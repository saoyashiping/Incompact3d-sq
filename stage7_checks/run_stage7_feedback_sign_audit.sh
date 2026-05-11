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

mkdir -p stage7_outputs
cmake -S . -B build_stage7 \
  -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft}"
cmake --build build_stage7 --target fibre_stage7_feedback_sign_audit
EXE="build_stage7/bin/fibre_stage7_feedback_sign_audit"; [[ -x "$EXE" ]] || EXE=$(find build_stage7 -type f -name fibre_stage7_feedback_sign_audit | head -n 1)
"$EXE"; cat stage7_outputs/fibre_stage7_feedback_sign_audit.dat
