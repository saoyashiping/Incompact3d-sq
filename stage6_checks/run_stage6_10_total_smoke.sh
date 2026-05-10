#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"
mkdir -p stage6_outputs
rm -f stage6_outputs/STAGE6_CLOSED.md
cmake -S . -B build_stage6
cmake --build build_stage6 --target fibre_stage6_total_smoke_check
EXE="build_stage6/bin/fibre_stage6_total_smoke_check"
if [[ ! -x "$EXE" ]]; then
  EXE=$(find build_stage6 -type f -name fibre_stage6_total_smoke_check | head -n 1)
fi
"$EXE"
cat stage6_outputs/fibre_stage6_total_smoke_check.dat
