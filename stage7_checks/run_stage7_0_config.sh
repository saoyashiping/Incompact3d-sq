#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"
mkdir -p stage7_outputs
cmake -S . -B build_stage7 -DFIBRE_STAGE_CHECKS_ONLY=ON
cmake --build build_stage7 --target fibre_stage7_config_check
EXE="build_stage7/bin/fibre_stage7_config_check"
if [[ ! -x "$EXE" ]]; then
  EXE=$(find build_stage7 -type f -name fibre_stage7_config_check | head -n 1)
fi
"$EXE"
cat stage7_outputs/fibre_stage7_config_check.dat
