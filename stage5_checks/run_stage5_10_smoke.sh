#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p stage5_outputs
cmake -S . -B build_stage5
cmake --build build_stage5 --target fibre_stage5_smoke_check
if [[ -x build_stage5/bin/fibre_stage5_smoke_check ]]; then
  build_stage5/bin/fibre_stage5_smoke_check
else
  "$(find build_stage5 -type f -name fibre_stage5_smoke_check | head -n 1)"
fi
cat stage5_outputs/fibre_stage5_smoke_check.dat
