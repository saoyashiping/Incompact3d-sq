#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p stage6_outputs
cmake -S . -B build_stage6
cmake --build build_stage6 --target fibre_stage6_projection_order_check
if [[ -x build_stage6/bin/fibre_stage6_projection_order_check ]]; then
  build_stage6/bin/fibre_stage6_projection_order_check
else
  "$(find build_stage6 -type f -name fibre_stage6_projection_order_check | head -n 1)"
fi
cat stage6_outputs/fibre_stage6_projection_order_check.dat
