#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(git rev-parse --show-toplevel)
cd "$ROOT_DIR"
mkdir -p stage5_outputs
cmake -S . -B build_stage5
cmake --build build_stage5 --target fibre_stage5_closed_loop_diagnostics_audit
build_stage5/bin/fibre_stage5_closed_loop_diagnostics_audit
cat stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat
