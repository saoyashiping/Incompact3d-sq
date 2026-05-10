#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(git rev-parse --show-toplevel 2>/dev/null || cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"
mkdir -p stage5_outputs
cmake -S . -B build_stage5 -DFIBRE_STAGE_CHECKS_ONLY=ON
cmake --build build_stage5 --target fibre_stage5_closed_loop_diagnostics_audit
EXE="build_stage5/bin/fibre_stage5_closed_loop_diagnostics_audit"; [[ -x "$EXE" ]] || EXE=$(find build_stage5 -type f -name fibre_stage5_closed_loop_diagnostics_audit | head -n 1)
"$EXE"; cat stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat
