#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(git rev-parse --show-toplevel 2>/dev/null || cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"
mkdir -p stage7_outputs
cmake -S . -B build_stage7
cmake --build build_stage7 --target fibre_stage7_rho_convention_audit
EXE="build_stage7/bin/fibre_stage7_rho_convention_audit"; [[ -x "$EXE" ]] || EXE=$(find build_stage7 -type f -name fibre_stage7_rho_convention_audit | head -n 1)
"$EXE"; cat stage7_outputs/fibre_stage7_rho_convention_audit.dat
