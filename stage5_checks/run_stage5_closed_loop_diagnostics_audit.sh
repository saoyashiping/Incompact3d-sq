#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p stage5_outputs

cmake -S . -B build_stage5
cmake --build build_stage5 --target fibre_stage5_closed_loop_diagnostics_audit

EXE="build_stage5/bin/fibre_stage5_closed_loop_diagnostics_audit"
if [[ ! -x "$EXE" ]]; then
  EXE="$(find build_stage5 -type f -name fibre_stage5_closed_loop_diagnostics_audit | head -n 1)"
fi

if [[ -z "${EXE:-}" || ! -x "$EXE" ]]; then
  echo "ERROR: executable fibre_stage5_closed_loop_diagnostics_audit not found"
  exit 1
fi

"$EXE"

if [[ ! -s stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat ]]; then
  echo "ERROR: stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat was not generated"
  exit 1
fi

cat stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat
