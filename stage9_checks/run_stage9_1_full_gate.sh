#!/usr/bin/env bash
set -u

echo "[INFO] Stage 9.1 full gate started"

fails=()

run_step() {
  local name="$1"; shift
  echo "[INFO] Running: $name"
  if "$@"; then
    echo "[PASS] $name"
  else
    echo "[FAIL] $name"
    fails+=("$name")
  fi
}

run_step "cmake configure" cmake -S . -B build_stage9
run_step "build xcompact3d" cmake --build build_stage9 --target xcompact3d -j
run_step "build fibre_stage9_dependency_gate_check" cmake --build build_stage9 --target fibre_stage9_dependency_gate_check -j
run_step "stage9.1 interface consistency" bash stage9_checks/run_stage9_1_interface_consistency.sh

if [ ${#fails[@]} -eq 0 ]; then
  echo "============================================================"
  echo "STAGE 9.1 FULL GATE VERDICT: PASS"
  echo "============================================================"
  exit 0
else
  echo "============================================================"
  echo "STAGE 9.1 FULL GATE VERDICT: FAIL"
  echo "Failed steps:"
  for item in "${fails[@]}"; do
    echo "  - $item"
  done
  echo "============================================================"
  exit 1
fi
