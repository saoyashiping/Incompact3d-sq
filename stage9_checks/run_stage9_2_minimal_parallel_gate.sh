#!/usr/bin/env bash
set -u

MPIEXEC=${MPIEXEC:-mpirun}
fails=()

pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1"; fails+=("$1"); }

run_step(){
  local name="$1"; shift
  echo "[INFO] Running: $name"
  if "$@"; then pass "$name"; else fail "$name"; fi
}

echo "[INFO] Stage 9.2 minimal parallel gate started"

if [ ! -d build_stage9 ]; then
  fail "build_stage9 directory missing"
else
  pass "build_stage9 directory exists"
fi

run_step "xcompact3d build" cmake --build build_stage9 --target xcompact3d -j
run_step "fibre_stage9_dependency_gate_check build" cmake --build build_stage9 --target fibre_stage9_dependency_gate_check -j
run_step "fibre_stage9_2_minimal_parallel_gate build" cmake --build build_stage9 --target fibre_stage9_2_minimal_parallel_gate -j
run_step "Stage 9.1 interface consistency" bash stage9_checks/run_stage9_1_interface_consistency.sh
run_step "np=1 minimal parallel gate" ${MPIEXEC} -np 1 ./build_stage9/bin/fibre_stage9_2_minimal_parallel_gate
run_step "np=2 minimal parallel gate" ${MPIEXEC} -np 2 ./build_stage9/bin/fibre_stage9_2_minimal_parallel_gate
run_step "np=4 minimal parallel gate" ${MPIEXEC} -np 4 ./build_stage9/bin/fibre_stage9_2_minimal_parallel_gate

if [ ${#fails[@]} -eq 0 ]; then
  echo "============================================================"
  echo "STAGE 9.2 FINAL VERDICT: PASS"
  echo "============================================================"
  exit 0
else
  echo "============================================================"
  echo "STAGE 9.2 FINAL VERDICT: FAIL"
  echo "Failed checks:"
  for f in "${fails[@]}"; do
    echo "  - $f"
  done
  echo "============================================================"
  exit 1
fi
