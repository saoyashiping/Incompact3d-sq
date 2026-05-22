#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
CHANNEL_INPUT=${CHANNEL_INPUT:-examples/Channel/input.i3d}

fails=()
pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1"; fails+=("$1"); }
run_step(){ local n="$1"; shift; echo "[INFO] Running: $n"; if "$@"; then pass "$n"; else fail "$n"; fi; }

if [ ! -d "${BUILD_DIR}" ]; then
  if [ -n "${DECOMP2D_ROOT}" ]; then
    run_step "cmake configure (${BUILD_DIR}) with DECOMP2D_ROOT" cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}"
  else
    run_step "cmake configure (${BUILD_DIR})" cmake -S . -B "${BUILD_DIR}"
  fi
fi

run_step "xcompact3d build" cmake --build "${BUILD_DIR}" --target xcompact3d -j
run_step "fibre_stage9_dependency_gate_check build" cmake --build "${BUILD_DIR}" --target fibre_stage9_dependency_gate_check -j
run_step "fibre_stage9_2_minimal_parallel_gate build" cmake --build "${BUILD_DIR}" --target fibre_stage9_2_minimal_parallel_gate -j
run_step "Stage 9.1 interface consistency" bash stage9_checks/run_stage9_1_interface_consistency.sh
run_step "Stage 9.2 minimal parallel gate" bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh

mkdir -p stage9_outputs
EXE="${BUILD_DIR}/bin/xcompact3d"

run_np(){
  local np="$1" log="stage9_outputs/stage9_3_channel_init_dryrun_np${1}.log"
  echo "[INFO] Running dryrun np=${np} -> ${log}"
  if X3D_STAGE9_3_CHANNEL_INIT_DRYRUN=1 "${MPIEXEC}" ${MPIEXEC_FLAGS} -np "${np}" "${EXE}" "${CHANNEL_INPUT}" >"${log}" 2>&1; then
    if grep -q "STAGE 9.3 CHANNEL INIT DRY-RUN VERDICT: PASS" "${log}"; then pass "np=${np} dryrun verdict"; else fail "np=${np} missing PASS verdict"; fi
  else
    fail "np=${np} dryrun execution failed"
  fi
}

if [ -x "${EXE}" ]; then
  run_np 1
  run_np 2
  run_np 4
else
  fail "xcompact3d executable missing: ${EXE}"
fi

if [ -f stage9_outputs/fibre_stage9_3_channel_init_dryrun.dat ] && grep -q "stage9_3_channel_init_dryrun_status[[:space:]]\+1" stage9_outputs/fibre_stage9_3_channel_init_dryrun.dat; then
  pass "stage9.3 diagnostic dat exists with PASS status"
else
  fail "stage9.3 diagnostic dat missing or status!=1"
fi

if [ ${#fails[@]} -eq 0 ]; then
  echo "STAGE 9.3 FINAL VERDICT: PASS"
  exit 0
else
  echo "STAGE 9.3 FINAL VERDICT: FAIL"
  echo "Failed checks:"
  for f in "${fails[@]}"; do echo "  - $f"; done
  exit 1
fi
