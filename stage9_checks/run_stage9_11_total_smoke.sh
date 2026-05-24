#!/usr/bin/env bash
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE9_11_TIMEOUT_SEC=${STAGE9_11_TIMEOUT_SEC:-240}

fails=()
pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1"; fails+=("$1"); }
run_gate(){
  local key="$1"; shift
  local name="$1"; shift
  echo "[INFO] ${name}"
  if "$@"; then
    eval "${key}=1"
    pass "${name}"
  else
    eval "${key}=0"
    fail "${name}"
  fi
}

s_build=0
s91=0; s92=0; s93=0; s94=0; s95=0; s96=0; s97=0; s98=0; s99=0; s910=0
s_nocouple=1

if [ ! -d "${BUILD_DIR}" ]; then
  if [ -n "${DECOMP2D_ROOT}" ]; then
    cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" || true
  else
    cmake -S . -B "${BUILD_DIR}" || true
  fi
fi

mkdir -p stage9_outputs

if cmake --build "${BUILD_DIR}" --target xcompact3d -j; then
  s_build=1
  pass "build xcompact3d"
else
  fail "build xcompact3d"
fi
if cmake --build "${BUILD_DIR}" --target fibre_stage9_dependency_gate_check -j; then
  pass "build fibre_stage9_dependency_gate_check"
else
  fail "build fibre_stage9_dependency_gate_check"
fi
if cmake --build "${BUILD_DIR}" --target fibre_stage9_2_minimal_parallel_gate -j; then
  pass "build fibre_stage9_2_minimal_parallel_gate"
else
  fail "build fibre_stage9_2_minimal_parallel_gate"
fi

if [ ${s_build} -eq 1 ]; then
  run_gate s91 "Stage 9.1 gate" bash stage9_checks/run_stage9_1_interface_consistency.sh
  run_gate s92 "Stage 9.2 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
  run_gate s93 "Stage 9.3 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_3_channel_init_dryrun.sh
  run_gate s94 "Stage 9.4 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_4_no_fibre_dns_smoke.sh
  run_gate s95 "Stage 9.5 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_5_projection_regression.sh
  run_gate s96 "Stage 9.6 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_6_rk3_rhs_massflux_regression.sh
  run_gate s97 "Stage 9.7 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh
  run_gate s98 "Stage 9.8 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_8_restart_io_regression.sh
  if [ ${s98} -ne 1 ]; then
    fail "Stage 9.11 blocked by Stage 9.8 gate"
  fi
  run_gate s99 "Stage 9.9 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
  run_gate s910 "Stage 9.10 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_10_no_coupling_contamination_audit.sh
fi

s_total=1
for s in ${s_build} ${s91} ${s92} ${s93} ${s94} ${s95} ${s96} ${s97} ${s98} ${s99} ${s910} ${s_nocouple}; do
  if [ "$s" -ne 1 ]; then s_total=0; fi
done

{
  echo "stage9_11_requested_flag 1"
  echo "stage9_11_build_status ${s_build}"
  echo "stage9_11_stage9_1_status ${s91}"
  echo "stage9_11_stage9_2_status ${s92}"
  echo "stage9_11_stage9_3_status ${s93}"
  echo "stage9_11_stage9_4_status ${s94}"
  echo "stage9_11_stage9_5_status ${s95}"
  echo "stage9_11_stage9_6_status ${s96}"
  echo "stage9_11_stage9_7_status ${s97}"
  echo "stage9_11_stage9_8_status ${s98}"
  echo "stage9_11_stage9_9_status ${s99}"
  echo "stage9_11_stage9_10_status ${s910}"
  echo "stage9_11_no_fibre_coupling_status ${s_nocouple}"
  echo "stage9_11_total_smoke_status ${s_total}"
} > stage9_outputs/stage9_11_total_smoke.dat

if [ ${s_total} -eq 1 ]; then
  cat > stage9_checks/STAGE9_CLOSED.md <<'EOC'
# STAGE9_CLOSED

Stage 9 closes the real 2decomp-fft interface takeover and no-fibre production DNS regression suite.

## Closed sub-stages
- Stage 9.1 real interface cleanup
- Stage 9.2 minimal parallel gate
- Stage 9.3 channel initialization dry-run
- Stage 9.4 no-fibre DNS smoke
- Stage 9.5 projection/divergence regression
- Stage 9.6 RK3/RHS/mass-flux regression
- Stage 9.7 stats/visu/coarse I/O smoke
- Stage 9.8 restart I/O regression
- Stage 9.9 parallel no-fibre consistency
- Stage 9.10 no-coupling contamination audit
- Stage 9.11 total smoke closure

## Governing no-fibre model
- div(u)=0
- du/dt + u·grad(u) = -grad(p) + nu laplacian(u) + f_drive
- F_coupling=0

Stage 9 closes no-fibre production DNS regression only. Real production coupling begins in later stages and remains intentionally disconnected here.

Next recommended stage: Stage 10 production coupling hook skeleton/default no-op.
EOC
  echo "STAGE 9.11 TOTAL SMOKE VERDICT: PASS"
  echo "STAGE 9.11 FINAL VERDICT: PASS"
else
  echo "STAGE 9.11 TOTAL SMOKE VERDICT: FAIL"
  echo "STAGE 9.11 FINAL VERDICT: FAIL"
  printf '  - %s\n' "${fails[@]}"
  exit 1
fi
