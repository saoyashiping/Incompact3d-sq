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
STAGE9_5_MAX_STEPS=${STAGE9_5_MAX_STEPS:-3}
STAGE9_5_DIV_MAX_TOL=${STAGE9_5_DIV_MAX_TOL:-1.0e-8}
STAGE9_5_DIV_MEAN_TOL=${STAGE9_5_DIV_MEAN_TOL:-1.0e-9}

fails=(); pass(){ echo "[PASS] $1"; }; fail(){ echo "[FAIL] $1"; fails+=("$1"); }
run(){ local n="$1"; shift; echo "[INFO] $n"; if "$@"; then pass "$n"; else fail "$n"; fi; }

if [ ! -d "${BUILD_DIR}" ]; then
  if [ -n "${DECOMP2D_ROOT}" ]; then run "cmake configure" cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}"; else run "cmake configure" cmake -S . -B "${BUILD_DIR}"; fi
fi
run "build xcompact3d" cmake --build "${BUILD_DIR}" --target xcompact3d -j
run "build fibre_stage9_dependency_gate_check" cmake --build "${BUILD_DIR}" --target fibre_stage9_dependency_gate_check -j
run "build fibre_stage9_2_minimal_parallel_gate" cmake --build "${BUILD_DIR}" --target fibre_stage9_2_minimal_parallel_gate -j
run "stage9.1 gate" bash stage9_checks/run_stage9_1_interface_consistency.sh
run "stage9.2 gate" bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
run "stage9.3 gate" bash stage9_checks/run_stage9_3_channel_init_dryrun.sh
run "stage9.4 gate" bash stage9_checks/run_stage9_4_no_fibre_dns_smoke.sh

# static guard against placeholder/divergence misuse
if rg -n "stage9_5_record_divergence_before_projection|stage9_5_record_divergence_after_projection" src/xcompact3d.f90 >/dev/null; then
  fail "xcompact3d contains forbidden placeholder stage9.5 divergence calls"
fi
if rg -n "public ::.*stage9_5_record_divergence_before_projection|public ::.*stage9_5_record_divergence_after_projection" src/fibre_stage9_5_projection_regression.f90 >/dev/null; then
  fail "stage9.5 module publicly exports placeholder divergence routines"
fi
for sub in pipe_bulk pipe_bulk_u pipe_bulk_phi pipe_volume_avg; do
  if awk "/subroutine ${sub}\(/,/end subroutine ${sub}/" src/navier.f90 | rg -n "fibre_stage9_5_projection_regression" >/dev/null; then
    fail "${sub} block imports stage9.5 diagnostics"
  fi
done

mkdir -p stage9_outputs
EXE="${BUILD_DIR}/bin/xcompact3d"
for np in 1 2 4; do
  log="stage9_outputs/stage9_5_projection_regression_np${np}.log"
  if X3D_STAGE9_5_PROJECTION_REGRESSION=1 X3D_STAGE9_5_MAX_STEPS="${STAGE9_5_MAX_STEPS}" X3D_STAGE9_5_DIV_MAX_TOL="${STAGE9_5_DIV_MAX_TOL}" X3D_STAGE9_5_DIV_MEAN_TOL="${STAGE9_5_DIV_MEAN_TOL}" "${MPIEXEC}" ${MPIEXEC_FLAGS} -np "${np}" "${EXE}" "${CHANNEL_INPUT}" >"${log}" 2>&1; then
    if grep -q "STAGE 9.5 PRESSURE PROJECTION REGRESSION VERDICT: PASS" "${log}"; then pass "np=${np} verdict"; else fail "np=${np} missing PASS verdict"; fi
  else
    fail "np=${np} run failed"
  fi
  dat="stage9_outputs/fibre_stage9_5_projection_regression_np${np}.dat"
  if [ -f stage9_outputs/fibre_stage9_5_projection_regression.dat ]; then cp stage9_outputs/fibre_stage9_5_projection_regression.dat "${dat}"; fi
  if [ -f "${dat}" ] && grep -q "stage9_5_projection_regression_status[[:space:]]\+1" "${dat}"; then pass "np=${np} dat PASS"; else fail "np=${np} dat missing/fail"; fi
done

if [ ${#fails[@]} -eq 0 ]; then
  echo "STAGE 9.5 FINAL VERDICT: PASS"; exit 0
else
  echo "STAGE 9.5 FINAL VERDICT: FAIL"; printf '  - %s\n' "${fails[@]}"; exit 1
fi
