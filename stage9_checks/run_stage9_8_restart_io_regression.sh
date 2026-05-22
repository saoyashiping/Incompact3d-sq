#!/usr/bin/env bash
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"
BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE9_SKIP_PREREQS=${STAGE9_SKIP_PREREQS:-0}
CHANNEL_INPUT=${CHANNEL_INPUT:-examples/Channel/input.i3d}
STAGE9_8_MAX_STEPS_BEFORE_RESTART=${STAGE9_8_MAX_STEPS_BEFORE_RESTART:-3}
STAGE9_8_MAX_STEPS_AFTER_RESTART=${STAGE9_8_MAX_STEPS_AFTER_RESTART:-3}
STAGE9_8_RESTART_SIGNATURE_TOL=${STAGE9_8_RESTART_SIGNATURE_TOL:-1.0e-8}
STAGE9_8_TIMEOUT_SEC=${STAGE9_8_TIMEOUT_SEC:-240}
fails=(); pass(){ echo "[PASS] $1"; }; fail(){ echo "[FAIL] $1"; fails+=("$1"); }
run(){ local n="$1"; shift; echo "[INFO] $n"; if "$@"; then pass "$n"; else fail "$n"; fi; }
if [ ! -d "${BUILD_DIR}" ]; then if [ -n "${DECOMP2D_ROOT}" ]; then run "cmake configure" cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}"; else run "cmake configure" cmake -S . -B "${BUILD_DIR}"; fi; fi
run "build xcompact3d" cmake --build "${BUILD_DIR}" --target xcompact3d -j
run "build fibre_stage9_dependency_gate_check" cmake --build "${BUILD_DIR}" --target fibre_stage9_dependency_gate_check -j
run "build fibre_stage9_2_minimal_parallel_gate" cmake --build "${BUILD_DIR}" --target fibre_stage9_2_minimal_parallel_gate -j
if [ "${STAGE9_SKIP_PREREQS}" != "1" ]; then
 run "stage9.1 gate" bash stage9_checks/run_stage9_1_interface_consistency.sh
 run "stage9.2 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
 run "stage9.3 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_3_channel_init_dryrun.sh
 run "stage9.4 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_4_no_fibre_dns_smoke.sh
 run "stage9.5 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_5_projection_regression.sh
 run "stage9.6 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_6_rk3_rhs_massflux_regression.sh
 run "stage9.7 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh
fi
mkdir -p stage9_outputs
EXE="${BUILD_DIR}/bin/xcompact3d"
for np in 1 2 4; do
 restart_input="stage9_outputs/stage9_8_input_restart_np${np}.i3d"
 awk '{ if ($0 ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 1"); print }' "${CHANNEL_INPUT}" > "${restart_input}"
 for phase in cold restart; do
  log="stage9_outputs/stage9_8_restart_io_regression_np${np}_${phase}.log"
  dat="stage9_outputs/fibre_stage9_8_restart_io_regression_np${np}_${phase}.dat"
  sig="stage9_outputs/stage9_8_restart_signature_np${np}.dat"
  input_file="${CHANNEL_INPUT}"
  if [ "${phase}" = "restart" ]; then
    input_file="${restart_input}"
    if ! find . -maxdepth 1 -type f -name 'restart*' -size +0c >/dev/null 2>&1; then
      fail "np=${np} restart files missing/empty before restart phase"
      continue
    fi
  fi
  timeout "${STAGE9_8_TIMEOUT_SEC}" env X3D_STAGE9_8_RESTART_IO_REGRESSION=1 X3D_STAGE9_8_PHASE="${phase}" X3D_STAGE9_8_MAX_STEPS_BEFORE_RESTART="${STAGE9_8_MAX_STEPS_BEFORE_RESTART}" X3D_STAGE9_8_MAX_STEPS_AFTER_RESTART="${STAGE9_8_MAX_STEPS_AFTER_RESTART}" X3D_STAGE9_8_RESTART_SIGNATURE_TOL="${STAGE9_8_RESTART_SIGNATURE_TOL}" X3D_STAGE9_8_SIGNATURE_FILE="${sig}" "${MPIEXEC}" ${MPIEXEC_FLAGS} -np "${np}" "${EXE}" "${input_file}" >"${log}" 2>&1
  rc=$?
  if [ ${rc} -ne 0 ]; then fail "np=${np} ${phase} failed/timeout"; tail -n 160 "${log}"; continue; fi
  grep -q "STAGE 9.8 RESTART IO REGRESSION VERDICT: PASS" "${log}" && pass "np=${np} ${phase} pass line" || fail "np=${np} ${phase} missing pass line"
  if [ -f stage9_outputs/fibre_stage9_8_restart_io_regression.dat ]; then cp stage9_outputs/fibre_stage9_8_restart_io_regression.dat "${dat}"; fi
  grep -q "stage9_8_restart_io_regression_status[[:space:]]\+1" "${dat}" && pass "np=${np} ${phase} dat pass" || fail "np=${np} ${phase} dat fail"
 done
 find . -maxdepth 1 -type f -name 'restart*' -size +0c >/dev/null 2>&1 && pass "np=${np} restart files exist/nonempty" || fail "np=${np} restart files missing/empty"
done
if [ ${#fails[@]} -eq 0 ]; then echo "STAGE 9.8 FINAL VERDICT: PASS"; else echo "STAGE 9.8 FINAL VERDICT: FAIL"; printf '  - %s\n' "${fails[@]}"; exit 1; fi
