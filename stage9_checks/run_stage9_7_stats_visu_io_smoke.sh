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
STAGE9_7_MAX_STEPS=${STAGE9_7_MAX_STEPS:-3}
STAGE9_7_REQUIRE_STATS=${STAGE9_7_REQUIRE_STATS:-1}
STAGE9_7_REQUIRE_VISU=${STAGE9_7_REQUIRE_VISU:-1}
STAGE9_7_REQUIRE_COARSE_IO=${STAGE9_7_REQUIRE_COARSE_IO:-1}
fails=(); pass(){ echo "[PASS] $1"; }; fail(){ echo "[FAIL] $1"; fails+=("$1"); }
run(){ local n="$1"; shift; echo "[INFO] $n"; if "$@"; then pass "$n"; else fail "$n"; fi; }
if [ ! -d "${BUILD_DIR}" ]; then if [ -n "${DECOMP2D_ROOT}" ]; then run "cmake configure" cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}"; else run "cmake configure" cmake -S . -B "${BUILD_DIR}"; fi; fi
run "build xcompact3d" cmake --build "${BUILD_DIR}" --target xcompact3d -j
run "build fibre_stage9_dependency_gate_check" cmake --build "${BUILD_DIR}" --target fibre_stage9_dependency_gate_check -j
run "build fibre_stage9_2_minimal_parallel_gate" cmake --build "${BUILD_DIR}" --target fibre_stage9_2_minimal_parallel_gate -j
if [ "${STAGE9_SKIP_PREREQS}" = "1" ]; then
  echo "[INFO] STAGE9_SKIP_PREREQS=1 -> skipping Stage 9.1-9.6 prerequisites"
else
  run "stage9.1 gate" bash stage9_checks/run_stage9_1_interface_consistency.sh
  run "stage9.2 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
  run "stage9.3 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_3_channel_init_dryrun.sh
  run "stage9.4 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_4_no_fibre_dns_smoke.sh
  run "stage9.5 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_5_projection_regression.sh
  run "stage9.6 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_6_rk3_rhs_massflux_regression.sh
fi
mkdir -p stage9_outputs
EXE="${BUILD_DIR}/bin/xcompact3d"
for np in 1 2 4; do
  echo "[INFO] stage9.7 run np=${np}"
  log="stage9_outputs/stage9_7_stats_visu_io_smoke_np${np}.log"
  dat="stage9_outputs/fibre_stage9_7_stats_visu_io_smoke_np${np}.dat"
  if X3D_STAGE9_7_STATS_VISU_IO_SMOKE=1 X3D_STAGE9_7_MAX_STEPS="${STAGE9_7_MAX_STEPS}" X3D_STAGE9_7_REQUIRE_STATS="${STAGE9_7_REQUIRE_STATS}" X3D_STAGE9_7_REQUIRE_VISU="${STAGE9_7_REQUIRE_VISU}" X3D_STAGE9_7_REQUIRE_COARSE_IO="${STAGE9_7_REQUIRE_COARSE_IO}" "${MPIEXEC}" ${MPIEXEC_FLAGS} -np "${np}" "${EXE}" "${CHANNEL_INPUT}" >"${log}" 2>&1; then
    grep -q "STAGE 9.7 STATS VISU IO SMOKE VERDICT: PASS" "${log}" && pass "np=${np} PASS line" || fail "np=${np} missing PASS line"
    if [ -f stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat ]; then
      cp stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat "${dat}"
      grep -q "stage9_7_stats_visu_io_smoke_status[[:space:]]\\+1" "${dat}" && pass "np=${np} dat status" || fail "np=${np} dat status!=1"
    else
      fail "np=${np} dat missing"
    fi
    find data statistics -type f -size +0c >/dev/null 2>&1 && pass "np=${np} output files non-empty" || fail "np=${np} expected output files missing/empty"
  else
    fail "np=${np} execution failed"
  fi
done
if [ ${#fails[@]} -eq 0 ]; then echo "STAGE 9.7 FINAL VERDICT: PASS"; else echo "STAGE 9.7 FINAL VERDICT: FAIL"; printf '  - %s\n' "${fails[@]}"; exit 1; fi
