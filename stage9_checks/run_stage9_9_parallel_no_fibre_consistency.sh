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
STAGE9_9_MAX_STEPS=${STAGE9_9_MAX_STEPS:-3}
STAGE9_9_SIGNATURE_TOL=${STAGE9_9_SIGNATURE_TOL:-1.0e-8}
STAGE9_9_FINAL_SIGNATURE_ABS_TOL=${STAGE9_9_FINAL_SIGNATURE_ABS_TOL:-1.0e-6}
STAGE9_9_FINAL_SIGNATURE_REL_TOL=${STAGE9_9_FINAL_SIGNATURE_REL_TOL:-1.0e-12}
STAGE9_9_DIVERGENCE_TOL=${STAGE9_9_DIVERGENCE_TOL:-1.0e-8}
STAGE9_9_MASSFLUX_TOL=${STAGE9_9_MASSFLUX_TOL:-1.0e-6}
STAGE9_9_TIMEOUT_SEC=${STAGE9_9_TIMEOUT_SEC:-240}

fails=()
pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1"; fails+=("$1"); }
run(){ local n="$1"; shift; echo "[INFO] $n"; if "$@"; then pass "$n"; else fail "$n"; fi; }

if [ ! -d "${BUILD_DIR}" ]; then
  if [ -n "${DECOMP2D_ROOT}" ]; then
    run "cmake configure" cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}"
  else
    run "cmake configure" cmake -S . -B "${BUILD_DIR}"
  fi
fi

if ! cmake --build "${BUILD_DIR}" --target xcompact3d -j; then
  echo "STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: FAIL"
  echo "STAGE 9.9 FINAL VERDICT: FAIL"
  exit 1
fi
pass "build xcompact3d"
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
  run "stage9.8 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_8_restart_io_regression.sh
fi

mkdir -p stage9_outputs
EXE="${BUILD_DIR}/bin/xcompact3d"
for np in 1 2 4; do
  input_file="stage9_outputs/stage9_9_input_np${np}.i3d"
  awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); print line }' "${CHANNEL_INPUT}" > "${input_file}"
  log="stage9_outputs/stage9_9_parallel_no_fibre_consistency_np${np}.log"
  dat="stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat"
  timeout "${STAGE9_9_TIMEOUT_SEC}" env \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS="${STAGE9_9_MAX_STEPS}" \
    X3D_STAGE9_9_SIGNATURE_TOL="${STAGE9_9_SIGNATURE_TOL}" \
    X3D_STAGE9_9_DIVERGENCE_TOL="${STAGE9_9_DIVERGENCE_TOL}" \
    X3D_STAGE9_9_MASSFLUX_TOL="${STAGE9_9_MASSFLUX_TOL}" \
    ${MPIEXEC} ${MPIEXEC_FLAGS} -np "${np}" "${EXE}" "${input_file}" >"${log}" 2>&1
  rc=$?
  if [ ${rc} -ne 0 ]; then
    fail "np=${np} failed/timeout"
    tail -n 160 "${log}"
    continue
  fi
  grep -q "STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS" "${log}" && pass "np=${np} pass line" || fail "np=${np} missing pass line"
  if [ -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat ]; then
    cp stage9_outputs/fibre_stage9_9_parallel_consistency.dat "${dat}"
  fi
  grep -q "stage9_9_parallel_consistency_local_status[[:space:]]\+1" "${dat}" && pass "np=${np} dat pass" || fail "np=${np} dat fail"
  grep -q "stage9_9_decomposition_invariant_initial_state_status[[:space:]]\+1" "${dat}" && pass "np=${np} invariant init status" || fail "np=${np} invariant init status!=1"
done

metric_check_abs() {
  local metric="$1"
  awk -v m="${metric}" -v tol="${STAGE9_9_SIGNATURE_TOL}" '
    BEGIN{ref="";ok=1}
    FNR==1{file=FILENAME}
    $1==m{
      v=$2+0
      if (ref=="") {ref=v; next}
      d=v-ref; if (d<0) d=-d
      if (d>tol) {print "[FAIL] metric",m,"file",file,"delta",d,"tol",tol; ok=0}
    }
    END{exit(ok?0:1)}' \
    stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
    stage9_outputs/fibre_stage9_9_parallel_consistency_np2.dat \
    stage9_outputs/fibre_stage9_9_parallel_consistency_np4.dat
}

metric_check_final() {
  local metric="$1"
  awk -v m="${metric}" -v abs_tol="${STAGE9_9_FINAL_SIGNATURE_ABS_TOL}" -v rel_tol="${STAGE9_9_FINAL_SIGNATURE_REL_TOL}" '
    BEGIN{ref_set=0;ok=1}
    FNR==1{file=FILENAME}
    $1==m{
      v=$2+0
      if (ref_set==0) {ref=v; ref_set=1; next}
      d=v-ref; if (d<0) d=-d
      aref=ref; if (aref<0) aref=-aref
      base=aref; if (base<1.0) base=1.0
      eff_tol=abs_tol
      rel_term=rel_tol*base
      if (rel_term>eff_tol) eff_tol=rel_term
      printf("[INFO] final metric=%s file=%s delta=%.12e abs_tol=%.12e rel_tol=%.12e reference=%.12e effective_tol=%.12e\n",m,file,d,abs_tol,rel_tol,ref,eff_tol)
      if (d>eff_tol) {ok=0}
    }
    END{exit(ok?0:1)}' \
    stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
    stage9_outputs/fibre_stage9_9_parallel_consistency_np2.dat \
    stage9_outputs/fibre_stage9_9_parallel_consistency_np4.dat
}

initial_metrics="stage9_9_initial_signature_sum_ux stage9_9_initial_signature_sum_uy stage9_9_initial_signature_sum_uz \
stage9_9_initial_signature_max_ux stage9_9_initial_signature_max_uy stage9_9_initial_signature_max_uz \
stage9_9_initial_signature_l2_ux stage9_9_initial_signature_l2_uy stage9_9_initial_signature_l2_uz"
for metric in ${initial_metrics}; do
  if metric_check_abs "${metric}"; then
    pass "initial metric ${metric}"
  else
    fail "Stage 9.9 deterministic initial state is not decomposition-invariant. (${metric})"
  fi
done

if [ ${#fails[@]} -eq 0 ]; then
  final_metrics="stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz"
  for metric in ${final_metrics}; do
    if metric_check_final "${metric}"; then
      pass "final metric ${metric}"
    else
      fail "Stage 9.9 parallel time-advance signatures differ across MPI decompositions. (${metric})"
    fi
  done
fi

if [ ${#fails[@]} -eq 0 ]; then
  echo "STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS"
  echo "STAGE 9.9 FINAL VERDICT: PASS"
else
  echo "STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: FAIL"
  echo "STAGE 9.9 FINAL VERDICT: FAIL"
  printf '  - %s\n' "${fails[@]}"
  exit 1
fi
