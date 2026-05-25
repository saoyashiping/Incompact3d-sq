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
STAGE9_10_TIMEOUT_SEC=${STAGE9_10_TIMEOUT_SEC:-240}
STAGE9_10_INVARIANCE_ABS_TOL=${STAGE9_10_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE9_10_INVARIANCE_REL_TOL=${STAGE9_10_INVARIANCE_REL_TOL:-1.0e-14}
STAGE9_9_MAX_STEPS=${STAGE9_9_MAX_STEPS:-3}

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
  echo "STAGE 9.10 NO-COUPLING CONTAMINATION AUDIT VERDICT: FAIL"
  echo "STAGE 9.10 FINAL VERDICT: FAIL"
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
  run "stage9.9 gate" env STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
fi

static_file_status() {
  local pattern="$1"
  shift
  local files="$*"
  if grep -E "^[[:space:]]*[^!].*${pattern}" ${files} >/dev/null 2>&1; then
    echo 0
  else
    echo 1
  fi
}

prod_files="src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90"
s_vel2fib=$(static_file_status "velocity_to_fibre|fibre_velocity_interpolation" ${prod_files})
s_feedback=$(static_file_status "feedback_force|fibre_feedback" ${prod_files})
s_twoway=$(static_file_status "two[_ ]*way|force_density" ${prod_files})
s_ibm_rhs=$(static_file_status "ibm.*spread|spread.*ibm" ${prod_files})
s_structure=$(static_file_status "fibre_structure_advance|advance_fibre" ${prod_files})
s_rhsinj=$(static_file_status "fibre_force|fsi_force|coupling_force" ${prod_files})
s_placeholders=$(static_file_status "stage9_5_record_divergence_before_projection|stage9_5_record_divergence_after_projection" src/xcompact3d.f90 src/navier.f90)

s_static=1
for s in $s_vel2fib $s_feedback $s_twoway $s_ibm_rhs $s_structure $s_rhsinj $s_placeholders; do
  if [ "$s" -ne 1 ]; then s_static=0; fi
done
if [ "$s_static" -eq 1 ]; then pass "stage9.10 static audit"; else fail "stage9.10 static audit"; fi

mkdir -p stage9_outputs
EXE="${BUILD_DIR}/bin/xcompact3d"
input_file="stage9_outputs/stage9_10_input_np1.i3d"
awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); print line }' "${CHANNEL_INPUT}" > "${input_file}"

run99() {
  local tag="$1"
  local extra_flag="$2"
  local log="stage9_outputs/stage9_10_${tag}.log"
  local dat="stage9_outputs/stage9_10_${tag}.dat"
  timeout "${STAGE9_10_TIMEOUT_SEC}" env \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS="${STAGE9_9_MAX_STEPS}" \
    ${extra_flag} \
    ${MPIEXEC} ${MPIEXEC_FLAGS} -np 1 "${EXE}" "${input_file}" >"${log}" 2>&1
  local rc=$?
  if [ ${rc} -ne 0 ]; then
    fail "stage9.10 ${tag} runtime failed/timeout"
    tail -n 160 "${log}"
    return 1
  fi
  grep -q "STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS" "${log}" || { fail "stage9.10 ${tag} missing stage9.9 PASS"; return 1; }
  if [ -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat ]; then
    cp stage9_outputs/fibre_stage9_9_parallel_consistency.dat "${dat}"
  else
    fail "stage9.10 ${tag} missing stage9.9 dat"
    return 1
  fi
  echo "${dat}"
  return 0
}

baseline_dat=$(run99 baseline "") || true
audit_dat=$(run99 audit "X3D_STAGE9_10_NO_COUPLING_AUDIT=1") || true

s_runtime_baseline=0
s_runtime_audit=0
[ -n "${baseline_dat:-}" ] && [ -f "${baseline_dat}" ] && s_runtime_baseline=1
[ -n "${audit_dat:-}" ] && [ -f "${audit_dat}" ] && s_runtime_audit=1

metric_cmp() {
  local metric="$1"
  awk -v m="${metric}" -v abs_tol="${STAGE9_10_INVARIANCE_ABS_TOL}" -v rel_tol="${STAGE9_10_INVARIANCE_REL_TOL}" '
    FNR==1{idx++}
    $1==m{val[idx]=$2+0}
    END{
      if (!(1 in val) || !(2 in val)) exit 2
      ref=val[1]; cur=val[2]
      d=cur-ref; if (d<0) d=-d
      aref=ref; if (aref<0) aref=-aref
      base=aref; if (base<1.0) base=1.0
      eff=abs_tol; rel=rel_tol*base; if (rel>eff) eff=rel
      printf("[INFO] stage9.10 invariance metric=%s delta=%.12e abs_tol=%.12e rel_tol=%.12e reference=%.12e effective_tol=%.12e\n",m,d,abs_tol,rel_tol,ref,eff)
      if (d<=eff) exit 0
      exit 1
    }' "${baseline_dat}" "${audit_dat}"
}

s_runtime_invariance=1
if [ ${s_runtime_baseline} -eq 1 ] && [ ${s_runtime_audit} -eq 1 ]; then
  for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
    if metric_cmp "${metric}"; then
      pass "stage9.10 invariance ${metric}"
    else
      fail "stage9.10 invariance ${metric}"
      s_runtime_invariance=0
    fi
  done
else
  s_runtime_invariance=0
fi

[ ${s_runtime_baseline} -eq 1 ] && pass "stage9.10 runtime baseline" || fail "stage9.10 runtime baseline"
[ ${s_runtime_audit} -eq 1 ] && pass "stage9.10 runtime audit" || fail "stage9.10 runtime audit"
[ ${s_runtime_invariance} -eq 1 ] && pass "stage9.10 runtime invariance" || fail "stage9.10 runtime invariance"

s_final=1
for s in ${s_static} ${s_runtime_baseline} ${s_runtime_audit} ${s_runtime_invariance}; do
  if [ "$s" -ne 1 ]; then s_final=0; fi
done

out="stage9_outputs/stage9_10_no_coupling_contamination_audit.dat"
{
  echo "stage9_10_requested_flag 1"
  echo "stage9_10_static_audit_status ${s_static}"
  echo "stage9_10_no_stage8_velocity_to_fibre_status ${s_vel2fib}"
  echo "stage9_10_no_stage8_feedback_status ${s_feedback}"
  echo "stage9_10_no_stage8_twoway_force_density_status ${s_twoway}"
  echo "stage9_10_no_ibm_spreading_rhs_status ${s_ibm_rhs}"
  echo "stage9_10_no_fibre_structure_advance_status ${s_structure}"
  echo "stage9_10_no_rhs_force_injection_status ${s_rhsinj}"
  echo "stage9_10_no_stage9_5_forbidden_placeholder_status ${s_placeholders}"
  echo "stage9_10_runtime_baseline_status ${s_runtime_baseline}"
  echo "stage9_10_runtime_audit_status ${s_runtime_audit}"
  echo "stage9_10_runtime_invariance_status ${s_runtime_invariance}"
  echo "stage9_10_no_coupling_contamination_status ${s_final}"
} > "${out}"

if [ ${s_final} -eq 1 ]; then
  echo "STAGE 9.10 NO-COUPLING CONTAMINATION AUDIT VERDICT: PASS"
else
  echo "STAGE 9.10 NO-COUPLING CONTAMINATION AUDIT VERDICT: FAIL"
fi

if [ ${#fails[@]} -eq 0 ] && [ ${s_final} -eq 1 ]; then
  echo "STAGE 9.10 FINAL VERDICT: PASS"
else
  echo "STAGE 9.10 FINAL VERDICT: FAIL"
  printf '  - %s\n' "${fails[@]}"
  exit 1
fi
