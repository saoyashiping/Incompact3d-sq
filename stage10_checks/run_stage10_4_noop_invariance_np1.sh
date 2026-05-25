#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}
STAGE10_4_TIMEOUT_SEC=${STAGE10_4_TIMEOUT_SEC:-240}
STAGE10_4_INVARIANCE_ABS_TOL=${STAGE10_4_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE10_4_INVARIANCE_REL_TOL=${STAGE10_4_INVARIANCE_REL_TOL:-1.0e-14}

ensure_build_dir() {
  if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    cmake -S . -B "$BUILD_DIR"
  fi
}

mkdir -p stage10_outputs stage9_outputs

failures=""
add_failure() {
  local msg="$1"
  if [ -z "$failures" ]; then
    failures="$msg"
  else
    failures="$failures\n$msg"
  fi
}

build_status=1
stage10_0_status=1
stage10_1_status=1
stage10_2_status=1
stage10_3_status=1
baseline_run_status=1
hook_run_status=1
hook_active_status=1
noop_safety_status=1
signature_invariance_status=1

build_target() {
  local target="$1"
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$target" -j; then
    build_status=0
    add_failure "build failed: $target"
  fi
}

build_target xcompact3d
build_target fibre_stage10_config_check
build_target fibre_stage10_noop_hook_check

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh || { stage10_0_status=0; add_failure "Stage 10.0 prerequisite failed"; }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_1_hook_build.sh || { stage10_1_status=0; add_failure "Stage 10.1 prerequisite failed"; }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_2_hook_site_audit.sh || { stage10_2_status=0; add_failure "Stage 10.2 prerequisite failed"; }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_3_main_noop_hook.sh || { stage10_3_status=0; add_failure "Stage 10.3 prerequisite failed"; }
fi

# Baseline deterministic no-fibre run (np=1 compared)
if ! env STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 timeout "$STAGE10_4_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh; then
  baseline_run_status=0
  add_failure "baseline stage9.9 run failed"
fi

if [ -s stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat ]; then
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat stage10_outputs/stage10_4_baseline_np1.dat
else
  baseline_run_status=0
  add_failure "missing baseline np1 dat"
fi

# Hook-enabled no-op run
if ! env X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 timeout "$STAGE10_4_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh; then
  hook_run_status=0
  add_failure "hook-enabled stage9.9 run failed"
fi

if [ -s stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat ]; then
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat stage10_outputs/stage10_4_hook_np1.dat
else
  hook_run_status=0
  add_failure "missing hook np1 dat"
fi

if [ ! -s stage10_outputs/fibre_stage10_3_main_noop_hook.dat ]; then
  hook_active_status=0
  add_failure "missing stage10.3 hook dat proof"
fi

base_dat="stage10_outputs/stage10_4_baseline_np1.dat"
hook_dat="stage10_outputs/stage10_4_hook_np1.dat"
main_hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"

getv() {
  awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1
}

check_key_one() {
  local key="$1"
  local file="$2"
  local v
  v=$(getv "$key" "$file")
  if [ "$v" != "1" ]; then
    add_failure "${file}: ${key} != 1"
    return 1
  fi
  return 0
}

check_key_zero() {
  local key="$1"
  local file="$2"
  local v
  v=$(getv "$key" "$file")
  if [ "$v" != "0" ]; then
    add_failure "${file}: ${key} != 0"
    return 1
  fi
  return 0
}

check_key_one stage9_9_parallel_consistency_local_status "$base_dat" || baseline_run_status=0
check_key_one stage9_9_decomposition_invariant_initial_state_status "$base_dat" || baseline_run_status=0
check_key_one stage9_9_parallel_consistency_local_status "$hook_dat" || hook_run_status=0
check_key_one stage9_9_decomposition_invariant_initial_state_status "$hook_dat" || hook_run_status=0

for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status stage10_3_hook_post_projection_status stage10_3_hook_post_step_status stage10_3_hook_finalize_status stage10_3_no_fibre_state_status stage10_3_no_force_status stage10_3_no_rhs_injection_status stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status stage10_3_main_noop_hook_status; do
  check_key_one "$k" "$main_hook_dat" || noop_safety_status=0
 done
check_key_zero stage10_3_field_modified_status "$main_hook_dat" || noop_safety_status=0

compare_metric() {
  local key="$1"
  local b h delta eff status
  b=$(getv "$key" "$base_dat")
  h=$(getv "$key" "$hook_dat")
  if [ -z "$b" ] || [ -z "$h" ]; then
    add_failure "missing metric ${key}"
    signature_invariance_status=0
    return
  fi
  delta=$(awk -v b="$b" -v h="$h" 'BEGIN{d=h-b; if(d<0)d=-d; printf "%.17g", d}')
  eff=$(awk -v a="$STAGE10_4_INVARIANCE_ABS_TOL" -v r="$STAGE10_4_INVARIANCE_REL_TOL" -v b="$b" 'BEGIN{ab=b; if(ab<0)ab=-ab; t=r*((ab>1)?ab:1); printf "%.17g", (a>t)?a:t}')
  status=$(awk -v d="$delta" -v e="$eff" 'BEGIN{if(d<=e)print "PASS"; else print "FAIL"}')
  echo "metric=${key} baseline=${b} hook=${h} delta=${delta} abs_tol=${STAGE10_4_INVARIANCE_ABS_TOL} rel_tol=${STAGE10_4_INVARIANCE_REL_TOL} effective_tol=${eff} status=${status}"
  if [ "$status" != "PASS" ]; then
    signature_invariance_status=0
    add_failure "metric ${key} failed invariance"
  fi
}

for m in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
  compare_metric "$m"
done

if [ "$hook_active_status" != "1" ] || [ "$noop_safety_status" != "1" ] || [ "$signature_invariance_status" != "1" ]; then
  :
fi

final_status=1
for v in $build_status $stage10_0_status $stage10_1_status $stage10_2_status $stage10_3_status $baseline_run_status $hook_run_status $hook_active_status $noop_safety_status $signature_invariance_status; do
  if [ "$v" -ne 1 ]; then
    final_status=0
  fi
done

cat > stage10_outputs/stage10_4_noop_invariance_np1.dat <<DAT
stage10_4_requested_flag 1
stage10_4_build_status $build_status
stage10_4_stage10_0_status $stage10_0_status
stage10_4_stage10_1_status $stage10_1_status
stage10_4_stage10_2_status $stage10_2_status
stage10_4_stage10_3_status $stage10_3_status
stage10_4_baseline_run_status $baseline_run_status
stage10_4_hook_run_status $hook_run_status
stage10_4_hook_active_status $hook_active_status
stage10_4_noop_safety_status $noop_safety_status
stage10_4_signature_invariance_status $signature_invariance_status
stage10_4_noop_invariance_np1_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 10.4 NOOP INVARIANCE NP1 VERDICT: PASS"
  echo "STAGE 10.4 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.4 NOOP INVARIANCE NP1 VERDICT: FAIL"
echo "STAGE 10.4 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
