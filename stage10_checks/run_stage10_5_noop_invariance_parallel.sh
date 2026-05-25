#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_5_TIMEOUT_SEC=${STAGE10_5_TIMEOUT_SEC:-240}
STAGE10_5_RUN_PREREQS=${STAGE10_5_RUN_PREREQS:-0}
STAGE10_5_INVARIANCE_ABS_TOL=${STAGE10_5_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE10_5_INVARIANCE_REL_TOL=${STAGE10_5_INVARIANCE_REL_TOL:-1.0e-14}

ensure_build_dir(){ if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then cmake -S . -B "$BUILD_DIR"; fi; }
mkdir -p stage10_outputs stage9_outputs

failures=""
add_failure(){ if [ -z "$failures" ]; then failures="$1"; else failures="$failures\n$1"; fi; }

build_status=1
build_target(){ ensure_build_dir; if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$1" -j; then build_status=0; add_failure "build failed: $1"; fi; }
build_target xcompact3d
build_target fibre_stage10_config_check
build_target fibre_stage10_noop_hook_check

if [ "$STAGE10_5_RUN_PREREQS" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_4_noop_invariance_np1.sh || add_failure "optional Stage 10.4 prerequisite failed"
fi

run_stage99_baseline(){
  env DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
    STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 \
    STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
    timeout "$STAGE10_5_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
}

run_stage99_hook(){
  env DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
    X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 \
    STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 \
    STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
    timeout "$STAGE10_5_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
}

baseline_run_status=1
hook_run_status=1
if ! run_stage99_baseline; then baseline_run_status=0; add_failure "baseline Stage 9.9 deterministic run failed"; fi
for np in 1 2 4; do
  if [ -s "stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat" ]; then
    cp "stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat" "stage10_outputs/stage10_5_baseline_np${np}.dat"
  else
    baseline_run_status=0; add_failure "missing baseline np=${np} dat"
  fi
done

rm -f stage10_outputs/fibre_stage10_3_main_noop_hook.dat
if ! run_stage99_hook; then hook_run_status=0; add_failure "hook-enabled Stage 9.9 deterministic run failed"; fi
for np in 1 2 4; do
  if [ -s "stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat" ]; then
    cp "stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat" "stage10_outputs/stage10_5_hook_np${np}.dat"
  else
    hook_run_status=0; add_failure "missing hook np=${np} dat"
  fi
done

hook_active_status=1
noop_safety_status=1
hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
getv(){ awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1; }
need_one(){ [ "$(getv "$1" "$2")" = "1" ]; }
need_zero(){ [ "$(getv "$1" "$2")" = "0" ]; }

if [ ! -s "$hook_dat" ]; then
  hook_active_status=0; noop_safety_status=0; add_failure "missing stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
else
  for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status stage10_3_hook_post_projection_status stage10_3_hook_post_step_status stage10_3_hook_finalize_status stage10_3_no_fibre_state_status stage10_3_no_force_status stage10_3_no_rhs_injection_status stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status stage10_3_main_noop_hook_status; do
    if ! need_one "$k" "$hook_dat"; then noop_safety_status=0; add_failure "hook dat: $k != 1"; fi
  done
  if ! need_zero stage10_3_field_modified_status "$hook_dat"; then noop_safety_status=0; add_failure "hook dat: stage10_3_field_modified_status != 0"; fi
fi

check_dat_status(){
  local file="$1" label="$2"
  for k in stage9_9_parallel_consistency_local_status stage9_9_decomposition_invariant_initial_state_status; do
    if ! need_one "$k" "$file"; then add_failure "$label: $k != 1"; return 1; fi
  done
  return 0
}
for np in 1 2 4; do
  check_dat_status "stage10_outputs/stage10_5_baseline_np${np}.dat" "baseline np=${np}" || baseline_run_status=0
  check_dat_status "stage10_outputs/stage10_5_hook_np${np}.dat" "hook np=${np}" || hook_run_status=0
done

np1_status=1; np2_status=1; np4_status=1
cmp_metric(){
  local np="$1" key="$2" b h delta eff status bfile hfile
  bfile="stage10_outputs/stage10_5_baseline_np${np}.dat"
  hfile="stage10_outputs/stage10_5_hook_np${np}.dat"
  b=$(getv "$key" "$bfile"); h=$(getv "$key" "$hfile")
  if [ -z "$b" ] || [ -z "$h" ]; then add_failure "np=${np} missing metric $key"; eval "np${np}_status=0"; return; fi
  delta=$(awk -v b="$b" -v h="$h" 'BEGIN{d=h-b; if(d<0)d=-d; printf "%.17g", d}')
  eff=$(awk -v a="$STAGE10_5_INVARIANCE_ABS_TOL" -v r="$STAGE10_5_INVARIANCE_REL_TOL" -v b="$b" 'BEGIN{ab=b; if(ab<0)ab=-ab; t=r*((ab>1)?ab:1); printf "%.17g", (a>t)?a:t}')
  status=$(awk -v d="$delta" -v e="$eff" 'BEGIN{if(d<=e)print "PASS"; else print "FAIL"}')
  echo "np=$np metric=$key baseline=$b hook=$h delta=$delta abs_tol=$STAGE10_5_INVARIANCE_ABS_TOL rel_tol=$STAGE10_5_INVARIANCE_REL_TOL effective_tol=$eff status=$status"
  if [ "$status" != "PASS" ]; then add_failure "np=${np} metric ${key} failed invariance"; eval "np${np}_status=0"; fi
}

for np in 1 2 4; do
  for m in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
    cmp_metric "$np" "$m"
  done
done

final_status=1
for v in $build_status $baseline_run_status $hook_run_status $hook_active_status $noop_safety_status $np1_status $np2_status $np4_status; do [ "$v" -eq 1 ] || final_status=0; done

cat > stage10_outputs/stage10_5_noop_invariance_parallel.dat <<DAT
stage10_5_requested_flag 1
stage10_5_build_status $build_status
stage10_5_baseline_run_status $baseline_run_status
stage10_5_hook_run_status $hook_run_status
stage10_5_hook_active_status $hook_active_status
stage10_5_noop_safety_status $noop_safety_status
stage10_5_np1_invariance_status $np1_status
stage10_5_np2_invariance_status $np2_status
stage10_5_np4_invariance_status $np4_status
stage10_5_parallel_noop_invariance_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 10.5 PARALLEL NOOP INVARIANCE VERDICT: PASS"
  echo "STAGE 10.5 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.5 PARALLEL NOOP INVARIANCE VERDICT: FAIL"
echo "STAGE 10.5 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
