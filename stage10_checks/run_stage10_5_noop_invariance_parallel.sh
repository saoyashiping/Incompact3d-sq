#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_5_RUN_PREREQS=${STAGE10_5_RUN_PREREQS:-0}
STAGE10_5_INVARIANCE_ABS_TOL=${STAGE10_5_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE10_5_INVARIANCE_REL_TOL=${STAGE10_5_INVARIANCE_REL_TOL:-1.0e-14}
STAGE10_5_TIMEOUT_SEC=${STAGE10_5_TIMEOUT_SEC:-240}

ensure_build_dir() {
  if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    cmake -S . -B "$BUILD_DIR"
  fi
}

mkdir -p stage10_outputs stage9_outputs

failures=""
add_failure() {
  local msg="$1"
  if [ -z "$failures" ]; then failures="$msg"; else failures="$failures\n$msg"; fi
}

build_status=1
baseline_run_status=1
hook_run_status=1
hook_active_status=1
noop_safety_status=1
np1_status=1
np2_status=1
np4_status=1

build_target() {
  local tgt="$1"
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$tgt" -j; then
    build_status=0
    add_failure "build failed: $tgt"
  fi
}

build_target xcompact3d
build_target fibre_stage10_config_check
build_target fibre_stage10_noop_hook_check

# Do not run Stage 10.2 or Stage 10.3 from Stage 10.5 by default.
# Optional prereq mode may run only the already higher-level 10.4 gate.
if [ "$STAGE10_5_RUN_PREREQS" = "1" ]; then
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
       bash stage10_checks/run_stage10_4_noop_invariance_np1.sh; then
    add_failure "optional Stage 10.4 prerequisite failed"
  fi
fi

run_stage99_baseline() {
  env DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
      STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
      X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 \
      STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
      timeout "$STAGE10_5_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
}

run_stage99_hook() {
  env DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
      X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 \
      STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
      X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 \
      STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
      timeout "$STAGE10_5_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
}

rm -f stage10_outputs/fibre_stage10_3_main_noop_hook.dat

if ! run_stage99_baseline; then
  baseline_run_status=0
  add_failure "baseline stage9.9 run failed"
fi

for np in 1 2 4; do
  src="stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat"
  dst="stage10_outputs/stage10_5_baseline_np${np}.dat"
  if [ -s "$src" ]; then
    cp "$src" "$dst"
  else
    baseline_run_status=0
    add_failure "missing baseline dat for np=${np}"
  fi
done

if ! run_stage99_hook; then
  hook_run_status=0
  add_failure "hook-enabled stage9.9 run failed"
fi

for np in 1 2 4; do
  src="stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat"
  dst="stage10_outputs/stage10_5_hook_np${np}.dat"
  if [ -s "$src" ]; then
    cp "$src" "$dst"
  else
    hook_run_status=0
    add_failure "missing hook dat for np=${np}"
  fi
done

if [ ! -s stage10_outputs/fibre_stage10_3_main_noop_hook.dat ]; then
  hook_active_status=0
  add_failure "missing stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
fi

getv() { awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1; }
need_one() { [ "$(getv "$1" "$2")" = "1" ]; }
need_zero() { [ "$(getv "$1" "$2")" = "0" ]; }

for np in 1 2 4; do
  b="stage10_outputs/stage10_5_baseline_np${np}.dat"
  h="stage10_outputs/stage10_5_hook_np${np}.dat"
  need_one stage9_9_parallel_consistency_local_status "$b" || { add_failure "baseline np=${np}: local status !=1"; eval "np${np}_status=0"; }
  need_one stage9_9_decomposition_invariant_initial_state_status "$b" || { add_failure "baseline np=${np}: invariant init status !=1"; eval "np${np}_status=0"; }
  need_one stage9_9_parallel_consistency_local_status "$h" || { add_failure "hook np=${np}: local status !=1"; eval "np${np}_status=0"; }
  need_one stage9_9_decomposition_invariant_initial_state_status "$h" || { add_failure "hook np=${np}: invariant init status !=1"; eval "np${np}_status=0"; }
done

hk="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status \
         stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status \
         stage10_3_hook_post_projection_status stage10_3_hook_post_step_status \
         stage10_3_hook_finalize_status stage10_3_no_fibre_state_status \
         stage10_3_no_force_status stage10_3_no_rhs_injection_status \
         stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status \
         stage10_3_main_noop_hook_status; do
  need_one "$k" "$hk" || { noop_safety_status=0; add_failure "hook dat: ${k} !=1"; }
done
need_zero stage10_3_field_modified_status "$hk" || { noop_safety_status=0; add_failure "hook dat: stage10_3_field_modified_status !=0"; }

cmp_metric() {
  local np="$1"
  local metric="$2"
  local bfile="stage10_outputs/stage10_5_baseline_np${np}.dat"
  local hfile="stage10_outputs/stage10_5_hook_np${np}.dat"
  local b h d eff pass
  b=$(getv "$metric" "$bfile")
  h=$(getv "$metric" "$hfile")
  if [ -z "$b" ] || [ -z "$h" ]; then
    add_failure "np=${np} missing metric ${metric}"
    eval "np${np}_status=0"
    return
  fi
  d=$(awk -v b="$b" -v h="$h" 'BEGIN{x=h-b;if(x<0)x=-x;printf "%.17g",x}')
  eff=$(awk -v a="$STAGE10_5_INVARIANCE_ABS_TOL" -v r="$STAGE10_5_INVARIANCE_REL_TOL" -v b="$b" 'BEGIN{ab=b;if(ab<0)ab=-ab;t=r*((ab>1)?ab:1);printf "%.17g",(a>t)?a:t}')
  pass=$(awk -v d="$d" -v e="$eff" 'BEGIN{if(d<=e)print "PASS"; else print "FAIL"}')
  echo "np=${np} metric=${metric} baseline=${b} hook=${h} delta=${d} abs_tol=${STAGE10_5_INVARIANCE_ABS_TOL} rel_tol=${STAGE10_5_INVARIANCE_REL_TOL} effective_tol=${eff} status=${pass}"
  if [ "$pass" != "PASS" ]; then
    add_failure "np=${np} metric ${metric} failed"
    eval "np${np}_status=0"
  fi
}

for np in 1 2 4; do
  for m in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
           stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
           stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
    cmp_metric "$np" "$m"
  done
done

final_status=1
for v in $build_status $baseline_run_status $hook_run_status $hook_active_status $noop_safety_status $np1_status $np2_status $np4_status; do
  [ "$v" -eq 1 ] || final_status=0
done

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
