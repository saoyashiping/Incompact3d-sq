#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}
STAGE10_3_TIMEOUT_SEC=${STAGE10_3_TIMEOUT_SEC:-240}
STAGE10_3_RUN_STAGE9_GUARD=${STAGE10_3_RUN_STAGE9_GUARD:-0}
STAGE10_3_RUN_SITE_AUDIT=${STAGE10_3_RUN_SITE_AUDIT:-0}

ensure_build_dir() {
  if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    cmake -S . -B "$BUILD_DIR"
  fi
}

mkdir -p stage10_outputs stage9_outputs

failures=""
add_failure(){ if [ -z "$failures" ]; then failures="$1"; else failures="$failures\n$1"; fi; }

build_status=1
build_target(){
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$1" -j; then
    build_status=0
    add_failure "build failed: $1"
    return 1
  fi
  return 0
}

build_target xcompact3d || true
build_target fibre_stage10_config_check || true
build_target fibre_stage10_noop_hook_check || true

s10_0=1; s10_1=1; s10_2=1
if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh || { s10_0=0; add_failure "stage10.0 prereq failed"; }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_1_hook_build.sh || { s10_1=0; add_failure "stage10.1 prereq failed"; }
  if [ "$STAGE10_3_RUN_SITE_AUDIT" = "1" ]; then
    DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 STAGE10_2_ALLOW_PRODUCTION_HOOKS=auto bash stage10_checks/run_stage10_2_hook_site_audit.sh || { s10_2=0; add_failure "stage10.2 site-audit prereq failed"; }
  fi
fi

run_log="stage10_outputs/stage10_3_main_noop_hook.log"
rm -f "$run_log" stage10_outputs/fibre_stage10_3_main_noop_hook.dat

# Do not bare-run xcompact3d here. Reuse the validated Stage 9.9 deterministic no-fibre path.
if ! env \
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
  X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 \
  STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 \
  timeout "$STAGE10_3_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh >"$run_log" 2>&1; then
  add_failure "Stage 9.9 deterministic no-fibre path failed under Stage 10 no-op hook"
fi

s99="stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat"
if [ ! -s "$s99" ]; then add_failure "missing Stage9.9 np1 dat output"; fi

hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
if [ ! -s "$hook_dat" ]; then add_failure "missing Stage 10.3 production hook dat output"; fi

getv(){ awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1; }
need_one(){ [ "$(getv "$1" "$2")" = "1" ]; }
need_zero(){ [ "$(getv "$1" "$2")" = "0" ]; }

hook_active_status=1
noop_safety_status=1
if [ -s "$hook_dat" ]; then
  for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status stage10_3_hook_post_projection_status stage10_3_hook_post_step_status stage10_3_hook_finalize_status stage10_3_no_fibre_state_status stage10_3_no_force_status stage10_3_no_rhs_injection_status stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status stage10_3_main_noop_hook_status; do
    if ! need_one "$k" "$hook_dat"; then noop_safety_status=0; add_failure "hook dat: $k != 1"; fi
  done
  if ! need_zero stage10_3_field_modified_status "$hook_dat"; then noop_safety_status=0; add_failure "hook dat: stage10_3_field_modified_status != 0"; fi
else
  hook_active_status=0
  noop_safety_status=0
fi

if [ "$STAGE10_3_RUN_STAGE9_GUARD" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh || add_failure "optional stage9.9 guard failed"
fi

main_status=1
for v in $build_status $s10_0 $s10_1 $s10_2 $hook_active_status $noop_safety_status; do [ "$v" -eq 1 ] || main_status=0; done

cat > stage10_outputs/stage10_3_main_noop_hook.dat <<DAT
stage10_3_requested_flag 1
stage10_3_build_status $build_status
stage10_3_stage10_0_status $s10_0
stage10_3_stage10_1_status $s10_1
stage10_3_stage10_2_status $s10_2
stage10_3_hook_active_status $hook_active_status
stage10_3_noop_safety_status $noop_safety_status
stage10_3_main_noop_hook_status $main_status
DAT

if [ "$main_status" -eq 1 ]; then
  echo "STAGE 10.3 MAIN NOOP HOOK VERDICT: PASS"
  echo "STAGE 10.3 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.3 MAIN NOOP HOOK VERDICT: FAIL"
echo "STAGE 10.3 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
