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
add_failure() { if [ -z "$failures" ]; then failures="$1"; else failures="$failures\n$1"; fi; }

build_status=1
stage10_0_status=1
stage10_1_status=1
stage10_2_status=1
run_status=1
hook_status=1

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

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" \
       bash stage10_checks/run_stage10_0_config.sh; then
    stage10_0_status=0
    add_failure "stage10.0 prereq failed"
  fi
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 \
       bash stage10_checks/run_stage10_1_hook_build.sh; then
    stage10_1_status=0
    add_failure "stage10.1 prereq failed"
  fi
fi

# Do not run the Stage 10.2 site audit by default after hooks have been connected.
# If explicitly requested, the corrected Stage 10.2 audit permits guarded xcompact3d hooks.
if [ "$STAGE10_3_RUN_SITE_AUDIT" = "1" ]; then
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 \
       STAGE10_2_ALLOW_PRODUCTION_HOOKS=1 bash stage10_checks/run_stage10_2_hook_site_audit.sh; then
    stage10_2_status=0
    add_failure "stage10.2 optional site-audit prereq failed"
  fi
fi

run_log="stage10_outputs/stage10_3_main_noop_hook.log"
rm -f "$run_log" stage10_outputs/fibre_stage10_3_main_noop_hook.dat

if ! env DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
        X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 \
        STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
        X3D_STAGE9_9_MAX_STEPS=3 STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 \
        STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
        timeout "$STAGE10_3_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh \
        > "$run_log" 2>&1; then
  run_status=0
  add_failure "Stage 9.9 deterministic no-fibre path failed under Stage 10 no-op hook"
fi

s99="stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat"
if [ ! -s "$s99" ]; then
  run_status=0
  add_failure "missing Stage9.9 np1 dat output"
fi

hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
if [ ! -s "$hook_dat" ]; then
  hook_status=0
  add_failure "missing Stage 10.3 production hook diagnostic dat"
fi

getv() { awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1; }
need_one() { [ "$(getv "$1" "$2")" = "1" ]; }
need_zero() { [ "$(getv "$1" "$2")" = "0" ]; }

if [ -s "$hook_dat" ]; then
  for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status \
           stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status \
           stage10_3_hook_post_projection_status stage10_3_hook_post_step_status \
           stage10_3_hook_finalize_status stage10_3_no_fibre_state_status \
           stage10_3_no_force_status stage10_3_no_rhs_injection_status \
           stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status \
           stage10_3_main_noop_hook_status; do
    need_one "$k" "$hook_dat" || { hook_status=0; add_failure "hook dat key not 1: $k"; }
  done
  need_zero stage10_3_field_modified_status "$hook_dat" || { hook_status=0; add_failure "hook dat key not 0: stage10_3_field_modified_status"; }
fi

if [ "$STAGE10_3_RUN_STAGE9_GUARD" = "1" ]; then
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE9_SKIP_PREREQS=1 \
       bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh; then
    add_failure "optional stage9.9 guard failed"
  fi
fi

final_status=1
for v in $build_status $stage10_0_status $stage10_1_status $stage10_2_status $run_status $hook_status; do
  [ "$v" -eq 1 ] || final_status=0
done

# Keep the canonical Stage 10.3 diagnostic dat produced by the production hook.
if [ "$final_status" -eq 1 ]; then
  echo "STAGE 10.3 MAIN NOOP HOOK VERDICT: PASS"
  echo "STAGE 10.3 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.3 MAIN NOOP HOOK VERDICT: FAIL"
echo "STAGE 10.3 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
