#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_7_RUN_PREREQS=${STAGE10_7_RUN_PREREQS:-0}

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
stage9_7_status=1
stage9_8_status=1
hook_active_status=1
noop_safety_status=1
no_restart_contamination_status=1
no_stats_visu_contamination_status=1

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

if [ "$STAGE10_7_RUN_PREREQS" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_5_noop_invariance_parallel.sh || add_failure "optional Stage 10.5 prerequisite failed"
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_6_rhs_contamination_audit.sh || add_failure "optional Stage 10.6 prerequisite failed"
fi

if ! env X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh; then
  stage9_7_status=0
  add_failure "Stage 9.7 under Stage 10 no-op hook failed"
fi

# Preserve Stage 9.7 artifacts (best effort)
for f in stage9_outputs/*stage9_7* stage9_outputs/*stats* stage9_outputs/*visu* stage9_outputs/*coarse*; do
  [ -e "$f" ] || continue
  cp "$f" stage10_outputs/ || true
done

# Canonical Stage 9.7 dat/log copies if present
for f in stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat stage9_outputs/stage9_7_stats_visu_io_smoke.log stage9_outputs/run_stage9_7_stats_visu_io_smoke.log; do
  [ -f "$f" ] && cp "$f" stage10_outputs/stage10_7_stage9_7_stats_visu_io_smoke.${f##*.} || true
done

if ! env X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_3_MAIN_NOOP_HOOK=1 STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_8_restart_io_regression.sh; then
  stage9_8_status=0
  add_failure "Stage 9.8 under Stage 10 no-op hook failed"
fi

# Preserve Stage 9.8 artifacts (best effort)
for f in stage9_outputs/*stage9_8* stage9_outputs/*restart*; do
  [ -e "$f" ] || continue
  cp "$f" stage10_outputs/ || true
done

# Canonical Stage 9.8 dat/log copies if present
for f in stage9_outputs/fibre_stage9_8_restart_io_regression_np1_restart.dat stage9_outputs/stage9_8_restart_io_regression.log stage9_outputs/run_stage9_8_restart_io_regression.log; do
  [ -f "$f" ] && cp "$f" stage10_outputs/stage10_7_stage9_8_restart_io_regression.${f##*.} || true
done

hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
if [ ! -s "$hook_dat" ]; then
  hook_active_status=0
  noop_safety_status=0
  add_failure "missing required hook diagnostic file: $hook_dat"
fi

getv() {
  awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1
}

need_one() {
  local k="$1"; local f="$2"; local v
  v=$(getv "$k" "$f")
  [ "$v" = "1" ]
}

need_zero() {
  local k="$1"; local f="$2"; local v
  v=$(getv "$k" "$f")
  [ "$v" = "0" ]
}

if [ -s "$hook_dat" ]; then
  for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status stage10_3_hook_post_projection_status stage10_3_hook_post_step_status stage10_3_hook_finalize_status stage10_3_no_fibre_state_status stage10_3_no_force_status stage10_3_no_rhs_injection_status stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status stage10_3_main_noop_hook_status; do
    need_one "$k" "$hook_dat" || { noop_safety_status=0; add_failure "hook dat key not 1: $k"; }
  done
  need_zero stage10_3_field_modified_status "$hook_dat" || { noop_safety_status=0; add_failure "hook dat key not 0: stage10_3_field_modified_status"; }
fi

# If 9.7/9.8 pass under hook env, contamination statuses pass
[ "$stage9_8_status" -eq 1 ] || no_restart_contamination_status=0
[ "$stage9_7_status" -eq 1 ] || no_stats_visu_contamination_status=0

final_status=1
for v in $build_status $stage9_7_status $stage9_8_status $hook_active_status $noop_safety_status $no_restart_contamination_status $no_stats_visu_contamination_status; do
  [ "$v" -eq 1 ] || final_status=0
done

cat > stage10_outputs/stage10_7_io_restart_stats_visu_noop.dat <<DAT
stage10_7_requested_flag 1
stage10_7_build_status $build_status
stage10_7_stage9_7_stats_visu_status $stage9_7_status
stage10_7_stage9_8_restart_status $stage9_8_status
stage10_7_hook_active_status $hook_active_status
stage10_7_noop_safety_status $noop_safety_status
stage10_7_no_restart_contamination_status $no_restart_contamination_status
stage10_7_no_stats_visu_contamination_status $no_stats_visu_contamination_status
stage10_7_io_restart_stats_visu_noop_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 10.7 IO RESTART STATS VISU NOOP VERDICT: PASS"
  echo "STAGE 10.7 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.7 IO RESTART STATS VISU NOOP VERDICT: FAIL"
echo "STAGE 10.7 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
