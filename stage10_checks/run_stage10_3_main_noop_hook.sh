#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}
STAGE10_3_TIMEOUT_SEC=${STAGE10_3_TIMEOUT_SEC:-240}
STAGE10_3_RUN_STAGE9_GUARD=${STAGE10_3_RUN_STAGE9_GUARD:-0}

ensure_build_dir() {
  if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    if [ -n "$DECOMP2D_ROOT" ]; then
      cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT"
    else
      cmake -S . -B "$BUILD_DIR"
    fi
  fi
}

mkdir -p stage10_outputs stage9_outputs

failures=""
add_failure() {
  local msg="$1"
  if [ -z "$failures" ]; then
    failures="$msg"
  else
    failures="$failures
$msg"
  fi
}

build_target() {
  local tgt="$1"
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$tgt" -j; then
    add_failure "build failed: $tgt"
    return 1
  fi
  return 0
}

getv() {
  local key="$1"
  local file="$2"
  awk -v k="$key" '$1==k {print $2; found=1} END{if(!found) exit 1}' "$file" 2>/dev/null
}

check_one() {
  local key="$1"
  local file="$2"
  local val
  val=$(getv "$key" "$file") || { add_failure "hook dat: missing $key"; return 1; }
  [ "$val" = "1" ] || { add_failure "hook dat: $key !=1"; return 1; }
  return 0
}

check_zero() {
  local key="$1"
  local file="$2"
  local val
  val=$(getv "$key" "$file") || { add_failure "hook dat: missing $key"; return 1; }
  [ "$val" = "0" ] || { add_failure "hook dat: $key !=0"; return 1; }
  return 0
}

build_target xcompact3d || true
build_target fibre_stage10_config_check || true
build_target fibre_stage10_noop_hook_check || true

s10_0=1
s10_1=1
s10_2=1
if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_0_config.sh || {
    s10_0=0
    add_failure "stage10.0 prereq failed"
  }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_1_hook_build.sh || {
    s10_1=0
    add_failure "stage10.1 prereq failed"
  }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 STAGE10_2_ALLOW_PRODUCTION_HOOKS=auto bash stage10_checks/run_stage10_2_hook_site_audit.sh || {
    s10_2=0
    add_failure "stage10.2 prereq failed"
  }
fi

run_log="stage10_outputs/stage10_3_main_noop_hook.log"
hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
rm -f "$hook_dat"

# Do not bare-run xcompact3d. Reuse the already validated Stage 9.9
# deterministic no-fibre run path, with Stage 10 production hook enabled.
if ! env \
  X3D_STAGE10_PRODUCTION_HOOK=1 \
  X3D_STAGE10_FORCE_NOOP=1 \
  X3D_STAGE10_3_MAIN_NOOP_HOOK=1 \
  STAGE9_SKIP_PREREQS=1 \
  X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
  X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
  X3D_STAGE9_9_MAX_STEPS=3 \
  STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 \
  STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
  timeout "$STAGE10_3_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh >"$run_log" 2>&1; then
  add_failure "stage9.9 deterministic hook smoke failed"
fi

if [ ! -s "stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat" ]; then
  add_failure "missing Stage9.9 np1 dat output"
fi

if [ ! -s "$hook_dat" ]; then
  add_failure "missing $hook_dat"
else
  for key in \
    stage10_3_requested_flag \
    stage10_3_noop_mode_status \
    stage10_3_hook_init_status \
    stage10_3_hook_pre_step_status \
    stage10_3_hook_pre_rhs_status \
    stage10_3_hook_post_projection_status \
    stage10_3_hook_post_step_status \
    stage10_3_hook_finalize_status \
    stage10_3_no_fibre_state_status \
    stage10_3_no_force_status \
    stage10_3_no_rhs_injection_status \
    stage10_3_no_ibm_call_status \
    stage10_3_no_structure_advance_status \
    stage10_3_main_noop_hook_status; do
    check_one "$key" "$hook_dat" || true
  done
  check_zero stage10_3_field_modified_status "$hook_dat" || true
fi

if [ "$STAGE10_3_RUN_STAGE9_GUARD" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh || add_failure "optional stage9.9 guard failed"
fi

main_status=1
[ -z "$failures" ] || main_status=0

cat > stage10_outputs/stage10_3_main_noop_hook_status.dat <<DAT
stage10_3_requested_flag 1
stage10_3_stage10_0_status $s10_0
stage10_3_stage10_1_status $s10_1
stage10_3_stage10_2_status $s10_2
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
