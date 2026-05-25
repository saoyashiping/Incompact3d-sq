#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}

mkdir -p stage10_outputs

failures=""

record_fail() {
  local reason="$1"
  if [ -z "$failures" ]; then
    failures="$reason"
  else
    failures="$failures\n$reason"
  fi
}

build_target() {
  local tgt="$1"
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$tgt" -j; then
    record_fail "build failed: $tgt"
    return 1
  fi
  return 0
}

echo "[INFO] Stage 10.1 build checks"
build_target xcompact3d || true
build_target fibre_stage10_config_check || true
build_target fibre_stage10_noop_hook_check || true

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  echo "[INFO] Stage 10.1 prerequisite: Stage 10.0"
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh; then
    record_fail "Stage 10.0 prerequisite failed"
  fi
fi

hook_log="stage10_outputs/stage10_1_hook_build.log"
dat_file="stage10_outputs/fibre_stage10_1_noop_hook.dat"

echo "[INFO] running fibre_stage10_noop_hook_check"
if ! env X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_MAX_STEPS=3 \
  "$BUILD_DIR/bin/fibre_stage10_noop_hook_check" >"$hook_log" 2>&1; then
  record_fail "fibre_stage10_noop_hook_check returned non-zero"
fi

if ! grep -Eq "STAGE 10\.1 HOOK BUILD VERDICT:[[:space:]]*PASS" "$hook_log"; then
  record_fail "missing Stage 10.1 PASS line"
fi

for key in \
  stage10_1_requested_flag \
  stage10_1_noop_mode_status \
  stage10_1_hook_init_status \
  stage10_1_hook_pre_step_status \
  stage10_1_hook_pre_rhs_status \
  stage10_1_hook_post_projection_status \
  stage10_1_hook_post_step_status \
  stage10_1_hook_finalize_status \
  stage10_1_no_fibre_state_status \
  stage10_1_no_force_status \
  stage10_1_no_rhs_injection_status \
  stage10_1_no_ibm_call_status \
  stage10_1_no_structure_advance_status \
  stage10_1_hook_build_status
 do
  if ! grep -Eq "^[[:space:]]*${key}[[:space:]]+1([[:space:]]*)$" "$dat_file"; then
    record_fail "dat check failed: ${key} != 1"
  fi
 done

if [ -n "$failures" ]; then
  echo "STAGE 10.1 FINAL VERDICT: FAIL"
  printf '%b\n' "$failures"
  exit 1
fi

echo "STAGE 10.1 FINAL VERDICT: PASS"
