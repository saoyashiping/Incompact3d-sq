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

ensure_build_dir() {
  if [ -d "$BUILD_DIR" ] && { [ -f "$BUILD_DIR/Makefile" ] || [ -f "$BUILD_DIR/build.ninja" ]; }; then
    return 0
  fi

  echo "[INFO] configure ${BUILD_DIR} for Stage 10.1"
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake -S . -B "$BUILD_DIR"; then
    record_fail "cmake configure failed: ${BUILD_DIR}"
    return 1
  fi
  return 0
}

build_target() {
  local tgt="$1"
  if ! ensure_build_dir; then
    return 1
  fi
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$tgt" -j; then
    record_fail "build failed: $tgt"
    return 1
  fi
  echo "[PASS] build $tgt"
  return 0
}

check_dat_key_one() {
  local key="$1"
  local dat_file="$2"
  if [ ! -s "$dat_file" ]; then
    record_fail "missing or empty dat file: ${dat_file}"
    return 1
  fi
  if ! grep -Eq "^[[:space:]]*${key}[[:space:]]+1([[:space:]]*)$" "$dat_file"; then
    record_fail "dat check failed: ${key} != 1"
    return 1
  fi
  echo "[PASS] dat ${key}[[:space:]]\\+1"
  return 0
}

echo "[INFO] Stage 10.1 prerequisite/build setup"
if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  echo "[INFO] Stage 10.1 prerequisite: Stage 10.0"
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh; then
    record_fail "Stage 10.0 prerequisite failed"
  fi
else
  echo "[INFO] STAGE10_SKIP_PREREQS=1 -> skipping Stage 10.0 prerequisite"
fi

echo "[INFO] Stage 10.1 build checks"
build_target xcompact3d || true
build_target fibre_stage10_config_check || true
build_target fibre_stage10_noop_hook_check || true

hook_log="stage10_outputs/stage10_1_hook_build.log"
dat_file="stage10_outputs/fibre_stage10_1_noop_hook.dat"
hook_exe="${BUILD_DIR}/bin/fibre_stage10_noop_hook_check"

rm -f "$hook_log" "$dat_file"

echo "[INFO] running fibre_stage10_noop_hook_check"
if [ ! -x "$hook_exe" ]; then
  record_fail "missing executable: ${hook_exe}"
else
  if ! env X3D_STAGE10_PRODUCTION_HOOK=1 X3D_STAGE10_FORCE_NOOP=1 X3D_STAGE10_MAX_STEPS=3 \
    "$hook_exe" >"$hook_log" 2>&1; then
    record_fail "fibre_stage10_noop_hook_check returned non-zero"
  fi
fi

if [ ! -s "$hook_log" ]; then
  record_fail "missing or empty Stage 10.1 hook log"
elif ! grep -Eq "STAGE 10\.1 HOOK BUILD VERDICT:[[:space:]]*PASS" "$hook_log"; then
  record_fail "missing Stage 10.1 PASS line"
else
  echo "[PASS] hook verdict pass line"
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
  check_dat_key_one "$key" "$dat_file" || true
 done

if [ -n "$failures" ]; then
  echo "STAGE 10.1 FINAL VERDICT: FAIL"
  printf '%b\n' "$failures"
  exit 1
fi

echo "STAGE 10.1 FINAL VERDICT: PASS"
