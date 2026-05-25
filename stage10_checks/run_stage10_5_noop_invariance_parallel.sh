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

getv() {
  local key="$1"
  local file="$2"
  awk -v k="$key" '$1==k {print $2; found=1} END{if(!found) exit 1}' "$file" 2>/dev/null
}

check_key_one() {
  local key="$1"
  local file="$2"
  local val
  val=$(getv "$key" "$file") || { add_failure "$file: missing $key"; return 1; }
  [ "$val" = "1" ] || { add_failure "$file: $key !=1"; return 1; }
  return 0
}

check_key_zero() {
  local key="$1"
  local file="$2"
  local val
  val=$(getv "$key" "$file") || { add_failure "$file: missing $key"; return 1; }
  [ "$val" = "0" ] || { add_failure "$file: $key !=0"; return 1; }
  return 0
}

build_target xcompact3d
build_target fibre_stage10_config_check
build_target fibre_stage10_noop_hook_check

# Do not run Stage 10.2 or Stage 10.3 here by default. That avoids
# re-entering any obsolete false-positive audit policy. Stage 10.5 validates
# runtime evidence directly.
if [ "$STAGE10_5_RUN_PREREQS" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_4_INVARIANCE_ABS_TOL="$STAGE10_5_INVARIANCE_ABS_TOL" STAGE10_4_INVARIANCE_REL_TOL="$STAGE10_5_INVARIANCE_REL_TOL" bash stage10_checks/run_stage10_4_noop_invariance_np1.sh || add_failure "optional Stage10.4 prerequisite failed"
fi

run_stage99_baseline() {
  env \
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 \
    STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
    timeout "$STAGE10_5_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
}

run_stage99_hook() {
  env \
    X3D_STAGE10_PRODUCTION_HOOK=1 \
    X3D_STAGE10_FORCE_NOOP=1 \
    X3D_STAGE10_3_MAIN_NOOP_HOOK=1 \
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    STAGE9_9_FINAL_SIGNATURE_ABS_TOL=1.0e-6 \
    STAGE9_9_FINAL_SIGNATURE_REL_TOL=1.0e-12 \
    timeout "$STAGE10_5_TIMEOUT_SEC" bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
}

if ! run_stage99_baseline; then
  baseline_run_status=0
  add_failure "baseline Stage9.9 deterministic run failed"
fi

for np in 1 2 4; do
  src="stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat"
  dst="stage10_outputs/stage10_5_baseline_np${np}.dat"
  if [ -s "$src" ]; then
    cp "$src" "$dst"
  else
    baseline_run_status=0
    add_failure "missing baseline np${np} dat"
  fi
done

rm -f stage10_outputs/fibre_stage10_3_main_noop_hook.dat

if ! run_stage99_hook; then
  hook_run_status=0
  add_failure "hook-enabled Stage9.9 deterministic run failed"
fi

for np in 1 2 4; do
  src="stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat"
  dst="stage10_outputs/stage10_5_hook_np${np}.dat"
  if [ -s "$src" ]; then
    cp "$src" "$dst"
  else
    hook_run_status=0
    add_failure "missing hook np${np} dat"
  fi
done

main_hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
if [ ! -s "$main_hook_dat" ]; then
  hook_active_status=0
  noop_safety_status=0
  add_failure "missing $main_hook_dat"
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
    check_key_one "$key" "$main_hook_dat" || { hook_active_status=0; noop_safety_status=0; }
  done
  check_key_zero stage10_3_field_modified_status "$main_hook_dat" || { hook_active_status=0; noop_safety_status=0; }
fi

sig_keys="stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz"

check_np_invariance() {
  local np="$1"
  local base="stage10_outputs/stage10_5_baseline_np${np}.dat"
  local hook="stage10_outputs/stage10_5_hook_np${np}.dat"
  local status=1

  check_key_one stage9_9_parallel_consistency_local_status "$base" || status=0
  check_key_one stage9_9_decomposition_invariant_initial_state_status "$base" || status=0
  check_key_one stage9_9_parallel_consistency_local_status "$hook" || status=0
  check_key_one stage9_9_decomposition_invariant_initial_state_status "$hook" || status=0

  for key in $sig_keys; do
    local b h cmp
    b=$(getv "$key" "$base") || { add_failure "np${np}: missing baseline $key"; status=0; continue; }
    h=$(getv "$key" "$hook") || { add_failure "np${np}: missing hook $key"; status=0; continue; }
    cmp=$(awk -v b="$b" -v h="$h" -v abs="$STAGE10_5_INVARIANCE_ABS_TOL" -v rel="$STAGE10_5_INVARIANCE_REL_TOL" -v key="$key" -v np="$np" '
      BEGIN {
        db=b+0.0; dh=h+0.0
        delta=db-dh; if (delta<0) delta=-delta
        ref=db; if (ref<0) ref=-ref
        scale=ref; if (scale<1.0) scale=1.0
        eff=rel*scale; if (eff<abs) eff=abs
        ok=(delta<=eff) ? 1 : 0
        printf("np=%s metric=%s baseline=%.17g hook=%.17g delta=%.17g abs_tol=%s rel_tol=%s effective_tol=%.17g status=%s\n", np, key, db, dh, delta, abs, rel, eff, ok ? "PASS" : "FAIL")
        exit(ok ? 0 : 1)
      }')
    echo "$cmp"
    if ! echo "$cmp" | grep -q 'status=PASS'; then
      status=0
      add_failure "np${np}: invariant signature mismatch for $key"
    fi
  done

  return $([ "$status" -eq 1 ] && echo 0 || echo 1)
}

check_np_invariance 1 || np1_status=0
check_np_invariance 2 || np2_status=0
check_np_invariance 4 || np4_status=0

final_status=1
for v in "$build_status" "$baseline_run_status" "$hook_run_status" "$hook_active_status" "$noop_safety_status" "$np1_status" "$np2_status" "$np4_status"; do
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
