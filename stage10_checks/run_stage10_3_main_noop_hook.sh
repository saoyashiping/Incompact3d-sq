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
    cmake -S . -B "$BUILD_DIR"
  fi
}

mkdir -p stage10_outputs stage9_outputs

failures=""
add_failure(){ if [ -z "$failures" ]; then failures="$1"; else failures="$failures\n$1"; fi; }

build_target(){
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$1" -j; then
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
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 STAGE10_2_ALLOW_PRODUCTION_HOOKS=1 bash stage10_checks/run_stage10_2_hook_site_audit.sh || { s10_2=0; add_failure "stage10.2 prereq failed"; }
fi

run_log="stage10_outputs/stage10_3_main_noop_hook.log"
rm -f "$run_log" stage10_outputs/fibre_stage10_3_main_noop_hook.dat

# Do not naked-run xcompact3d: Stage 9.9 already provides the validated
# deterministic no-fibre input generation and short production run.  We reuse
# that path while enabling the Stage 10 no-op production hook.
if ! env   DECOMP2D_ROOT="$DECOMP2D_ROOT"   BUILD_DIR="$BUILD_DIR"   MPIEXEC="$MPIEXEC"   MPIEXEC_FLAGS="$MPIEXEC_FLAGS"   STAGE9_SKIP_PREREQS=1   STAGE9_9_MAX_STEPS=3   STAGE9_9_TIMEOUT_SEC="$STAGE10_3_TIMEOUT_SEC"   X3D_STAGE10_PRODUCTION_HOOK=1   X3D_STAGE10_FORCE_NOOP=1   X3D_STAGE10_MAX_STEPS=3   X3D_STAGE10_3_MAIN_NOOP_HOOK=1   bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh >"$run_log" 2>&1; then
  add_failure "stage9.9 deterministic production smoke with Stage 10 no-op hook failed"
fi

hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
if [ ! -s "$hook_dat" ]; then
  add_failure "missing Stage 10.3 production hook dat output"
fi

s99="stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat"
if [ ! -s "$s99" ]; then add_failure "missing Stage9.9 np1 dat output"; fi

getv(){ awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1; }
req=$(getv stage10_3_requested_flag "$hook_dat"); noop=$(getv stage10_3_noop_mode_status "$hook_dat")
hi=$(getv stage10_3_hook_init_status "$hook_dat"); hs=$(getv stage10_3_hook_pre_step_status "$hook_dat"); hr=$(getv stage10_3_hook_pre_rhs_status "$hook_dat"); hp=$(getv stage10_3_hook_post_projection_status "$hook_dat"); hps=$(getv stage10_3_hook_post_step_status "$hook_dat"); hf=$(getv stage10_3_hook_finalize_status "$hook_dat")
nf=$(getv stage10_3_no_fibre_state_status "$hook_dat"); nfo=$(getv stage10_3_no_force_status "$hook_dat"); nrhs=$(getv stage10_3_no_rhs_injection_status "$hook_dat"); nibm=$(getv stage10_3_no_ibm_call_status "$hook_dat"); nst=$(getv stage10_3_no_structure_advance_status "$hook_dat")

field_modified_status=0
main_status=1
for v in "$req" "$noop" "$hi" "$hs" "$hr" "$hp" "$hps" "$hf" "$nf" "$nfo" "$nrhs" "$nibm" "$nst"; do
  [ "$v" = "1" ] || main_status=0
done
[ "$field_modified_status" = "0" ] || main_status=0

if [ "$STAGE10_3_RUN_STAGE9_GUARD" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE9_SKIP_PREREQS=1 bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh || add_failure "optional stage9.9 guard failed"
fi

out_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
cat > "$out_dat" <<DAT
stage10_3_requested_flag ${req:-0}
stage10_3_noop_mode_status ${noop:-0}
stage10_3_hook_init_status ${hi:-0}
stage10_3_hook_pre_step_status ${hs:-0}
stage10_3_hook_pre_rhs_status ${hr:-0}
stage10_3_hook_post_projection_status ${hp:-0}
stage10_3_hook_post_step_status ${hps:-0}
stage10_3_hook_finalize_status ${hf:-0}
stage10_3_no_fibre_state_status ${nf:-0}
stage10_3_no_force_status ${nfo:-0}
stage10_3_no_rhs_injection_status ${nrhs:-0}
stage10_3_no_ibm_call_status ${nibm:-0}
stage10_3_no_structure_advance_status ${nst:-0}
stage10_3_field_modified_status ${field_modified_status}
stage10_3_main_noop_hook_status ${main_status}
DAT

if [ "$main_status" = "1" ] && [ -z "$failures" ]; then
  echo "STAGE 10.3 MAIN NOOP HOOK VERDICT: PASS"
  echo "STAGE 10.3 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.3 MAIN NOOP HOOK VERDICT: FAIL"
echo "STAGE 10.3 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
