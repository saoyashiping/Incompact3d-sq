#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-1}

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
stage10_4_status=1
stage10_5_status=1
stage10_6_status=1
stage10_7_status=1
hook_active_status=1
noop_safety_status=1
no_fibre_state_status=1
no_force_status=1
no_rhs_injection_status=1
no_ibm_call_status=1
no_structure_advance_status=1

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

# Stage 10.8 runs only late stable gates by default (10.4-10.7), with forced skip/prereq knobs.
if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_4_noop_invariance_np1.sh; then
  stage10_4_status=0
  add_failure "Stage 10.4 gate failed"
fi
if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_5_RUN_PREREQS=0 bash stage10_checks/run_stage10_5_noop_invariance_parallel.sh; then
  stage10_5_status=0
  add_failure "Stage 10.5 gate failed"
fi
if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_6_RUN_PREREQS=0 bash stage10_checks/run_stage10_6_rhs_contamination_audit.sh; then
  stage10_6_status=0
  add_failure "Stage 10.6 gate failed"
fi
if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_7_RUN_PREREQS=0 bash stage10_checks/run_stage10_7_io_restart_stats_visu_noop.sh; then
  stage10_7_status=0
  add_failure "Stage 10.7 gate failed"
fi

hook_dat="stage10_outputs/fibre_stage10_3_main_noop_hook.dat"
if [ ! -s "$hook_dat" ]; then
  hook_active_status=0
  noop_safety_status=0
  add_failure "missing required hook dat: $hook_dat"
fi

getv() {
  awk -v k="$1" '$1==k{print $2}' "$2" 2>/dev/null | tail -n1
}

need_one() {
  local key="$1"
  local file="$2"
  [ "$(getv "$key" "$file")" = "1" ]
}

need_zero() {
  local key="$1"
  local file="$2"
  [ "$(getv "$key" "$file")" = "0" ]
}

if [ -s "$hook_dat" ]; then
  for k in stage10_3_requested_flag stage10_3_noop_mode_status stage10_3_hook_init_status stage10_3_hook_pre_step_status stage10_3_hook_pre_rhs_status stage10_3_hook_post_projection_status stage10_3_hook_post_step_status stage10_3_hook_finalize_status stage10_3_no_fibre_state_status stage10_3_no_force_status stage10_3_no_rhs_injection_status stage10_3_no_ibm_call_status stage10_3_no_structure_advance_status stage10_3_main_noop_hook_status; do
    if ! need_one "$k" "$hook_dat"; then
      noop_safety_status=0
      add_failure "hook dat key not 1: $k"
    fi
  done
  if ! need_zero stage10_3_field_modified_status "$hook_dat"; then
    noop_safety_status=0
    add_failure "hook dat key not 0: stage10_3_field_modified_status"
  fi

  if ! need_one stage10_3_no_fibre_state_status "$hook_dat"; then no_fibre_state_status=0; fi
  if ! need_one stage10_3_no_force_status "$hook_dat"; then no_force_status=0; fi
  if ! need_one stage10_3_no_rhs_injection_status "$hook_dat"; then no_rhs_injection_status=0; fi
  if ! need_one stage10_3_no_ibm_call_status "$hook_dat"; then no_ibm_call_status=0; fi
  if ! need_one stage10_3_no_structure_advance_status "$hook_dat"; then no_structure_advance_status=0; fi
fi

final_status=1
for v in $build_status $stage10_4_status $stage10_5_status $stage10_6_status $stage10_7_status $hook_active_status $noop_safety_status $no_fibre_state_status $no_force_status $no_rhs_injection_status $no_ibm_call_status $no_structure_advance_status; do
  [ "$v" -eq 1 ] || final_status=0
done

cat > stage10_outputs/stage10_8_total_smoke.dat <<DAT
stage10_8_requested_flag 1
stage10_8_build_status $build_status
stage10_8_stage10_4_status $stage10_4_status
stage10_8_stage10_5_status $stage10_5_status
stage10_8_stage10_6_status $stage10_6_status
stage10_8_stage10_7_status $stage10_7_status
stage10_8_hook_active_status $hook_active_status
stage10_8_noop_safety_status $noop_safety_status
stage10_8_no_fibre_state_status $no_fibre_state_status
stage10_8_no_force_status $no_force_status
stage10_8_no_rhs_injection_status $no_rhs_injection_status
stage10_8_no_ibm_call_status $no_ibm_call_status
stage10_8_no_structure_advance_status $no_structure_advance_status
stage10_8_total_smoke_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  cat > stage10_checks/STAGE10_CLOSED.md <<'MD'
# STAGE10_CLOSED

## Stage 10 purpose

Production coupling hook skeleton / default no-op closure.

## Closed sub-stages

- Stage 10.0 config and global switches
- Stage 10.1 no-op hook build
- Stage 10.2 hook placement audit
- Stage 10.3 main-loop no-op hook connection
- Stage 10.4 np=1 no-op invariance
- Stage 10.5 parallel no-op invariance
- Stage 10.6 RHS contamination audit
- Stage 10.7 restart/stats/visu/coarse I/O no-op validation
- Stage 10.8 total smoke closure

## Governing no-op model

`div(u)=0`

`du/dt + u·grad(u) = -grad(p) + nu laplacian(u) + f_drive + f_fsi`

`f_fsi=0`

`RHS_stage10=RHS_stage9`

## Explicit statement

Stage 10 closes production hook skeleton/default no-op only.
No real IBM, no real fibre force, no feedback force, and no two-way coupling are activated.
Real production one-way fluid-to-fibre hook begins in Stage 11.

## Next recommended stage

Stage 11 production one-way fluid-to-fibre hook.
MD

  echo "STAGE 10.8 TOTAL SMOKE VERDICT: PASS"
  echo "STAGE 10.8 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.8 TOTAL SMOKE VERDICT: FAIL"
echo "STAGE 10.8 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
