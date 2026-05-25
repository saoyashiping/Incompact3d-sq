#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}


ensure_build_dir() {
  if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    cmake -S . -B "$BUILD_DIR"
  fi
}

mkdir -p stage10_outputs

dat_file="stage10_outputs/stage10_2_hook_site_audit.dat"
log_file="stage10_outputs/stage10_2_hook_site_audit.log"
: > "$log_file"

failures=""

add_failure() {
  local msg="$1"
  if [ -z "$failures" ]; then
    failures="$msg"
  else
    failures="$failures\n$msg"
  fi
}

status_from_bool() {
  if [ "$1" -eq 0 ]; then
    echo 1
  else
    echo 0
  fi
}

build_status=1
stage10_0_status=1
stage10_1_status=1

run_build_target() {
  local target="$1"
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$target" -j >>"$log_file" 2>&1; then
    build_status=0
    add_failure "build failed for target: $target"
  fi
}

run_build_target xcompact3d
run_build_target fibre_stage10_config_check
run_build_target fibre_stage10_noop_hook_check

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh >>"$log_file" 2>&1; then
    stage10_0_status=0
    add_failure "Stage 10.0 prerequisite failed"
  fi
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_1_hook_build.sh >>"$log_file" 2>&1; then
    stage10_1_status=0
    add_failure "Stage 10.1 prerequisite failed"
  fi
fi

prod_files="src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90"

# Candidate site discovery (static, informational but required as status)
hook_init_site_status=0
hook_pre_step_site_status=0
hook_pre_rhs_site_status=0
hook_post_projection_site_status=0
hook_post_step_site_status=0
hook_finalize_site_status=0

if grep -n "call init_xcompact3d" src/xcompact3d.f90 >>"$log_file" 2>&1; then
  hook_init_site_status=1
else
  add_failure "hook_init candidate not identified"
fi

if grep -n "do itime" src/xcompact3d.f90 >>"$log_file" 2>&1; then
  hook_pre_step_site_status=1
else
  add_failure "hook_pre_step candidate not identified"
fi

if grep -n "call momentum" src/xcompact3d.f90 >>"$log_file" 2>&1; then
  hook_pre_rhs_site_status=1
else
  add_failure "hook_pre_rhs candidate not identified"
fi

if grep -n "call pre_correc" src/xcompact3d.f90 >>"$log_file" 2>&1; then
  hook_post_projection_site_status=1
else
  add_failure "hook_post_projection candidate not identified"
fi

if grep -n "call test_speed_min_max" src/xcompact3d.f90 >>"$log_file" 2>&1; then
  hook_post_step_site_status=1
else
  add_failure "hook_post_step candidate not identified"
fi

if grep -n "call finalise_xcompact3d" src/xcompact3d.f90 >>"$log_file" 2>&1; then
  hook_finalize_site_status=1
else
  add_failure "hook_finalize candidate not identified"
fi

no_production_hook_call_status=1
for sym in stage10_hook_init stage10_hook_pre_step stage10_hook_pre_rhs stage10_hook_post_projection stage10_hook_post_step stage10_hook_finalize; do
  if grep -n "$sym" $prod_files >>"$log_file" 2>&1; then
    no_production_hook_call_status=0
    add_failure "forbidden production hook symbol found in production files: $sym"
  fi
done

no_rhs_modification_status=1
if grep -n "f_fsi\|fsi_force\|two_way\|twoway\|feedback_force" src/navier.f90 src/time_integrators.f90 >>"$log_file" 2>&1; then
  no_rhs_modification_status=0
  add_failure "possible RHS-side coupling symbol found"
fi

no_poisson_modification_status=1
if grep -n "stage10_hook_" src/poisson.f90 >>"$log_file" 2>&1; then
  no_poisson_modification_status=0
  add_failure "forbidden hook reference found in poisson solver"
fi

no_projection_modification_status=1
if grep -n "stage10_hook_" src/time_integrators.f90 src/derive.f90 >>"$log_file" 2>&1; then
  no_projection_modification_status=0
  add_failure "forbidden hook reference found in projection/derivative path"
fi

no_restart_logic_modification_status=1
if grep -n "stage10_hook_" src/xcompact3d.f90 | grep -n "restart" >>"$log_file" 2>&1; then
  no_restart_logic_modification_status=0
  add_failure "forbidden hook reference found near restart logic"
fi

no_stage9_logic_modification_status=1
if grep -n "stage10_hook_" stage9_checks/*.sh stage9_checks/*.md >>"$log_file" 2>&1; then
  no_stage9_logic_modification_status=0
  add_failure "Stage 9 files reference Stage 10 hooks"
fi

no_ibm_call_status=1
if grep -n "ibm\|spread\|interpol" src/fibre_stage10_noop_hook.f90 | grep -v "no_ibm_call_status" >>"$log_file" 2>&1; then
  no_ibm_call_status=0
  add_failure "possible IBM-related activity found in Stage 10 hook module"
fi

no_structure_advance_status=1
if grep -n "structure\|fibre_" src/fibre_stage10_noop_hook.f90 | grep -v "no_structure_advance_status" >>"$log_file" 2>&1; then
  no_structure_advance_status=0
  add_failure "possible fibre structure activity found in Stage 10 hook module"
fi

hook_site_audit_status=1
for v in \
  $build_status $stage10_0_status $stage10_1_status \
  $hook_init_site_status $hook_pre_step_site_status $hook_pre_rhs_site_status \
  $hook_post_projection_site_status $hook_post_step_site_status $hook_finalize_site_status \
  $no_production_hook_call_status $no_rhs_modification_status $no_poisson_modification_status \
  $no_projection_modification_status $no_restart_logic_modification_status $no_stage9_logic_modification_status \
  $no_ibm_call_status $no_structure_advance_status; do
  if [ "$v" -ne 1 ]; then
    hook_site_audit_status=0
  fi
done

{
  echo "stage10_2_requested_flag 1"
  echo "stage10_2_build_status $build_status"
  echo "stage10_2_stage10_0_status $stage10_0_status"
  echo "stage10_2_stage10_1_status $stage10_1_status"
  echo "stage10_2_hook_init_site_status $hook_init_site_status"
  echo "stage10_2_hook_pre_step_site_status $hook_pre_step_site_status"
  echo "stage10_2_hook_pre_rhs_site_status $hook_pre_rhs_site_status"
  echo "stage10_2_hook_post_projection_site_status $hook_post_projection_site_status"
  echo "stage10_2_hook_post_step_site_status $hook_post_step_site_status"
  echo "stage10_2_hook_finalize_site_status $hook_finalize_site_status"
  echo "stage10_2_no_production_hook_call_status $no_production_hook_call_status"
  echo "stage10_2_no_rhs_modification_status $no_rhs_modification_status"
  echo "stage10_2_no_poisson_modification_status $no_poisson_modification_status"
  echo "stage10_2_no_projection_modification_status $no_projection_modification_status"
  echo "stage10_2_no_restart_logic_modification_status $no_restart_logic_modification_status"
  echo "stage10_2_no_stage9_logic_modification_status $no_stage9_logic_modification_status"
  echo "stage10_2_no_ibm_call_status $no_ibm_call_status"
  echo "stage10_2_no_structure_advance_status $no_structure_advance_status"
  echo "stage10_2_hook_site_audit_status $hook_site_audit_status"
} > "$dat_file"

if [ "$hook_site_audit_status" -eq 1 ]; then
  echo "STAGE 10.2 HOOK SITE AUDIT VERDICT: PASS"
  echo "STAGE 10.2 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.2 HOOK SITE AUDIT VERDICT: FAIL"
echo "STAGE 10.2 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
