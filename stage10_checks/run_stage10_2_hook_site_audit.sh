#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}

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

ensure_build_dir() {
  if [ -f "$BUILD_DIR/CMakeCache.txt" ]; then
    return 0
  fi

  echo "[INFO] configuring ${BUILD_DIR}" | tee -a "$log_file" >/dev/null
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake -S . -B "$BUILD_DIR" >>"$log_file" 2>&1; then
    add_failure "cmake configure failed for ${BUILD_DIR}"
    return 1
  fi
  return 0
}

build_status=1
stage10_0_status=1
stage10_1_status=1

run_build_target() {
  local target="$1"
  if ! ensure_build_dir; then
    build_status=0
    add_failure "build setup failed before target: $target"
    return 1
  fi

  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$target" -j >>"$log_file" 2>&1; then
    build_status=0
    add_failure "build failed for target: $target"
    return 1
  fi
  return 0
}

# Configure once up front so fresh-unzip runs do not misclassify missing build
# directories as compiler failures.
ensure_build_dir || build_status=0

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  echo "[INFO] Stage 10.2 prerequisite: Stage 10.0" | tee -a "$log_file" >/dev/null
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh >>"$log_file" 2>&1; then
    stage10_0_status=0
    add_failure "Stage 10.0 prerequisite failed"
  fi

  echo "[INFO] Stage 10.2 prerequisite: Stage 10.1" | tee -a "$log_file" >/dev/null
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_1_hook_build.sh >>"$log_file" 2>&1; then
    stage10_1_status=0
    add_failure "Stage 10.1 prerequisite failed"
  fi
fi

# Rebuild the required targets after prerequisites. This keeps Stage 10.2
# self-contained and makes build failures attributable to the current tree.
run_build_target xcompact3d || true
run_build_target fibre_stage10_config_check || true
run_build_target fibre_stage10_noop_hook_check || true

prod_files="src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90"

# Helper: strip simple Fortran comments before searching for active statements.
active_text() {
  sed 's/!.*$//' "$@"
}

# Candidate site discovery (static audit only; no source modification).
hook_init_site_status=0
hook_pre_step_site_status=0
hook_pre_rhs_site_status=0
hook_post_projection_site_status=0
hook_post_step_site_status=0
hook_finalize_site_status=0

if active_text src/xcompact3d.f90 | grep -n "call[[:space:]]\+init_xcompact3d" >>"$log_file" 2>&1; then
  hook_init_site_status=1
else
  add_failure "hook_init candidate not identified"
fi

if active_text src/xcompact3d.f90 | grep -n "do[[:space:]]\+itime" >>"$log_file" 2>&1; then
  hook_pre_step_site_status=1
else
  add_failure "hook_pre_step candidate not identified"
fi

if active_text src/xcompact3d.f90 | grep -n "call[[:space:]]\+calculate_transeq_rhs" >>"$log_file" 2>&1; then
  hook_pre_rhs_site_status=1
else
  add_failure "hook_pre_rhs candidate not identified"
fi

if active_text src/xcompact3d.f90 | grep -n "call[[:space:]]\+cor_vel" >>"$log_file" 2>&1; then
  hook_post_projection_site_status=1
else
  add_failure "hook_post_projection candidate not identified"
fi

if active_text src/xcompact3d.f90 | grep -n "call[[:space:]]\+test_speed_min_max" >>"$log_file" 2>&1; then
  hook_post_step_site_status=1
else
  add_failure "hook_post_step candidate not identified"
fi

if active_text src/xcompact3d.f90 | grep -n "call[[:space:]]\+finalise_xcompact3d" >>"$log_file" 2>&1; then
  hook_finalize_site_status=1
else
  add_failure "hook_finalize candidate not identified"
fi

no_production_hook_call_status=1
for sym in stage10_hook_init stage10_hook_pre_step stage10_hook_pre_rhs stage10_hook_post_projection stage10_hook_post_step stage10_hook_finalize; do
  if active_text $prod_files | grep -n "call[[:space:]]\+${sym}\>" >>"$log_file" 2>&1; then
    no_production_hook_call_status=0
    add_failure "forbidden production hook call found in production files: $sym"
  fi
  if active_text $prod_files | grep -n "use[[:space:]].*fibre_stage10_noop_hook" >>"$log_file" 2>&1; then
    no_production_hook_call_status=0
    add_failure "forbidden production import of fibre_stage10_noop_hook found"
    break
  fi
done

no_rhs_modification_status=1
if active_text src/navier.f90 src/time_integrators.f90 | grep -nE "f_fsi|fsi_force|two_way|twoway|feedback_force|stage10_hook_|fibre_ibm_spreading|fibre_stage8_.*force" >>"$log_file" 2>&1; then
  no_rhs_modification_status=0
  add_failure "possible active RHS-side coupling symbol found"
fi

no_poisson_modification_status=1
if active_text src/poisson.f90 | grep -n "stage10_hook_" >>"$log_file" 2>&1; then
  no_poisson_modification_status=0
  add_failure "forbidden hook reference found in poisson solver"
fi

no_projection_modification_status=1
if active_text src/time_integrators.f90 src/derive.f90 | grep -n "stage10_hook_" >>"$log_file" 2>&1; then
  no_projection_modification_status=0
  add_failure "forbidden hook reference found in projection/derivative path"
fi

no_restart_logic_modification_status=1
if active_text src/xcompact3d.f90 | grep -n "stage10_hook_" | grep -n "restart" >>"$log_file" 2>&1; then
  no_restart_logic_modification_status=0
  add_failure "forbidden hook reference found near restart logic"
fi

no_stage9_logic_modification_status=1
if active_text stage9_checks/*.sh stage9_checks/*.md | grep -n "stage10_hook_" >>"$log_file" 2>&1; then
  no_stage9_logic_modification_status=0
  add_failure "Stage 9 files reference Stage 10 hooks"
fi

no_ibm_call_status=1
if active_text src/fibre_stage10_noop_hook.f90 | grep -nE "^[[:space:]]*use[[:space:]].*fibre_ibm|^[[:space:]]*call[[:space:]].*(ibm|spread|interpol)" >>"$log_file" 2>&1; then
  no_ibm_call_status=0
  add_failure "possible active IBM call/import found in Stage 10 hook module"
fi

no_structure_advance_status=1
if active_text src/fibre_stage10_noop_hook.f90 | grep -nE "^[[:space:]]*use[[:space:]].*fibre_(structure|types|tension|implicit|bending)|^[[:space:]]*call[[:space:]].*(structure|advance|tension|bending)" >>"$log_file" 2>&1; then
  no_structure_advance_status=0
  add_failure "possible active fibre-structure call/import found in Stage 10 hook module"
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
