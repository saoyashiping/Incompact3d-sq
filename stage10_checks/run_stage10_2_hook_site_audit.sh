#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}
STAGE10_2_ALLOW_PRODUCTION_HOOKS=${STAGE10_2_ALLOW_PRODUCTION_HOOKS:-auto}

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
  if [ -z "$failures" ]; then failures="$msg"; else failures="$failures\n$msg"; fi
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

strip_comments() {
  # Stage 10 audit is intentionally conservative: strip comments before checking active code.
  sed 's/!.*$//' "$1"
}

tmpdir="stage10_outputs/stage10_2_active_sources"
rm -rf "$tmpdir"
mkdir -p "$tmpdir"
for f in src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90 src/fibre_stage10_noop_hook.f90; do
  strip_comments "$f" > "$tmpdir/$(basename "$f").active"
done

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

# Candidate site discovery. These are placement evidence only; Stage 10.3+ may already have valid guarded calls.
hook_init_site_status=0
hook_pre_step_site_status=0
hook_pre_rhs_site_status=0
hook_post_projection_site_status=0
hook_post_step_site_status=0
hook_finalize_site_status=0

grep -n "call init_xcompact3d" src/xcompact3d.f90 >>"$log_file" 2>&1 && hook_init_site_status=1 || add_failure "hook_init candidate not identified"
grep -n "do itime" src/xcompact3d.f90 >>"$log_file" 2>&1 && hook_pre_step_site_status=1 || add_failure "hook_pre_step candidate not identified"
grep -n "call momentum\|call calculate_transeq_rhs" src/xcompact3d.f90 >>"$log_file" 2>&1 && hook_pre_rhs_site_status=1 || add_failure "hook_pre_rhs candidate not identified"
grep -n "call pre_correc\|call cor_vel" src/xcompact3d.f90 >>"$log_file" 2>&1 && hook_post_projection_site_status=1 || add_failure "hook_post_projection candidate not identified"
grep -n "call test_speed_min_max\|call postprocessing" src/xcompact3d.f90 >>"$log_file" 2>&1 && hook_post_step_site_status=1 || add_failure "hook_post_step candidate not identified"
grep -n "call finalise_xcompact3d" src/xcompact3d.f90 >>"$log_file" 2>&1 && hook_finalize_site_status=1 || add_failure "hook_finalize candidate not identified"

no_production_hook_call_status=1
x3d_active="$tmpdir/xcompact3d.f90.active"

# After Stage 10.3, guarded hook calls in xcompact3d.f90 are valid and must not be rejected.
allow_prod_hooks="$STAGE10_2_ALLOW_PRODUCTION_HOOKS"
if [ "$allow_prod_hooks" = "auto" ]; then
  if grep -Eiq '^[[:space:]]*if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]*stage10_hook_' "$x3d_active"; then
    allow_prod_hooks=1
  else
    allow_prod_hooks=0
  fi
fi

if [ "$allow_prod_hooks" = "1" ]; then
  # Fail only on unguarded hook calls in xcompact3d.f90.
  if grep -Ein '^[[:space:]]*call[[:space:]]+stage10_hook_|^[[:space:]]*if[[:space:]]*\([^)]*\)[[:space:]]*call[[:space:]]*stage10_hook_' "$x3d_active" \
      | grep -Eiv '^[0-9]+:[[:space:]]*if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]*stage10_hook_' >>"$log_file" 2>&1; then
    no_production_hook_call_status=0
    add_failure "unguarded or incorrectly guarded Stage 10 hook call found in xcompact3d.f90"
  fi
else
  if grep -Eiq 'stage10_hook_(init|pre_step|pre_rhs|post_projection|post_step|finalize)' "$x3d_active"; then
    no_production_hook_call_status=0
    add_failure "production hook call found before Stage 10.3 connection mode"
  fi
fi

# Hook calls are forbidden in RHS/projection/Poisson/derivative/case production files.
for f in navier.f90 time_integrators.f90 derive.f90 poisson.f90 Case-Channel.f90; do
  if grep -Eiq 'stage10_hook_(init|pre_step|pre_rhs|post_projection|post_step|finalize)' "$tmpdir/$f.active"; then
    no_production_hook_call_status=0
    add_failure "forbidden Stage 10 hook call/reference found in $f"
  fi
done

no_rhs_modification_status=1
for f in navier.f90 time_integrators.f90 Case-Channel.f90; do
  if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_stage8_(twoway_force_density|oneway_forcing|feedback_candidate)|^[[:space:]]*use[[:space:]]+fibre_ibm_(spreading|feedback|force_buffer)' "$tmpdir/$f.active"; then
    no_rhs_modification_status=0
    add_failure "$f: active coupling module import in RHS-related path"
  fi
  if grep -Eiq '^[[:space:]]*call[[:space:]].*(twoway|feedback|spreading|force_density|rhs_injection|fibre.*force|ibm.*force)' "$tmpdir/$f.active"; then
    no_rhs_modification_status=0
    add_failure "$f: active coupling call in RHS-related path"
  fi
  if grep -Eiq '(f_fsi|f_ibm|f_feedback|force_density|fibre_force|twoway_force|rhs_fibre|rhs_ibm)' "$tmpdir/$f.active"; then
    no_rhs_modification_status=0
    add_failure "$f: active RHS coupling assignment/symbol found"
  fi
done

no_poisson_modification_status=1
if grep -Eiq 'stage10_hook_|^[[:space:]]*use[[:space:]]+fibre_.*(force|ibm)' "$tmpdir/poisson.f90.active"; then
  no_poisson_modification_status=0
  add_failure "forbidden hook/coupling reference found in Poisson solver"
fi

no_projection_modification_status=1
if grep -Eiq 'stage10_hook_|^[[:space:]]*use[[:space:]]+fibre_.*(force|ibm)' "$tmpdir/time_integrators.f90.active" "$tmpdir/derive.f90.active"; then
  no_projection_modification_status=0
  add_failure "forbidden hook/coupling reference found in projection/derivative path"
fi

no_restart_logic_modification_status=1
# Valid guarded hook calls in xcompact3d are allowed. Only reject a hook line that explicitly references restart text.
if grep -Ein 'stage10_hook_.*restart|restart.*stage10_hook_' "$x3d_active" >>"$log_file" 2>&1; then
  no_restart_logic_modification_status=0
  add_failure "Stage 10 hook reference found directly in restart logic"
fi

no_stage9_logic_modification_status=1
if grep -n "stage10_hook_" stage9_checks/*.sh stage9_checks/*.md >>"$log_file" 2>&1; then
  no_stage9_logic_modification_status=0
  add_failure "Stage 9 files reference Stage 10 hooks"
fi

no_ibm_call_status=1
hook_active="$tmpdir/fibre_stage10_noop_hook.f90.active"
if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_ibm|^[[:space:]]*call[[:space:]].*(ibm|spread|interpol)' "$hook_active"; then
  no_ibm_call_status=0
  add_failure "active IBM use/call found in Stage 10 no-op hook module"
fi

no_structure_advance_status=1
if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_(structure|tension|bending)|^[[:space:]]*call[[:space:]].*(structure_advance|advance_structure|tension|bending)' "$hook_active"; then
  no_structure_advance_status=0
  add_failure "active fibre-structure use/call found in Stage 10 no-op hook module"
fi

hook_site_audit_status=1
for v in \
  $build_status $stage10_0_status $stage10_1_status \
  $hook_init_site_status $hook_pre_step_site_status $hook_pre_rhs_site_status \
  $hook_post_projection_site_status $hook_post_step_site_status $hook_finalize_site_status \
  $no_production_hook_call_status $no_rhs_modification_status $no_poisson_modification_status \
  $no_projection_modification_status $no_restart_logic_modification_status $no_stage9_logic_modification_status \
  $no_ibm_call_status $no_structure_advance_status; do
  [ "$v" -eq 1 ] || hook_site_audit_status=0
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
