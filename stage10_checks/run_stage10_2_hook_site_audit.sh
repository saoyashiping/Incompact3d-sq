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
  if [ -z "$failures" ]; then
    failures="$msg"
  else
    failures="$failures\n$msg"
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

strip_fortran_comments() {
  # Remove text after ! for static active-code checks. This intentionally treats
  # documentation/comments as non-active and keeps the audit portable.
  sed 's/!.*$//' "$1"
}

active_grep() {
  local pattern="$1"
  shift
  for f in "$@"; do
    [ -f "$f" ] || continue
    if strip_fortran_comments "$f" | grep -inE "$pattern" >>"$log_file" 2>&1; then
      return 0
    fi
  done
  return 1
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
non_main_prod_files="src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90"

# Candidate site discovery (static, informational but required as status)
hook_init_site_status=0
hook_pre_step_site_status=0
hook_pre_rhs_site_status=0
hook_post_projection_site_status=0
hook_post_step_site_status=0
hook_finalize_site_status=0

if active_grep "call[[:space:]]+init_xcompact3d" src/xcompact3d.f90; then
  hook_init_site_status=1
else
  add_failure "hook_init candidate not identified"
fi

if active_grep "do[[:space:]]+itime" src/xcompact3d.f90; then
  hook_pre_step_site_status=1
else
  add_failure "hook_pre_step candidate not identified"
fi

if active_grep "call[[:space:]]+momentum" src/xcompact3d.f90; then
  hook_pre_rhs_site_status=1
else
  add_failure "hook_pre_rhs candidate not identified"
fi

if active_grep "call[[:space:]]+pre_correc" src/xcompact3d.f90; then
  hook_post_projection_site_status=1
else
  add_failure "hook_post_projection candidate not identified"
fi

if active_grep "call[[:space:]]+test_speed_min_max" src/xcompact3d.f90; then
  hook_post_step_site_status=1
else
  add_failure "hook_post_step candidate not identified"
fi

if active_grep "call[[:space:]]+finalise_xcompact3d" src/xcompact3d.f90; then
  hook_finalize_site_status=1
else
  add_failure "hook_finalize candidate not identified"
fi

# Stage 10.2 was originally the pre-connection site audit.  After Stage 10.3,
# xcompact3d.f90 legitimately contains guarded hook calls.  Therefore this audit
# has two modes:
#   - pre-10.3 mode: no production hook call is allowed;
#   - 10.3+ mode: only guarded hook calls in xcompact3d.f90 are allowed.
allow_production_hooks="$STAGE10_2_ALLOW_PRODUCTION_HOOKS"
if [ "$allow_production_hooks" = "auto" ]; then
  if active_grep "if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]+stage10_hook_" src/xcompact3d.f90; then
    allow_production_hooks=1
  else
    allow_production_hooks=0
  fi
fi

no_production_hook_call_status=1
for sym in stage10_hook_init stage10_hook_pre_step stage10_hook_pre_rhs stage10_hook_post_projection stage10_hook_post_step stage10_hook_finalize; do
  if [ "$allow_production_hooks" = "1" ]; then
    # Hook symbols are allowed only in xcompact3d.f90 and only as guarded calls.
    if active_grep "$sym" $non_main_prod_files; then
      no_production_hook_call_status=0
      add_failure "forbidden Stage 10 hook symbol found outside xcompact3d.f90: $sym"
    fi
    if active_grep "call[[:space:]]+$sym" src/xcompact3d.f90; then
      if ! strip_fortran_comments src/xcompact3d.f90 | grep -inE "if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]+$sym" >>"$log_file" 2>&1; then
        no_production_hook_call_status=0
        add_failure "unguarded Stage 10 hook call found in xcompact3d.f90: $sym"
      fi
    fi
  else
    if active_grep "$sym" $prod_files; then
      no_production_hook_call_status=0
      add_failure "forbidden production hook symbol found in production files: $sym"
    fi
  fi
done

no_rhs_modification_status=1
if active_grep "f_fsi|fsi_force|two_way|twoway|feedback_force" src/navier.f90 src/time_integrators.f90; then
  no_rhs_modification_status=0
  add_failure "possible RHS-side coupling symbol found"
fi

no_poisson_modification_status=1
if active_grep "stage10_hook_" src/poisson.f90; then
  no_poisson_modification_status=0
  add_failure "forbidden hook reference found in poisson solver"
fi

no_projection_modification_status=1
if active_grep "stage10_hook_" src/time_integrators.f90 src/derive.f90; then
  no_projection_modification_status=0
  add_failure "forbidden hook reference found in projection/derivative path"
fi

no_restart_logic_modification_status=1
# Do not do a loose line-near-restart grep here; Stage 10 hooks may appear in
# xcompact3d.f90 after Stage 10.3.  Instead, reject only active hook calls on
# lines that also mention restart/checkpoint.
if strip_fortran_comments src/xcompact3d.f90 | grep -inE "stage10_hook_.*(restart|checkpoint)|(restart|checkpoint).*stage10_hook_" >>"$log_file" 2>&1; then
  no_restart_logic_modification_status=0
  add_failure "forbidden hook reference found directly in restart/checkpoint logic"
fi

no_stage9_logic_modification_status=1
if grep -inE "stage10_hook_" stage9_checks/*.sh stage9_checks/*.md >>"$log_file" 2>&1; then
  no_stage9_logic_modification_status=0
  add_failure "Stage 9 files reference Stage 10 hooks"
fi

no_ibm_call_status=1
# Only active imports/calls are contamination.  Negative status names such as
# no_ibm_call_status are expected and must not be treated as IBM activity.
if active_grep "^[[:space:]]*use[[:space:]]+fibre_ibm|^[[:space:]]*call[[:space:]].*(ibm|spread|interpol)" src/fibre_stage10_noop_hook.f90; then
  no_ibm_call_status=0
  add_failure "active IBM-related import/call found in Stage 10 hook module"
fi

no_structure_advance_status=1
# Only active structure/tension/bending imports/calls are contamination.  Module
# names and negative diagnostic fields containing fibre/structure words are OK.
if active_grep "^[[:space:]]*use[[:space:]]+(fibre_structure|fibre_tension|fibre_bending)|^[[:space:]]*call[[:space:]].*(structure_advance|fibre_structure|tension|bending)" src/fibre_stage10_noop_hook.f90; then
  no_structure_advance_status=0
  add_failure "active fibre-structure import/call found in Stage 10 hook module"
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
