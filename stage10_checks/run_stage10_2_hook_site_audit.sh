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
    if [ -n "$DECOMP2D_ROOT" ]; then
      cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT"
    else
      cmake -S . -B "$BUILD_DIR"
    fi
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
    failures="$failures
$msg"
  fi
}

run_build_target() {
  local target="$1"
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$target" -j >>"$log_file" 2>&1; then
    build_status=0
    add_failure "build failed for target: $target"
  fi
}

strip_comments() {
  sed 's/!.*$//' "$1"
}

active_grep() {
  local pattern="$1"
  shift
  grep -Eiq "$pattern" "$@"
}

build_status=1
stage10_0_status=1
stage10_1_status=1
hook_init_site_status=1
hook_pre_step_site_status=1
hook_pre_rhs_site_status=1
hook_post_projection_site_status=1
hook_post_step_site_status=1
hook_finalize_site_status=1
no_production_hook_call_status=1
no_rhs_modification_status=1
no_poisson_modification_status=1
no_projection_modification_status=1
no_restart_logic_modification_status=1
no_stage9_logic_modification_status=1
no_ibm_call_status=1
no_structure_advance_status=1

run_build_target xcompact3d
run_build_target fibre_stage10_config_check
run_build_target fibre_stage10_noop_hook_check

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_0_config.sh >>"$log_file" 2>&1 || {
    stage10_0_status=0
    add_failure "Stage 10.0 prerequisite failed"
  }
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 bash stage10_checks/run_stage10_1_hook_build.sh >>"$log_file" 2>&1 || {
    stage10_1_status=0
    add_failure "Stage 10.1 prerequisite failed"
  }
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

prod_files="src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90"
for f in $prod_files src/fibre_stage10_noop_hook.f90; do
  strip_comments "$f" > "$tmpdir/$(basename "$f").active"
done

# Candidate future hook sites. These checks are intentionally descriptive:
# Stage 10.2 documents where hooks should be placed; Stage 10.3+ may already
# have inserted valid guarded calls in xcompact3d.f90.
if ! active_grep 'call[[:space:]]+init_xcompact3d' "$tmpdir/xcompact3d.f90.active"; then
  hook_init_site_status=0
  add_failure "hook_init candidate site not found near init_xcompact3d"
fi
if ! active_grep 'do[[:space:]]+itime' "$tmpdir/xcompact3d.f90.active"; then
  hook_pre_step_site_status=0
  add_failure "hook_pre_step candidate outer time-loop site not found"
fi
if ! active_grep 'call[[:space:]]+momentum' "$tmpdir/xcompact3d.f90.active"; then
  hook_pre_rhs_site_status=0
  add_failure "hook_pre_rhs candidate site near momentum/RHS not found"
fi
if ! active_grep 'call[[:space:]]+pre_correc' "$tmpdir/xcompact3d.f90.active"; then
  hook_post_projection_site_status=0
  add_failure "hook_post_projection candidate site near projection/correction not found"
fi
if ! active_grep 'call[[:space:]]+test_speed_min_max' "$tmpdir/xcompact3d.f90.active"; then
  hook_post_step_site_status=0
  add_failure "hook_post_step candidate diagnostic site not found"
fi
if ! active_grep 'call[[:space:]]+finalise_xcompact3d' "$tmpdir/xcompact3d.f90.active"; then
  hook_finalize_site_status=0
  add_failure "hook_finalize candidate finalise site not found"
fi

# Stage 10.2 originally forbade production hook calls. After Stage 10.3,
# guarded calls in xcompact3d.f90 are valid and must not be treated as
# contamination. The audit still forbids hook calls in RHS/projection/Poisson/
# derivative/channel internals.
allow_hooks="$STAGE10_2_ALLOW_PRODUCTION_HOOKS"
if [ "$allow_hooks" = "auto" ]; then
  if active_grep '^[[:space:]]*if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]+stage10_hook_' "$tmpdir/xcompact3d.f90.active"; then
    allow_hooks=1
  else
    allow_hooks=0
  fi
fi

if active_grep '^[[:space:]]*call[[:space:]]+stage10_hook_' "$tmpdir/xcompact3d.f90.active"; then
  if [ "$allow_hooks" != "1" ]; then
    no_production_hook_call_status=0
    add_failure "unallowed Stage 10 hook call found in xcompact3d.f90"
  else
    total_calls=$(grep -Ei '^[[:space:]]*(if[[:space:]]*\([^)]*\)[[:space:]]*)?call[[:space:]]+stage10_hook_' "$tmpdir/xcompact3d.f90.active" | wc -l | awk '{print $1}')
    guarded_calls=$(grep -Ei '^[[:space:]]*if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]+stage10_hook_' "$tmpdir/xcompact3d.f90.active" | wc -l | awk '{print $1}')
    if [ "$total_calls" -ne "$guarded_calls" ]; then
      no_production_hook_call_status=0
      add_failure "Stage 10 hook calls in xcompact3d.f90 are not all guarded by if (stage10_reg)"
    fi
  fi
fi

for f in navier.f90 time_integrators.f90 derive.f90 poisson.f90 Case-Channel.f90; do
  if active_grep '^[[:space:]]*call[[:space:]]+stage10_hook_' "$tmpdir/$f.active"; then
    no_production_hook_call_status=0
    add_failure "$f: Stage 10 hook call found outside xcompact3d main loop"
  fi
done

# No RHS/projection/Poisson/restart/Stage9 logic contamination.
if active_grep '(f_fsi|fsi_force|two_way|twoway|feedback_force)' "$tmpdir/navier.f90.active" "$tmpdir/time_integrators.f90.active"; then
  no_rhs_modification_status=0
  add_failure "possible RHS force-injection token found in active RHS files"
fi

if active_grep 'stage10_hook_' "$tmpdir/poisson.f90.active"; then
  no_poisson_modification_status=0
  add_failure "stage10 hook found in poisson internals"
fi

if active_grep 'stage10_hook_' "$tmpdir/time_integrators.f90.active" "$tmpdir/derive.f90.active"; then
  no_projection_modification_status=0
  add_failure "stage10 hook found in time_integrators/derive internals"
fi

if active_grep 'stage10_hook_.*restart|restart.*stage10_hook_' "$tmpdir/xcompact3d.f90.active"; then
  no_restart_logic_modification_status=0
  add_failure "stage10 hook appears entangled with restart logic"
fi

if grep -Eiq 'stage10_hook_' stage9_checks/*.sh stage9_checks/*.md 2>/dev/null; then
  no_stage9_logic_modification_status=0
  add_failure "Stage 9 script/documentation contains Stage 10 hook reference"
fi

# Precise active-code checks only. Do not broad-grep ibm/structure/fibre_ and
# do not classify negative diagnostic fields as activity.
if active_grep '^[[:space:]]*use[[:space:]]+fibre_ibm|^[[:space:]]*call[[:space:]].*(ibm|spread|interpol)' "$tmpdir/fibre_stage10_noop_hook.f90.active"; then
  no_ibm_call_status=0
  add_failure "active IBM use/call found in Stage 10 no-op hook module"
fi

if active_grep '^[[:space:]]*use[[:space:]]+fibre_(structure|tension|bending)|^[[:space:]]*call[[:space:]].*(structure_advance|tension|bending)' "$tmpdir/fibre_stage10_noop_hook.f90.active"; then
  no_structure_advance_status=0
  add_failure "active fibre-structure use/call found in Stage 10 no-op hook module"
fi

final_status=1
for v in "$build_status" "$stage10_0_status" "$stage10_1_status" \
         "$hook_init_site_status" "$hook_pre_step_site_status" "$hook_pre_rhs_site_status" \
         "$hook_post_projection_site_status" "$hook_post_step_site_status" "$hook_finalize_site_status" \
         "$no_production_hook_call_status" "$no_rhs_modification_status" \
         "$no_poisson_modification_status" "$no_projection_modification_status" \
         "$no_restart_logic_modification_status" "$no_stage9_logic_modification_status" \
         "$no_ibm_call_status" "$no_structure_advance_status"; do
  [ "$v" -eq 1 ] || final_status=0
done

cat > "$dat_file" <<DAT
stage10_2_requested_flag 1
stage10_2_build_status $build_status
stage10_2_stage10_0_status $stage10_0_status
stage10_2_stage10_1_status $stage10_1_status
stage10_2_hook_init_site_status $hook_init_site_status
stage10_2_hook_pre_step_site_status $hook_pre_step_site_status
stage10_2_hook_pre_rhs_site_status $hook_pre_rhs_site_status
stage10_2_hook_post_projection_site_status $hook_post_projection_site_status
stage10_2_hook_post_step_site_status $hook_post_step_site_status
stage10_2_hook_finalize_site_status $hook_finalize_site_status
stage10_2_no_production_hook_call_status $no_production_hook_call_status
stage10_2_no_rhs_modification_status $no_rhs_modification_status
stage10_2_no_poisson_modification_status $no_poisson_modification_status
stage10_2_no_projection_modification_status $no_projection_modification_status
stage10_2_no_restart_logic_modification_status $no_restart_logic_modification_status
stage10_2_no_stage9_logic_modification_status $no_stage9_logic_modification_status
stage10_2_no_ibm_call_status $no_ibm_call_status
stage10_2_no_structure_advance_status $no_structure_advance_status
stage10_2_hook_site_audit_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 10.2 HOOK SITE AUDIT VERDICT: PASS"
  echo "STAGE 10.2 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.2 HOOK SITE AUDIT VERDICT: FAIL"
echo "STAGE 10.2 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
