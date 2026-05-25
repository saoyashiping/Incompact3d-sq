#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_SKIP_PREREQS=${STAGE10_SKIP_PREREQS:-0}
STAGE10_2_ALLOW_PRODUCTION_HOOKS=${STAGE10_2_ALLOW_PRODUCTION_HOOKS:-auto}

ensure_build_dir() {
  if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    cmake -S . -B "$BUILD_DIR"
  fi
}

mkdir -p stage10_outputs
log_file="stage10_outputs/stage10_2_hook_site_audit.log"
: > "$log_file"

failures=""
add_failure() {
  local msg="$1"
  if [ -z "$failures" ]; then failures="$msg"; else failures="$failures\n$msg"; fi
  echo "[FAIL] $msg" >> "$log_file"
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

build_target() {
  local tgt="$1"
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$tgt" -j >> "$log_file" 2>&1; then
    build_status=0
    add_failure "build failed: $tgt"
  fi
}

build_target xcompact3d
build_target fibre_stage10_config_check
build_target fibre_stage10_noop_hook_check

if [ "$STAGE10_SKIP_PREREQS" != "1" ]; then
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_0_config.sh >> "$log_file" 2>&1; then
    stage10_0_status=0
    add_failure "Stage 10.0 prerequisite failed"
  fi
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" STAGE10_SKIP_PREREQS=1 \
       bash stage10_checks/run_stage10_1_hook_build.sh >> "$log_file" 2>&1; then
    stage10_1_status=0
    add_failure "Stage 10.1 prerequisite failed"
  fi
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
strip_comments() { sed 's/!.*$//' "$1"; }

prod_files="src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90"
for f in $prod_files; do
  strip_comments "$f" > "$tmpdir/$(basename "$f").active"
done
strip_comments src/fibre_stage10_noop_hook.f90 > "$tmpdir/fibre_stage10_noop_hook.f90.active"

# Candidate site checks are intentionally conservative static anchors. They document safe future/current hook regions.
grep -Eiq '^[[:space:]]*call[[:space:]]+init_xcompact3d\b' "$tmpdir/xcompact3d.f90.active" || { hook_init_site_status=0; add_failure "hook_init candidate anchor not found"; }
grep -Eiq '^[[:space:]]*do[[:space:]]+itime\b|^[[:space:]]*do[[:space:]]+while' "$tmpdir/xcompact3d.f90.active" || { hook_pre_step_site_status=0; add_failure "hook_pre_step candidate loop anchor not found"; }
grep -Eiq '^[[:space:]]*call[[:space:]]+(momentum|intt|solve)' "$tmpdir/xcompact3d.f90.active" || { hook_pre_rhs_site_status=0; add_failure "hook_pre_rhs candidate time-advance/RHS anchor not found"; }
grep -Eiq '^[[:space:]]*call[[:space:]]+(pre_correc|correc|test_speed_min_max)' "$tmpdir/xcompact3d.f90.active" || { hook_post_projection_site_status=0; add_failure "hook_post_projection candidate anchor not found"; }
grep -Eiq '^[[:space:]]*call[[:space:]]+(test_speed_min_max|statistics|visu|restart)' "$tmpdir/xcompact3d.f90.active" || { hook_post_step_site_status=0; add_failure "hook_post_step candidate anchor not found"; }
grep -Eiq '^[[:space:]]*call[[:space:]]+finalise_xcompact3d\b|^[[:space:]]*end[[:space:]]+program' "$tmpdir/xcompact3d.f90.active" || { hook_finalize_site_status=0; add_failure "hook_finalize candidate anchor not found"; }

x3d_active="$tmpdir/xcompact3d.f90.active"
guarded_hook_count=$(grep -Eic '^[[:space:]]*if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]+stage10_hook_' "$x3d_active" || true)
any_x3d_hook_count=$(grep -Eic '^[[:space:]]*call[[:space:]]+stage10_hook_|^[[:space:]]*if[[:space:]]*\([^)]*\)[[:space:]]*call[[:space:]]+stage10_hook_' "$x3d_active" || true)

allow_hooks="$STAGE10_2_ALLOW_PRODUCTION_HOOKS"
if [ "$allow_hooks" = "auto" ]; then
  if [ "$guarded_hook_count" -gt 0 ]; then allow_hooks=1; else allow_hooks=0; fi
fi

if [ "$allow_hooks" = "1" ]; then
  if [ "$any_x3d_hook_count" -ne "$guarded_hook_count" ]; then
    no_production_hook_call_status=0
    add_failure "xcompact3d.f90 contains unguarded or non-stage10_reg Stage 10 hook calls"
  fi
else
  if grep -Eiq 'stage10_hook_' "$x3d_active"; then
    no_production_hook_call_status=0
    add_failure "Stage 10 hook call found before hook-connection stage"
  fi
fi

# Hook calls are never allowed in RHS/projection/poisson/channel implementation files.
for f in navier.f90 time_integrators.f90 derive.f90 poisson.f90 Case-Channel.f90; do
  if grep -Eiq '^[[:space:]]*(if[[:space:]]*\([^)]*\)[[:space:]]*)?call[[:space:]]+stage10_hook_' "$tmpdir/$f.active"; then
    no_production_hook_call_status=0
    add_failure "$f contains Stage 10 hook call; hooks must remain outside RHS/projection internals"
  fi
done

# Do not use broad ibm/structure/fibre keyword greps. Only active use/call statements count.
hook_active="$tmpdir/fibre_stage10_noop_hook.f90.active"
if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_ibm|^[[:space:]]*call[[:space:]].*(ibm|spread|spreading|interpol|interpolation)' "$hook_active"; then
  no_ibm_call_status=0
  add_failure "active IBM/interpolation/spreading use/call found in Stage 10 no-op hook"
fi
if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_(structure|tension|bending)|^[[:space:]]*call[[:space:]].*(structure_advance|advance_structure|tension|bending)' "$hook_active"; then
  no_structure_advance_status=0
  add_failure "active fibre-structure/tension/bending use/call found in Stage 10 no-op hook"
fi

# Active coupling symbols must not appear in RHS/projection/poisson production files.
for f in navier.f90 time_integrators.f90 Case-Channel.f90; do
  af="$tmpdir/$f.active"
  if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_stage8_(twoway_force_density|oneway_forcing|feedback_candidate)|^[[:space:]]*use[[:space:]]+fibre_ibm_(spreading|feedback|force_buffer)|^[[:space:]]*call[[:space:]].*(twoway|feedback|spreading|force_density|rhs_injection|fibre.*force)' "$af"; then
    no_rhs_modification_status=0
    add_failure "$f contains active RHS coupling contamination pattern"
  fi
done

if grep -Eiq '^[[:space:]]*(use|call)[[:space:]].*(fibre|ibm|stage10_hook)' "$tmpdir/poisson.f90.active"; then
  no_poisson_modification_status=0
  add_failure "poisson.f90 contains active fibre/IBM/hook contamination pattern"
fi
if grep -Eiq '^[[:space:]]*(use|call)[[:space:]].*(fibre|ibm|stage10_hook)' "$tmpdir/derive.f90.active"; then
  no_projection_modification_status=0
  add_failure "derive.f90 contains active fibre/IBM/hook contamination pattern"
fi
if grep -Eiq 'stage10_hook_.*restart|restart.*stage10_hook_' "$x3d_active"; then
  no_restart_logic_modification_status=0
  add_failure "Stage 10 hook appears to be mixed with restart logic"
fi
if grep -Eiq 'stage10_hook_' stage9_checks/*.sh stage9_checks/*.md 2>/dev/null; then
  no_stage9_logic_modification_status=0
  add_failure "Stage 9 files contain Stage 10 hook reference"
fi

hook_site_audit_status=1
for v in $build_status $stage10_0_status $stage10_1_status \
         $hook_init_site_status $hook_pre_step_site_status $hook_pre_rhs_site_status \
         $hook_post_projection_site_status $hook_post_step_site_status $hook_finalize_site_status \
         $no_production_hook_call_status $no_rhs_modification_status $no_poisson_modification_status \
         $no_projection_modification_status $no_restart_logic_modification_status \
         $no_stage9_logic_modification_status $no_ibm_call_status $no_structure_advance_status; do
  [ "$v" -eq 1 ] || hook_site_audit_status=0
done

cat > stage10_outputs/stage10_2_hook_site_audit.dat <<DAT
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
stage10_2_hook_site_audit_status $hook_site_audit_status
DAT

if [ "$hook_site_audit_status" -eq 1 ]; then
  echo "STAGE 10.2 HOOK SITE AUDIT VERDICT: PASS"
  echo "STAGE 10.2 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.2 HOOK SITE AUDIT VERDICT: FAIL"
echo "STAGE 10.2 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
