#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE10_6_RUN_PREREQS=${STAGE10_6_RUN_PREREQS:-0}

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

mkdir -p stage10_outputs

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
no_stage10_hook_in_rhs_files_status=1
xcompact3d_guarded_hook_status=1
no_ibm_spreading_rhs_status=1
no_feedback_rhs_status=1
no_twoway_force_density_rhs_status=1
no_fibre_force_rhs_status=1
no_rhs_assignment_contamination_status=1
no_poisson_contamination_status=1
no_projection_contamination_status=1
no_rk3_contamination_status=1
rhs_equivalence_status=1

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

if [ "$STAGE10_6_RUN_PREREQS" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_5_noop_invariance_parallel.sh || add_failure "optional Stage10.5 prerequisite failed"
fi

strip_comments() {
  sed 's/!.*$//' "$1"
}

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

for f in src/xcompact3d.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/poisson.f90 src/Case-Channel.f90; do
  strip_comments "$f" > "$tmpdir/$(basename "$f").active"
done

# Guarded hook calls in xcompact3d are allowed; require at least one guarded hook call
if ! grep -Eiq '^[[:space:]]*if[[:space:]]*\([[:space:]]*stage10_reg[[:space:]]*\)[[:space:]]*call[[:space:]]*stage10_hook_' "$tmpdir/xcompact3d.f90.active"; then
  xcompact3d_guarded_hook_status=0
  add_failure "xcompact3d: missing guarded stage10 hook call pattern"
fi

# No hook calls allowed in RHS/projection/poisson/channel files
for f in navier.f90 time_integrators.f90 derive.f90 poisson.f90 Case-Channel.f90; do
  if grep -Eiq '\bstage10_hook_(init|pre_step|pre_rhs|post_projection|post_step|finalize)\b' "$tmpdir/$f.active"; then
    no_stage10_hook_in_rhs_files_status=0
    add_failure "$f: stage10_hook_* call/reference found in active code"
  fi
done

# Contamination imports/calls checks in rhs-relevant files only
rhs_files="navier.f90 time_integrators.f90 Case-Channel.f90"
for f in $rhs_files; do
  af="$tmpdir/$f.active"
  if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_ibm_spreading\b|^[[:space:]]*use[[:space:]]+fibre_ibm_feedback\b|^[[:space:]]*use[[:space:]]+fibre_ibm_force_buffer\b|^[[:space:]]*call[[:space:]].*spreading|^[[:space:]]*call[[:space:]].*ibm.*feedback' "$af"; then
    no_ibm_spreading_rhs_status=0
    add_failure "$f: IBM spreading/feedback contamination pattern found"
  fi
  if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_stage8_feedback_candidate\b|^[[:space:]]*call[[:space:]].*feedback' "$af"; then
    no_feedback_rhs_status=0
    add_failure "$f: feedback contamination pattern found"
  fi
  if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_stage8_twoway_force_density\b|^[[:space:]]*call[[:space:]].*twoway|^[[:space:]]*call[[:space:]].*force_density' "$af"; then
    no_twoway_force_density_rhs_status=0
    add_failure "$f: two-way/force-density contamination pattern found"
  fi
  if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_stage8_oneway_forcing\b|^[[:space:]]*call[[:space:]].*rhs_injection|^[[:space:]]*call[[:space:]].*fibre.*force' "$af"; then
    no_fibre_force_rhs_status=0
    add_failure "$f: fibre-force contamination pattern found"
  fi
done

# RHS assignment contamination (active code only) in navier/time_integrators
for f in navier.f90 time_integrators.f90; do
  af="$tmpdir/$f.active"
  if grep -Eiq '=[^\n]*(f_fsi|f_ibm|f_feedback|force_density|fibre_force|twoway_force|rhs_fibre|rhs_ibm)\b' "$af"; then
    no_rhs_assignment_contamination_status=0
    add_failure "$f: suspicious RHS assignment contamination token found"
  fi
done

# Poisson / projection / RK3 protections
if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_.*(force|ibm)|\bstage10_hook_\b' "$tmpdir/poisson.f90.active"; then
  no_poisson_contamination_status=0
  add_failure "poisson.f90: contamination token found"
fi
if grep -Eiq '\bstage10_hook_\b|^[[:space:]]*use[[:space:]]+fibre_.*(force|ibm)' "$tmpdir/derive.f90.active"; then
  no_projection_contamination_status=0
  add_failure "derive.f90: contamination token found"
fi
if grep -Eiq '^[[:space:]]*use[[:space:]]+fibre_.*(force|ibm)|^[[:space:]]*call[[:space:]].*(twoway|feedback|spreading|force_density|rhs_injection)' "$tmpdir/time_integrators.f90.active"; then
  no_rk3_contamination_status=0
  add_failure "time_integrators.f90: contamination token found"
fi

for v in $no_stage10_hook_in_rhs_files_status $xcompact3d_guarded_hook_status $no_ibm_spreading_rhs_status $no_feedback_rhs_status $no_twoway_force_density_rhs_status $no_fibre_force_rhs_status $no_rhs_assignment_contamination_status $no_poisson_contamination_status $no_projection_contamination_status $no_rk3_contamination_status; do
  [ "$v" -eq 1 ] || rhs_equivalence_status=0
done

final_status=1
for v in $build_status $rhs_equivalence_status; do
  [ "$v" -eq 1 ] || final_status=0
done

cat > stage10_outputs/stage10_6_rhs_contamination_audit.dat <<DAT
stage10_6_requested_flag 1
stage10_6_build_status $build_status
stage10_6_no_stage10_hook_in_rhs_files_status $no_stage10_hook_in_rhs_files_status
stage10_6_xcompact3d_guarded_hook_status $xcompact3d_guarded_hook_status
stage10_6_no_ibm_spreading_rhs_status $no_ibm_spreading_rhs_status
stage10_6_no_feedback_rhs_status $no_feedback_rhs_status
stage10_6_no_twoway_force_density_rhs_status $no_twoway_force_density_rhs_status
stage10_6_no_fibre_force_rhs_status $no_fibre_force_rhs_status
stage10_6_no_rhs_assignment_contamination_status $no_rhs_assignment_contamination_status
stage10_6_no_poisson_contamination_status $no_poisson_contamination_status
stage10_6_no_projection_contamination_status $no_projection_contamination_status
stage10_6_no_rk3_contamination_status $no_rk3_contamination_status
stage10_6_rhs_equivalence_status $rhs_equivalence_status
stage10_6_rhs_contamination_audit_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 10.6 RHS CONTAMINATION AUDIT VERDICT: PASS"
  echo "STAGE 10.6 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 10.6 RHS CONTAMINATION AUDIT VERDICT: FAIL"
echo "STAGE 10.6 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
