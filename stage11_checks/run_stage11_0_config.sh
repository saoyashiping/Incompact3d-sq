#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE11_0_RUN_STAGE10_CLOSURE=${STAGE11_0_RUN_STAGE10_CLOSURE:-0}

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

mkdir -p stage11_outputs

failures=""
add_failure(){ if [ -z "$failures" ]; then failures="$1"; else failures="$failures\n$1"; fi; }

build_status=1
config_check_status=1
no_lagrangian_state_status=1
no_velocity_sampling_status=1
no_fluid_field_modification_status=1
no_rhs_injection_status=1
no_ibm_spreading_status=1
no_feedback_force_status=1
no_twoway_force_status=1
no_structure_advance_status=1

build_target(){
  ensure_build_dir
  if ! DECOMP2D_ROOT="$DECOMP2D_ROOT" cmake --build "$BUILD_DIR" --target "$1" -j; then
    build_status=0
    add_failure "build failed: $1"
  fi
}

build_target xcompact3d
build_target fibre_stage10_config_check
build_target fibre_stage10_noop_hook_check
build_target fibre_stage11_config_check

if [ "$STAGE11_0_RUN_STAGE10_CLOSURE" = "1" ]; then
  DECOMP2D_ROOT="$DECOMP2D_ROOT" BUILD_DIR="$BUILD_DIR" bash stage10_checks/run_stage10_8_total_smoke.sh || add_failure "optional Stage10 closure failed"
fi

log="stage11_outputs/stage11_0_config.log"
if ! env X3D_STAGE11_ONEWAY_HOOK=1 X3D_STAGE11_FORCE_READONLY=1 X3D_STAGE11_MAX_POINTS=8 X3D_STAGE11_MAX_STEPS=3 "$BUILD_DIR/bin/fibre_stage11_config_check" >"$log" 2>&1; then
  config_check_status=0
  add_failure "run fibre_stage11_config_check failed"
fi
if ! grep -Eq 'STAGE 11\.0 CONFIG VERDICT:[[:space:]]*PASS' "$log"; then
  config_check_status=0
  add_failure "missing STAGE 11.0 CONFIG VERDICT: PASS"
fi

dat="stage11_outputs/fibre_stage11_0_config.dat"
getv(){ awk -v k="$1" '$1==k{print $2}' "$dat" 2>/dev/null | tail -n1; }

for k in stage11_0_requested_flag stage11_0_readonly_mode_status stage11_0_no_lagrangian_state_status stage11_0_no_velocity_sampling_status stage11_0_no_fluid_field_access_status stage11_0_no_fluid_field_modification_status stage11_0_no_rhs_injection_status stage11_0_no_ibm_spreading_status stage11_0_no_feedback_force_status stage11_0_no_twoway_force_status stage11_0_no_structure_advance_status stage11_0_config_status; do
  v=$(getv "$k")
  if [ "$v" != "1" ]; then
    config_check_status=0
    add_failure "dat key not 1: $k"
  fi
done

[ "$(getv stage11_0_no_lagrangian_state_status)" = "1" ] || no_lagrangian_state_status=0
[ "$(getv stage11_0_no_velocity_sampling_status)" = "1" ] || no_velocity_sampling_status=0
[ "$(getv stage11_0_no_fluid_field_modification_status)" = "1" ] || no_fluid_field_modification_status=0
[ "$(getv stage11_0_no_rhs_injection_status)" = "1" ] || no_rhs_injection_status=0
[ "$(getv stage11_0_no_ibm_spreading_status)" = "1" ] || no_ibm_spreading_status=0
[ "$(getv stage11_0_no_feedback_force_status)" = "1" ] || no_feedback_force_status=0
[ "$(getv stage11_0_no_twoway_force_status)" = "1" ] || no_twoway_force_status=0
[ "$(getv stage11_0_no_structure_advance_status)" = "1" ] || no_structure_advance_status=0

final_status=1
for v in $build_status $config_check_status $no_lagrangian_state_status $no_velocity_sampling_status $no_fluid_field_modification_status $no_rhs_injection_status $no_ibm_spreading_status $no_feedback_force_status $no_twoway_force_status $no_structure_advance_status; do
  [ "$v" -eq 1 ] || final_status=0
done

cat > stage11_outputs/stage11_0_config_gate.dat <<DAT
stage11_0_gate_requested_flag 1
stage11_0_gate_build_status $build_status
stage11_0_gate_config_check_status $config_check_status
stage11_0_gate_no_lagrangian_state_status $no_lagrangian_state_status
stage11_0_gate_no_velocity_sampling_status $no_velocity_sampling_status
stage11_0_gate_no_fluid_field_modification_status $no_fluid_field_modification_status
stage11_0_gate_no_rhs_injection_status $no_rhs_injection_status
stage11_0_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage11_0_gate_no_feedback_force_status $no_feedback_force_status
stage11_0_gate_no_twoway_force_status $no_twoway_force_status
stage11_0_gate_no_structure_advance_status $no_structure_advance_status
stage11_0_gate_status $final_status
DAT

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 11.0 FINAL VERDICT: PASS"
  exit 0
fi

echo "STAGE 11.0 FINAL VERDICT: FAIL"
printf '%b\n' "$failures"
exit 1
