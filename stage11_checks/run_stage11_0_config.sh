#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_0_RUN_STAGE10_CLOSURE=${STAGE11_0_RUN_STAGE10_CLOSURE:-0}
mkdir -p stage11_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_0_RUN_STAGE10_CLOSURE" = "1" ]; then
    sh stage10_checks/run_stage10_8_total_smoke.sh
fi

build_status=1
config_check_status=1
requested_flag=0
no_lagrangian_state_status=0
no_velocity_sampling_status=0
no_fluid_field_modification_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_force_status=0
no_twoway_force_status=0
no_structure_advance_status=0
gate_status=1
reasons=""

if ! cmake --build "$BUILD_DIR" --target xcompact3d -j; then
    build_status=0
fi
if ! cmake --build "$BUILD_DIR" --target fibre_stage10_config_check -j; then
    build_status=0
fi
if ! cmake --build "$BUILD_DIR" --target fibre_stage10_noop_hook_check -j; then
    build_status=0
fi
if ! cmake --build "$BUILD_DIR" --target fibre_stage11_config_check -j; then
    build_status=0
fi

LOG_FILE=stage11_outputs/stage11_0_config.log
if [ "$build_status" -eq 1 ]; then
    if ! X3D_STAGE11_ONEWAY_HOOK=1 \
         X3D_STAGE11_FORCE_READONLY=1 \
         X3D_STAGE11_MAX_POINTS=8 \
         X3D_STAGE11_MAX_STEPS=3 \
         "$BUILD_DIR/src/fibre_stage11_config_check" > "$LOG_FILE" 2>&1; then
        config_check_status=0
    fi
else
    config_check_status=0
fi

if [ "$config_check_status" -eq 1 ]; then
    if ! grep -q "STAGE 11.0 CONFIG VERDICT: PASS" "$LOG_FILE"; then
        config_check_status=0
    fi
fi

DAT_FILE=stage11_outputs/fibre_stage11_0_config.dat
if [ "$config_check_status" -eq 1 ] && [ -f "$DAT_FILE" ]; then
    get_val() {
        key=$1
        awk -v k="$key" '$1==k {print $2}' "$DAT_FILE"
    }
    requested_flag=$(get_val stage11_0_requested_flag)
    no_lagrangian_state_status=$(get_val stage11_0_no_lagrangian_state_status)
    no_velocity_sampling_status=$(get_val stage11_0_no_velocity_sampling_status)
    no_fluid_field_modification_status=$(get_val stage11_0_no_fluid_field_modification_status)
    no_rhs_injection_status=$(get_val stage11_0_no_rhs_injection_status)
    no_ibm_spreading_status=$(get_val stage11_0_no_ibm_spreading_status)
    no_feedback_force_status=$(get_val stage11_0_no_feedback_force_status)
    no_twoway_force_status=$(get_val stage11_0_no_twoway_force_status)
    no_structure_advance_status=$(get_val stage11_0_no_structure_advance_status)
    if [ "$requested_flag" != "1" ] || [ "$no_lagrangian_state_status" != "1" ] || \
       [ "$no_velocity_sampling_status" != "1" ] || [ "$no_fluid_field_modification_status" != "1" ] || \
       [ "$no_rhs_injection_status" != "1" ] || [ "$no_ibm_spreading_status" != "1" ] || \
       [ "$no_feedback_force_status" != "1" ] || [ "$no_twoway_force_status" != "1" ] || \
       [ "$no_structure_advance_status" != "1" ]; then
        gate_status=0
    fi
else
    gate_status=0
fi

if [ "$build_status" -ne 1 ]; then reasons="$reasons build_failed;"; gate_status=0; fi
if [ "$config_check_status" -ne 1 ]; then reasons="$reasons config_check_failed;"; gate_status=0; fi

cat > stage11_outputs/stage11_0_config_gate.dat <<EOD
stage11_0_gate_requested_flag $requested_flag
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
stage11_0_gate_status $gate_status
EOD

if [ "$gate_status" -eq 1 ]; then
    echo "STAGE 11.0 FINAL VERDICT: PASS"
else
    echo "STAGE 11.0 FINAL VERDICT: FAIL"
    echo "Reasons:$reasons"
    exit 1
fi
