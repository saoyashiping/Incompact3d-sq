#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_1_RUN_STAGE11_0=${STAGE11_1_RUN_STAGE11_0:-0}
mkdir -p stage11_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_1_RUN_STAGE11_0" = "1" ]; then
    sh stage11_checks/run_stage11_0_config.sh
fi

build_status=1
state_check_status=1
gate_status=1
reasons=""

for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check; do
    if ! cmake --build "$BUILD_DIR" --target "$tgt" -j; then
        build_status=0
    fi
done

LOG_FILE=stage11_outputs/stage11_1_lagrangian_state.log
EXE_BIN="$BUILD_DIR/bin/fibre_stage11_lagrangian_state_check"
EXE_SRC="$BUILD_DIR/src/fibre_stage11_lagrangian_state_check"
EXE_PATH="$EXE_BIN"
if [ ! -x "$EXE_PATH" ]; then
    EXE_PATH="$EXE_SRC"
fi

if [ "$build_status" -eq 1 ]; then
    if ! X3D_STAGE11_ONEWAY_HOOK=1 \
         X3D_STAGE11_FORCE_READONLY=1 \
         X3D_STAGE11_MAX_POINTS=8 \
         X3D_STAGE11_MAX_STEPS=3 \
         "$EXE_PATH" > "$LOG_FILE" 2>&1; then
        state_check_status=0
    fi
else
    state_check_status=0
fi

if [ "$state_check_status" -eq 1 ]; then
    if ! grep -q "STAGE 11.1 LAGRANGIAN STATE VERDICT: PASS" "$LOG_FILE"; then
        state_check_status=0
    fi
fi

requested_flag=0
allocated_status=0
point_count_status=0
coordinates_finite_status=0
velocity_placeholder_status=0
valid_flag_status=0
no_fluid_field_access_status=0
no_velocity_sampling_status=0
no_fluid_field_modification_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_force_status=0
no_twoway_force_status=0
no_structure_advance_status=0
lagrangian_state_status=0

DAT_FILE=stage11_outputs/fibre_stage11_1_lagrangian_state.dat
if [ "$state_check_status" -eq 1 ] && [ -f "$DAT_FILE" ]; then
    get_val() {
        awk -v k="$1" '$1==k {print $2}' "$DAT_FILE"
    }
    requested_flag=$(get_val stage11_1_requested_flag)
    readonly_mode_status=$(get_val stage11_1_readonly_mode_status)
    allocated_status=$(get_val stage11_1_allocated_status)
    point_count_status=$(get_val stage11_1_point_count_status)
    coordinates_finite_status=$(get_val stage11_1_coordinates_finite_status)
    velocity_placeholder_status=$(get_val stage11_1_velocity_placeholder_status)
    valid_flag_status=$(get_val stage11_1_valid_flag_status)
    no_fluid_field_access_status=$(get_val stage11_1_no_fluid_field_access_status)
    no_velocity_sampling_status=$(get_val stage11_1_no_velocity_sampling_status)
    no_fluid_field_modification_status=$(get_val stage11_1_no_fluid_field_modification_status)
    no_rhs_injection_status=$(get_val stage11_1_no_rhs_injection_status)
    no_ibm_spreading_status=$(get_val stage11_1_no_ibm_spreading_status)
    no_feedback_force_status=$(get_val stage11_1_no_feedback_force_status)
    no_twoway_force_status=$(get_val stage11_1_no_twoway_force_status)
    no_structure_advance_status=$(get_val stage11_1_no_structure_advance_status)
    lagrangian_state_status=$(get_val stage11_1_lagrangian_state_status)
    if [ "$requested_flag" != "1" ] || [ "$readonly_mode_status" != "1" ] || [ "$allocated_status" != "1" ] || \
       [ "$point_count_status" != "1" ] || [ "$coordinates_finite_status" != "1" ] || \
       [ "$velocity_placeholder_status" != "1" ] || [ "$valid_flag_status" != "1" ] || \
       [ "$no_fluid_field_access_status" != "1" ] || [ "$no_velocity_sampling_status" != "1" ] || \
       [ "$no_fluid_field_modification_status" != "1" ] || [ "$no_rhs_injection_status" != "1" ] || \
       [ "$no_ibm_spreading_status" != "1" ] || [ "$no_feedback_force_status" != "1" ] || \
       [ "$no_twoway_force_status" != "1" ] || [ "$no_structure_advance_status" != "1" ] || \
       [ "$lagrangian_state_status" != "1" ]; then
        gate_status=0
    fi
else
    gate_status=0
fi

if [ "$build_status" -ne 1 ]; then reasons="$reasons build_failed;"; gate_status=0; fi
if [ "$state_check_status" -ne 1 ]; then reasons="$reasons state_check_failed;"; gate_status=0; fi

cat > stage11_outputs/stage11_1_lagrangian_state_gate.dat <<EOD
stage11_1_gate_requested_flag $requested_flag
stage11_1_gate_build_status $build_status
stage11_1_gate_state_check_status $state_check_status
stage11_1_gate_allocated_status $allocated_status
stage11_1_gate_point_count_status $point_count_status
stage11_1_gate_coordinates_finite_status $coordinates_finite_status
stage11_1_gate_velocity_placeholder_status $velocity_placeholder_status
stage11_1_gate_valid_flag_status $valid_flag_status
stage11_1_gate_no_fluid_field_access_status $no_fluid_field_access_status
stage11_1_gate_no_velocity_sampling_status $no_velocity_sampling_status
stage11_1_gate_no_fluid_field_modification_status $no_fluid_field_modification_status
stage11_1_gate_no_rhs_injection_status $no_rhs_injection_status
stage11_1_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage11_1_gate_no_feedback_force_status $no_feedback_force_status
stage11_1_gate_no_twoway_force_status $no_twoway_force_status
stage11_1_gate_no_structure_advance_status $no_structure_advance_status
stage11_1_gate_status $gate_status
EOD

if [ "$gate_status" -eq 1 ]; then
    echo "STAGE 11.1 FINAL VERDICT: PASS"
else
    echo "STAGE 11.1 FINAL VERDICT: FAIL"
    echo "Reasons:$reasons"
    exit 1
fi
