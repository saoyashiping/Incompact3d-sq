#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_2_RUN_STAGE11_1=${STAGE11_2_RUN_STAGE11_1:-0}
mkdir -p stage11_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_2_RUN_STAGE11_1" = "1" ]; then
    sh stage11_checks/run_stage11_1_lagrangian_state.sh
fi

build_status=1
grid_metadata_check_status=1
gate_status=1
reasons=""

for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check; do
    if ! cmake --build "$BUILD_DIR" --target "$tgt" -j; then
        build_status=0
    fi
done

LOG_FILE=stage11_outputs/stage11_2_grid_metadata.log
EXE_PATH="$BUILD_DIR/bin/fibre_stage11_grid_metadata_check"
if [ ! -x "$EXE_PATH" ]; then
    EXE_PATH="$BUILD_DIR/src/fibre_stage11_grid_metadata_check"
fi

if [ "$build_status" -eq 1 ]; then
    if ! X3D_STAGE11_ONEWAY_HOOK=1 \
         X3D_STAGE11_FORCE_READONLY=1 \
         X3D_STAGE11_MAX_POINTS=8 \
         X3D_STAGE11_MAX_STEPS=3 \
         "$EXE_PATH" > "$LOG_FILE" 2>&1; then
        grid_metadata_check_status=0
    fi
else
    grid_metadata_check_status=0
fi

if [ "$grid_metadata_check_status" -eq 1 ]; then
    if ! grep -q "STAGE 11.2 GRID METADATA VERDICT: PASS" "$LOG_FILE"; then
        grid_metadata_check_status=0
    fi
fi

DAT_FILE=stage11_outputs/fibre_stage11_2_grid_metadata.dat
requested_flag=0
lagrangian_state_available_status=0
grid_initialized_status=0
global_size_status=0
local_bounds_status=0
extents_finite_status=0
spacing_positive_status=0
periodic_flags_status=0
staggered_layout_status=0
nonuniform_y_policy_status=0
domain_safety_status=0
no_fluid_field_access_status=0
no_velocity_sampling_status=0
no_fluid_field_modification_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_force_status=0
no_twoway_force_status=0
no_structure_advance_status=0
grid_metadata_status=0

if [ "$grid_metadata_check_status" -eq 1 ] && [ -f "$DAT_FILE" ]; then
    get_val() { awk -v k="$1" '$1==k {print $2}' "$DAT_FILE"; }
    requested_flag=$(get_val stage11_2_requested_flag)
    readonly_mode_status=$(get_val stage11_2_readonly_mode_status)
    lagrangian_state_available_status=$(get_val stage11_2_lagrangian_state_available_status)
    grid_initialized_status=$(get_val stage11_2_grid_initialized_status)
    global_size_status=$(get_val stage11_2_global_size_status)
    local_bounds_status=$(get_val stage11_2_local_bounds_status)
    extents_finite_status=$(get_val stage11_2_extents_finite_status)
    spacing_positive_status=$(get_val stage11_2_spacing_positive_status)
    periodic_flags_status=$(get_val stage11_2_periodic_flags_status)
    staggered_layout_status=$(get_val stage11_2_staggered_layout_status)
    nonuniform_y_policy_status=$(get_val stage11_2_nonuniform_y_policy_status)
    domain_safety_status=$(get_val stage11_2_domain_safety_status)
    no_fluid_field_access_status=$(get_val stage11_2_no_fluid_field_access_status)
    no_velocity_sampling_status=$(get_val stage11_2_no_velocity_sampling_status)
    no_fluid_field_modification_status=$(get_val stage11_2_no_fluid_field_modification_status)
    no_rhs_injection_status=$(get_val stage11_2_no_rhs_injection_status)
    no_ibm_spreading_status=$(get_val stage11_2_no_ibm_spreading_status)
    no_feedback_force_status=$(get_val stage11_2_no_feedback_force_status)
    no_twoway_force_status=$(get_val stage11_2_no_twoway_force_status)
    no_structure_advance_status=$(get_val stage11_2_no_structure_advance_status)
    grid_metadata_status=$(get_val stage11_2_grid_metadata_status)
    if [ "$requested_flag" != "1" ] || [ "$readonly_mode_status" != "1" ] || [ "$lagrangian_state_available_status" != "1" ] || \
       [ "$grid_initialized_status" != "1" ] || [ "$global_size_status" != "1" ] || [ "$local_bounds_status" != "1" ] || \
       [ "$extents_finite_status" != "1" ] || [ "$spacing_positive_status" != "1" ] || [ "$periodic_flags_status" != "1" ] || \
       [ "$staggered_layout_status" != "1" ] || [ "$nonuniform_y_policy_status" != "1" ] || [ "$domain_safety_status" != "1" ] || \
       [ "$no_fluid_field_access_status" != "1" ] || [ "$no_velocity_sampling_status" != "1" ] || [ "$no_fluid_field_modification_status" != "1" ] || \
       [ "$no_rhs_injection_status" != "1" ] || [ "$no_ibm_spreading_status" != "1" ] || [ "$no_feedback_force_status" != "1" ] || \
       [ "$no_twoway_force_status" != "1" ] || [ "$no_structure_advance_status" != "1" ] || [ "$grid_metadata_status" != "1" ]; then
        gate_status=0
    fi
else
    gate_status=0
fi

if [ "$build_status" -ne 1 ]; then reasons="$reasons build_failed;"; gate_status=0; fi
if [ "$grid_metadata_check_status" -ne 1 ]; then reasons="$reasons grid_metadata_check_failed;"; gate_status=0; fi

cat > stage11_outputs/stage11_2_grid_metadata_gate.dat <<EOD
stage11_2_gate_requested_flag $requested_flag
stage11_2_gate_build_status $build_status
stage11_2_gate_grid_metadata_check_status $grid_metadata_check_status
stage11_2_gate_lagrangian_state_available_status $lagrangian_state_available_status
stage11_2_gate_grid_initialized_status $grid_initialized_status
stage11_2_gate_global_size_status $global_size_status
stage11_2_gate_local_bounds_status $local_bounds_status
stage11_2_gate_extents_finite_status $extents_finite_status
stage11_2_gate_spacing_positive_status $spacing_positive_status
stage11_2_gate_periodic_flags_status $periodic_flags_status
stage11_2_gate_staggered_layout_status $staggered_layout_status
stage11_2_gate_nonuniform_y_policy_status $nonuniform_y_policy_status
stage11_2_gate_domain_safety_status $domain_safety_status
stage11_2_gate_no_fluid_field_access_status $no_fluid_field_access_status
stage11_2_gate_no_velocity_sampling_status $no_velocity_sampling_status
stage11_2_gate_no_fluid_field_modification_status $no_fluid_field_modification_status
stage11_2_gate_no_rhs_injection_status $no_rhs_injection_status
stage11_2_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage11_2_gate_no_feedback_force_status $no_feedback_force_status
stage11_2_gate_no_twoway_force_status $no_twoway_force_status
stage11_2_gate_no_structure_advance_status $no_structure_advance_status
stage11_2_gate_status $gate_status
EOD

if [ "$gate_status" -eq 1 ]; then
    echo "STAGE 11.2 FINAL VERDICT: PASS"
else
    echo "STAGE 11.2 FINAL VERDICT: FAIL"
    echo "Reasons:$reasons"
    exit 1
fi
