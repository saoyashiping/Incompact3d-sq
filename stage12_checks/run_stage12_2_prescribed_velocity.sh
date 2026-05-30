#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE12_2_RUN_STAGE12_1=${STAGE12_2_RUN_STAGE12_1:-0}

mkdir -p stage12_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE12_2_RUN_STAGE12_1" = "1" ]; then
    bash stage12_checks/run_stage12_1_force_buffer.sh
fi

build_status=1
prescribed_velocity_check_status=0
requested_flag=0
allocated_status=0
point_count_status=0
zero_velocity_status=0
constant_velocity_status=0
velocity_norm_finite_status=0
velocity_valid_flag_status=0
clear_status=0
no_force_computation_status=0
no_eulerian_force_density_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_fluid_field_access_status=0
no_fluid_field_modification_status=0
final_status=1
reasons="init"

for tgt in xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check; do
    cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done

exe="$BUILD_DIR/bin/fibre_stage12_prescribed_velocity_check"
if [ ! -x "$exe" ]; then
    exe="$BUILD_DIR/src/fibre_stage12_prescribed_velocity_check"
fi

log_file=stage12_outputs/stage12_2_prescribed_velocity.log
if [ "$build_status" -eq 1 ] && [ -x "$exe" ]; then
    X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
    X3D_STAGE12_FORCE_READONLY=1 \
    X3D_STAGE12_FEEDBACK_GAIN=1.0 \
    X3D_STAGE12_FORCE_SIGN=1 \
    X3D_STAGE12_MAX_POINTS=8 \
    "$exe" > "$log_file" 2>&1 || prescribed_velocity_check_status=0
    grep -q "STAGE 12.2 PRESCRIBED VELOCITY VERDICT: PASS" "$log_file" && prescribed_velocity_check_status=1 || prescribed_velocity_check_status=0
else
    prescribed_velocity_check_status=0
    [ "$reasons" = "init" ] && reasons="missing_or_unbuilt_fibre_stage12_prescribed_velocity_check"
fi

dat=stage12_outputs/fibre_stage12_2_prescribed_velocity.dat
get_val(){ awk -v k="$1" '$1==k{print $2}' "$2"; }

if [ -f "$dat" ]; then
    for k in \
      stage12_2_requested_flag \
      stage12_2_readonly_mode_status \
      stage12_2_allocated_status \
      stage12_2_point_count_status \
      stage12_2_zero_velocity_status \
      stage12_2_constant_velocity_status \
      stage12_2_velocity_norm_finite_status \
      stage12_2_velocity_valid_flag_status \
      stage12_2_clear_status \
      stage12_2_no_force_computation_status \
      stage12_2_no_eulerian_force_density_status \
      stage12_2_no_rhs_injection_status \
      stage12_2_no_ibm_spreading_status \
      stage12_2_no_feedback_application_status \
      stage12_2_no_twoway_force_status \
      stage12_2_no_structure_advance_status \
      stage12_2_no_fluid_field_access_status \
      stage12_2_no_fluid_field_modification_status \
      stage12_2_prescribed_velocity_status; do
        v=$(get_val "$k" "$dat")
        if [ "$v" != "1" ]; then
            final_status=0
            [ "$reasons" = "init" ] && reasons="$k"
        fi
    done

    requested_flag=$(get_val stage12_2_requested_flag "$dat")
    allocated_status=$(get_val stage12_2_allocated_status "$dat")
    point_count_status=$(get_val stage12_2_point_count_status "$dat")
    zero_velocity_status=$(get_val stage12_2_zero_velocity_status "$dat")
    constant_velocity_status=$(get_val stage12_2_constant_velocity_status "$dat")
    velocity_norm_finite_status=$(get_val stage12_2_velocity_norm_finite_status "$dat")
    velocity_valid_flag_status=$(get_val stage12_2_velocity_valid_flag_status "$dat")
    clear_status=$(get_val stage12_2_clear_status "$dat")
    no_force_computation_status=$(get_val stage12_2_no_force_computation_status "$dat")
    no_eulerian_force_density_status=$(get_val stage12_2_no_eulerian_force_density_status "$dat")
    no_rhs_injection_status=$(get_val stage12_2_no_rhs_injection_status "$dat")
    no_ibm_spreading_status=$(get_val stage12_2_no_ibm_spreading_status "$dat")
    no_feedback_application_status=$(get_val stage12_2_no_feedback_application_status "$dat")
    no_twoway_force_status=$(get_val stage12_2_no_twoway_force_status "$dat")
    no_structure_advance_status=$(get_val stage12_2_no_structure_advance_status "$dat")
    no_fluid_field_access_status=$(get_val stage12_2_no_fluid_field_access_status "$dat")
    no_fluid_field_modification_status=$(get_val stage12_2_no_fluid_field_modification_status "$dat")
else
    final_status=0
    [ "$reasons" = "init" ] && reasons="missing_fibre_stage12_2_prescribed_velocity_dat"
fi

if [ "$build_status" -ne 1 ] || [ "$prescribed_velocity_check_status" -ne 1 ]; then
    final_status=0
    [ "$reasons" = "init" ] && reasons="build_or_prescribed_velocity_check_failed"
fi

cat > stage12_outputs/stage12_2_prescribed_velocity_gate.dat <<EOD
stage12_2_gate_requested_flag $requested_flag
stage12_2_gate_build_status $build_status
stage12_2_gate_prescribed_velocity_check_status $prescribed_velocity_check_status
stage12_2_gate_allocated_status $allocated_status
stage12_2_gate_point_count_status $point_count_status
stage12_2_gate_zero_velocity_status $zero_velocity_status
stage12_2_gate_constant_velocity_status $constant_velocity_status
stage12_2_gate_velocity_norm_finite_status $velocity_norm_finite_status
stage12_2_gate_velocity_valid_flag_status $velocity_valid_flag_status
stage12_2_gate_clear_status $clear_status
stage12_2_gate_no_force_computation_status $no_force_computation_status
stage12_2_gate_no_eulerian_force_density_status $no_eulerian_force_density_status
stage12_2_gate_no_rhs_injection_status $no_rhs_injection_status
stage12_2_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage12_2_gate_no_feedback_application_status $no_feedback_application_status
stage12_2_gate_no_twoway_force_status $no_twoway_force_status
stage12_2_gate_no_structure_advance_status $no_structure_advance_status
stage12_2_gate_no_fluid_field_access_status $no_fluid_field_access_status
stage12_2_gate_no_fluid_field_modification_status $no_fluid_field_modification_status
stage12_2_gate_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
    echo "STAGE 12.2 FINAL VERDICT: PASS"
else
    [ "$reasons" = "init" ] && reasons="unknown_failure"
    echo "STAGE 12.2 FINAL VERDICT: FAIL"
    echo "Reasons:$reasons"
    exit 1
fi
