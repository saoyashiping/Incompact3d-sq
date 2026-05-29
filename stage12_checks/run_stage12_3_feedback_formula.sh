#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE12_3_RUN_STAGE12_2=${STAGE12_3_RUN_STAGE12_2:-0}

mkdir -p stage12_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE12_3_RUN_STAGE12_2" = "1" ]; then
    bash stage12_checks/run_stage12_2_prescribed_velocity.sh
fi

build_status=1
feedback_formula_check_status=0
requested_flag=0
zero_slip_status=0
positive_slip_sign_status=0
negative_slip_sign_status=0
force_sign_minus_status=0
multicomponent_slip_status=0
gain_scaling_status=0
finite_force_status=0
force_norm_finite_status=0
no_eulerian_force_density_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_fluid_field_access_status=0
no_fluid_field_modification_status=0
final_status=1
reasons=""

add_reason() {
    if [ -z "$reasons" ]; then
        reasons="$1"
    else
        reasons="$reasons,$1"
    fi
}

for tgt in xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check; do
    cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done

exe="$BUILD_DIR/bin/fibre_stage12_feedback_formula_check"
if [ ! -x "$exe" ]; then
    exe="$BUILD_DIR/src/fibre_stage12_feedback_formula_check"
fi

log_file=stage12_outputs/stage12_3_feedback_formula.log
if [ "$build_status" -eq 1 ] && [ -x "$exe" ]; then
    X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
    X3D_STAGE12_FORCE_READONLY=1 \
    X3D_STAGE12_FEEDBACK_GAIN=1.0 \
    X3D_STAGE12_FORCE_SIGN=1 \
    X3D_STAGE12_MAX_POINTS=8 \
    "$exe" > "$log_file" 2>&1 || feedback_formula_check_status=0
    grep -q "STAGE 12.3 FEEDBACK FORMULA VERDICT: PASS" "$log_file" && feedback_formula_check_status=1 || feedback_formula_check_status=0
else
    feedback_formula_check_status=0
    add_reason "missing_or_unbuilt_fibre_stage12_feedback_formula_check"
fi

dat=stage12_outputs/fibre_stage12_3_feedback_formula.dat
get_val(){ awk -v k="$1" '$1==k{print $2}' "$2"; }
check_error_le(){ awk -v v="$1" -v limit="$2" 'BEGIN{if ((v+0) <= (limit+0)) exit 0; exit 1}'; }

if [ -f "$dat" ]; then
    for k in \
      stage12_3_requested_flag \
      stage12_3_readonly_mode_status \
      stage12_3_initialized_status \
      stage12_3_zero_slip_status \
      stage12_3_positive_slip_sign_status \
      stage12_3_negative_slip_sign_status \
      stage12_3_force_sign_minus_status \
      stage12_3_multicomponent_slip_status \
      stage12_3_gain_scaling_status \
      stage12_3_finite_force_status \
      stage12_3_force_norm_finite_status \
      stage12_3_no_eulerian_force_density_status \
      stage12_3_no_rhs_injection_status \
      stage12_3_no_ibm_spreading_status \
      stage12_3_no_feedback_application_status \
      stage12_3_no_twoway_force_status \
      stage12_3_no_structure_advance_status \
      stage12_3_no_fluid_field_access_status \
      stage12_3_no_fluid_field_modification_status \
      stage12_3_feedback_formula_status; do
        v=$(get_val "$k" "$dat")
        if [ "$v" != "1" ]; then
            final_status=0
            add_reason "$k"
        fi
    done

    for k in \
      stage12_3_zero_slip_max_error \
      stage12_3_positive_slip_max_error \
      stage12_3_negative_slip_max_error \
      stage12_3_force_sign_minus_max_error \
      stage12_3_multicomponent_max_error \
      stage12_3_gain_scaling_max_error \
      stage12_3_force_norm_max_error; do
        v=$(get_val "$k" "$dat")
        if [ -z "$v" ] || ! check_error_le "$v" "1.0e-12"; then
            final_status=0
            add_reason "$k"
        fi
    done

    requested_flag=$(get_val stage12_3_requested_flag "$dat")
    zero_slip_status=$(get_val stage12_3_zero_slip_status "$dat")
    positive_slip_sign_status=$(get_val stage12_3_positive_slip_sign_status "$dat")
    negative_slip_sign_status=$(get_val stage12_3_negative_slip_sign_status "$dat")
    force_sign_minus_status=$(get_val stage12_3_force_sign_minus_status "$dat")
    multicomponent_slip_status=$(get_val stage12_3_multicomponent_slip_status "$dat")
    gain_scaling_status=$(get_val stage12_3_gain_scaling_status "$dat")
    finite_force_status=$(get_val stage12_3_finite_force_status "$dat")
    force_norm_finite_status=$(get_val stage12_3_force_norm_finite_status "$dat")
    no_eulerian_force_density_status=$(get_val stage12_3_no_eulerian_force_density_status "$dat")
    no_rhs_injection_status=$(get_val stage12_3_no_rhs_injection_status "$dat")
    no_ibm_spreading_status=$(get_val stage12_3_no_ibm_spreading_status "$dat")
    no_feedback_application_status=$(get_val stage12_3_no_feedback_application_status "$dat")
    no_twoway_force_status=$(get_val stage12_3_no_twoway_force_status "$dat")
    no_structure_advance_status=$(get_val stage12_3_no_structure_advance_status "$dat")
    no_fluid_field_access_status=$(get_val stage12_3_no_fluid_field_access_status "$dat")
    no_fluid_field_modification_status=$(get_val stage12_3_no_fluid_field_modification_status "$dat")
else
    final_status=0
    add_reason "missing_fibre_stage12_3_feedback_formula_dat"
fi

if [ "$build_status" -ne 1 ] || [ "$feedback_formula_check_status" -ne 1 ]; then
    final_status=0
    add_reason "build_or_feedback_formula_check_failed"
fi

cat > stage12_outputs/stage12_3_feedback_formula_gate.dat <<EOD
stage12_3_gate_requested_flag $requested_flag
stage12_3_gate_build_status $build_status
stage12_3_gate_feedback_formula_check_status $feedback_formula_check_status
stage12_3_gate_zero_slip_status $zero_slip_status
stage12_3_gate_positive_slip_sign_status $positive_slip_sign_status
stage12_3_gate_negative_slip_sign_status $negative_slip_sign_status
stage12_3_gate_force_sign_minus_status $force_sign_minus_status
stage12_3_gate_multicomponent_slip_status $multicomponent_slip_status
stage12_3_gate_gain_scaling_status $gain_scaling_status
stage12_3_gate_finite_force_status $finite_force_status
stage12_3_gate_force_norm_finite_status $force_norm_finite_status
stage12_3_gate_no_eulerian_force_density_status $no_eulerian_force_density_status
stage12_3_gate_no_rhs_injection_status $no_rhs_injection_status
stage12_3_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage12_3_gate_no_feedback_application_status $no_feedback_application_status
stage12_3_gate_no_twoway_force_status $no_twoway_force_status
stage12_3_gate_no_structure_advance_status $no_structure_advance_status
stage12_3_gate_no_fluid_field_access_status $no_fluid_field_access_status
stage12_3_gate_no_fluid_field_modification_status $no_fluid_field_modification_status
stage12_3_gate_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
    echo "STAGE 12.3 FINAL VERDICT: PASS"
else
    [ -z "$reasons" ] && reasons="unknown_failure"
    echo "STAGE 12.3 FINAL VERDICT: FAIL"
    echo "Reasons:$reasons"
    exit 1
fi
