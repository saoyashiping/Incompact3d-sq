#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
LOG_FILE="$OUTPUT_DIR/stage13_5_conservation_sign_audit.log"
DAT_FILE="$OUTPUT_DIR/fibre_stage13_5_conservation_sign_audit.dat"
GATE_FILE="$OUTPUT_DIR/stage13_5_conservation_sign_audit_gate.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_5_conservation_sign_audit_reasons.tmp"

mkdir -p "$OUTPUT_DIR"
: > "$REASONS_FILE"

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

add_reason() {
    echo "$1" >> "$REASONS_FILE"
}

get_dat_value() {
    key=$1
    awk -v key="$key" '$1 == key { value=$2 } END { if (value != "") print value }' "$DAT_FILE"
}

check_status_key() {
    key=$1
    value=$(get_dat_value "$key")
    if [ "$value" = "1" ]; then
        return 0
    fi
    add_reason "$key expected 1 but found ${value:-missing}"
    return 1
}

check_scalar_le() {
    key=$1
    limit=$2
    value=$(get_dat_value "$key")
    if [ -z "$value" ]; then
        add_reason "$key missing"
        return 1
    fi
    awk -v value="$value" -v limit="$limit" 'BEGIN { exit !(value <= limit) }'
    if [ $? -eq 0 ]; then
        return 0
    fi
    add_reason "$key expected <= $limit but found $value"
    return 1
}

check_scalar_ge() {
    key=$1
    limit=$2
    value=$(get_dat_value "$key")
    if [ -z "$value" ]; then
        add_reason "$key missing"
        return 1
    fi
    awk -v value="$value" -v limit="$limit" 'BEGIN { exit !(value >= limit) }'
    if [ $? -eq 0 ]; then
        return 0
    fi
    add_reason "$key expected >= $limit but found $value"
    return 1
}

write_gate_file() {
    gate_status=$1
    {
        echo "stage13_5_gate_requested_flag $requested_flag"
        echo "stage13_5_gate_build_status $build_status"
        echo "stage13_5_gate_conservation_sign_audit_check_status $check_status"
        echo "stage13_5_gate_fluid_to_fibre_sign_status $fluid_to_fibre_sign_status"
        echo "stage13_5_gate_fibre_to_fluid_sign_status $fibre_to_fluid_sign_status"
        echo "stage13_5_gate_action_reaction_status $action_reaction_status"
        echo "stage13_5_gate_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_5_gate_integrated_force_fibre_to_fluid_status $integrated_force_fibre_to_fluid_status"
        echo "stage13_5_gate_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_5_gate_component_sign_conservation_status $component_sign_conservation_status"
        echo "stage13_5_gate_multipoint_sign_conservation_status $multipoint_sign_conservation_status"
        echo "stage13_5_gate_nonuniform_volume_sign_conservation_status $nonuniform_volume_sign_conservation_status"
        echo "stage13_5_gate_boundary_sign_conservation_status $boundary_sign_conservation_status"
        echo "stage13_5_gate_finite_force_density_status $finite_force_density_status"
        echo "stage13_5_gate_clear_status $clear_status"
        echo "stage13_5_gate_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_5_gate_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage13_5_gate_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_5_gate_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_5_gate_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_5_gate_no_fluid_field_access_status $no_fluid_field_access_status"
        echo "stage13_5_gate_no_fluid_field_modification_status $no_fluid_field_modification_status"
        echo "stage13_5_gate_status $gate_status"
    } > "$GATE_FILE"
}

requested_flag=1
build_status=1
check_status=0
fluid_to_fibre_sign_status=0
fibre_to_fluid_sign_status=0
action_reaction_status=0
spreading_input_sign_status=0
integrated_force_fibre_to_fluid_status=0
wrong_sign_rejection_status=0
component_sign_conservation_status=0
multipoint_sign_conservation_status=0
nonuniform_volume_sign_conservation_status=0
boundary_sign_conservation_status=0
finite_force_density_status=0
clear_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_fluid_field_access_status=0
no_fluid_field_modification_status=0

STAGE13_5_RUN_STAGE13_4=${STAGE13_5_RUN_STAGE13_4:-0}
if [ "$STAGE13_5_RUN_STAGE13_4" = "1" ]; then
    if ! sh stage13_checks/run_stage13_4_volume_normalization_audit.sh; then
        add_reason "optional Stage 13.4 prerequisite failed"
    fi
fi

ensure_build_dir

targets="xcompact3d \
fibre_stage11_config_check \
fibre_stage11_lagrangian_state_check \
fibre_stage11_grid_metadata_check \
fibre_stage11_oneway_interpolation_check \
fibre_stage11_controlled_interpolation_check \
fibre_stage11_production_oneway_hook_check \
fibre_stage12_config_check \
fibre_stage12_force_buffer_check \
fibre_stage12_prescribed_velocity_check \
fibre_stage12_feedback_formula_check \
fibre_stage12_sign_convention_audit_check \
fibre_stage12_power_diagnostics_check \
fibre_stage12_production_feedback_candidate_check \
fibre_stage13_config_check \
fibre_stage13_force_density_buffer_check \
fibre_stage13_spreading_kernel_check \
fibre_stage13_volume_normalization_audit_check \
fibre_stage13_conservation_sign_audit_check"

for target in $targets; do
    if ! cmake --build "$BUILD_DIR" --target "$target" -j; then
        build_status=0
        add_reason "build failed for $target"
    fi
done

if [ -x "$BUILD_DIR/bin/fibre_stage13_conservation_sign_audit_check" ]; then
    EXE="$BUILD_DIR/bin/fibre_stage13_conservation_sign_audit_check"
elif [ -x "$BUILD_DIR/src/fibre_stage13_conservation_sign_audit_check" ]; then
    EXE="$BUILD_DIR/src/fibre_stage13_conservation_sign_audit_check"
else
    EXE=""
    add_reason "fibre_stage13_conservation_sign_audit_check executable not found"
fi

if [ -n "$EXE" ]; then
    X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1 \
    X3D_STAGE13_FORCE_READONLY=1 \
    X3D_STAGE13_SPREADING_READONLY=1 \
    X3D_STAGE13_MAX_POINTS=8 \
    X3D_STAGE13_MAX_EULERIAN_POINTS=64 \
    X3D_STAGE13_SPREADING_NORMALIZATION=conservative \
        "$EXE" > "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        add_reason "fibre_stage13_conservation_sign_audit_check execution failed"
    fi
else
    : > "$LOG_FILE"
fi

if grep -q "STAGE 13.5 CONSERVATION SIGN AUDIT VERDICT: PASS" "$LOG_FILE"; then
    check_status=1
else
    add_reason "PASS verdict missing from Stage 13.5 log"
fi

if [ ! -f "$DAT_FILE" ]; then
    add_reason "$DAT_FILE missing"
else
    check_status_key stage13_5_requested_flag >/dev/null && requested_flag=1
    check_status_key stage13_5_readonly_mode_status >/dev/null || true
    check_status_key stage13_5_spreading_readonly_status >/dev/null || true
    check_status_key stage13_5_initialized_status >/dev/null || true
    check_status_key stage13_5_fluid_to_fibre_sign_status >/dev/null && fluid_to_fibre_sign_status=1
    check_status_key stage13_5_fibre_to_fluid_sign_status >/dev/null && fibre_to_fluid_sign_status=1
    check_status_key stage13_5_action_reaction_status >/dev/null && action_reaction_status=1
    check_status_key stage13_5_spreading_input_sign_status >/dev/null && spreading_input_sign_status=1
    check_status_key stage13_5_integrated_force_fibre_to_fluid_status >/dev/null && integrated_force_fibre_to_fluid_status=1
    check_status_key stage13_5_wrong_sign_rejection_status >/dev/null && wrong_sign_rejection_status=1
    check_status_key stage13_5_component_sign_conservation_status >/dev/null && component_sign_conservation_status=1
    check_status_key stage13_5_multipoint_sign_conservation_status >/dev/null && multipoint_sign_conservation_status=1
    check_status_key stage13_5_nonuniform_volume_sign_conservation_status >/dev/null && nonuniform_volume_sign_conservation_status=1
    check_status_key stage13_5_boundary_sign_conservation_status >/dev/null && boundary_sign_conservation_status=1
    check_status_key stage13_5_finite_force_density_status >/dev/null && finite_force_density_status=1
    check_status_key stage13_5_clear_status >/dev/null && clear_status=1
    check_status_key stage13_5_no_rhs_injection_status >/dev/null && no_rhs_injection_status=1
    check_status_key stage13_5_no_ibm_spreading_status >/dev/null && no_ibm_spreading_status=1
    check_status_key stage13_5_no_feedback_application_status >/dev/null && no_feedback_application_status=1
    check_status_key stage13_5_no_twoway_force_status >/dev/null && no_twoway_force_status=1
    check_status_key stage13_5_no_structure_advance_status >/dev/null && no_structure_advance_status=1
    check_status_key stage13_5_no_fluid_field_access_status >/dev/null && no_fluid_field_access_status=1
    check_status_key stage13_5_no_fluid_field_modification_status >/dev/null && no_fluid_field_modification_status=1
    check_status_key stage13_5_conservation_sign_audit_status >/dev/null || true

    check_scalar_le stage13_5_action_reaction_max_error 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_integrated_force_fibre_to_fluid_error_l2 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_component_error_x 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_component_error_y 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_component_error_z 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_nonuniform_sign_error_l2 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_boundary_sign_error_l2 1.0e-12 >/dev/null || true
    check_scalar_le stage13_5_max_abs_force_density_norm_after_clear 1.0e-12 >/dev/null || true
    check_scalar_ge stage13_5_wrong_sign_separation 1.0e-8 >/dev/null || true
fi

if [ ! -s "$REASONS_FILE" ] && [ "$build_status" = "1" ] && [ "$check_status" = "1" ]; then
    write_gate_file 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.5 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.5 gate failed without a captured reason"
fi
write_gate_file 0
echo "STAGE 13.5 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
