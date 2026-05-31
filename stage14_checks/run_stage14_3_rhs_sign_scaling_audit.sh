#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_3_RUN_STAGE14_2=${STAGE14_3_RUN_STAGE14_2:-0}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_3_rhs_sign_scaling_audit_gate.dat
REASONS_FILE=$OUTPUT_DIR/stage14_3_rhs_sign_scaling_audit_reasons.tmp
BUILD_STATUS=1
CHECK_STATUS=0

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

get_value() {
    key=$1
    file=$2
    awk -v key="$key" '$1 == key { print $2; found = 1; exit } END { if (!found) exit 1 }' "$file"
}

require_key_value() {
    file=$1
    key=$2
    expected=$3
    value=$(get_value "$key" "$file" 2>/dev/null) || {
        add_reason "missing_${key}"
        return 1
    }
    if [ "$value" != "$expected" ]; then
        add_reason "${key}_expected_${expected}_got_${value}"
        return 1
    fi
    return 0
}

require_real_le() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$key" "$file" 2>/dev/null) || {
        add_reason "missing_${key}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value + 0.0) <= (limit + 0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_le_${limit}_got_${value}"
        return 1
    }
    return 0
}

require_real_ge() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$key" "$file" 2>/dev/null) || {
        add_reason "missing_${key}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value + 0.0) >= (limit + 0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_ge_${limit}_got_${value}"
        return 1
    }
    return 0
}

find_check_exe() {
    if [ -x "$BUILD_DIR/bin/fibre_stage14_rhs_sign_scaling_audit_check" ]; then
        echo "$BUILD_DIR/bin/fibre_stage14_rhs_sign_scaling_audit_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/src/fibre_stage14_rhs_sign_scaling_audit_check" ]; then
        echo "$BUILD_DIR/src/fibre_stage14_rhs_sign_scaling_audit_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/fibre_stage14_rhs_sign_scaling_audit_check" ]; then
        echo "$BUILD_DIR/fibre_stage14_rhs_sign_scaling_audit_check"
        return 0
    fi
    return 1
}

build_target() {
    target=$1
    ensure_build_dir || return 1
    cmake --build "$BUILD_DIR" --target "$target" -j
}

write_gate_dat() {
    gate_status=$1
    cat > "$GATE_DAT" <<EOF_GATE
stage14_3_gate_requested_flag 1
stage14_3_gate_build_status $BUILD_STATUS
stage14_3_gate_rhs_sign_scaling_audit_check_status $CHECK_STATUS
stage14_3_gate_fluid_to_fibre_sign_status $GATE_FLUID_TO_FIBRE_SIGN_STATUS
stage14_3_gate_fibre_to_fluid_sign_status $GATE_FIBRE_TO_FLUID_SIGN_STATUS
stage14_3_gate_action_reaction_status $GATE_ACTION_REACTION_STATUS
stage14_3_gate_rhs_increment_uses_fibre_to_fluid_status $GATE_RHS_INCREMENT_USES_FIBRE_TO_FLUID_STATUS
stage14_3_gate_wrong_sign_rejection_status $GATE_WRONG_SIGN_REJECTION_STATUS
stage14_3_gate_lambda_zero_scaling_status $GATE_LAMBDA_ZERO_SCALING_STATUS
stage14_3_gate_lambda_one_scaling_status $GATE_LAMBDA_ONE_SCALING_STATUS
stage14_3_gate_lambda_fractional_scaling_status $GATE_LAMBDA_FRACTIONAL_SCALING_STATUS
stage14_3_gate_lambda_negative_scaling_status $GATE_LAMBDA_NEGATIVE_SCALING_STATUS
stage14_3_gate_component_scaling_status $GATE_COMPONENT_SCALING_STATUS
stage14_3_gate_integrated_rhs_sign_status $GATE_INTEGRATED_RHS_SIGN_STATUS
stage14_3_gate_integrated_rhs_scaling_status $GATE_INTEGRATED_RHS_SCALING_STATUS
stage14_3_gate_finite_rhs_increment_status $GATE_FINITE_RHS_INCREMENT_STATUS
stage14_3_gate_no_production_rhs_modification_status $GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS
stage14_3_gate_no_xcompact3d_hook_status $GATE_NO_XCOMPACT3D_HOOK_STATUS
stage14_3_gate_no_fluid_field_access_status $GATE_NO_FLUID_FIELD_ACCESS_STATUS
stage14_3_gate_no_fluid_field_modification_status $GATE_NO_FLUID_FIELD_MODIFICATION_STATUS
stage14_3_gate_no_pressure_modification_status $GATE_NO_PRESSURE_MODIFICATION_STATUS
stage14_3_gate_no_projection_modification_status $GATE_NO_PROJECTION_MODIFICATION_STATUS
stage14_3_gate_no_rk3_modification_status $GATE_NO_RK3_MODIFICATION_STATUS
stage14_3_gate_no_channel_forcing_modification_status $GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_3_gate_no_production_ibm_forcing_status $GATE_NO_PRODUCTION_IBM_FORCING_STATUS
stage14_3_gate_no_feedback_application_status $GATE_NO_FEEDBACK_APPLICATION_STATUS
stage14_3_gate_no_twoway_force_status $GATE_NO_TWOWAY_FORCE_STATUS
stage14_3_gate_no_structure_advance_status $GATE_NO_STRUCTURE_ADVANCE_STATUS
stage14_3_gate_status $gate_status
EOF_GATE
}

run_audit_check() {
    log_file=$OUTPUT_DIR/stage14_3_rhs_sign_scaling_audit.log
    dat_file=$OUTPUT_DIR/fibre_stage14_3_rhs_sign_scaling_audit.dat
    rm -f "$dat_file"

    check_exe=$(find_check_exe) || {
        add_reason "missing_fibre_stage14_rhs_sign_scaling_audit_check_executable"
        return 1
    }

    (
        X3D_STAGE14_RHS_INJECTION=1
        X3D_STAGE14_INJECTION_GAIN=0.0
        X3D_STAGE14_MAX_STEPS=3
        X3D_STAGE14_REQUIRE_STAGE13=1
        X3D_STAGE14_DIAGNOSTIC_ONLY=1
        export X3D_STAGE14_RHS_INJECTION X3D_STAGE14_INJECTION_GAIN X3D_STAGE14_MAX_STEPS
        export X3D_STAGE14_REQUIRE_STAGE13 X3D_STAGE14_DIAGNOSTIC_ONLY
        "$check_exe"
    ) > "$log_file" 2>&1 || {
        add_reason "fibre_stage14_rhs_sign_scaling_audit_check_failed"
        cat "$log_file"
        return 1
    }

    grep 'STAGE 14.3 RHS SIGN SCALING AUDIT VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage14_3_pass_verdict"
        return 1
    }

    if [ ! -f "$dat_file" ]; then
        add_reason "missing_fibre_stage14_3_rhs_sign_scaling_audit_dat"
        return 1
    fi

    status=0
    require_key_value "$dat_file" stage14_3_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_3_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_3_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_3_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_3_fluid_to_fibre_sign_status 1 || status=1
    require_key_value "$dat_file" stage14_3_fibre_to_fluid_sign_status 1 || status=1
    require_key_value "$dat_file" stage14_3_action_reaction_status 1 || status=1
    require_key_value "$dat_file" stage14_3_rhs_increment_uses_fibre_to_fluid_status 1 || status=1
    require_key_value "$dat_file" stage14_3_wrong_sign_rejection_status 1 || status=1
    require_key_value "$dat_file" stage14_3_lambda_zero_scaling_status 1 || status=1
    require_key_value "$dat_file" stage14_3_lambda_one_scaling_status 1 || status=1
    require_key_value "$dat_file" stage14_3_lambda_fractional_scaling_status 1 || status=1
    require_key_value "$dat_file" stage14_3_lambda_negative_scaling_status 1 || status=1
    require_key_value "$dat_file" stage14_3_component_scaling_status 1 || status=1
    require_key_value "$dat_file" stage14_3_integrated_rhs_sign_status 1 || status=1
    require_key_value "$dat_file" stage14_3_integrated_rhs_scaling_status 1 || status=1
    require_key_value "$dat_file" stage14_3_finite_rhs_increment_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_production_rhs_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_xcompact3d_hook_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_fluid_field_access_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_fluid_field_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_pressure_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_projection_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_rk3_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_feedback_application_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_twoway_force_status 1 || status=1
    require_key_value "$dat_file" stage14_3_no_structure_advance_status 1 || status=1
    require_key_value "$dat_file" stage14_3_rhs_sign_scaling_audit_status 1 || status=1

    require_real_le "$dat_file" stage14_3_action_reaction_max_error 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_lambda_zero_integrated_error_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_lambda_one_integrated_error_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_lambda_fractional_integrated_error_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_lambda_negative_integrated_error_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_component_error_x 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_component_error_y 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_3_component_error_z 1.0e-12 || status=1
    require_real_ge "$dat_file" stage14_3_wrong_sign_separation 1.0e-8 || status=1

    return $status
}

GATE_FLUID_TO_FIBRE_SIGN_STATUS=0
GATE_FIBRE_TO_FLUID_SIGN_STATUS=0
GATE_ACTION_REACTION_STATUS=0
GATE_RHS_INCREMENT_USES_FIBRE_TO_FLUID_STATUS=0
GATE_WRONG_SIGN_REJECTION_STATUS=0
GATE_LAMBDA_ZERO_SCALING_STATUS=0
GATE_LAMBDA_ONE_SCALING_STATUS=0
GATE_LAMBDA_FRACTIONAL_SCALING_STATUS=0
GATE_LAMBDA_NEGATIVE_SCALING_STATUS=0
GATE_COMPONENT_SCALING_STATUS=0
GATE_INTEGRATED_RHS_SIGN_STATUS=0
GATE_INTEGRATED_RHS_SCALING_STATUS=0
GATE_FINITE_RHS_INCREMENT_STATUS=0
GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS=0
GATE_NO_XCOMPACT3D_HOOK_STATUS=0
GATE_NO_FLUID_FIELD_ACCESS_STATUS=0
GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=0
GATE_NO_PRESSURE_MODIFICATION_STATUS=0
GATE_NO_PROJECTION_MODIFICATION_STATUS=0
GATE_NO_RK3_MODIFICATION_STATUS=0
GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
GATE_NO_PRODUCTION_IBM_FORCING_STATUS=0
GATE_NO_FEEDBACK_APPLICATION_STATUS=0
GATE_NO_TWOWAY_FORCE_STATUS=0
GATE_NO_STRUCTURE_ADVANCE_STATUS=0

if [ "$STAGE14_3_RUN_STAGE14_2" = "1" ]; then
    bash stage14_checks/run_stage14_2_rhs_addition_formula.sh || {
        BUILD_STATUS=0
        add_reason "stage14_2_rhs_addition_formula_failed"
    }
fi

targets="xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check fibre_stage12_sign_convention_audit_check fibre_stage12_power_diagnostics_check fibre_stage12_production_feedback_candidate_check fibre_stage13_config_check fibre_stage13_force_density_buffer_check fibre_stage13_spreading_kernel_check fibre_stage13_volume_normalization_audit_check fibre_stage13_conservation_sign_audit_check fibre_stage13_production_force_density_candidate_check fibre_stage14_config_check fibre_stage14_rhs_accumulator_check fibre_stage14_rhs_addition_formula_check fibre_stage14_rhs_sign_scaling_audit_check"

for target in $targets; do
    build_target "$target" || {
        BUILD_STATUS=0
        add_reason "build_failed_${target}"
        break
    }
done

if [ "$BUILD_STATUS" = "1" ]; then
    run_audit_check && CHECK_STATUS=1
fi

if [ "$CHECK_STATUS" = "1" ]; then
    dat_file=$OUTPUT_DIR/fibre_stage14_3_rhs_sign_scaling_audit.dat
    GATE_FLUID_TO_FIBRE_SIGN_STATUS=$(get_value stage14_3_fluid_to_fibre_sign_status "$dat_file" 2>/dev/null || echo 0)
    GATE_FIBRE_TO_FLUID_SIGN_STATUS=$(get_value stage14_3_fibre_to_fluid_sign_status "$dat_file" 2>/dev/null || echo 0)
    GATE_ACTION_REACTION_STATUS=$(get_value stage14_3_action_reaction_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_INCREMENT_USES_FIBRE_TO_FLUID_STATUS=$(get_value stage14_3_rhs_increment_uses_fibre_to_fluid_status "$dat_file" 2>/dev/null || echo 0)
    GATE_WRONG_SIGN_REJECTION_STATUS=$(get_value stage14_3_wrong_sign_rejection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_ZERO_SCALING_STATUS=$(get_value stage14_3_lambda_zero_scaling_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_ONE_SCALING_STATUS=$(get_value stage14_3_lambda_one_scaling_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_FRACTIONAL_SCALING_STATUS=$(get_value stage14_3_lambda_fractional_scaling_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_NEGATIVE_SCALING_STATUS=$(get_value stage14_3_lambda_negative_scaling_status "$dat_file" 2>/dev/null || echo 0)
    GATE_COMPONENT_SCALING_STATUS=$(get_value stage14_3_component_scaling_status "$dat_file" 2>/dev/null || echo 0)
    GATE_INTEGRATED_RHS_SIGN_STATUS=$(get_value stage14_3_integrated_rhs_sign_status "$dat_file" 2>/dev/null || echo 0)
    GATE_INTEGRATED_RHS_SCALING_STATUS=$(get_value stage14_3_integrated_rhs_scaling_status "$dat_file" 2>/dev/null || echo 0)
    GATE_FINITE_RHS_INCREMENT_STATUS=$(get_value stage14_3_finite_rhs_increment_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS=$(get_value stage14_3_no_production_rhs_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_XCOMPACT3D_HOOK_STATUS=$(get_value stage14_3_no_xcompact3d_hook_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_ACCESS_STATUS=$(get_value stage14_3_no_fluid_field_access_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=$(get_value stage14_3_no_fluid_field_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_3_no_pressure_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_3_no_projection_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_RK3_MODIFICATION_STATUS=$(get_value stage14_3_no_rk3_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_3_no_channel_forcing_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_3_no_production_ibm_forcing_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_3_no_feedback_application_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_TWOWAY_FORCE_STATUS=$(get_value stage14_3_no_twoway_force_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_3_no_structure_advance_status "$dat_file" 2>/dev/null || echo 0)
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$CHECK_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.3 RHS SIGN SCALING AUDIT VERDICT: PASS'
    echo 'STAGE 14.3 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.3 RHS SIGN SCALING AUDIT VERDICT: FAIL'
echo 'STAGE 14.3 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_3_failure'
fi
exit 1
