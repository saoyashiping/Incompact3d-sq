#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_2_RUN_STAGE14_1=${STAGE14_2_RUN_STAGE14_1:-0}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_2_rhs_addition_formula_gate.dat
REASONS_FILE=$OUTPUT_DIR/stage14_2_rhs_addition_formula_reasons.tmp
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

find_check_exe() {
    if [ -x "$BUILD_DIR/bin/fibre_stage14_rhs_addition_formula_check" ]; then
        echo "$BUILD_DIR/bin/fibre_stage14_rhs_addition_formula_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/src/fibre_stage14_rhs_addition_formula_check" ]; then
        echo "$BUILD_DIR/src/fibre_stage14_rhs_addition_formula_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/fibre_stage14_rhs_addition_formula_check" ]; then
        echo "$BUILD_DIR/fibre_stage14_rhs_addition_formula_check"
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
stage14_2_gate_requested_flag 1
stage14_2_gate_build_status $BUILD_STATUS
stage14_2_gate_rhs_addition_formula_check_status $CHECK_STATUS
stage14_2_gate_shape_status $GATE_SHAPE_STATUS
stage14_2_gate_lambda_zero_invariance_status $GATE_LAMBDA_ZERO_INVARIANCE_STATUS
stage14_2_gate_lambda_one_addition_status $GATE_LAMBDA_ONE_ADDITION_STATUS
stage14_2_gate_lambda_fractional_addition_status $GATE_LAMBDA_FRACTIONAL_ADDITION_STATUS
stage14_2_gate_component_addition_status $GATE_COMPONENT_ADDITION_STATUS
stage14_2_gate_additive_not_overwrite_status $GATE_ADDITIVE_NOT_OVERWRITE_STATUS
stage14_2_gate_rhs_old_preserved_status $GATE_RHS_OLD_PRESERVED_STATUS
stage14_2_gate_finite_rhs_old_status $GATE_FINITE_RHS_OLD_STATUS
stage14_2_gate_finite_rhs_increment_status $GATE_FINITE_RHS_INCREMENT_STATUS
stage14_2_gate_finite_rhs_new_status $GATE_FINITE_RHS_NEW_STATUS
stage14_2_gate_no_production_rhs_modification_status $GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS
stage14_2_gate_no_xcompact3d_hook_status $GATE_NO_XCOMPACT3D_HOOK_STATUS
stage14_2_gate_no_fluid_field_access_status $GATE_NO_FLUID_FIELD_ACCESS_STATUS
stage14_2_gate_no_fluid_field_modification_status $GATE_NO_FLUID_FIELD_MODIFICATION_STATUS
stage14_2_gate_no_pressure_modification_status $GATE_NO_PRESSURE_MODIFICATION_STATUS
stage14_2_gate_no_projection_modification_status $GATE_NO_PROJECTION_MODIFICATION_STATUS
stage14_2_gate_no_rk3_modification_status $GATE_NO_RK3_MODIFICATION_STATUS
stage14_2_gate_no_channel_forcing_modification_status $GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_2_gate_no_production_ibm_forcing_status $GATE_NO_PRODUCTION_IBM_FORCING_STATUS
stage14_2_gate_no_feedback_application_status $GATE_NO_FEEDBACK_APPLICATION_STATUS
stage14_2_gate_no_twoway_force_status $GATE_NO_TWOWAY_FORCE_STATUS
stage14_2_gate_no_structure_advance_status $GATE_NO_STRUCTURE_ADVANCE_STATUS
stage14_2_gate_status $gate_status
EOF_GATE
}

run_formula_check() {
    log_file=$OUTPUT_DIR/stage14_2_rhs_addition_formula.log
    dat_file=$OUTPUT_DIR/fibre_stage14_2_rhs_addition_formula.dat
    rm -f "$dat_file"

    check_exe=$(find_check_exe) || {
        add_reason "missing_fibre_stage14_rhs_addition_formula_check_executable"
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
        add_reason "fibre_stage14_rhs_addition_formula_check_failed"
        cat "$log_file"
        return 1
    }

    grep 'STAGE 14.2 RHS ADDITION FORMULA VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage14_2_pass_verdict"
        return 1
    }

    if [ ! -f "$dat_file" ]; then
        add_reason "missing_fibre_stage14_2_rhs_addition_formula_dat"
        return 1
    fi

    status=0
    require_key_value "$dat_file" stage14_2_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_2_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_2_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_2_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_2_shape_status 1 || status=1
    require_key_value "$dat_file" stage14_2_lambda_zero_invariance_status 1 || status=1
    require_key_value "$dat_file" stage14_2_lambda_one_addition_status 1 || status=1
    require_key_value "$dat_file" stage14_2_lambda_fractional_addition_status 1 || status=1
    require_key_value "$dat_file" stage14_2_component_addition_status 1 || status=1
    require_key_value "$dat_file" stage14_2_additive_not_overwrite_status 1 || status=1
    require_key_value "$dat_file" stage14_2_rhs_old_preserved_status 1 || status=1
    require_key_value "$dat_file" stage14_2_finite_rhs_old_status 1 || status=1
    require_key_value "$dat_file" stage14_2_finite_rhs_increment_status 1 || status=1
    require_key_value "$dat_file" stage14_2_finite_rhs_new_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_production_rhs_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_xcompact3d_hook_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_fluid_field_access_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_fluid_field_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_pressure_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_projection_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_rk3_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_feedback_application_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_twoway_force_status 1 || status=1
    require_key_value "$dat_file" stage14_2_no_structure_advance_status 1 || status=1
    require_key_value "$dat_file" stage14_2_rhs_addition_formula_status 1 || status=1

    require_real_le "$dat_file" stage14_2_lambda_zero_max_abs 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_2_lambda_one_max_error 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_2_lambda_fractional_max_error 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_2_component_error_x 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_2_component_error_y 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_2_component_error_z 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_2_rhs_old_preservation_error 1.0e-12 || status=1

    return $status
}

GATE_SHAPE_STATUS=0
GATE_LAMBDA_ZERO_INVARIANCE_STATUS=0
GATE_LAMBDA_ONE_ADDITION_STATUS=0
GATE_LAMBDA_FRACTIONAL_ADDITION_STATUS=0
GATE_COMPONENT_ADDITION_STATUS=0
GATE_ADDITIVE_NOT_OVERWRITE_STATUS=0
GATE_RHS_OLD_PRESERVED_STATUS=0
GATE_FINITE_RHS_OLD_STATUS=0
GATE_FINITE_RHS_INCREMENT_STATUS=0
GATE_FINITE_RHS_NEW_STATUS=0
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

if [ "$STAGE14_2_RUN_STAGE14_1" = "1" ]; then
    bash stage14_checks/run_stage14_1_rhs_accumulator.sh || add_reason "stage14_1_rhs_accumulator_failed"
fi

targets="xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check fibre_stage12_sign_convention_audit_check fibre_stage12_power_diagnostics_check fibre_stage12_production_feedback_candidate_check fibre_stage13_config_check fibre_stage13_force_density_buffer_check fibre_stage13_spreading_kernel_check fibre_stage13_volume_normalization_audit_check fibre_stage13_conservation_sign_audit_check fibre_stage13_production_force_density_candidate_check fibre_stage14_config_check fibre_stage14_rhs_accumulator_check fibre_stage14_rhs_addition_formula_check"

for target in $targets; do
    build_target "$target" || {
        BUILD_STATUS=0
        add_reason "build_failed_${target}"
        break
    }
done

if [ "$BUILD_STATUS" = "1" ]; then
    run_formula_check && CHECK_STATUS=1
fi

if [ "$CHECK_STATUS" = "1" ]; then
    dat_file=$OUTPUT_DIR/fibre_stage14_2_rhs_addition_formula.dat
    GATE_SHAPE_STATUS=$(get_value stage14_2_shape_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_ZERO_INVARIANCE_STATUS=$(get_value stage14_2_lambda_zero_invariance_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_ONE_ADDITION_STATUS=$(get_value stage14_2_lambda_one_addition_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_FRACTIONAL_ADDITION_STATUS=$(get_value stage14_2_lambda_fractional_addition_status "$dat_file" 2>/dev/null || echo 0)
    GATE_COMPONENT_ADDITION_STATUS=$(get_value stage14_2_component_addition_status "$dat_file" 2>/dev/null || echo 0)
    GATE_ADDITIVE_NOT_OVERWRITE_STATUS=$(get_value stage14_2_additive_not_overwrite_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_OLD_PRESERVED_STATUS=$(get_value stage14_2_rhs_old_preserved_status "$dat_file" 2>/dev/null || echo 0)
    GATE_FINITE_RHS_OLD_STATUS=$(get_value stage14_2_finite_rhs_old_status "$dat_file" 2>/dev/null || echo 0)
    GATE_FINITE_RHS_INCREMENT_STATUS=$(get_value stage14_2_finite_rhs_increment_status "$dat_file" 2>/dev/null || echo 0)
    GATE_FINITE_RHS_NEW_STATUS=$(get_value stage14_2_finite_rhs_new_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS=$(get_value stage14_2_no_production_rhs_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_XCOMPACT3D_HOOK_STATUS=$(get_value stage14_2_no_xcompact3d_hook_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_ACCESS_STATUS=$(get_value stage14_2_no_fluid_field_access_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=$(get_value stage14_2_no_fluid_field_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_2_no_pressure_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_2_no_projection_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_RK3_MODIFICATION_STATUS=$(get_value stage14_2_no_rk3_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_2_no_channel_forcing_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_2_no_production_ibm_forcing_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_2_no_feedback_application_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_TWOWAY_FORCE_STATUS=$(get_value stage14_2_no_twoway_force_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_2_no_structure_advance_status "$dat_file" 2>/dev/null || echo 0)
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$CHECK_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.2 RHS ADDITION FORMULA VERDICT: PASS'
    echo 'STAGE 14.2 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.2 RHS ADDITION FORMULA VERDICT: FAIL'
echo 'STAGE 14.2 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_2_failure'
fi
exit 1
