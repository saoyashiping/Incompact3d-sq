#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_4_RUN_STAGE14_3=${STAGE14_4_RUN_STAGE14_3:-0}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_4_rk_timing_audit_gate.dat
REASONS_FILE=$OUTPUT_DIR/stage14_4_rk_timing_audit_reasons.tmp
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
    if [ -x "$BUILD_DIR/bin/fibre_stage14_rk_timing_audit_check" ]; then
        echo "$BUILD_DIR/bin/fibre_stage14_rk_timing_audit_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/src/fibre_stage14_rk_timing_audit_check" ]; then
        echo "$BUILD_DIR/src/fibre_stage14_rk_timing_audit_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/fibre_stage14_rk_timing_audit_check" ]; then
        echo "$BUILD_DIR/fibre_stage14_rk_timing_audit_check"
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
stage14_4_gate_requested_flag 1
stage14_4_gate_build_status $BUILD_STATUS
stage14_4_gate_rk_timing_audit_check_status $CHECK_STATUS
stage14_4_gate_rk_substep_count_status $GATE_RK_SUBSTEP_COUNT_STATUS
stage14_4_gate_stage13_candidate_before_increment_status $GATE_STAGE13_CANDIDATE_BEFORE_INCREMENT_STATUS
stage14_4_gate_increment_before_rhs_addition_status $GATE_INCREMENT_BEFORE_RHS_ADDITION_STATUS
stage14_4_gate_rhs_addition_before_projection_status $GATE_RHS_ADDITION_BEFORE_PROJECTION_STATUS
stage14_4_gate_once_per_substep_status $GATE_ONCE_PER_SUBSTEP_STATUS
stage14_4_gate_no_duplicate_injection_status $GATE_NO_DUPLICATE_INJECTION_STATUS
stage14_4_gate_no_skipped_substep_status $GATE_NO_SKIPPED_SUBSTEP_STATUS
stage14_4_gate_duplicate_detection_status $GATE_DUPLICATE_DETECTION_STATUS
stage14_4_gate_skipped_substep_detection_status $GATE_SKIPPED_SUBSTEP_DETECTION_STATUS
stage14_4_gate_wrong_order_detection_status $GATE_WRONG_ORDER_DETECTION_STATUS
stage14_4_gate_invalid_event_rejection_status $GATE_INVALID_EVENT_REJECTION_STATUS
stage14_4_gate_invalid_substep_rejection_status $GATE_INVALID_SUBSTEP_REJECTION_STATUS
stage14_4_gate_lambda_zero_increment_status $GATE_LAMBDA_ZERO_INCREMENT_STATUS
stage14_4_gate_timing_diagnostics_finite_status $GATE_TIMING_DIAGNOSTICS_FINITE_STATUS
stage14_4_gate_no_production_rhs_modification_status $GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS
stage14_4_gate_no_xcompact3d_hook_status $GATE_NO_XCOMPACT3D_HOOK_STATUS
stage14_4_gate_no_fluid_field_access_status $GATE_NO_FLUID_FIELD_ACCESS_STATUS
stage14_4_gate_no_fluid_field_modification_status $GATE_NO_FLUID_FIELD_MODIFICATION_STATUS
stage14_4_gate_no_pressure_modification_status $GATE_NO_PRESSURE_MODIFICATION_STATUS
stage14_4_gate_no_projection_modification_status $GATE_NO_PROJECTION_MODIFICATION_STATUS
stage14_4_gate_no_poisson_modification_status $GATE_NO_POISSON_MODIFICATION_STATUS
stage14_4_gate_no_rk3_modification_status $GATE_NO_RK3_MODIFICATION_STATUS
stage14_4_gate_no_channel_forcing_modification_status $GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_4_gate_no_production_ibm_forcing_status $GATE_NO_PRODUCTION_IBM_FORCING_STATUS
stage14_4_gate_no_feedback_application_status $GATE_NO_FEEDBACK_APPLICATION_STATUS
stage14_4_gate_no_twoway_force_status $GATE_NO_TWOWAY_FORCE_STATUS
stage14_4_gate_no_structure_advance_status $GATE_NO_STRUCTURE_ADVANCE_STATUS
stage14_4_gate_status $gate_status
EOF_GATE
}

run_timing_check() {
    log_file=$OUTPUT_DIR/stage14_4_rk_timing_audit.log
    dat_file=$OUTPUT_DIR/fibre_stage14_4_rk_timing_audit.dat
    rm -f "$dat_file"

    check_exe=$(find_check_exe) || {
        add_reason "missing_fibre_stage14_rk_timing_audit_check_executable"
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
        add_reason "fibre_stage14_rk_timing_audit_check_failed"
        cat "$log_file"
        return 1
    }

    grep 'STAGE 14.4 RK TIMING AUDIT VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage14_4_pass_verdict"
        return 1
    }

    if [ ! -f "$dat_file" ]; then
        add_reason "missing_fibre_stage14_4_rk_timing_audit_dat"
        return 1
    fi

    status=0
    require_key_value "$dat_file" stage14_4_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_4_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_4_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_4_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_4_rk_substep_count_status 1 || status=1
    require_key_value "$dat_file" stage14_4_stage13_candidate_before_increment_status 1 || status=1
    require_key_value "$dat_file" stage14_4_increment_before_rhs_addition_status 1 || status=1
    require_key_value "$dat_file" stage14_4_rhs_addition_before_projection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_once_per_substep_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_duplicate_injection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_skipped_substep_status 1 || status=1
    require_key_value "$dat_file" stage14_4_duplicate_detection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_skipped_substep_detection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_wrong_order_detection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_invalid_event_rejection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_invalid_substep_rejection_status 1 || status=1
    require_key_value "$dat_file" stage14_4_lambda_zero_increment_status 1 || status=1
    require_key_value "$dat_file" stage14_4_timing_diagnostics_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_production_rhs_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_xcompact3d_hook_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_fluid_field_access_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_fluid_field_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_pressure_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_projection_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_poisson_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_rk3_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_feedback_application_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_twoway_force_status 1 || status=1
    require_key_value "$dat_file" stage14_4_no_structure_advance_status 1 || status=1
    require_key_value "$dat_file" stage14_4_rk_timing_audit_status 1 || status=1

    require_key_value "$dat_file" stage14_4_expected_nsubsteps 3 || status=1
    require_key_value "$dat_file" stage14_4_recorded_rhs_addition_count 3 || status=1
    require_key_value "$dat_file" stage14_4_duplicate_injection_count 0 || status=1
    require_key_value "$dat_file" stage14_4_skipped_substep_count 0 || status=1
    require_key_value "$dat_file" stage14_4_wrong_order_count 0 || status=1
    require_real_le "$dat_file" stage14_4_lambda_zero_max_abs_rhs_increment 1.0e-12 || status=1

    return $status
}

GATE_RK_SUBSTEP_COUNT_STATUS=0
GATE_STAGE13_CANDIDATE_BEFORE_INCREMENT_STATUS=0
GATE_INCREMENT_BEFORE_RHS_ADDITION_STATUS=0
GATE_RHS_ADDITION_BEFORE_PROJECTION_STATUS=0
GATE_ONCE_PER_SUBSTEP_STATUS=0
GATE_NO_DUPLICATE_INJECTION_STATUS=0
GATE_NO_SKIPPED_SUBSTEP_STATUS=0
GATE_DUPLICATE_DETECTION_STATUS=0
GATE_SKIPPED_SUBSTEP_DETECTION_STATUS=0
GATE_WRONG_ORDER_DETECTION_STATUS=0
GATE_INVALID_EVENT_REJECTION_STATUS=0
GATE_INVALID_SUBSTEP_REJECTION_STATUS=0
GATE_LAMBDA_ZERO_INCREMENT_STATUS=0
GATE_TIMING_DIAGNOSTICS_FINITE_STATUS=0
GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS=0
GATE_NO_XCOMPACT3D_HOOK_STATUS=0
GATE_NO_FLUID_FIELD_ACCESS_STATUS=0
GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=0
GATE_NO_PRESSURE_MODIFICATION_STATUS=0
GATE_NO_PROJECTION_MODIFICATION_STATUS=0
GATE_NO_POISSON_MODIFICATION_STATUS=0
GATE_NO_RK3_MODIFICATION_STATUS=0
GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
GATE_NO_PRODUCTION_IBM_FORCING_STATUS=0
GATE_NO_FEEDBACK_APPLICATION_STATUS=0
GATE_NO_TWOWAY_FORCE_STATUS=0
GATE_NO_STRUCTURE_ADVANCE_STATUS=0

if [ "$STAGE14_4_RUN_STAGE14_3" = "1" ]; then
    bash stage14_checks/run_stage14_3_rhs_sign_scaling_audit.sh || {
        BUILD_STATUS=0
        add_reason "stage14_3_rhs_sign_scaling_audit_failed"
    }
fi

targets="xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check fibre_stage12_sign_convention_audit_check fibre_stage12_power_diagnostics_check fibre_stage12_production_feedback_candidate_check fibre_stage13_config_check fibre_stage13_force_density_buffer_check fibre_stage13_spreading_kernel_check fibre_stage13_volume_normalization_audit_check fibre_stage13_conservation_sign_audit_check fibre_stage13_production_force_density_candidate_check fibre_stage14_config_check fibre_stage14_rhs_accumulator_check fibre_stage14_rhs_addition_formula_check fibre_stage14_rhs_sign_scaling_audit_check fibre_stage14_rk_timing_audit_check"

for target in $targets; do
    build_target "$target" || {
        BUILD_STATUS=0
        add_reason "build_failed_${target}"
        break
    }
done

if [ "$BUILD_STATUS" = "1" ]; then
    run_timing_check && CHECK_STATUS=1
fi

if [ "$CHECK_STATUS" = "1" ]; then
    dat_file=$OUTPUT_DIR/fibre_stage14_4_rk_timing_audit.dat
    GATE_RK_SUBSTEP_COUNT_STATUS=$(get_value stage14_4_rk_substep_count_status "$dat_file" 2>/dev/null || echo 0)
    GATE_STAGE13_CANDIDATE_BEFORE_INCREMENT_STATUS=$(get_value stage14_4_stage13_candidate_before_increment_status "$dat_file" 2>/dev/null || echo 0)
    GATE_INCREMENT_BEFORE_RHS_ADDITION_STATUS=$(get_value stage14_4_increment_before_rhs_addition_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_ADDITION_BEFORE_PROJECTION_STATUS=$(get_value stage14_4_rhs_addition_before_projection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_ONCE_PER_SUBSTEP_STATUS=$(get_value stage14_4_once_per_substep_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_DUPLICATE_INJECTION_STATUS=$(get_value stage14_4_no_duplicate_injection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_SKIPPED_SUBSTEP_STATUS=$(get_value stage14_4_no_skipped_substep_status "$dat_file" 2>/dev/null || echo 0)
    GATE_DUPLICATE_DETECTION_STATUS=$(get_value stage14_4_duplicate_detection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_SKIPPED_SUBSTEP_DETECTION_STATUS=$(get_value stage14_4_skipped_substep_detection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_WRONG_ORDER_DETECTION_STATUS=$(get_value stage14_4_wrong_order_detection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_INVALID_EVENT_REJECTION_STATUS=$(get_value stage14_4_invalid_event_rejection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_INVALID_SUBSTEP_REJECTION_STATUS=$(get_value stage14_4_invalid_substep_rejection_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_ZERO_INCREMENT_STATUS=$(get_value stage14_4_lambda_zero_increment_status "$dat_file" 2>/dev/null || echo 0)
    GATE_TIMING_DIAGNOSTICS_FINITE_STATUS=$(get_value stage14_4_timing_diagnostics_finite_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_RHS_MODIFICATION_STATUS=$(get_value stage14_4_no_production_rhs_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_XCOMPACT3D_HOOK_STATUS=$(get_value stage14_4_no_xcompact3d_hook_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_ACCESS_STATUS=$(get_value stage14_4_no_fluid_field_access_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=$(get_value stage14_4_no_fluid_field_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_4_no_pressure_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_4_no_projection_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_POISSON_MODIFICATION_STATUS=$(get_value stage14_4_no_poisson_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_RK3_MODIFICATION_STATUS=$(get_value stage14_4_no_rk3_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_4_no_channel_forcing_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_4_no_production_ibm_forcing_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_4_no_feedback_application_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_TWOWAY_FORCE_STATUS=$(get_value stage14_4_no_twoway_force_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_4_no_structure_advance_status "$dat_file" 2>/dev/null || echo 0)
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$CHECK_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.4 RK TIMING AUDIT VERDICT: PASS'
    echo 'STAGE 14.4 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.4 RK TIMING AUDIT VERDICT: FAIL'
echo 'STAGE 14.4 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_4_failure'
fi
exit 1
