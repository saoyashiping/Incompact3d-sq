#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_0_RUN_STAGE13_CLOSURE=${STAGE14_0_RUN_STAGE13_CLOSURE:-0}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_0_config_gate.dat
REASONS_FILE=$OUTPUT_DIR/stage14_0_config_reasons.tmp
BUILD_STATUS=1
DEFAULT_STATUS=0
ZERO_GAIN_STATUS=0
INVALID_STATUS=0

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
    phase=$4
    value=$(get_value "$key" "$file" 2>/dev/null) || {
        add_reason "${phase}:missing_${key}"
        return 1
    }
    if [ "$value" != "$expected" ]; then
        add_reason "${phase}:${key}_expected_${expected}_got_${value}"
        return 1
    fi
    return 0
}

require_real_zero() {
    file=$1
    key=$2
    phase=$3
    value=$(get_value "$key" "$file" 2>/dev/null) || {
        add_reason "${phase}:missing_${key}"
        return 1
    }
    awk -v value="$value" 'BEGIN { if (value + 0.0 == 0.0) exit 0; exit 1 }' || {
        add_reason "${phase}:${key}_expected_zero_got_${value}"
        return 1
    }
    return 0
}

find_check_exe() {
    if [ -x "$BUILD_DIR/bin/fibre_stage14_config_check" ]; then
        echo "$BUILD_DIR/bin/fibre_stage14_config_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/src/fibre_stage14_config_check" ]; then
        echo "$BUILD_DIR/src/fibre_stage14_config_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/fibre_stage14_config_check" ]; then
        echo "$BUILD_DIR/fibre_stage14_config_check"
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
stage14_0_gate_requested_flag 1
stage14_0_gate_build_status $BUILD_STATUS
stage14_0_gate_default_disabled_status $DEFAULT_STATUS
stage14_0_gate_requested_zero_gain_status $ZERO_GAIN_STATUS
stage14_0_gate_invalid_fallback_status $INVALID_STATUS
stage14_0_gate_gain_parse_status $GATE_GAIN_PARSE_STATUS
stage14_0_gate_gain_default_zero_status $GATE_GAIN_DEFAULT_ZERO_STATUS
stage14_0_gate_gain_finite_status $GATE_GAIN_FINITE_STATUS
stage14_0_gate_safe_fallback_status $GATE_SAFE_FALLBACK_STATUS
stage14_0_gate_max_steps_parse_status $GATE_MAX_STEPS_PARSE_STATUS
stage14_0_gate_require_stage13_status $GATE_REQUIRE_STAGE13_STATUS
stage14_0_gate_diagnostic_only_status $GATE_DIAGNOSTIC_ONLY_STATUS
stage14_0_gate_no_rhs_buffer_allocation_status $GATE_NO_RHS_BUFFER_ALLOCATION_STATUS
stage14_0_gate_no_rhs_injection_status $GATE_NO_RHS_INJECTION_STATUS
stage14_0_gate_no_production_hook_status $GATE_NO_PRODUCTION_HOOK_STATUS
stage14_0_gate_no_fluid_field_modification_status $GATE_NO_FLUID_FIELD_MODIFICATION_STATUS
stage14_0_gate_no_pressure_modification_status $GATE_NO_PRESSURE_MODIFICATION_STATUS
stage14_0_gate_no_projection_modification_status $GATE_NO_PROJECTION_MODIFICATION_STATUS
stage14_0_gate_no_rk3_modification_status $GATE_NO_RK3_MODIFICATION_STATUS
stage14_0_gate_no_channel_forcing_modification_status $GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_0_gate_no_production_ibm_forcing_status $GATE_NO_PRODUCTION_IBM_FORCING_STATUS
stage14_0_gate_no_feedback_application_status $GATE_NO_FEEDBACK_APPLICATION_STATUS
stage14_0_gate_no_twoway_force_status $GATE_NO_TWOWAY_FORCE_STATUS
stage14_0_gate_no_structure_advance_status $GATE_NO_STRUCTURE_ADVANCE_STATUS
stage14_0_gate_status $gate_status
EOF_GATE
}

run_check_case() {
    phase=$1
    log_file=$2
    dat_copy=$3
    mode=$4

    rm -f "$OUTPUT_DIR/fibre_stage14_0_config.dat"
    check_exe=$(find_check_exe) || {
        add_reason "${phase}:missing_fibre_stage14_config_check_executable"
        return 1
    }

    if [ "$mode" = "default" ]; then
        (
            unset X3D_STAGE14_RHS_INJECTION
            unset X3D_STAGE14_INJECTION_GAIN
            unset X3D_STAGE14_MAX_STEPS
            unset X3D_STAGE14_REQUIRE_STAGE13
            unset X3D_STAGE14_DIAGNOSTIC_ONLY
            "$check_exe"
        ) > "$log_file" 2>&1
    elif [ "$mode" = "zero_gain" ]; then
        (
            X3D_STAGE14_RHS_INJECTION=1
            X3D_STAGE14_INJECTION_GAIN=0.0
            X3D_STAGE14_MAX_STEPS=3
            X3D_STAGE14_REQUIRE_STAGE13=1
            X3D_STAGE14_DIAGNOSTIC_ONLY=1
            export X3D_STAGE14_RHS_INJECTION X3D_STAGE14_INJECTION_GAIN X3D_STAGE14_MAX_STEPS
            export X3D_STAGE14_REQUIRE_STAGE13 X3D_STAGE14_DIAGNOSTIC_ONLY
            "$check_exe"
        ) > "$log_file" 2>&1
    else
        (
            X3D_STAGE14_RHS_INJECTION=1
            X3D_STAGE14_INJECTION_GAIN=not_a_number
            X3D_STAGE14_MAX_STEPS=-5
            X3D_STAGE14_REQUIRE_STAGE13=invalid_value
            X3D_STAGE14_DIAGNOSTIC_ONLY=0
            export X3D_STAGE14_RHS_INJECTION X3D_STAGE14_INJECTION_GAIN X3D_STAGE14_MAX_STEPS
            export X3D_STAGE14_REQUIRE_STAGE13 X3D_STAGE14_DIAGNOSTIC_ONLY
            "$check_exe"
        ) > "$log_file" 2>&1
    fi

    if [ $? -ne 0 ]; then
        add_reason "${phase}:fibre_stage14_config_check_failed"
        cat "$log_file"
        return 1
    fi

    grep 'STAGE 14.0 CONFIG VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "${phase}:missing_pass_verdict"
        return 1
    }

    if [ ! -f "$OUTPUT_DIR/fibre_stage14_0_config.dat" ]; then
        add_reason "${phase}:missing_fibre_stage14_0_config_dat"
        return 1
    fi

    cp "$OUTPUT_DIR/fibre_stage14_0_config.dat" "$dat_copy"
    return 0
}

run_default_disabled() {
    log_file=$OUTPUT_DIR/stage14_0_config_default_disabled.log
    dat_file=$OUTPUT_DIR/fibre_stage14_0_config_default_disabled.dat
    run_check_case default_disabled "$log_file" "$dat_file" default || return 1

    status=0
    require_key_value "$dat_file" stage14_0_requested_flag 0 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_rhs_injection_enabled_flag 0 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_disabled_by_default_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_gain_default_zero_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_gain_finite_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_require_stage13_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_diagnostic_only_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_no_rhs_buffer_allocation_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_no_rhs_injection_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_no_production_hook_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_no_fluid_field_modification_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_no_structure_advance_status 1 default_disabled || status=1
    require_key_value "$dat_file" stage14_0_config_status 1 default_disabled || status=1
    return $status
}

run_requested_zero_gain() {
    log_file=$OUTPUT_DIR/stage14_0_config_requested_zero_gain.log
    dat_file=$OUTPUT_DIR/fibre_stage14_0_config_requested_zero_gain.dat
    run_check_case requested_zero_gain "$log_file" "$dat_file" zero_gain || return 1

    status=0
    require_key_value "$dat_file" stage14_0_requested_flag 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_rhs_injection_enabled_flag 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_gain_parse_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_gain_default_zero_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_gain_finite_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_max_steps_parse_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_require_stage13_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_require_stage13_parse_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_diagnostic_only_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_diagnostic_only_forced_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_rhs_buffer_allocation_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_rhs_injection_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_production_hook_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_fluid_field_access_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_fluid_field_modification_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_pressure_modification_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_projection_modification_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_rk3_modification_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_channel_forcing_modification_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_production_ibm_forcing_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_feedback_application_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_twoway_force_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_no_structure_advance_status 1 requested_zero_gain || status=1
    require_key_value "$dat_file" stage14_0_config_status 1 requested_zero_gain || status=1
    return $status
}

run_invalid_fallback() {
    log_file=$OUTPUT_DIR/stage14_0_config_invalid_fallback.log
    dat_file=$OUTPUT_DIR/fibre_stage14_0_config_invalid_fallback.dat
    run_check_case invalid_fallback "$log_file" "$dat_file" invalid || return 1

    status=0
    require_key_value "$dat_file" stage14_0_gain_parse_status 0 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_safe_fallback_status 1 invalid_fallback || status=1
    require_real_zero "$dat_file" stage14_0_injection_gain invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_max_steps_parse_status 0 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_require_stage13_parse_status 0 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_require_stage13_status 1 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_diagnostic_only_forced_status 1 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_diagnostic_only_status 1 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_no_rhs_injection_status 1 invalid_fallback || status=1
    require_key_value "$dat_file" stage14_0_config_status 1 invalid_fallback || status=1
    return $status
}

GATE_GAIN_PARSE_STATUS=0
GATE_GAIN_DEFAULT_ZERO_STATUS=0
GATE_GAIN_FINITE_STATUS=0
GATE_SAFE_FALLBACK_STATUS=0
GATE_MAX_STEPS_PARSE_STATUS=0
GATE_REQUIRE_STAGE13_STATUS=0
GATE_DIAGNOSTIC_ONLY_STATUS=0
GATE_NO_RHS_BUFFER_ALLOCATION_STATUS=0
GATE_NO_RHS_INJECTION_STATUS=0
GATE_NO_PRODUCTION_HOOK_STATUS=0
GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=0
GATE_NO_PRESSURE_MODIFICATION_STATUS=0
GATE_NO_PROJECTION_MODIFICATION_STATUS=0
GATE_NO_RK3_MODIFICATION_STATUS=0
GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
GATE_NO_PRODUCTION_IBM_FORCING_STATUS=0
GATE_NO_FEEDBACK_APPLICATION_STATUS=0
GATE_NO_TWOWAY_FORCE_STATUS=0
GATE_NO_STRUCTURE_ADVANCE_STATUS=0

if [ "$STAGE14_0_RUN_STAGE13_CLOSURE" = "1" ]; then
    if [ -x stage13_checks/run_stage13_11_total_smoke.sh ]; then
        bash stage13_checks/run_stage13_11_total_smoke.sh || add_reason "stage13_11_total_smoke_failed"
    elif [ -x stage13_checks/run_stage13_11_total_closure.sh ]; then
        bash stage13_checks/run_stage13_11_total_closure.sh || add_reason "stage13_11_total_closure_failed"
    else
        add_reason "stage13_11_closure_script_missing"
    fi
fi

targets="xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check fibre_stage12_sign_convention_audit_check fibre_stage12_power_diagnostics_check fibre_stage12_production_feedback_candidate_check fibre_stage13_config_check fibre_stage13_force_density_buffer_check fibre_stage13_spreading_kernel_check fibre_stage13_volume_normalization_audit_check fibre_stage13_conservation_sign_audit_check fibre_stage13_production_force_density_candidate_check fibre_stage14_config_check"

for target in $targets; do
    build_target "$target" || {
        BUILD_STATUS=0
        add_reason "build_failed_${target}"
        break
    }
done

if [ "$BUILD_STATUS" = "1" ]; then
    run_default_disabled && DEFAULT_STATUS=1
    run_requested_zero_gain && ZERO_GAIN_STATUS=1
    run_invalid_fallback && INVALID_STATUS=1
fi

if [ "$ZERO_GAIN_STATUS" = "1" ]; then
    zero_dat=$OUTPUT_DIR/fibre_stage14_0_config_requested_zero_gain.dat
    GATE_GAIN_PARSE_STATUS=$(get_value stage14_0_gain_parse_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_GAIN_DEFAULT_ZERO_STATUS=$(get_value stage14_0_gain_default_zero_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_GAIN_FINITE_STATUS=$(get_value stage14_0_gain_finite_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_MAX_STEPS_PARSE_STATUS=$(get_value stage14_0_max_steps_parse_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_REQUIRE_STAGE13_STATUS=$(get_value stage14_0_require_stage13_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_DIAGNOSTIC_ONLY_STATUS=$(get_value stage14_0_diagnostic_only_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_RHS_BUFFER_ALLOCATION_STATUS=$(get_value stage14_0_no_rhs_buffer_allocation_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_RHS_INJECTION_STATUS=$(get_value stage14_0_no_rhs_injection_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_HOOK_STATUS=$(get_value stage14_0_no_production_hook_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_FLUID_FIELD_MODIFICATION_STATUS=$(get_value stage14_0_no_fluid_field_modification_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_0_no_pressure_modification_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_0_no_projection_modification_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_RK3_MODIFICATION_STATUS=$(get_value stage14_0_no_rk3_modification_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_0_no_channel_forcing_modification_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_0_no_production_ibm_forcing_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_0_no_feedback_application_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_TWOWAY_FORCE_STATUS=$(get_value stage14_0_no_twoway_force_status "$zero_dat" 2>/dev/null || echo 0)
    GATE_NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_0_no_structure_advance_status "$zero_dat" 2>/dev/null || echo 0)
fi

if [ "$INVALID_STATUS" = "1" ]; then
    invalid_dat=$OUTPUT_DIR/fibre_stage14_0_config_invalid_fallback.dat
    GATE_SAFE_FALLBACK_STATUS=$(get_value stage14_0_safe_fallback_status "$invalid_dat" 2>/dev/null || echo 0)
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$DEFAULT_STATUS" = "1" ] && [ "$ZERO_GAIN_STATUS" = "1" ] && [ "$INVALID_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.0 CONFIG VERDICT: PASS'
    echo 'STAGE 14.0 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.0 CONFIG VERDICT: FAIL'
echo 'STAGE 14.0 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_0_failure'
fi
exit 1
