#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_5_RUN_STAGE14_4=${STAGE14_5_RUN_STAGE14_4:-0}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_5_production_rhs_hook_gate.dat
REASONS_FILE=$OUTPUT_DIR/stage14_5_production_rhs_hook_reasons.tmp
BUILD_STATUS=1
STANDALONE_CHECK_STATUS=0
PRODUCTION_SMOKE_STATUS=0

mkdir -p stage14_outputs stage13_outputs stage12_outputs stage11_outputs stage9_outputs
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
    if [ -x "$BUILD_DIR/bin/fibre_stage14_production_rhs_injection_check" ]; then
        echo "$BUILD_DIR/bin/fibre_stage14_production_rhs_injection_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/src/fibre_stage14_production_rhs_injection_check" ]; then
        echo "$BUILD_DIR/src/fibre_stage14_production_rhs_injection_check"
        return 0
    fi
    if [ -x "$BUILD_DIR/fibre_stage14_production_rhs_injection_check" ]; then
        echo "$BUILD_DIR/fibre_stage14_production_rhs_injection_check"
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
stage14_5_gate_requested_flag 1
stage14_5_gate_build_status $BUILD_STATUS
stage14_5_gate_standalone_check_status $STANDALONE_CHECK_STATUS
stage14_5_gate_production_smoke_status $PRODUCTION_SMOKE_STATUS
stage14_5_gate_hook_active_status $GATE_HOOK_ACTIVE_STATUS
stage14_5_gate_lambda_zero_status $GATE_LAMBDA_ZERO_STATUS
stage14_5_gate_nonzero_lambda_blocked_status $GATE_NONZERO_LAMBDA_BLOCKED_STATUS
stage14_5_gate_stage13_dependency_status $GATE_STAGE13_DEPENDENCY_STATUS
stage14_5_gate_rhs_arrays_available_status $GATE_RHS_ARRAYS_AVAILABLE_STATUS
stage14_5_gate_rhs_increment_computed_status $GATE_RHS_INCREMENT_COMPUTED_STATUS
stage14_5_gate_rhs_increment_zero_status $GATE_RHS_INCREMENT_ZERO_STATUS
stage14_5_gate_rhs_unchanged_status $GATE_RHS_UNCHANGED_STATUS
stage14_5_gate_no_pressure_modification_status $GATE_NO_PRESSURE_MODIFICATION_STATUS
stage14_5_gate_no_projection_modification_status $GATE_NO_PROJECTION_MODIFICATION_STATUS
stage14_5_gate_no_poisson_modification_status $GATE_NO_POISSON_MODIFICATION_STATUS
stage14_5_gate_no_rk3_modification_status $GATE_NO_RK3_MODIFICATION_STATUS
stage14_5_gate_no_channel_forcing_modification_status $GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_5_gate_no_production_ibm_forcing_status $GATE_NO_PRODUCTION_IBM_FORCING_STATUS
stage14_5_gate_no_feedback_application_status $GATE_NO_FEEDBACK_APPLICATION_STATUS
stage14_5_gate_no_twoway_force_status $GATE_NO_TWOWAY_FORCE_STATUS
stage14_5_gate_no_structure_advance_status $GATE_NO_STRUCTURE_ADVANCE_STATUS
stage14_5_gate_status $gate_status
EOF_GATE
}

run_standalone_check() {
    log_file=$OUTPUT_DIR/stage14_5_production_rhs_hook_check.log
    dat_file=$OUTPUT_DIR/fibre_stage14_5_production_rhs_hook_check.dat
    rm -f "$dat_file"

    check_exe=$(find_check_exe) || {
        add_reason "missing_fibre_stage14_production_rhs_injection_check_executable"
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
        add_reason "fibre_stage14_production_rhs_injection_check_failed"
        cat "$log_file"
        return 1
    }

    grep 'STAGE 14.5 PRODUCTION RHS HOOK CHECK VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage14_5_check_pass_verdict"
        return 1
    }

    if [ ! -f "$dat_file" ]; then
        add_reason "missing_fibre_stage14_5_production_rhs_hook_check_dat"
        return 1
    fi

    status=0
    require_key_value "$dat_file" stage14_5_check_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_check_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_check_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_lambda_zero_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_hook_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_hook_apply_called_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_rhs_unchanged_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_rhs_increment_zero_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_pressure_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_projection_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_poisson_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_rk3_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_feedback_application_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_twoway_force_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_no_structure_advance_status 1 || status=1
    require_key_value "$dat_file" stage14_5_check_production_rhs_hook_status 1 || status=1
    require_real_le "$dat_file" stage14_5_check_rhs_signature_delta_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_5_check_rhs_increment_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_5_check_rhs_increment_max_abs 1.0e-12 || status=1
    return $status
}

run_production_smoke() {
    log_file=$OUTPUT_DIR/stage14_5_stage9_9_production_smoke.log
    (
        X3D_STAGE11_ONEWAY_HOOK=1
        X3D_STAGE11_FORCE_READONLY=1
        X3D_STAGE11_MAX_POINTS=8
        X3D_STAGE11_MAX_STEPS=3
        X3D_STAGE12_FEEDBACK_CANDIDATE=1
        X3D_STAGE12_FORCE_READONLY=1
        X3D_STAGE12_FEEDBACK_GAIN=1.0
        X3D_STAGE12_FORCE_SIGN=1
        X3D_STAGE12_MAX_POINTS=8
        X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1
        X3D_STAGE13_FORCE_READONLY=1
        X3D_STAGE13_SPREADING_READONLY=1
        X3D_STAGE13_MAX_POINTS=8
        X3D_STAGE13_MAX_EULERIAN_POINTS=64
        X3D_STAGE13_SPREADING_NORMALIZATION=conservative
        X3D_STAGE14_RHS_INJECTION=1
        X3D_STAGE14_INJECTION_GAIN=0.0
        X3D_STAGE14_MAX_STEPS=3
        X3D_STAGE14_REQUIRE_STAGE13=1
        X3D_STAGE14_DIAGNOSTIC_ONLY=1
        STAGE9_SKIP_PREREQS=1
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1
        X3D_STAGE9_9_DETERMINISTIC_INIT=1
        X3D_STAGE9_9_MAX_STEPS=3
        export X3D_STAGE11_ONEWAY_HOOK X3D_STAGE11_FORCE_READONLY X3D_STAGE11_MAX_POINTS X3D_STAGE11_MAX_STEPS
        export X3D_STAGE12_FEEDBACK_CANDIDATE X3D_STAGE12_FORCE_READONLY X3D_STAGE12_FEEDBACK_GAIN
        export X3D_STAGE12_FORCE_SIGN X3D_STAGE12_MAX_POINTS X3D_STAGE13_FORCE_DENSITY_CANDIDATE
        export X3D_STAGE13_FORCE_READONLY X3D_STAGE13_SPREADING_READONLY X3D_STAGE13_MAX_POINTS
        export X3D_STAGE13_MAX_EULERIAN_POINTS X3D_STAGE13_SPREADING_NORMALIZATION X3D_STAGE14_RHS_INJECTION
        export X3D_STAGE14_INJECTION_GAIN X3D_STAGE14_MAX_STEPS X3D_STAGE14_REQUIRE_STAGE13
        export X3D_STAGE14_DIAGNOSTIC_ONLY STAGE9_SKIP_PREREQS X3D_STAGE9_9_PARALLEL_CONSISTENCY
        export X3D_STAGE9_9_DETERMINISTIC_INIT X3D_STAGE9_9_MAX_STEPS
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
    ) > "$log_file" 2>&1 || {
        add_reason "stage9_9_production_smoke_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_pass_verdict"
        return 1
    }
    return 0
}

verify_production_dat() {
    dat_file=$OUTPUT_DIR/fibre_stage14_5_production_rhs_hook.dat
    if [ ! -f "$dat_file" ]; then
        add_reason "missing_fibre_stage14_5_production_rhs_hook_dat"
        return 1
    fi
    status=0
    require_key_value "$dat_file" stage14_5_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_5_lambda_zero_status 1 || status=1
    require_key_value "$dat_file" stage14_5_nonzero_lambda_blocked_status 1 || status=1
    require_key_value "$dat_file" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$dat_file" stage14_5_stage13_dependency_status 1 || status=1
    require_key_value "$dat_file" stage14_5_stage13_candidate_required_status 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_arrays_available_status 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_increment_computed_status 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_increment_zero_status 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_unchanged_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_feedback_application_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_twoway_force_status 1 || status=1
    require_key_value "$dat_file" stage14_5_no_structure_advance_status 1 || status=1
    require_key_value "$dat_file" stage14_5_production_rhs_hook_status 1 || status=1
    require_real_le "$dat_file" stage14_5_rhs_signature_delta_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_5_rhs_increment_l2 1.0e-12 || status=1
    require_real_le "$dat_file" stage14_5_rhs_increment_max_abs 1.0e-12 || status=1
    require_real_ge "$dat_file" stage14_5_apply_call_count 1 || status=1
    return $status
}

GATE_HOOK_ACTIVE_STATUS=0
GATE_LAMBDA_ZERO_STATUS=0
GATE_NONZERO_LAMBDA_BLOCKED_STATUS=0
GATE_STAGE13_DEPENDENCY_STATUS=0
GATE_RHS_ARRAYS_AVAILABLE_STATUS=0
GATE_RHS_INCREMENT_COMPUTED_STATUS=0
GATE_RHS_INCREMENT_ZERO_STATUS=0
GATE_RHS_UNCHANGED_STATUS=0
GATE_NO_PRESSURE_MODIFICATION_STATUS=0
GATE_NO_PROJECTION_MODIFICATION_STATUS=0
GATE_NO_POISSON_MODIFICATION_STATUS=0
GATE_NO_RK3_MODIFICATION_STATUS=0
GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
GATE_NO_PRODUCTION_IBM_FORCING_STATUS=0
GATE_NO_FEEDBACK_APPLICATION_STATUS=0
GATE_NO_TWOWAY_FORCE_STATUS=0
GATE_NO_STRUCTURE_ADVANCE_STATUS=0

if [ "$STAGE14_5_RUN_STAGE14_4" = "1" ]; then
    bash stage14_checks/run_stage14_4_rk_timing_audit.sh || {
        BUILD_STATUS=0
        add_reason "stage14_4_rk_timing_audit_failed"
    }
fi

targets="xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check fibre_stage12_sign_convention_audit_check fibre_stage12_power_diagnostics_check fibre_stage12_production_feedback_candidate_check fibre_stage13_config_check fibre_stage13_force_density_buffer_check fibre_stage13_spreading_kernel_check fibre_stage13_volume_normalization_audit_check fibre_stage13_conservation_sign_audit_check fibre_stage13_production_force_density_candidate_check fibre_stage14_config_check fibre_stage14_rhs_accumulator_check fibre_stage14_rhs_addition_formula_check fibre_stage14_rhs_sign_scaling_audit_check fibre_stage14_rk_timing_audit_check fibre_stage14_production_rhs_injection_check"

for target in $targets; do
    build_target "$target" || {
        BUILD_STATUS=0
        add_reason "build_failed_${target}"
        break
    }
done

if [ "$BUILD_STATUS" = "1" ]; then
    run_standalone_check && STANDALONE_CHECK_STATUS=1
fi

if [ "$BUILD_STATUS" = "1" ]; then
    run_production_smoke && PRODUCTION_SMOKE_STATUS=1
fi

if [ "$PRODUCTION_SMOKE_STATUS" = "1" ]; then
    verify_production_dat || PRODUCTION_SMOKE_STATUS=0
fi

if [ "$PRODUCTION_SMOKE_STATUS" = "1" ]; then
    dat_file=$OUTPUT_DIR/fibre_stage14_5_production_rhs_hook.dat
    GATE_HOOK_ACTIVE_STATUS=$(get_value stage14_5_hook_apply_called_status "$dat_file" 2>/dev/null || echo 0)
    GATE_LAMBDA_ZERO_STATUS=$(get_value stage14_5_lambda_zero_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NONZERO_LAMBDA_BLOCKED_STATUS=$(get_value stage14_5_nonzero_lambda_blocked_status "$dat_file" 2>/dev/null || echo 0)
    GATE_STAGE13_DEPENDENCY_STATUS=$(get_value stage14_5_stage13_dependency_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_ARRAYS_AVAILABLE_STATUS=$(get_value stage14_5_rhs_arrays_available_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_INCREMENT_COMPUTED_STATUS=$(get_value stage14_5_rhs_increment_computed_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_INCREMENT_ZERO_STATUS=$(get_value stage14_5_rhs_increment_zero_status "$dat_file" 2>/dev/null || echo 0)
    GATE_RHS_UNCHANGED_STATUS=$(get_value stage14_5_rhs_unchanged_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_5_no_pressure_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_5_no_projection_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_POISSON_MODIFICATION_STATUS=$(get_value stage14_5_no_poisson_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_RK3_MODIFICATION_STATUS=$(get_value stage14_5_no_rk3_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_5_no_channel_forcing_modification_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_5_no_production_ibm_forcing_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_5_no_feedback_application_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_TWOWAY_FORCE_STATUS=$(get_value stage14_5_no_twoway_force_status "$dat_file" 2>/dev/null || echo 0)
    GATE_NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_5_no_structure_advance_status "$dat_file" 2>/dev/null || echo 0)
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STANDALONE_CHECK_STATUS" = "1" ] && [ "$PRODUCTION_SMOKE_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.5 PRODUCTION RHS HOOK VERDICT: PASS'
    echo 'STAGE 14.5 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.5 PRODUCTION RHS HOOK VERDICT: FAIL'
echo 'STAGE 14.5 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_5_failure'
fi
exit 1
