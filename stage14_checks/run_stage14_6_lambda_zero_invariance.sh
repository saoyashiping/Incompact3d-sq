#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_6_RUN_STAGE14_5=${STAGE14_6_RUN_STAGE14_5:-0}
STAGE14_6_INVARIANCE_ABS_TOL=${STAGE14_6_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE14_6_INVARIANCE_REL_TOL=${STAGE14_6_INVARIANCE_REL_TOL:-1.0e-14}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_6_lambda_zero_invariance.dat
REASONS_FILE=$OUTPUT_DIR/stage14_6_lambda_zero_invariance_reasons.tmp
BUILD_STATUS=1
BASELINE_RUN_STATUS=0
CANDIDATE_RUN_STATUS=0
NP1_SIGNATURE_INVARIANCE_STATUS=0
HOOK_ACTIVE_STATUS=0
LAMBDA_ZERO_STATUS=0
RHS_INCREMENT_ZERO_STATUS=0
RHS_UNCHANGED_STATUS=0
NO_PRESSURE_MODIFICATION_STATUS=0
NO_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_MODIFICATION_STATUS=0
NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
NO_PRODUCTION_IBM_FORCING_STATUS=0
NO_FEEDBACK_APPLICATION_STATUS=0
NO_TWOWAY_FORCE_STATUS=0
NO_STRUCTURE_ADVANCE_STATUS=0

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

build_target() {
    target=$1
    ensure_build_dir || return 1
    cmake --build "$BUILD_DIR" --target "$target" -j
}

write_gate_dat() {
    gate_status=$1
    cat > "$GATE_DAT" <<EOF_GATE
stage14_6_requested_flag 1
stage14_6_build_status $BUILD_STATUS
stage14_6_baseline_run_status $BASELINE_RUN_STATUS
stage14_6_candidate_run_status $CANDIDATE_RUN_STATUS
stage14_6_hook_active_status $HOOK_ACTIVE_STATUS
stage14_6_lambda_zero_status $LAMBDA_ZERO_STATUS
stage14_6_rhs_increment_zero_status $RHS_INCREMENT_ZERO_STATUS
stage14_6_rhs_unchanged_status $RHS_UNCHANGED_STATUS
stage14_6_np1_signature_invariance_status $NP1_SIGNATURE_INVARIANCE_STATUS
stage14_6_no_pressure_modification_status $NO_PRESSURE_MODIFICATION_STATUS
stage14_6_no_projection_modification_status $NO_PROJECTION_MODIFICATION_STATUS
stage14_6_no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
stage14_6_no_rk3_modification_status $NO_RK3_MODIFICATION_STATUS
stage14_6_no_channel_forcing_modification_status $NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_6_no_production_ibm_forcing_status $NO_PRODUCTION_IBM_FORCING_STATUS
stage14_6_no_feedback_application_status $NO_FEEDBACK_APPLICATION_STATUS
stage14_6_no_twoway_force_status $NO_TWOWAY_FORCE_STATUS
stage14_6_no_structure_advance_status $NO_STRUCTURE_ADVANCE_STATUS
stage14_6_lambda_zero_invariance_status $gate_status
EOF_GATE
}

verify_np1_dat() {
    dat_file=$1
    status=0
    if [ ! -f "$dat_file" ]; then
        add_reason "missing_${dat_file}"
        return 1
    fi
    require_key_value "$dat_file" stage9_9_parallel_consistency_local_status 1 || status=1
    require_key_value "$dat_file" stage9_9_decomposition_invariant_initial_state_status 1 || status=1
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        get_value "$metric" "$dat_file" >/dev/null 2>&1 || {
            add_reason "missing_${metric}_in_${dat_file}"
            status=1
        }
    done
    return $status
}

run_baseline() {
    log_file=$OUTPUT_DIR/stage14_6_baseline_stage13.log
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat
    (
        unset X3D_STAGE14_RHS_INJECTION
        unset X3D_STAGE14_INJECTION_GAIN
        unset X3D_STAGE14_MAX_STEPS
        unset X3D_STAGE14_REQUIRE_STAGE13
        unset X3D_STAGE14_DIAGNOSTIC_ONLY
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
        STAGE9_SKIP_PREREQS=1
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1
        X3D_STAGE9_9_DETERMINISTIC_INIT=1
        X3D_STAGE9_9_MAX_STEPS=3
        export X3D_STAGE11_ONEWAY_HOOK X3D_STAGE11_FORCE_READONLY X3D_STAGE11_MAX_POINTS X3D_STAGE11_MAX_STEPS
        export X3D_STAGE12_FEEDBACK_CANDIDATE X3D_STAGE12_FORCE_READONLY X3D_STAGE12_FEEDBACK_GAIN
        export X3D_STAGE12_FORCE_SIGN X3D_STAGE12_MAX_POINTS X3D_STAGE13_FORCE_DENSITY_CANDIDATE
        export X3D_STAGE13_FORCE_READONLY X3D_STAGE13_SPREADING_READONLY X3D_STAGE13_MAX_POINTS
        export X3D_STAGE13_MAX_EULERIAN_POINTS X3D_STAGE13_SPREADING_NORMALIZATION STAGE9_SKIP_PREREQS
        export X3D_STAGE9_9_PARALLEL_CONSISTENCY X3D_STAGE9_9_DETERMINISTIC_INIT X3D_STAGE9_9_MAX_STEPS
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
    ) > "$log_file" 2>&1 || {
        add_reason "stage14_6_baseline_stage9_9_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_baseline_pass_verdict"
        return 1
    }
    verify_np1_dat stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat || return 1
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
       stage14_outputs/stage14_6_baseline_stage13_np1.dat
    return 0
}

run_candidate() {
    log_file=$OUTPUT_DIR/stage14_6_candidate_lambda0.log
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat
    rm -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
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
        add_reason "stage14_6_candidate_stage9_9_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_candidate_pass_verdict"
        return 1
    }
    verify_np1_dat stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat || return 1
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
       stage14_outputs/stage14_6_candidate_lambda0_np1.dat
    return 0
}

compare_signatures() {
    baseline_dat=stage14_outputs/stage14_6_baseline_stage13_np1.dat
    candidate_dat=stage14_outputs/stage14_6_candidate_lambda0_np1.dat
    status=0
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        baseline=$(get_value "$metric" "$baseline_dat" 2>/dev/null) || {
            add_reason "missing_baseline_${metric}"
            status=1
            continue
        }
        candidate=$(get_value "$metric" "$candidate_dat" 2>/dev/null) || {
            add_reason "missing_candidate_${metric}"
            status=1
            continue
        }
        comparison=$(awk -v metric="$metric" -v baseline="$baseline" -v candidate="$candidate" \
                         -v abs_tol="$STAGE14_6_INVARIANCE_ABS_TOL" \
                         -v rel_tol="$STAGE14_6_INVARIANCE_REL_TOL" '
            function absval(x) { return x < 0.0 ? -x : x }
            function maxval(a,b) { return a > b ? a : b }
            BEGIN {
              b = baseline + 0.0
              c = candidate + 0.0
              delta = c - b
              eff_tol = maxval(abs_tol + 0.0, (rel_tol + 0.0) * maxval(1.0, absval(b)))
              verdict = (absval(delta) <= eff_tol) ? "PASS" : "FAIL"
              printf "%s %.17e %.17e %.17e %.17e %.17e %.17e %s", metric, b, c, delta, abs_tol + 0.0, rel_tol + 0.0, eff_tol, verdict
              exit verdict == "PASS" ? 0 : 1
            }')
        compare_status=$?
        set -- $comparison
        printf 'metric %s baseline %s candidate %s delta %s abs_tol %s rel_tol %s effective_tolerance %s %s\n' \
               "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
        if [ "$compare_status" -ne 0 ]; then
            add_reason "signature_invariance_failed_${metric}"
            status=1
        fi
    done
    return $status
}

verify_stage14_diagnostics() {
    dat_file=stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
    status=0
    if [ ! -f "$dat_file" ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics"
        return 1
    fi
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
    if [ "$status" = "1" ]; then
        return 1
    fi
    HOOK_ACTIVE_STATUS=$(get_value stage14_5_hook_apply_called_status "$dat_file" 2>/dev/null || echo 0)
    LAMBDA_ZERO_STATUS=$(get_value stage14_5_lambda_zero_status "$dat_file" 2>/dev/null || echo 0)
    RHS_INCREMENT_ZERO_STATUS=$(get_value stage14_5_rhs_increment_zero_status "$dat_file" 2>/dev/null || echo 0)
    RHS_UNCHANGED_STATUS=$(get_value stage14_5_rhs_unchanged_status "$dat_file" 2>/dev/null || echo 0)
    NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_5_no_pressure_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_5_no_projection_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_POISSON_MODIFICATION_STATUS=$(get_value stage14_5_no_poisson_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_RK3_MODIFICATION_STATUS=$(get_value stage14_5_no_rk3_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_5_no_channel_forcing_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_5_no_production_ibm_forcing_status "$dat_file" 2>/dev/null || echo 0)
    NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_5_no_feedback_application_status "$dat_file" 2>/dev/null || echo 0)
    NO_TWOWAY_FORCE_STATUS=$(get_value stage14_5_no_twoway_force_status "$dat_file" 2>/dev/null || echo 0)
    NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_5_no_structure_advance_status "$dat_file" 2>/dev/null || echo 0)
    return 0
}

if [ "$STAGE14_6_RUN_STAGE14_5" = "1" ]; then
    bash stage14_checks/run_stage14_5_production_rhs_hook.sh || {
        BUILD_STATUS=0
        add_reason "stage14_5_production_rhs_hook_failed"
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
    run_baseline && BASELINE_RUN_STATUS=1
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$BASELINE_RUN_STATUS" = "1" ]; then
    run_candidate && CANDIDATE_RUN_STATUS=1
fi

if [ "$BASELINE_RUN_STATUS" = "1" ] && [ "$CANDIDATE_RUN_STATUS" = "1" ]; then
    compare_signatures && NP1_SIGNATURE_INVARIANCE_STATUS=1
fi

if [ "$CANDIDATE_RUN_STATUS" = "1" ]; then
    verify_stage14_diagnostics || CANDIDATE_RUN_STATUS=0
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$BASELINE_RUN_STATUS" = "1" ] && \
   [ "$CANDIDATE_RUN_STATUS" = "1" ] && [ "$NP1_SIGNATURE_INVARIANCE_STATUS" = "1" ] && \
   [ "$HOOK_ACTIVE_STATUS" = "1" ] && [ "$LAMBDA_ZERO_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_ZERO_STATUS" = "1" ] && [ "$RHS_UNCHANGED_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.6 LAMBDA ZERO INVARIANCE VERDICT: PASS'
    echo 'STAGE 14.6 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.6 LAMBDA ZERO INVARIANCE VERDICT: FAIL'
echo 'STAGE 14.6 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_6_failure'
fi
exit 1
