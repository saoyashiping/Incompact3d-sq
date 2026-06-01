#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE14_7_RUN_STAGE14_6=${STAGE14_7_RUN_STAGE14_6:-0}
STAGE14_7_SMALL_LAMBDA=${STAGE14_7_SMALL_LAMBDA:-1.0e-8}
STAGE14_7_RESPONSE_ABS_TOL=${STAGE14_7_RESPONSE_ABS_TOL:-1.0e-14}
STAGE14_7_RESPONSE_REL_TOL=${STAGE14_7_RESPONSE_REL_TOL:-1.0e-8}
STAGE14_7_MAX_FLUID_DELTA=${STAGE14_7_MAX_FLUID_DELTA:-1.0e-4}
STAGE14_7_MAX_RHS_INCREMENT=${STAGE14_7_MAX_RHS_INCREMENT:-1.0e-4}
OUTPUT_DIR=stage14_outputs
GATE_DAT=$OUTPUT_DIR/stage14_7_small_lambda_response_np1.dat
REASONS_FILE=$OUTPUT_DIR/stage14_7_small_lambda_response_np1_reasons.tmp
BUILD_STATUS=1
LAMBDA0_RUN_STATUS=0
SMALL_LAMBDA_RUN_STATUS=0
HOOK_ACTIVE_STATUS=0
LAMBDA_NONZERO_STATUS=0
RHS_INCREMENT_NONZERO_STATUS=0
RHS_INCREMENT_FINITE_STATUS=0
RHS_INCREMENT_BOUNDED_STATUS=0
FINAL_FLUID_DELTA_FINITE_STATUS=0
FINAL_FLUID_DELTA_BOUNDED_STATUS=0
NO_PRESSURE_MODIFICATION_STATUS=0
NO_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_MODIFICATION_STATUS=0
NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
NO_PRODUCTION_IBM_FORCING_STATUS=0
NO_FEEDBACK_APPLICATION_STATUS=0
NO_TWOWAY_FORCE_STATUS=0
NO_STRUCTURE_ADVANCE_STATUS=0
RHS_INCREMENT_L2=0.0
RHS_INCREMENT_MAX_ABS=0.0
MAX_FLUID_SIGNATURE_DELTA=0.0
SUM_UX_DELTA=0.0
SUM_UY_DELTA=0.0
SUM_UZ_DELTA=0.0
L2_UX_DELTA=0.0
L2_UY_DELTA=0.0
L2_UZ_DELTA=0.0

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

is_finite_number() {
    value=$1
    awk -v value="$value" 'BEGIN { v = value + 0.0; if (v == v && v < 1.0e300 && v > -1.0e300) exit 0; exit 1 }'
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

require_real_gt() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$key" "$file" 2>/dev/null) || {
        add_reason "missing_${key}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value + 0.0) > (limit + 0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_gt_${limit}_got_${value}"
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
stage14_7_requested_flag 1
stage14_7_build_status $BUILD_STATUS
stage14_7_lambda0_run_status $LAMBDA0_RUN_STATUS
stage14_7_small_lambda_run_status $SMALL_LAMBDA_RUN_STATUS
stage14_7_hook_active_status $HOOK_ACTIVE_STATUS
stage14_7_lambda_nonzero_status $LAMBDA_NONZERO_STATUS
stage14_7_rhs_increment_nonzero_status $RHS_INCREMENT_NONZERO_STATUS
stage14_7_rhs_increment_finite_status $RHS_INCREMENT_FINITE_STATUS
stage14_7_rhs_increment_bounded_status $RHS_INCREMENT_BOUNDED_STATUS
stage14_7_final_fluid_delta_finite_status $FINAL_FLUID_DELTA_FINITE_STATUS
stage14_7_final_fluid_delta_bounded_status $FINAL_FLUID_DELTA_BOUNDED_STATUS
stage14_7_no_pressure_modification_status $NO_PRESSURE_MODIFICATION_STATUS
stage14_7_no_projection_modification_status $NO_PROJECTION_MODIFICATION_STATUS
stage14_7_no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
stage14_7_no_rk3_modification_status $NO_RK3_MODIFICATION_STATUS
stage14_7_no_channel_forcing_modification_status $NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_7_no_production_ibm_forcing_status $NO_PRODUCTION_IBM_FORCING_STATUS
stage14_7_no_feedback_application_status $NO_FEEDBACK_APPLICATION_STATUS
stage14_7_no_twoway_force_status $NO_TWOWAY_FORCE_STATUS
stage14_7_no_structure_advance_status $NO_STRUCTURE_ADVANCE_STATUS
stage14_7_small_lambda_response_np1_status $gate_status
stage14_7_small_lambda $STAGE14_7_SMALL_LAMBDA
stage14_7_rhs_increment_l2 $RHS_INCREMENT_L2
stage14_7_rhs_increment_max_abs $RHS_INCREMENT_MAX_ABS
stage14_7_max_fluid_signature_delta $MAX_FLUID_SIGNATURE_DELTA
stage14_7_sum_ux_delta $SUM_UX_DELTA
stage14_7_sum_uy_delta $SUM_UY_DELTA
stage14_7_sum_uz_delta $SUM_UZ_DELTA
stage14_7_l2_ux_delta $L2_UX_DELTA
stage14_7_l2_uy_delta $L2_UY_DELTA
stage14_7_l2_uz_delta $L2_UZ_DELTA
EOF_GATE
}

verify_small_lambda_value() {
    is_finite_number "$STAGE14_7_SMALL_LAMBDA" || {
        add_reason "stage14_7_small_lambda_not_finite"
        return 1
    }
    awk -v value="$STAGE14_7_SMALL_LAMBDA" 'BEGIN { v = value + 0.0; if (v > 0.0 && v <= 1.0e-4) exit 0; exit 1 }' || {
        add_reason "stage14_7_small_lambda_not_small_positive"
        return 1
    }
    return 0
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

run_lambda0() {
    log_file=$OUTPUT_DIR/stage14_7_lambda0.log
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
        add_reason "stage14_7_lambda0_stage9_9_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_lambda0_pass_verdict"
        return 1
    }
    verify_np1_dat stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat || return 1
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
       stage14_outputs/stage14_7_lambda0_np1.dat
    if [ ! -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics"
        return 1
    fi
    cp stage14_outputs/fibre_stage14_5_production_rhs_hook.dat \
       stage14_outputs/stage14_7_lambda0_rhs_hook.dat
    return 0
}

run_small_lambda() {
    log_file=$OUTPUT_DIR/stage14_7_small_lambda.log
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
        X3D_STAGE14_INJECTION_GAIN=$STAGE14_7_SMALL_LAMBDA
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
        add_reason "stage14_7_small_lambda_stage9_9_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_small_lambda_pass_verdict"
        return 1
    }
    verify_np1_dat stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat || return 1
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
       stage14_outputs/stage14_7_small_lambda_np1.dat
    if [ ! -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics"
        add_reason "nonzero_lambda_still_blocked_by_stage14_5"
        return 1
    fi
    cp stage14_outputs/fibre_stage14_5_production_rhs_hook.dat \
       stage14_outputs/stage14_7_small_lambda_rhs_hook.dat
    return 0
}

verify_lambda0_diagnostics() {
    dat_file=stage14_outputs/stage14_7_lambda0_rhs_hook.dat
    status=0
    if [ ! -f "$dat_file" ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics"
        return 1
    fi
    require_key_value "$dat_file" stage14_5_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_5_lambda_zero_status 1 || status=1
    require_key_value "$dat_file" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_5_hook_apply_called_status 1 || status=1
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
    return $status
}

verify_small_lambda_diagnostics() {
    dat_file=stage14_outputs/stage14_7_small_lambda_rhs_hook.dat
    status=0
    if [ ! -f "$dat_file" ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics"
        add_reason "nonzero_lambda_still_blocked_by_stage14_5"
        return 1
    fi
    require_key_value "$dat_file" stage14_5_requested_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$dat_file" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$dat_file" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$dat_file" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$dat_file" stage14_5_stage13_dependency_status 1 || status=1
    require_key_value "$dat_file" stage14_5_stage13_candidate_required_status 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_arrays_available_status 1 || status=1
    require_key_value "$dat_file" stage14_5_rhs_increment_computed_status 1 || status=1
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

    lambda_zero=$(get_value stage14_5_lambda_zero_status "$dat_file" 2>/dev/null || echo missing)
    if [ "$lambda_zero" != "0" ]; then
        add_reason "nonzero_lambda_still_blocked_by_stage14_5"
        status=1
    else
        LAMBDA_NONZERO_STATUS=1
    fi

    injection_gain=$(get_value stage14_5_injection_gain "$dat_file" 2>/dev/null || echo missing)
    if ! awk -v observed="$injection_gain" -v expected="$STAGE14_7_SMALL_LAMBDA" \
            -v abs_tol="$STAGE14_7_RESPONSE_ABS_TOL" -v rel_tol="$STAGE14_7_RESPONSE_REL_TOL" '
        function absval(x) { return x < 0.0 ? -x : x }
        function maxval(a,b) { return a > b ? a : b }
        BEGIN {
          o = observed + 0.0
          e = expected + 0.0
          delta = o - e
          tol = maxval(abs_tol + 0.0, (rel_tol + 0.0) * maxval(1.0, absval(e)))
          if (delta == delta && absval(delta) <= tol) exit 0
          exit 1
        }'; then
        add_reason "stage14_7_injection_gain_not_scaled_to_small_lambda"
        status=1
    fi

    RHS_INCREMENT_L2=$(get_value stage14_5_rhs_increment_l2 "$dat_file" 2>/dev/null || echo missing)
    RHS_INCREMENT_MAX_ABS=$(get_value stage14_5_rhs_increment_max_abs "$dat_file" 2>/dev/null || echo missing)
    if is_finite_number "$RHS_INCREMENT_L2" && is_finite_number "$RHS_INCREMENT_MAX_ABS"; then
        RHS_INCREMENT_FINITE_STATUS=1
    else
        add_reason "stage14_7_rhs_increment_not_finite"
        status=1
    fi
    if require_real_gt "$dat_file" stage14_5_rhs_increment_l2 0.0 && \
       require_real_gt "$dat_file" stage14_5_rhs_increment_max_abs 0.0; then
        RHS_INCREMENT_NONZERO_STATUS=1
    else
        add_reason "nonzero_lambda_still_blocked_by_stage14_5"
        status=1
    fi
    if require_real_le "$dat_file" stage14_5_rhs_increment_max_abs "$STAGE14_7_MAX_RHS_INCREMENT"; then
        RHS_INCREMENT_BOUNDED_STATUS=1
    else
        status=1
    fi

    HOOK_ACTIVE_STATUS=$(get_value stage14_5_hook_apply_called_status "$dat_file" 2>/dev/null || echo 0)
    NO_PRESSURE_MODIFICATION_STATUS=$(get_value stage14_5_no_pressure_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_PROJECTION_MODIFICATION_STATUS=$(get_value stage14_5_no_projection_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_POISSON_MODIFICATION_STATUS=$(get_value stage14_5_no_poisson_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_RK3_MODIFICATION_STATUS=$(get_value stage14_5_no_rk3_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value stage14_5_no_channel_forcing_modification_status "$dat_file" 2>/dev/null || echo 0)
    NO_PRODUCTION_IBM_FORCING_STATUS=$(get_value stage14_5_no_production_ibm_forcing_status "$dat_file" 2>/dev/null || echo 0)
    NO_FEEDBACK_APPLICATION_STATUS=$(get_value stage14_5_no_feedback_application_status "$dat_file" 2>/dev/null || echo 0)
    NO_TWOWAY_FORCE_STATUS=$(get_value stage14_5_no_twoway_force_status "$dat_file" 2>/dev/null || echo 0)
    NO_STRUCTURE_ADVANCE_STATUS=$(get_value stage14_5_no_structure_advance_status "$dat_file" 2>/dev/null || echo 0)
    return $status
}

compare_final_signatures() {
    lambda0_dat=stage14_outputs/stage14_7_lambda0_np1.dat
    small_dat=stage14_outputs/stage14_7_small_lambda_np1.dat
    status=0
    max_delta=0.0
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        lambda0=$(get_value "$metric" "$lambda0_dat" 2>/dev/null) || {
            add_reason "missing_lambda0_${metric}"
            status=1
            continue
        }
        small=$(get_value "$metric" "$small_dat" 2>/dev/null) || {
            add_reason "missing_small_lambda_${metric}"
            status=1
            continue
        }
        comparison=$(awk -v metric="$metric" -v lambda0="$lambda0" -v small="$small" \
                         -v max_delta_limit="$STAGE14_7_MAX_FLUID_DELTA" '
            function absval(x) { return x < 0.0 ? -x : x }
            BEGIN {
              z = lambda0 + 0.0
              s = small + 0.0
              delta = s - z
              abs_delta = absval(delta)
              finite = (delta == delta && abs_delta < 1.0e300)
              verdict = (finite && abs_delta <= (max_delta_limit + 0.0)) ? "PASS" : "FAIL"
              printf "%s %.17e %.17e %.17e %.17e %.17e %s", metric, z, s, delta, abs_delta, max_delta_limit + 0.0, verdict
              exit verdict == "PASS" ? 0 : 1
            }')
        compare_status=$?
        set -- $comparison
        printf 'metric %s lambda0 %s small_lambda %s delta %s abs_delta %s bounded_tolerance %s %s\n' \
               "$1" "$2" "$3" "$4" "$5" "$6" "$7"
        if ! is_finite_number "$4"; then
            status=1
        fi
        max_delta=$(awk -v current="$max_delta" -v candidate="$5" 'BEGIN { print ((candidate + 0.0) > (current + 0.0)) ? candidate : current }')
        case "$metric" in
            stage9_9_signature_sum_ux) SUM_UX_DELTA=$4 ;;
            stage9_9_signature_sum_uy) SUM_UY_DELTA=$4 ;;
            stage9_9_signature_sum_uz) SUM_UZ_DELTA=$4 ;;
            stage9_9_signature_l2_ux) L2_UX_DELTA=$4 ;;
            stage9_9_signature_l2_uy) L2_UY_DELTA=$4 ;;
            stage9_9_signature_l2_uz) L2_UZ_DELTA=$4 ;;
        esac
        if [ "$compare_status" -ne 0 ]; then
            add_reason "final_fluid_delta_bounded_failed_${metric}"
            status=1
        fi
    done
    MAX_FLUID_SIGNATURE_DELTA=$max_delta
    if [ "$status" = "0" ]; then
        FINAL_FLUID_DELTA_FINITE_STATUS=1
        FINAL_FLUID_DELTA_BOUNDED_STATUS=1
    else
        is_finite_number "$MAX_FLUID_SIGNATURE_DELTA" && FINAL_FLUID_DELTA_FINITE_STATUS=1
    fi
    return $status
}

if [ "$STAGE14_7_RUN_STAGE14_6" = "1" ]; then
    bash stage14_checks/run_stage14_6_lambda_zero_invariance.sh || {
        BUILD_STATUS=0
        add_reason "stage14_6_lambda_zero_invariance_failed"
    }
fi

verify_small_lambda_value || BUILD_STATUS=0

targets="xcompact3d fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check fibre_stage12_config_check fibre_stage12_force_buffer_check fibre_stage12_prescribed_velocity_check fibre_stage12_feedback_formula_check fibre_stage12_sign_convention_audit_check fibre_stage12_power_diagnostics_check fibre_stage12_production_feedback_candidate_check fibre_stage13_config_check fibre_stage13_force_density_buffer_check fibre_stage13_spreading_kernel_check fibre_stage13_volume_normalization_audit_check fibre_stage13_conservation_sign_audit_check fibre_stage13_production_force_density_candidate_check fibre_stage14_config_check fibre_stage14_rhs_accumulator_check fibre_stage14_rhs_addition_formula_check fibre_stage14_rhs_sign_scaling_audit_check fibre_stage14_rk_timing_audit_check fibre_stage14_production_rhs_injection_check"

if [ "$BUILD_STATUS" = "1" ]; then
    for target in $targets; do
        build_target "$target" || {
            BUILD_STATUS=0
            add_reason "build_failed_${target}"
            break
        }
    done
fi

if [ "$BUILD_STATUS" = "1" ]; then
    run_lambda0 && LAMBDA0_RUN_STATUS=1
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$LAMBDA0_RUN_STATUS" = "1" ]; then
    run_small_lambda && SMALL_LAMBDA_RUN_STATUS=1
fi

if [ "$LAMBDA0_RUN_STATUS" = "1" ]; then
    verify_lambda0_diagnostics || LAMBDA0_RUN_STATUS=0
fi

if [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ]; then
    verify_small_lambda_diagnostics || SMALL_LAMBDA_RUN_STATUS=0
fi

if [ "$LAMBDA0_RUN_STATUS" = "1" ] && [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ]; then
    compare_final_signatures || SMALL_LAMBDA_RUN_STATUS=0
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$LAMBDA0_RUN_STATUS" = "1" ] && \
   [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ] && [ "$HOOK_ACTIVE_STATUS" = "1" ] && \
   [ "$LAMBDA_NONZERO_STATUS" = "1" ] && [ "$RHS_INCREMENT_NONZERO_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_FINITE_STATUS" = "1" ] && [ "$RHS_INCREMENT_BOUNDED_STATUS" = "1" ] && \
   [ "$FINAL_FLUID_DELTA_FINITE_STATUS" = "1" ] && [ "$FINAL_FLUID_DELTA_BOUNDED_STATUS" = "1" ]; then
    write_gate_dat 1
    echo 'STAGE 14.7 SMALL LAMBDA RESPONSE NP1 VERDICT: PASS'
    echo 'STAGE 14.7 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_gate_dat 0
echo 'STAGE 14.7 SMALL LAMBDA RESPONSE NP1 VERDICT: FAIL'
echo 'STAGE 14.7 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_7_failure'
fi
exit 1
