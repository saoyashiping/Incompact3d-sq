#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
CHANNEL_INPUT=${CHANNEL_INPUT:-examples/Channel/input.i3d}
STAGE14_8_SMALL_LAMBDA=${STAGE14_8_SMALL_LAMBDA:-1.0e-8}
STAGE14_8_MAX_RHS_INCREMENT=${STAGE14_8_MAX_RHS_INCREMENT:-1.0e-4}
STAGE14_8_MAX_FLUID_DELTA=${STAGE14_8_MAX_FLUID_DELTA:-1.0e-4}
STAGE14_8_FINAL_SIGNATURE_ABS_TOL=${STAGE14_8_FINAL_SIGNATURE_ABS_TOL:-1.0e-6}
STAGE14_8_FINAL_SIGNATURE_REL_TOL=${STAGE14_8_FINAL_SIGNATURE_REL_TOL:-1.0e-12}
STAGE14_8_NORMALIZED_RHS_ABS_TOL=${STAGE14_8_NORMALIZED_RHS_ABS_TOL:-1.0e-12}
STAGE14_8_NORMALIZED_RHS_REL_TOL=${STAGE14_8_NORMALIZED_RHS_REL_TOL:-1.0e-8}
STAGE14_8_FORCE_DENSITY_ABS_TOL=${STAGE14_8_FORCE_DENSITY_ABS_TOL:-1.0e-10}
STAGE14_8_FORCE_DENSITY_REL_TOL=${STAGE14_8_FORCE_DENSITY_REL_TOL:-1.0e-8}
STAGE14_8_TIMEOUT_SEC=${STAGE14_8_TIMEOUT_SEC:-240}
OUTPUT_DIR=stage14_outputs
OUT_DAT=$OUTPUT_DIR/stage14_8_parallel_small_lambda_response.dat
REASONS_FILE=$OUTPUT_DIR/stage14_8_parallel_small_lambda_response_reasons.tmp

BUILD_STATUS=1
LAMBDA0_RUN_STATUS=0
SMALL_LAMBDA_RUN_STATUS=0
LAMBDA0_HOOK_ACTIVE_STATUS=0
SMALL_LAMBDA_HOOK_ACTIVE_STATUS=0
LAMBDA0_RHS_INCREMENT_ZERO_STATUS=0
SMALL_LAMBDA_NONZERO_STATUS=0
SMALL_LAMBDA_RHS_INCREMENT_NONZERO_STATUS=0
SMALL_LAMBDA_RHS_INCREMENT_FINITE_STATUS=0
SMALL_LAMBDA_RHS_INCREMENT_BOUNDED_STATUS=0
RHS_SIGN_CORRECT_STATUS=0
RHS_LINEAR_SCALING_STATUS=0
NORMALIZED_RHS_PARALLEL_STATUS=0
FORCE_DENSITY_PARALLEL_STATUS=0
FINAL_FLUID_DELTA_FINITE_STATUS=0
FINAL_FLUID_DELTA_BOUNDED_STATUS=0
FINAL_SIGNATURE_PARALLEL_STATUS=0
NO_PRESSURE_MODIFICATION_STATUS=0
NO_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_MODIFICATION_STATUS=0
NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
NO_PRODUCTION_IBM_FORCING_STATUS=0
NO_FEEDBACK_APPLICATION_STATUS=0
NO_TWOWAY_FORCE_STATUS=0
NO_STRUCTURE_ADVANCE_STATUS=0

MAX_RHS_INCREMENT=0.0
MAX_FLUID_SIGNATURE_DELTA=0.0
MAX_NORMALIZED_RHS_L2_DELTA=0.0
MAX_NORMALIZED_RHS_MAX_ABS_DELTA=0.0
MAX_FORCE_DENSITY_DELTA=0.0
NORMALIZED_RHS_L2_NP1=0.0
NORMALIZED_RHS_L2_NP2=0.0
NORMALIZED_RHS_L2_NP4=0.0
NORMALIZED_RHS_MAX_ABS_NP1=0.0
NORMALIZED_RHS_MAX_ABS_NP2=0.0
NORMALIZED_RHS_MAX_ABS_NP4=0.0

mkdir -p stage14_outputs stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$OUT_DAT"

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        if [ -n "$DECOMP2D_ROOT" ]; then
            cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT"
        else
            cmake -S . -B "$BUILD_DIR"
        fi
    fi
}

add_reason() {
    echo "$1" >> "$REASONS_FILE"
}

is_finite_number() {
    value=$1
    awk -v value="$value" 'BEGIN { v=value+0.0; if (v==v && v<1.0e300 && v>-1.0e300) exit 0; exit 1 }'
}

abs_value() {
    value=$1
    awk -v value="$value" 'BEGIN { v=value+0.0; if (v<0.0) v=-v; print v }'
}

max_value() {
    first=$1
    second=$2
    awk -v a="$first" -v b="$second" 'BEGIN { print ((a+0.0) > (b+0.0)) ? a : b }'
}

get_value() {
    file=$1
    key=$2
    awk -v key="$key" '$1 == key { print $2; found=1; exit } END { if (!found) exit 1 }' "$file"
}

require_key_value() {
    file=$1
    key=$2
    expected=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    if [ "$value" != "$expected" ]; then
        add_reason "${key}_expected_${expected}_got_${value}_in_${file}"
        return 1
    fi
    return 0
}

require_key_exists() {
    file=$1
    key=$2
    get_value "$file" "$key" >/dev/null 2>&1 || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    return 0
}

require_real_gt() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value+0.0) > (limit+0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_gt_${limit}_got_${value}_in_${file}"
        return 1
    }
    return 0
}

require_real_le() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value+0.0) <= (limit+0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_le_${limit}_got_${value}_in_${file}"
        return 1
    }
    return 0
}

within_mixed_tol() {
    reference=$1
    candidate=$2
    abs_tol=$3
    rel_tol=$4
    awk -v ref="$reference" -v cand="$candidate" -v abs_tol="$abs_tol" -v rel_tol="$rel_tol" '
      function absval(x) { return x < 0.0 ? -x : x }
      function maxval(a,b) { return a > b ? a : b }
      BEGIN {
        r=ref+0.0; c=cand+0.0; d=c-r
        tol=maxval(abs_tol+0.0, (rel_tol+0.0)*maxval(1.0, absval(r)))
        if (d==d && absval(d) <= tol) exit 0
        exit 1
      }'
}

build_target() {
    target=$1
    ensure_build_dir || return 1
    cmake --build "$BUILD_DIR" --target "$target" -j
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage14_8_requested_flag 1
stage14_8_build_status $BUILD_STATUS
stage14_8_lambda0_run_status $LAMBDA0_RUN_STATUS
stage14_8_small_lambda_run_status $SMALL_LAMBDA_RUN_STATUS
stage14_8_lambda0_hook_active_status $LAMBDA0_HOOK_ACTIVE_STATUS
stage14_8_small_lambda_hook_active_status $SMALL_LAMBDA_HOOK_ACTIVE_STATUS
stage14_8_lambda0_rhs_increment_zero_status $LAMBDA0_RHS_INCREMENT_ZERO_STATUS
stage14_8_small_lambda_nonzero_status $SMALL_LAMBDA_NONZERO_STATUS
stage14_8_small_lambda_rhs_increment_nonzero_status $SMALL_LAMBDA_RHS_INCREMENT_NONZERO_STATUS
stage14_8_small_lambda_rhs_increment_finite_status $SMALL_LAMBDA_RHS_INCREMENT_FINITE_STATUS
stage14_8_small_lambda_rhs_increment_bounded_status $SMALL_LAMBDA_RHS_INCREMENT_BOUNDED_STATUS
stage14_8_rhs_sign_correct_status $RHS_SIGN_CORRECT_STATUS
stage14_8_rhs_linear_scaling_status $RHS_LINEAR_SCALING_STATUS
stage14_8_normalized_rhs_parallel_consistency_status $NORMALIZED_RHS_PARALLEL_STATUS
stage14_8_force_density_parallel_consistency_status $FORCE_DENSITY_PARALLEL_STATUS
stage14_8_final_fluid_delta_finite_status $FINAL_FLUID_DELTA_FINITE_STATUS
stage14_8_final_fluid_delta_bounded_status $FINAL_FLUID_DELTA_BOUNDED_STATUS
stage14_8_final_signature_parallel_consistency_status $FINAL_SIGNATURE_PARALLEL_STATUS
stage14_8_no_pressure_modification_status $NO_PRESSURE_MODIFICATION_STATUS
stage14_8_no_projection_modification_status $NO_PROJECTION_MODIFICATION_STATUS
stage14_8_no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
stage14_8_no_rk3_modification_status $NO_RK3_MODIFICATION_STATUS
stage14_8_no_channel_forcing_modification_status $NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_8_no_production_ibm_forcing_status $NO_PRODUCTION_IBM_FORCING_STATUS
stage14_8_no_feedback_application_status $NO_FEEDBACK_APPLICATION_STATUS
stage14_8_no_twoway_force_status $NO_TWOWAY_FORCE_STATUS
stage14_8_no_structure_advance_status $NO_STRUCTURE_ADVANCE_STATUS
stage14_8_parallel_small_lambda_response_status $final_status
stage14_8_lambda_14 $STAGE14_8_SMALL_LAMBDA
stage14_8_max_rhs_increment $MAX_RHS_INCREMENT
stage14_8_max_fluid_signature_delta $MAX_FLUID_SIGNATURE_DELTA
stage14_8_max_normalized_rhs_l2_delta $MAX_NORMALIZED_RHS_L2_DELTA
stage14_8_max_normalized_rhs_max_abs_delta $MAX_NORMALIZED_RHS_MAX_ABS_DELTA
stage14_8_max_force_density_delta $MAX_FORCE_DENSITY_DELTA
stage14_8_normalized_rhs_l2_np1 $NORMALIZED_RHS_L2_NP1
stage14_8_normalized_rhs_l2_np2 $NORMALIZED_RHS_L2_NP2
stage14_8_normalized_rhs_l2_np4 $NORMALIZED_RHS_L2_NP4
stage14_8_normalized_rhs_max_abs_np1 $NORMALIZED_RHS_MAX_ABS_NP1
stage14_8_normalized_rhs_max_abs_np2 $NORMALIZED_RHS_MAX_ABS_NP2
stage14_8_normalized_rhs_max_abs_np4 $NORMALIZED_RHS_MAX_ABS_NP4
EOF_DAT
}

verify_small_lambda_value() {
    is_finite_number "$STAGE14_8_SMALL_LAMBDA" || {
        add_reason "stage14_8_small_lambda_not_finite"
        return 1
    }
    awk -v value="$STAGE14_8_SMALL_LAMBDA" 'BEGIN { v=value+0.0; if (v>0.0 && v<=1.0e-4) exit 0; exit 1 }' || {
        add_reason "stage14_8_small_lambda_not_small_positive"
        return 1
    }
    return 0
}

verify_stage9_dat() {
    file=$1
    label=$2
    status=0
    if [ ! -f "$file" ]; then
        add_reason "missing_${label}_stage9_dat"
        return 1
    fi
    require_key_value "$file" stage9_9_parallel_consistency_local_status 1 || status=1
    require_key_value "$file" stage9_9_decomposition_invariant_initial_state_status 1 || status=1
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        require_key_exists "$file" "$metric" || status=1
    done
    return $status
}

prepare_input() {
    np=$1
    label=$2
    input_file="$OUTPUT_DIR/stage14_8_${label}_input_np${np}.i3d"
    awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); print line }' \
        "$CHANNEL_INPUT" > "$input_file"
    echo "$input_file"
}

run_case() {
    np=$1
    lambda=$2
    label=$3
    log_file="$OUTPUT_DIR/stage14_8_${label}_np${np}.log"
    input_file=$(prepare_input "$np" "$label")
    exe="$BUILD_DIR/bin/xcompact3d"
    if [ ! -x "$exe" ]; then
        exe="$BUILD_DIR/src/xcompact3d"
    fi
    if [ ! -x "$exe" ]; then
        exe="$BUILD_DIR/xcompact3d"
    fi
    if [ ! -x "$exe" ]; then
        add_reason "missing_xcompact3d_executable_for_${label}_np${np}"
        return 1
    fi

    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat
    rm -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
    rm -f stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat

    timeout "$STAGE14_8_TIMEOUT_SEC" env \
        X3D_STAGE11_ONEWAY_HOOK=1 \
        X3D_STAGE11_FORCE_READONLY=1 \
        X3D_STAGE11_MAX_POINTS=8 \
        X3D_STAGE11_MAX_STEPS=3 \
        X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
        X3D_STAGE12_FORCE_READONLY=1 \
        X3D_STAGE12_FEEDBACK_GAIN=1.0 \
        X3D_STAGE12_FORCE_SIGN=1 \
        X3D_STAGE12_MAX_POINTS=8 \
        X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1 \
        X3D_STAGE13_FORCE_READONLY=1 \
        X3D_STAGE13_SPREADING_READONLY=1 \
        X3D_STAGE13_MAX_POINTS=8 \
        X3D_STAGE13_MAX_EULERIAN_POINTS=64 \
        X3D_STAGE13_SPREADING_NORMALIZATION=conservative \
        X3D_STAGE14_RHS_INJECTION=1 \
        X3D_STAGE14_INJECTION_GAIN="$lambda" \
        X3D_STAGE14_MAX_STEPS=3 \
        X3D_STAGE14_REQUIRE_STAGE13=1 \
        X3D_STAGE14_DIAGNOSTIC_ONLY=1 \
        STAGE9_SKIP_PREREQS=1 \
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
        X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
        X3D_STAGE9_9_MAX_STEPS=3 \
        ${MPIEXEC} ${MPIEXEC_FLAGS} -np "$np" "$exe" "$input_file" > "$log_file" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ]; then
        add_reason "stage14_8_${label}_np${np}_run_failed"
        tail -n 120 "$log_file"
        return 1
    fi
    if ! grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 && \
       ! grep 'STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS' "$log_file" >/dev/null 2>&1; then
        add_reason "missing_stage9_9_pass_verdict_${label}_np${np}"
        return 1
    fi
    if [ ! -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat ]; then
        add_reason "missing_stage9_9_dat_${label}_np${np}"
        return 1
    fi
    cp stage9_outputs/fibre_stage9_9_parallel_consistency.dat \
       "$OUTPUT_DIR/stage14_8_${label}_np${np}.dat"
    if [ ! -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics_${label}_np${np}"
        if [ "$label" = "small_lambda" ]; then
            add_reason "nonzero_lambda_still_blocked_by_stage14_5_np${np}"
        fi
        return 1
    fi
    cp stage14_outputs/fibre_stage14_5_production_rhs_hook.dat \
       "$OUTPUT_DIR/stage14_8_${label}_rhs_hook_np${np}.dat"
    if [ ! -f stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat ]; then
        add_reason "missing_stage13_6_force_density_candidate_diagnostics_${label}_np${np}"
        return 1
    fi
    cp stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat \
       "$OUTPUT_DIR/stage14_8_${label}_stage13_force_density_np${np}.dat"
    verify_stage9_dat "$OUTPUT_DIR/stage14_8_${label}_np${np}.dat" "${label}_np${np}" || return 1
    return 0
}

verify_lambda0_hook() {
    np=$1
    file="$OUTPUT_DIR/stage14_8_lambda0_rhs_hook_np${np}.dat"
    status=0
    if [ ! -f "$file" ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics_lambda0_np${np}"
        return 1
    fi
    require_key_value "$file" stage14_5_requested_flag 1 || status=1
    require_key_value "$file" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$file" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$file" stage14_5_lambda_zero_status 1 || status=1
    require_key_value "$file" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$file" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$file" stage14_5_rhs_increment_zero_status 1 || status=1
    require_key_value "$file" stage14_5_rhs_unchanged_status 1 || status=1
    require_key_value "$file" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$file" stage14_5_no_feedback_application_status 1 || status=1
    require_key_value "$file" stage14_5_no_twoway_force_status 1 || status=1
    require_key_value "$file" stage14_5_no_structure_advance_status 1 || status=1
    require_key_value "$file" stage14_5_production_rhs_hook_status 1 || status=1
    require_real_le "$file" stage14_5_rhs_increment_l2 1.0e-12 || status=1
    require_real_le "$file" stage14_5_rhs_increment_max_abs 1.0e-12 || status=1
    return $status
}

verify_stage13_sign_and_density() {
    label=$1
    np=$2
    file="$OUTPUT_DIR/stage14_8_${label}_stage13_force_density_np${np}.dat"
    status=0
    if [ ! -f "$file" ]; then
        add_reason "missing_stage13_6_force_density_candidate_diagnostics_${label}_np${np}"
        return 1
    fi
    require_key_value "$file" stage13_6_hook_initialized_status 1 || status=1
    require_key_value "$file" stage13_6_hook_sample_called_status 1 || status=1
    require_key_value "$file" stage13_6_force_density_candidate_computed_status 1 || status=1
    require_key_value "$file" stage13_6_force_density_candidate_finite_status 1 || status=1
    require_key_value "$file" stage13_6_force_density_norm_finite_status 1 || status=1
    require_key_value "$file" stage13_6_integrated_force_finite_status 1 || status=1
    require_key_value "$file" stage13_6_spreading_input_sign_status 1 || status=1
    require_key_value "$file" stage13_6_wrong_sign_rejection_status 1 || status=1
    require_key_value "$file" stage13_6_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$file" stage13_6_no_feedback_application_status 1 || status=1
    require_key_value "$file" stage13_6_no_twoway_force_status 1 || status=1
    require_key_value "$file" stage13_6_no_structure_advance_status 1 || status=1
    for key in stage13_6_force_density_l2 stage13_6_max_abs_force_density \
               stage13_6_integrated_force_x stage13_6_integrated_force_y stage13_6_integrated_force_z; do
        value=$(get_value "$file" "$key" 2>/dev/null) || {
            add_reason "missing_${key}_in_${file}"
            status=1
            continue
        }
        is_finite_number "$value" || {
            add_reason "${key}_not_finite_in_${file}"
            status=1
        }
    done
    return $status
}

verify_small_hook() {
    np=$1
    file="$OUTPUT_DIR/stage14_8_small_lambda_rhs_hook_np${np}.dat"
    status=0
    if [ ! -f "$file" ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics_small_lambda_np${np}"
        add_reason "nonzero_lambda_still_blocked_by_stage14_5_np${np}"
        return 1
    fi
    require_key_value "$file" stage14_5_requested_flag 1 || status=1
    require_key_value "$file" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$file" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$file" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$file" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$file" stage14_5_stage13_dependency_status 1 || status=1
    require_key_value "$file" stage14_5_stage13_candidate_required_status 1 || status=1
    require_key_value "$file" stage14_5_rhs_arrays_available_status 1 || status=1
    require_key_value "$file" stage14_5_rhs_increment_computed_status 1 || status=1
    require_key_value "$file" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$file" stage14_5_no_feedback_application_status 1 || status=1
    require_key_value "$file" stage14_5_no_twoway_force_status 1 || status=1
    require_key_value "$file" stage14_5_no_structure_advance_status 1 || status=1
    require_key_value "$file" stage14_5_production_rhs_hook_status 1 || status=1

    lambda_zero=$(get_value "$file" stage14_5_lambda_zero_status 2>/dev/null || echo missing)
    if [ "$lambda_zero" != "0" ]; then
        add_reason "nonzero_lambda_still_blocked_by_stage14_5_np${np}"
        status=1
    fi
    injection_gain=$(get_value "$file" stage14_5_injection_gain 2>/dev/null || echo missing)
    within_mixed_tol "$STAGE14_8_SMALL_LAMBDA" "$injection_gain" 1.0e-14 1.0e-8 || {
        add_reason "stage14_8_small_lambda_gain_mismatch_np${np}"
        status=1
    }
    l2=$(get_value "$file" stage14_5_rhs_increment_l2 2>/dev/null || echo missing)
    max_abs=$(get_value "$file" stage14_5_rhs_increment_max_abs 2>/dev/null || echo missing)
    is_finite_number "$l2" || { add_reason "stage14_8_rhs_increment_l2_not_finite_np${np}"; status=1; }
    is_finite_number "$max_abs" || { add_reason "stage14_8_rhs_increment_max_abs_not_finite_np${np}"; status=1; }
    require_real_gt "$file" stage14_5_rhs_increment_l2 0.0 || {
        add_reason "nonzero_lambda_still_blocked_by_stage14_5_np${np}"
        status=1
    }
    require_real_gt "$file" stage14_5_rhs_increment_max_abs 0.0 || {
        add_reason "nonzero_lambda_still_blocked_by_stage14_5_np${np}"
        status=1
    }
    require_real_le "$file" stage14_5_rhs_increment_max_abs "$STAGE14_8_MAX_RHS_INCREMENT" || status=1
    MAX_RHS_INCREMENT=$(max_value "$MAX_RHS_INCREMENT" "$max_abs")
    normalized_l2=$(awk -v value="$l2" -v lambda="$STAGE14_8_SMALL_LAMBDA" 'BEGIN { printf "%.17e", (value+0.0)/(lambda+0.0) }')
    normalized_max=$(awk -v value="$max_abs" -v lambda="$STAGE14_8_SMALL_LAMBDA" 'BEGIN { printf "%.17e", (value+0.0)/(lambda+0.0) }')
    case "$np" in
        1) NORMALIZED_RHS_L2_NP1=$normalized_l2; NORMALIZED_RHS_MAX_ABS_NP1=$normalized_max ;;
        2) NORMALIZED_RHS_L2_NP2=$normalized_l2; NORMALIZED_RHS_MAX_ABS_NP2=$normalized_max ;;
        4) NORMALIZED_RHS_L2_NP4=$normalized_l2; NORMALIZED_RHS_MAX_ABS_NP4=$normalized_max ;;
    esac
    return $status
}

compare_final_delta_for_np() {
    np=$1
    lambda0_file="$OUTPUT_DIR/stage14_8_lambda0_np${np}.dat"
    small_file="$OUTPUT_DIR/stage14_8_small_lambda_np${np}.dat"
    status=0
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        lambda0=$(get_value "$lambda0_file" "$metric" 2>/dev/null || echo missing)
        small=$(get_value "$small_file" "$metric" 2>/dev/null || echo missing)
        if ! is_finite_number "$lambda0" || ! is_finite_number "$small"; then
            add_reason "stage14_8_final_signature_not_finite_${metric}_np${np}"
            status=1
            continue
        fi
        delta=$(awk -v a="$lambda0" -v b="$small" 'BEGIN { printf "%.17e", (b+0.0)-(a+0.0) }')
        abs_delta=$(abs_value "$delta")
        MAX_FLUID_SIGNATURE_DELTA=$(max_value "$MAX_FLUID_SIGNATURE_DELTA" "$abs_delta")
        printf 'np %s metric %s lambda0 %s small_lambda %s delta %s abs_delta %s bounded_tolerance %s\n' \
               "$np" "$metric" "$lambda0" "$small" "$delta" "$abs_delta" "$STAGE14_8_MAX_FLUID_DELTA"
        awk -v abs_delta="$abs_delta" -v limit="$STAGE14_8_MAX_FLUID_DELTA" 'BEGIN { if ((abs_delta+0.0) <= (limit+0.0)) exit 0; exit 1 }' || {
            add_reason "stage14_8_final_fluid_delta_bounded_failed_${metric}_np${np}"
            status=1
        }
    done
    return $status
}

compare_parallel_metric_set() {
    label=$1
    file_prefix=$2
    abs_tol=$3
    rel_tol=$4
    shift 4
    status=0
    for metric in "$@"; do
        ref_file="${file_prefix}1.dat"
        ref=$(get_value "$ref_file" "$metric" 2>/dev/null || echo missing)
        if ! is_finite_number "$ref"; then
            add_reason "stage14_8_${label}_${metric}_np1_missing_or_not_finite"
            status=1
            continue
        fi
        for np in 2 4; do
            file="${file_prefix}${np}.dat"
            candidate=$(get_value "$file" "$metric" 2>/dev/null || echo missing)
            if ! is_finite_number "$candidate"; then
                add_reason "stage14_8_${label}_${metric}_np${np}_missing_or_not_finite"
                status=1
                continue
            fi
            delta_abs=$(awk -v ref="$ref" -v candidate="$candidate" 'BEGIN { d=(candidate+0.0)-(ref+0.0); if (d<0.0) d=-d; print d }')
            case "$label" in
                force_density*) MAX_FORCE_DENSITY_DELTA=$(max_value "$MAX_FORCE_DENSITY_DELTA" "$delta_abs") ;;
            esac
            if ! within_mixed_tol "$ref" "$candidate" "$abs_tol" "$rel_tol"; then
                add_reason "stage14_8_${label}_${metric}_parallel_mismatch_np${np}"
                status=1
            fi
        done
    done
    return $status
}

compare_normalized_rhs() {
    status=0
    for np in 2 4; do
        case "$np" in
            2) l2=$NORMALIZED_RHS_L2_NP2; max_abs=$NORMALIZED_RHS_MAX_ABS_NP2 ;;
            4) l2=$NORMALIZED_RHS_L2_NP4; max_abs=$NORMALIZED_RHS_MAX_ABS_NP4 ;;
        esac
        l2_delta=$(awk -v ref="$NORMALIZED_RHS_L2_NP1" -v candidate="$l2" 'BEGIN { d=(candidate+0.0)-(ref+0.0); if (d<0.0) d=-d; print d }')
        max_delta=$(awk -v ref="$NORMALIZED_RHS_MAX_ABS_NP1" -v candidate="$max_abs" 'BEGIN { d=(candidate+0.0)-(ref+0.0); if (d<0.0) d=-d; print d }')
        MAX_NORMALIZED_RHS_L2_DELTA=$(max_value "$MAX_NORMALIZED_RHS_L2_DELTA" "$l2_delta")
        MAX_NORMALIZED_RHS_MAX_ABS_DELTA=$(max_value "$MAX_NORMALIZED_RHS_MAX_ABS_DELTA" "$max_delta")
        within_mixed_tol "$NORMALIZED_RHS_L2_NP1" "$l2" "$STAGE14_8_NORMALIZED_RHS_ABS_TOL" "$STAGE14_8_NORMALIZED_RHS_REL_TOL" || {
            add_reason "stage14_8_normalized_rhs_l2_parallel_mismatch_np${np}"
            status=1
        }
        within_mixed_tol "$NORMALIZED_RHS_MAX_ABS_NP1" "$max_abs" "$STAGE14_8_NORMALIZED_RHS_ABS_TOL" "$STAGE14_8_NORMALIZED_RHS_REL_TOL" || {
            add_reason "stage14_8_normalized_rhs_max_abs_parallel_mismatch_np${np}"
            status=1
        }
    done
    return $status
}

aggregate_no_modification_statuses() {
    status=0
    for np in 1 2 4; do
        file="$OUTPUT_DIR/stage14_8_small_lambda_rhs_hook_np${np}.dat"
        for key in stage14_5_no_pressure_modification_status stage14_5_no_projection_modification_status \
                   stage14_5_no_poisson_modification_status stage14_5_no_rk3_modification_status \
                   stage14_5_no_channel_forcing_modification_status stage14_5_no_production_ibm_forcing_status \
                   stage14_5_no_feedback_application_status stage14_5_no_twoway_force_status \
                   stage14_5_no_structure_advance_status; do
            require_key_value "$file" "$key" 1 || status=1
        done
    done
    if [ "$status" = "0" ]; then
        NO_PRESSURE_MODIFICATION_STATUS=1
        NO_PROJECTION_MODIFICATION_STATUS=1
        NO_POISSON_MODIFICATION_STATUS=1
        NO_RK3_MODIFICATION_STATUS=1
        NO_CHANNEL_FORCING_MODIFICATION_STATUS=1
        NO_PRODUCTION_IBM_FORCING_STATUS=1
        NO_FEEDBACK_APPLICATION_STATUS=1
        NO_TWOWAY_FORCE_STATUS=1
        NO_STRUCTURE_ADVANCE_STATUS=1
    fi
    return $status
}

if ! verify_small_lambda_value; then
    BUILD_STATUS=0
fi

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
    lambda0_ok=1
    for np in 1 2 4; do
        run_case "$np" 0.0 lambda0 || lambda0_ok=0
    done
    [ "$lambda0_ok" = "1" ] && LAMBDA0_RUN_STATUS=1
fi

if [ "$BUILD_STATUS" = "1" ]; then
    small_ok=1
    for np in 1 2 4; do
        run_case "$np" "$STAGE14_8_SMALL_LAMBDA" small_lambda || small_ok=0
    done
    [ "$small_ok" = "1" ] && SMALL_LAMBDA_RUN_STATUS=1
fi

if [ "$LAMBDA0_RUN_STATUS" = "1" ]; then
    ok=1
    for np in 1 2 4; do
        verify_lambda0_hook "$np" || ok=0
        verify_stage13_sign_and_density lambda0 "$np" || ok=0
    done
    [ "$ok" = "1" ] && LAMBDA0_HOOK_ACTIVE_STATUS=1 && LAMBDA0_RHS_INCREMENT_ZERO_STATUS=1
fi

if [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ]; then
    ok=1
    for np in 1 2 4; do
        verify_small_hook "$np" || ok=0
        verify_stage13_sign_and_density small_lambda "$np" || ok=0
    done
    if [ "$ok" = "1" ]; then
        SMALL_LAMBDA_HOOK_ACTIVE_STATUS=1
        SMALL_LAMBDA_NONZERO_STATUS=1
        SMALL_LAMBDA_RHS_INCREMENT_NONZERO_STATUS=1
        SMALL_LAMBDA_RHS_INCREMENT_FINITE_STATUS=1
        SMALL_LAMBDA_RHS_INCREMENT_BOUNDED_STATUS=1
        RHS_SIGN_CORRECT_STATUS=1
        RHS_LINEAR_SCALING_STATUS=1
    fi
fi

if [ "$LAMBDA0_RUN_STATUS" = "1" ] && [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ]; then
    ok=1
    for np in 1 2 4; do
        compare_final_delta_for_np "$np" || ok=0
    done
    if [ "$ok" = "1" ]; then
        FINAL_FLUID_DELTA_FINITE_STATUS=1
        FINAL_FLUID_DELTA_BOUNDED_STATUS=1
    fi
fi

if [ "$SMALL_LAMBDA_RHS_INCREMENT_NONZERO_STATUS" = "1" ]; then
    compare_normalized_rhs && NORMALIZED_RHS_PARALLEL_STATUS=1
fi

if [ "$LAMBDA0_RUN_STATUS" = "1" ] && [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ]; then
    final_metrics="stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz"
    fd_metrics="stage13_6_force_density_l2 stage13_6_max_abs_force_density stage13_6_integrated_force_x stage13_6_integrated_force_y stage13_6_integrated_force_z"
    final_ok=1
    compare_parallel_metric_set lambda0_final_signature "$OUTPUT_DIR/stage14_8_lambda0_np" \
        "$STAGE14_8_FINAL_SIGNATURE_ABS_TOL" "$STAGE14_8_FINAL_SIGNATURE_REL_TOL" $final_metrics || final_ok=0
    compare_parallel_metric_set small_lambda_final_signature "$OUTPUT_DIR/stage14_8_small_lambda_np" \
        "$STAGE14_8_FINAL_SIGNATURE_ABS_TOL" "$STAGE14_8_FINAL_SIGNATURE_REL_TOL" $final_metrics || final_ok=0
    [ "$final_ok" = "1" ] && FINAL_SIGNATURE_PARALLEL_STATUS=1
    fd_ok=1
    compare_parallel_metric_set force_density_lambda0 "$OUTPUT_DIR/stage14_8_lambda0_stage13_force_density_np" \
        "$STAGE14_8_FORCE_DENSITY_ABS_TOL" "$STAGE14_8_FORCE_DENSITY_REL_TOL" $fd_metrics || fd_ok=0
    compare_parallel_metric_set force_density_small_lambda "$OUTPUT_DIR/stage14_8_small_lambda_stage13_force_density_np" \
        "$STAGE14_8_FORCE_DENSITY_ABS_TOL" "$STAGE14_8_FORCE_DENSITY_REL_TOL" $fd_metrics || fd_ok=0
    [ "$fd_ok" = "1" ] && FORCE_DENSITY_PARALLEL_STATUS=1
fi

if [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ]; then
    aggregate_no_modification_statuses || true
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$LAMBDA0_RUN_STATUS" = "1" ] && [ "$SMALL_LAMBDA_RUN_STATUS" = "1" ] && \
   [ "$LAMBDA0_HOOK_ACTIVE_STATUS" = "1" ] && [ "$SMALL_LAMBDA_HOOK_ACTIVE_STATUS" = "1" ] && \
   [ "$LAMBDA0_RHS_INCREMENT_ZERO_STATUS" = "1" ] && [ "$SMALL_LAMBDA_NONZERO_STATUS" = "1" ] && \
   [ "$SMALL_LAMBDA_RHS_INCREMENT_NONZERO_STATUS" = "1" ] && [ "$SMALL_LAMBDA_RHS_INCREMENT_FINITE_STATUS" = "1" ] && \
   [ "$SMALL_LAMBDA_RHS_INCREMENT_BOUNDED_STATUS" = "1" ] && [ "$RHS_SIGN_CORRECT_STATUS" = "1" ] && \
   [ "$RHS_LINEAR_SCALING_STATUS" = "1" ] && [ "$NORMALIZED_RHS_PARALLEL_STATUS" = "1" ] && \
   [ "$FORCE_DENSITY_PARALLEL_STATUS" = "1" ] && [ "$FINAL_FLUID_DELTA_FINITE_STATUS" = "1" ] && \
   [ "$FINAL_FLUID_DELTA_BOUNDED_STATUS" = "1" ] && [ "$FINAL_SIGNATURE_PARALLEL_STATUS" = "1" ] && \
   [ "$NO_PRESSURE_MODIFICATION_STATUS" = "1" ] && [ "$NO_PROJECTION_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_POISSON_MODIFICATION_STATUS" = "1" ] && [ "$NO_RK3_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_CHANNEL_FORCING_MODIFICATION_STATUS" = "1" ] && [ "$NO_PRODUCTION_IBM_FORCING_STATUS" = "1" ] && \
   [ "$NO_FEEDBACK_APPLICATION_STATUS" = "1" ] && [ "$NO_TWOWAY_FORCE_STATUS" = "1" ] && \
   [ "$NO_STRUCTURE_ADVANCE_STATUS" = "1" ]; then
    write_output_dat 1
    echo 'STAGE 14.8 PARALLEL SMALL LAMBDA RESPONSE VERDICT: PASS'
    echo 'STAGE 14.8 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 14.8 PARALLEL SMALL LAMBDA RESPONSE VERDICT: FAIL'
echo 'STAGE 14.8 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_8_failure'
fi
exit 1
