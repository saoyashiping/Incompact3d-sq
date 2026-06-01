#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
STAGE13_DAT="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate.dat"
OUT_DAT="$OUTPUT_DIR/stage13_8_force_density_parallel_consistency.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_8_force_density_parallel_consistency_reasons.tmp"
BASELINE_LOG="$OUTPUT_DIR/stage13_8_baseline_stage12.log"
CANDIDATE_LOG="$OUTPUT_DIR/stage13_8_candidate.log"
ABS_TOL=${STAGE13_8_INVARIANCE_ABS_TOL:-1.0e-12}
REL_TOL=${STAGE13_8_INVARIANCE_REL_TOL:-1.0e-14}

mkdir -p stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$OUT_DAT" "$STAGE13_DAT"
rm -f "$OUTPUT_DIR"/stage13_8_baseline_stage12_np1.dat \
      "$OUTPUT_DIR"/stage13_8_baseline_stage12_np2.dat \
      "$OUTPUT_DIR"/stage13_8_baseline_stage12_np4.dat \
      "$OUTPUT_DIR"/stage13_8_candidate_np1.dat \
      "$OUTPUT_DIR"/stage13_8_candidate_np2.dat \
      "$OUTPUT_DIR"/stage13_8_candidate_np4.dat

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

add_reason() {
    echo "$1" >> "$REASONS_FILE"
}

get_dat_value() {
    file=$1
    key=$2
    awk -v key="$key" '$1 == key { value=$2 } END { if (value != "") print value }' "$file"
}

check_key_equals() {
    file=$1
    key=$2
    expected=$3
    value=$(get_dat_value "$file" "$key")
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$key expected $expected but found ${value:-missing} in $file"
    return 1
}

check_key_exists() {
    file=$1
    key=$2
    value=$(get_dat_value "$file" "$key")
    if [ -n "$value" ]; then
        return 0
    fi
    add_reason "$key missing from $file"
    return 1
}

copy_stage9_outputs() {
    prefix=$1
    status=1
    for np in 1 2 4; do
        src="stage9_outputs/fibre_stage9_9_parallel_consistency_np${np}.dat"
        dst="$OUTPUT_DIR/stage13_8_${prefix}_np${np}.dat"
        if [ -f "$src" ]; then
            cp "$src" "$dst"
        else
            add_reason "$prefix np=$np Stage 9.9 output missing"
            status=0
        fi
    done
    return "$status"
}

validate_stage9_dat() {
    file=$1
    np=$2
    label=$3
    status=1
    if [ ! -f "$file" ]; then
        add_reason "$label np=$np dat file missing: $file"
        return 1
    fi
    check_key_equals "$file" stage9_9_parallel_consistency_local_status 1 >/dev/null || status=0
    check_key_equals "$file" stage9_9_decomposition_invariant_initial_state_status 1 >/dev/null || status=0
    for key in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
               stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
               stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        check_key_exists "$file" "$key" >/dev/null || status=0
    done
    return "$status"
}

compare_metric() {
    np=$1
    metric=$2
    baseline_file="$OUTPUT_DIR/stage13_8_baseline_stage12_np${np}.dat"
    candidate_file="$OUTPUT_DIR/stage13_8_candidate_np${np}.dat"
    baseline_value=$(get_dat_value "$baseline_file" "$metric")
    candidate_value=$(get_dat_value "$candidate_file" "$metric")
    if [ -z "$baseline_value" ] || [ -z "$candidate_value" ]; then
        echo "np value $np"
        echo "metric name $metric"
        echo "baseline value ${baseline_value:-missing}"
        echo "candidate value ${candidate_value:-missing}"
        echo "delta missing"
        echo "abs_tol $ABS_TOL"
        echo "rel_tol $REL_TOL"
        echo "effective tolerance missing"
        echo "FAIL"
        add_reason "$metric missing for np=$np signature comparison"
        return 1
    fi
    awk -v np="$np" -v metric="$metric" -v baseline="$baseline_value" \
        -v candidate="$candidate_value" -v abs_tol="$ABS_TOL" -v rel_tol="$REL_TOL" '
      BEGIN {
        b = baseline + 0.0
        c = candidate + 0.0
        delta = c - b
        abs_delta = delta
        if (abs_delta < 0.0) abs_delta = -abs_delta
        ref_abs = b
        if (ref_abs < 0.0) ref_abs = -ref_abs
        scale = ref_abs
        if (scale < 1.0) scale = 1.0
        eff_tol = abs_tol + 0.0
        rel_eff = (rel_tol + 0.0) * scale
        if (rel_eff > eff_tol) eff_tol = rel_eff
        status = "PASS"
        ok = 1
        if (abs_delta > eff_tol) { status = "FAIL"; ok = 0 }
        printf("np value %s\n", np)
        printf("metric name %s\n", metric)
        printf("baseline value %.16e\n", b)
        printf("candidate value %.16e\n", c)
        printf("delta %.16e\n", delta)
        printf("abs_tol %.16e\n", abs_tol + 0.0)
        printf("rel_tol %.16e\n", rel_tol + 0.0)
        printf("effective tolerance %.16e\n", eff_tol)
        printf("%s\n", status)
        exit(ok ? 0 : 1)
      }'
    if [ $? -ne 0 ]; then
        add_reason "$metric np=$np signature invariance failed"
        return 1
    fi
    return 0
}

write_output_dat() {
    final_status=$1
    {
        echo "stage13_8_requested_flag 1"
        echo "stage13_8_build_status $build_status"
        echo "stage13_8_baseline_run_status $baseline_run_status"
        echo "stage13_8_candidate_run_status $candidate_run_status"
        echo "stage13_8_hook_active_status $hook_active_status"
        echo "stage13_8_force_density_candidate_computed_status $candidate_computed_status"
        echo "stage13_8_force_density_candidate_finite_status $candidate_finite_status"
        echo "stage13_8_force_density_norm_finite_status $norm_finite_status"
        echo "stage13_8_integrated_force_finite_status $integrated_finite_status"
        echo "stage13_8_integrated_force_conservation_status $integrated_conservation_status"
        echo "stage13_8_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_8_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_8_no_field_modification_status $no_field_modification_status"
        echo "stage13_8_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage13_8_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_8_no_production_ibm_forcing_status $no_production_ibm_forcing_status"
        echo "stage13_8_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_8_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_8_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_8_np1_signature_invariance_status $np1_signature_invariance_status"
        echo "stage13_8_np2_signature_invariance_status $np2_signature_invariance_status"
        echo "stage13_8_np4_signature_invariance_status $np4_signature_invariance_status"
        echo "stage13_8_parallel_consistency_status $parallel_consistency_status"
        echo "stage13_8_force_density_parallel_consistency_status $final_status"
    } > "$OUT_DAT"
}

build_status=1
baseline_run_status=0
candidate_run_status=0
hook_active_status=0
candidate_computed_status=0
candidate_finite_status=0
norm_finite_status=0
integrated_finite_status=0
integrated_conservation_status=0
spreading_input_sign_status=0
wrong_sign_rejection_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_rhs_injection_status=0
no_production_ibm_forcing_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
np1_signature_invariance_status=0
np2_signature_invariance_status=0
np4_signature_invariance_status=0
parallel_consistency_status=0

STAGE13_8_RUN_STAGE13_7=${STAGE13_8_RUN_STAGE13_7:-0}
if [ "$STAGE13_8_RUN_STAGE13_7" = "1" ]; then
    if ! sh stage13_checks/run_stage13_7_force_density_invariance_np1.sh; then
        add_reason "optional Stage 13.7 prerequisite failed"
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
fibre_stage13_conservation_sign_audit_check \
fibre_stage13_production_force_density_candidate_check"

for target in $targets; do
    if ! cmake --build "$BUILD_DIR" --target "$target" -j; then
        build_status=0
        add_reason "build failed for $target"
    fi
done

if [ "$build_status" = "1" ]; then
    rm -f "$STAGE13_DAT"
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
          stage9_outputs/fibre_stage9_9_parallel_consistency_np2.dat \
          stage9_outputs/fibre_stage9_9_parallel_consistency_np4.dat
    (
        unset X3D_STAGE13_FORCE_DENSITY_CANDIDATE
        unset X3D_STAGE13_FORCE_READONLY
        unset X3D_STAGE13_SPREADING_READONLY
        unset X3D_STAGE13_MAX_POINTS
        unset X3D_STAGE13_MAX_EULERIAN_POINTS
        unset X3D_STAGE13_SPREADING_NORMALIZATION
        X3D_STAGE11_ONEWAY_HOOK=1 \
        X3D_STAGE11_FORCE_READONLY=1 \
        X3D_STAGE11_MAX_POINTS=8 \
        X3D_STAGE11_MAX_STEPS=3 \
        X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
        X3D_STAGE12_FORCE_READONLY=1 \
        X3D_STAGE12_FEEDBACK_GAIN=1.0 \
        X3D_STAGE12_FORCE_SIGN=1 \
        X3D_STAGE12_MAX_POINTS=8 \
        STAGE9_SKIP_PREREQS=1 \
        STAGE9_9_MAX_STEPS=3 \
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
        X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
        X3D_STAGE9_9_MAX_STEPS=3 \
            bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
    ) > "$BASELINE_LOG" 2>&1
    if [ $? -eq 0 ]; then
        baseline_run_status=1
    else
        add_reason "baseline Stage 9.9 deterministic np=1/2/4 run failed"
    fi
    copy_stage9_outputs baseline_stage12 >/dev/null || true

    rm -f "$STAGE13_DAT"
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat \
          stage9_outputs/fibre_stage9_9_parallel_consistency_np2.dat \
          stage9_outputs/fibre_stage9_9_parallel_consistency_np4.dat
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
    STAGE9_SKIP_PREREQS=1 \
    STAGE9_9_MAX_STEPS=3 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$CANDIDATE_LOG" 2>&1
    if [ $? -eq 0 ]; then
        candidate_run_status=1
    else
        add_reason "Stage 13 candidate-enabled Stage 9.9 deterministic np=1/2/4 run failed"
    fi
    copy_stage9_outputs candidate >/dev/null || true
    if grep 'STAGE 9.9 FINAL VERDICT: PASS' "$CANDIDATE_LOG" >/dev/null 2>&1; then
        parallel_consistency_status=1
    else
        add_reason "candidate Stage 9.9 final verdict missing or not PASS"
    fi
else
    add_reason "run phase skipped because required build failed"
fi

for np in 1 2 4; do
    baseline_file="$OUTPUT_DIR/stage13_8_baseline_stage12_np${np}.dat"
    candidate_file="$OUTPUT_DIR/stage13_8_candidate_np${np}.dat"
    validate_stage9_dat "$baseline_file" "$np" baseline >/dev/null || true
    validate_stage9_dat "$candidate_file" "$np" candidate >/dev/null || true
done

if [ -f "$STAGE13_DAT" ]; then
    check_key_equals "$STAGE13_DAT" stage13_6_requested_flag 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_readonly_mode_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_spreading_readonly_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_hook_initialized_status 1 >/dev/null && hook_active_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_hook_sample_called_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_sampled_velocity_available_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_force_density_candidate_computed_status 1 >/dev/null && candidate_computed_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_force_density_candidate_finite_status 1 >/dev/null && candidate_finite_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_force_density_norm_finite_status 1 >/dev/null && norm_finite_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_integrated_force_finite_status 1 >/dev/null && integrated_finite_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_integrated_force_conservation_status 1 >/dev/null && integrated_conservation_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_spreading_input_sign_status 1 >/dev/null && spreading_input_sign_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_wrong_sign_rejection_status 1 >/dev/null && wrong_sign_rejection_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_field_modified_status 0 >/dev/null && no_field_modification_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_rhs_modified_status 0 >/dev/null && no_rhs_modification_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_rhs_injection_status 1 >/dev/null && no_rhs_injection_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_production_ibm_forcing_status 1 >/dev/null && no_production_ibm_forcing_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_feedback_application_status 1 >/dev/null && no_feedback_application_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_twoway_force_status 1 >/dev/null && no_twoway_force_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_structure_advance_status 1 >/dev/null && no_structure_advance_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_production_force_density_candidate_status 1 >/dev/null || true
else
    add_reason "missing_stage13_6_force_density_candidate_diagnostics"
fi

for np in 1 2 4; do
    signature_status=1
    baseline_file="$OUTPUT_DIR/stage13_8_baseline_stage12_np${np}.dat"
    candidate_file="$OUTPUT_DIR/stage13_8_candidate_np${np}.dat"
    if [ -f "$baseline_file" ] && [ -f "$candidate_file" ]; then
        for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                      stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                      stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
            if ! compare_metric "$np" "$metric"; then
                signature_status=0
            fi
        done
    else
        signature_status=0
        add_reason "np=$np signature comparison skipped because baseline or candidate dat is missing"
    fi
    if [ "$np" = "1" ]; then
        np1_signature_invariance_status=$signature_status
    elif [ "$np" = "2" ]; then
        np2_signature_invariance_status=$signature_status
    else
        np4_signature_invariance_status=$signature_status
    fi
done

if [ ! -s "$REASONS_FILE" ] && [ "$build_status" = "1" ] && [ "$baseline_run_status" = "1" ] && \
   [ "$candidate_run_status" = "1" ] && [ "$hook_active_status" = "1" ] && \
   [ "$candidate_computed_status" = "1" ] && [ "$candidate_finite_status" = "1" ] && \
   [ "$norm_finite_status" = "1" ] && [ "$integrated_finite_status" = "1" ] && \
   [ "$integrated_conservation_status" = "1" ] && [ "$spreading_input_sign_status" = "1" ] && \
   [ "$wrong_sign_rejection_status" = "1" ] && [ "$no_field_modification_status" = "1" ] && \
   [ "$no_rhs_modification_status" = "1" ] && [ "$no_rhs_injection_status" = "1" ] && \
   [ "$no_production_ibm_forcing_status" = "1" ] && [ "$no_feedback_application_status" = "1" ] && \
   [ "$no_twoway_force_status" = "1" ] && [ "$no_structure_advance_status" = "1" ] && \
   [ "$np1_signature_invariance_status" = "1" ] && [ "$np2_signature_invariance_status" = "1" ] && \
   [ "$np4_signature_invariance_status" = "1" ] && [ "$parallel_consistency_status" = "1" ]; then
    write_output_dat 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.8 FORCE DENSITY PARALLEL CONSISTENCY VERDICT: PASS"
    echo "STAGE 13.8 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.8 gate failed without a captured reason"
fi
write_output_dat 0
echo "STAGE 13.8 FORCE DENSITY PARALLEL CONSISTENCY VERDICT: FAIL"
echo "STAGE 13.8 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
