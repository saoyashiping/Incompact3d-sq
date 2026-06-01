#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
BASELINE_COPY="$OUTPUT_DIR/stage13_7_baseline_stage12_np1.dat"
CANDIDATE_COPY="$OUTPUT_DIR/stage13_7_candidate_np1.dat"
STAGE13_DAT="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate.dat"
OUT_DAT="$OUTPUT_DIR/stage13_7_force_density_invariance_np1.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_7_force_density_invariance_np1_reasons.tmp"
BASELINE_LOG="$OUTPUT_DIR/stage13_7_baseline_stage12_np1.log"
CANDIDATE_LOG="$OUTPUT_DIR/stage13_7_candidate_np1.log"
ABS_TOL=${STAGE13_7_INVARIANCE_ABS_TOL:-1.0e-12}
REL_TOL=${STAGE13_7_INVARIANCE_REL_TOL:-1.0e-14}

mkdir -p stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$BASELINE_COPY" "$CANDIDATE_COPY" "$STAGE13_DAT" "$OUT_DAT"

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
    add_reason "$key expected $expected but found ${value:-missing}"
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

write_output_dat() {
    final_status=$1
    {
        echo "stage13_7_requested_flag 1"
        echo "stage13_7_build_status $build_status"
        echo "stage13_7_baseline_run_status $baseline_run_status"
        echo "stage13_7_candidate_run_status $candidate_run_status"
        echo "stage13_7_hook_active_status $hook_active_status"
        echo "stage13_7_force_density_candidate_computed_status $candidate_computed_status"
        echo "stage13_7_force_density_candidate_finite_status $candidate_finite_status"
        echo "stage13_7_force_density_norm_finite_status $norm_finite_status"
        echo "stage13_7_integrated_force_finite_status $integrated_finite_status"
        echo "stage13_7_integrated_force_conservation_status $integrated_conservation_status"
        echo "stage13_7_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_7_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_7_no_field_modification_status $no_field_modification_status"
        echo "stage13_7_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage13_7_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_7_no_production_ibm_forcing_status $no_production_ibm_forcing_status"
        echo "stage13_7_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_7_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_7_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_7_np1_signature_invariance_status $np1_signature_invariance_status"
        echo "stage13_7_force_density_invariance_np1_status $final_status"
    } > "$OUT_DAT"
}

compare_metric() {
    metric=$1
    baseline_value=$(get_dat_value "$BASELINE_COPY" "$metric")
    candidate_value=$(get_dat_value "$CANDIDATE_COPY" "$metric")
    if [ -z "$baseline_value" ] || [ -z "$candidate_value" ]; then
        echo "metric $metric"
        echo "baseline value ${baseline_value:-missing}"
        echo "candidate value ${candidate_value:-missing}"
        echo "delta missing"
        echo "abs_tol $ABS_TOL"
        echo "rel_tol $REL_TOL"
        echo "effective tolerance missing"
        echo "FAIL"
        add_reason "$metric missing for np=1 signature comparison"
        return 1
    fi
    awk -v metric="$metric" -v baseline="$baseline_value" -v candidate="$candidate_value" \
        -v abs_tol="$ABS_TOL" -v rel_tol="$REL_TOL" '
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
        printf("metric %s\n", metric)
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
        add_reason "$metric np=1 signature invariance failed"
        return 1
    fi
    return 0
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

STAGE13_7_RUN_STAGE13_6=${STAGE13_7_RUN_STAGE13_6:-0}
if [ "$STAGE13_7_RUN_STAGE13_6" = "1" ]; then
    if ! sh stage13_checks/run_stage13_6_production_force_density_candidate.sh; then
        add_reason "optional Stage 13.6 prerequisite failed"
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
    rm -f "$BASELINE_COPY" "$CANDIDATE_COPY" "$STAGE13_DAT"
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
        add_reason "baseline Stage 9.9 np=1 deterministic run failed"
    fi
    if [ -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat ]; then
        cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$BASELINE_COPY"
    else
        add_reason "baseline np=1 Stage 9.9 output missing"
    fi

    rm -f "$STAGE13_DAT"
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
        add_reason "Stage 13 candidate-enabled Stage 9.9 np=1 deterministic run failed"
    fi
    if [ -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat ]; then
        cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$CANDIDATE_COPY"
    else
        add_reason "candidate np=1 Stage 9.9 output missing"
    fi
else
    add_reason "run phase skipped because required build failed"
fi

for file in "$BASELINE_COPY" "$CANDIDATE_COPY"; do
    if [ -f "$file" ]; then
        check_key_equals "$file" stage9_9_parallel_consistency_local_status 1 >/dev/null || true
        check_key_equals "$file" stage9_9_decomposition_invariant_initial_state_status 1 >/dev/null || true
        for key in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                   stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                   stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
            check_key_exists "$file" "$key" >/dev/null || true
        done
    fi
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

if [ -f "$BASELINE_COPY" ] && [ -f "$CANDIDATE_COPY" ]; then
    metric_status=1
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        if ! compare_metric "$metric"; then
            metric_status=0
        fi
    done
    np1_signature_invariance_status=$metric_status
fi

if [ ! -s "$REASONS_FILE" ] && [ "$build_status" = "1" ] && [ "$baseline_run_status" = "1" ] && \
   [ "$candidate_run_status" = "1" ] && [ "$hook_active_status" = "1" ] && \
   [ "$np1_signature_invariance_status" = "1" ]; then
    write_output_dat 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.7 FORCE DENSITY INVARIANCE NP1 VERDICT: PASS"
    echo "STAGE 13.7 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.7 gate failed without a captured reason"
fi
write_output_dat 0
echo "STAGE 13.7 FORCE DENSITY INVARIANCE NP1 VERDICT: FAIL"
echo "STAGE 13.7 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
