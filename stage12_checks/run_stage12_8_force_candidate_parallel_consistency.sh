#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage12_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
ABS_TOL=${STAGE12_8_INVARIANCE_ABS_TOL:-1.0e-12}
REL_TOL=${STAGE12_8_INVARIANCE_REL_TOL:-1.0e-14}
STAGE12_DIAG="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate.dat"
GATE_DAT="$OUTPUT_DIR/stage12_8_force_candidate_parallel_consistency.dat"
BASELINE_LOG="$OUTPUT_DIR/stage12_8_baseline_stage11_parallel.log"
CANDIDATE_LOG="$OUTPUT_DIR/stage12_8_candidate_parallel.log"
BUILD_LOG="$OUTPUT_DIR/stage12_8_build.log"
REASONS=""

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

add_reason() {
    if [ -z "$REASONS" ]; then
        REASONS="$1"
    else
        REASONS="$REASONS; $1"
    fi
}

get_value() {
    awk -v key="$1" '$1 == key { print $2; found=1 } END { if (!found) exit 1 }' "$2"
}

check_key() {
    key=$1
    expected=$2
    file=$3
    value=$(get_value "$key" "$file" 2>/dev/null)
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$key"
    return 1
}

check_key_status() {
    key=$1
    expected=$2
    file=$3
    value=$(get_value "$key" "$file" 2>/dev/null)
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$key"
    return 1
}

write_gate() {
    {
        echo "stage12_8_requested_flag $requested_flag"
        echo "stage12_8_build_status $build_status"
        echo "stage12_8_baseline_run_status $baseline_run_status"
        echo "stage12_8_candidate_run_status $candidate_run_status"
        echo "stage12_8_hook_active_status $hook_active_status"
        echo "stage12_8_force_candidate_computed_status $force_candidate_computed_status"
        echo "stage12_8_force_candidate_finite_status $force_candidate_finite_status"
        echo "stage12_8_force_norm_finite_status $force_norm_finite_status"
        echo "stage12_8_power_diagnostics_finite_status $power_diagnostics_finite_status"
        echo "stage12_8_no_field_modification_status $no_field_modification_status"
        echo "stage12_8_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage12_8_no_eulerian_force_density_status $no_eulerian_force_density_status"
        echo "stage12_8_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage12_8_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage12_8_no_feedback_application_status $no_feedback_application_status"
        echo "stage12_8_no_twoway_force_status $no_twoway_force_status"
        echo "stage12_8_no_structure_advance_status $no_structure_advance_status"
        echo "stage12_8_np1_signature_invariance_status $np1_signature_invariance_status"
        echo "stage12_8_np2_signature_invariance_status $np2_signature_invariance_status"
        echo "stage12_8_np4_signature_invariance_status $np4_signature_invariance_status"
        echo "stage12_8_parallel_consistency_status $parallel_consistency_status"
        echo "stage12_8_force_candidate_parallel_consistency_status $force_candidate_parallel_consistency_status"
    } > "$GATE_DAT"
}

baseline_dat_for_np() {
    echo "$OUTPUT_DIR/stage12_8_baseline_stage11_np$1.dat"
}

candidate_dat_for_np() {
    echo "$OUTPUT_DIR/stage12_8_candidate_np$1.dat"
}

copy_stage9_dat() {
    prefix=$1
    failed=0
    for np in 1 2 4; do
        src="$STAGE9_OUTPUT_DIR/fibre_stage9_9_parallel_consistency_np${np}.dat"
        if [ "$prefix" = "baseline" ]; then
            dst=$(baseline_dat_for_np "$np")
        else
            dst=$(candidate_dat_for_np "$np")
        fi
        if [ -f "$src" ]; then
            cp "$src" "$dst"
        else
            add_reason "${prefix}_np${np}_stage9_9_dat_missing"
            failed=1
        fi
    done
    return $failed
}

check_stage9_dat_for_np() {
    np=$1
    file=$2
    prefix=$3
    failed=0
    if [ ! -f "$file" ]; then
        add_reason "${prefix}_np${np}_dat_missing"
        return 1
    fi
    if ! check_key stage9_9_parallel_consistency_local_status 1 "$file" >/dev/null; then
        failed=1
    fi
    if ! check_key stage9_9_decomposition_invariant_initial_state_status 1 "$file" >/dev/null; then
        failed=1
    fi
    for metric in \
        stage9_9_signature_sum_ux \
        stage9_9_signature_sum_uy \
        stage9_9_signature_sum_uz \
        stage9_9_signature_max_ux \
        stage9_9_signature_max_uy \
        stage9_9_signature_max_uz \
        stage9_9_signature_l2_ux \
        stage9_9_signature_l2_uy \
        stage9_9_signature_l2_uz
    do
        if ! get_value "$metric" "$file" >/dev/null 2>&1; then
            add_reason "${prefix}_np${np}_${metric}"
            failed=1
        fi
    done
    return $failed
}

compare_metric() {
    np=$1
    metric=$2
    baseline_file=$(baseline_dat_for_np "$np")
    candidate_file=$(candidate_dat_for_np "$np")
    reference=$(get_value "$metric" "$baseline_file" 2>/dev/null)
    candidate=$(get_value "$metric" "$candidate_file" 2>/dev/null)
    if [ -z "$reference" ] || [ -z "$candidate" ]; then
        echo "np=$np metric=$metric baseline=${reference:-missing} candidate=${candidate:-missing} delta=missing abs_tol=$ABS_TOL rel_tol=$REL_TOL effective_tol=missing FAIL"
        add_reason "np${np}_${metric}"
        return 1
    fi
    awk -v np="$np" -v metric="$metric" -v ref="$reference" -v cand="$candidate" -v abs_tol="$ABS_TOL" -v rel_tol="$REL_TOL" '
      BEGIN {
        delta = cand - ref
        if (delta < 0.0) delta = -delta
        aref = ref
        if (aref < 0.0) aref = -aref
        base = aref
        if (base < 1.0) base = 1.0
        eff = abs_tol
        rel_eff = rel_tol * base
        if (rel_eff > eff) eff = rel_eff
        status = "PASS"
        ok = 1
        if (delta > eff) { status = "FAIL"; ok = 0 }
        printf("np=%s metric=%s baseline=%.16e candidate=%.16e delta=%.16e abs_tol=%.16e rel_tol=%.16e effective_tol=%.16e %s\n", np, metric, ref, cand, delta, abs_tol, rel_tol, eff, status)
        exit(ok ? 0 : 1)
      }'
    rc=$?
    if [ $rc -ne 0 ]; then
        add_reason "np${np}_${metric}"
        return 1
    fi
    return 0
}

check_np_invariance() {
    np=$1
    failed=0
    for metric in \
        stage9_9_signature_sum_ux \
        stage9_9_signature_sum_uy \
        stage9_9_signature_sum_uz \
        stage9_9_signature_max_ux \
        stage9_9_signature_max_uy \
        stage9_9_signature_max_uz \
        stage9_9_signature_l2_ux \
        stage9_9_signature_l2_uy \
        stage9_9_signature_l2_uz
    do
        if ! compare_metric "$np" "$metric"; then
            failed=1
        fi
    done
    return $failed
}

mkdir -p "$OUTPUT_DIR" "$STAGE9_OUTPUT_DIR" "$STAGE11_OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$BASELINE_LOG"
: > "$CANDIDATE_LOG"

requested_flag=1
build_status=0
baseline_run_status=0
candidate_run_status=0
hook_active_status=0
force_candidate_computed_status=0
force_candidate_finite_status=0
force_norm_finite_status=0
power_diagnostics_finite_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_eulerian_force_density_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
np1_signature_invariance_status=0
np2_signature_invariance_status=0
np4_signature_invariance_status=0
parallel_consistency_status=0
force_candidate_parallel_consistency_status=0

if [ "${STAGE12_8_RUN_STAGE12_7:-0}" = "1" ]; then
    if ! bash stage12_checks/run_stage12_7_force_candidate_invariance_np1.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage12_7_prerequisite_failed"
    fi
fi

if ensure_build_dir >> "$BUILD_LOG" 2>&1; then
    build_ok=1
else
    build_ok=0
    add_reason "cmake_configure_failed"
fi

if [ "$build_ok" -eq 1 ]; then
    for target in \
        xcompact3d \
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
        fibre_stage12_production_feedback_candidate_check
    do
        if ! cmake --build "$BUILD_DIR" --target "$target" -j >> "$BUILD_LOG" 2>&1; then
            build_ok=0
            add_reason "build_failed_$target"
        fi
    done
fi
build_status=$build_ok

if [ "$build_ok" -eq 1 ]; then
    rm -f "$OUTPUT_DIR"/stage12_8_baseline_stage11_np*.dat
    (
        unset X3D_STAGE12_FEEDBACK_CANDIDATE
        unset X3D_STAGE12_FORCE_READONLY
        unset X3D_STAGE12_FEEDBACK_GAIN
        unset X3D_STAGE12_FORCE_SIGN
        unset X3D_STAGE12_MAX_POINTS
        X3D_STAGE11_ONEWAY_HOOK=1 \
        X3D_STAGE11_FORCE_READONLY=1 \
        X3D_STAGE11_MAX_POINTS=8 \
        X3D_STAGE11_MAX_STEPS=3 \
        STAGE9_SKIP_PREREQS=1 \
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
        X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
        X3D_STAGE9_9_MAX_STEPS=3 \
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
    ) >> "$BASELINE_LOG" 2>&1
    if [ $? -eq 0 ]; then
        copy_stage9_dat baseline
        if [ $? -eq 0 ]; then
            baseline_run_status=1
        else
            add_reason "baseline_stage11_parallel_dat_copy_failed"
        fi
    else
        add_reason "baseline_stage11_parallel_run_failed"
    fi
else
    add_reason "baseline_stage11_parallel_run_skipped"
fi

if [ "$build_ok" -eq 1 ]; then
    rm -f "$OUTPUT_DIR"/stage12_8_candidate_np*.dat "$STAGE12_DIAG"
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
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh >> "$CANDIDATE_LOG" 2>&1
    if [ $? -eq 0 ]; then
        copy_stage9_dat candidate
        if [ $? -eq 0 ]; then
            candidate_run_status=1
        else
            add_reason "candidate_stage12_parallel_dat_copy_failed"
        fi
    else
        add_reason "candidate_stage12_parallel_run_failed"
    fi
else
    add_reason "candidate_stage12_parallel_run_skipped"
fi

baseline_dat_status=1
candidate_dat_status=1
for np in 1 2 4; do
    if ! check_stage9_dat_for_np "$np" "$(baseline_dat_for_np "$np")" baseline >/dev/null; then
        baseline_dat_status=0
    fi
    if ! check_stage9_dat_for_np "$np" "$(candidate_dat_for_np "$np")" candidate >/dev/null; then
        candidate_dat_status=0
    fi
done
if [ "$baseline_dat_status" -ne 1 ]; then
    baseline_run_status=0
fi
if [ "$candidate_dat_status" -ne 1 ]; then
    candidate_run_status=0
fi

if [ -f "$STAGE12_DIAG" ]; then
    requested_ok=0
    readonly_ok=0
    initialized_ok=0
    sample_called_ok=0
    check_key_status stage12_6_requested_flag 1 "$STAGE12_DIAG" && requested_ok=1
    check_key_status stage12_6_readonly_mode_status 1 "$STAGE12_DIAG" && readonly_ok=1
    check_key_status stage12_6_hook_initialized_status 1 "$STAGE12_DIAG" && initialized_ok=1
    check_key_status stage12_6_hook_sample_called_status 1 "$STAGE12_DIAG" && sample_called_ok=1
    if [ "$requested_ok" -eq 1 ] && [ "$readonly_ok" -eq 1 ] && [ "$initialized_ok" -eq 1 ] && [ "$sample_called_ok" -eq 1 ]; then
        hook_active_status=1
    fi
    check_key_status stage12_6_sampled_velocity_available_status 1 "$STAGE12_DIAG" >/dev/null
    check_key_status stage12_6_force_candidate_computed_status 1 "$STAGE12_DIAG" && force_candidate_computed_status=1
    check_key_status stage12_6_force_candidate_finite_status 1 "$STAGE12_DIAG" && force_candidate_finite_status=1
    check_key_status stage12_6_force_norm_finite_status 1 "$STAGE12_DIAG" && force_norm_finite_status=1
    check_key_status stage12_6_power_diagnostics_finite_status 1 "$STAGE12_DIAG" && power_diagnostics_finite_status=1
    check_key_status stage12_6_action_reaction_status 1 "$STAGE12_DIAG" >/dev/null
    check_key_status stage12_6_pair_power_consistency_status 1 "$STAGE12_DIAG" >/dev/null
    check_key_status stage12_6_field_modified_status 0 "$STAGE12_DIAG" && no_field_modification_status=1
    check_key_status stage12_6_rhs_modified_status 0 "$STAGE12_DIAG" && no_rhs_modification_status=1
    check_key_status stage12_6_no_eulerian_force_density_status 1 "$STAGE12_DIAG" && no_eulerian_force_density_status=1
    check_key_status stage12_6_no_rhs_injection_status 1 "$STAGE12_DIAG" && no_rhs_injection_status=1
    check_key_status stage12_6_no_ibm_spreading_status 1 "$STAGE12_DIAG" && no_ibm_spreading_status=1
    check_key_status stage12_6_no_feedback_application_status 1 "$STAGE12_DIAG" && no_feedback_application_status=1
    check_key_status stage12_6_no_twoway_force_status 1 "$STAGE12_DIAG" && no_twoway_force_status=1
    check_key_status stage12_6_no_structure_advance_status 1 "$STAGE12_DIAG" && no_structure_advance_status=1
    check_key_status stage12_6_production_feedback_candidate_status 1 "$STAGE12_DIAG" >/dev/null
else
    add_reason "missing_stage12_6_feedback_candidate_diagnostics"
fi

if [ "$candidate_run_status" -eq 1 ] && grep -q 'STAGE 9.9 FINAL VERDICT: PASS' "$CANDIDATE_LOG"; then
    parallel_consistency_status=1
else
    add_reason "stage9_9_candidate_final_verdict"
fi

if check_np_invariance 1; then
    np1_signature_invariance_status=1
fi
if check_np_invariance 2; then
    np2_signature_invariance_status=1
fi
if check_np_invariance 4; then
    np4_signature_invariance_status=1
fi

if [ "$requested_flag" -eq 1 ] && \
   [ "$build_status" -eq 1 ] && \
   [ "$baseline_run_status" -eq 1 ] && \
   [ "$candidate_run_status" -eq 1 ] && \
   [ "$hook_active_status" -eq 1 ] && \
   [ "$force_candidate_computed_status" -eq 1 ] && \
   [ "$force_candidate_finite_status" -eq 1 ] && \
   [ "$force_norm_finite_status" -eq 1 ] && \
   [ "$power_diagnostics_finite_status" -eq 1 ] && \
   [ "$no_field_modification_status" -eq 1 ] && \
   [ "$no_rhs_modification_status" -eq 1 ] && \
   [ "$no_eulerian_force_density_status" -eq 1 ] && \
   [ "$no_rhs_injection_status" -eq 1 ] && \
   [ "$no_ibm_spreading_status" -eq 1 ] && \
   [ "$no_feedback_application_status" -eq 1 ] && \
   [ "$no_twoway_force_status" -eq 1 ] && \
   [ "$no_structure_advance_status" -eq 1 ] && \
   [ "$np1_signature_invariance_status" -eq 1 ] && \
   [ "$np2_signature_invariance_status" -eq 1 ] && \
   [ "$np4_signature_invariance_status" -eq 1 ] && \
   [ "$parallel_consistency_status" -eq 1 ]; then
    force_candidate_parallel_consistency_status=1
fi

write_gate

if [ "$force_candidate_parallel_consistency_status" -eq 1 ]; then
    echo 'STAGE 12.8 FORCE CANDIDATE PARALLEL CONSISTENCY VERDICT: PASS'
    echo 'STAGE 12.8 FINAL VERDICT: PASS'
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage12_8_failure"
fi

echo 'STAGE 12.8 FORCE CANDIDATE PARALLEL CONSISTENCY VERDICT: FAIL'
echo 'STAGE 12.8 FINAL VERDICT: FAIL'
echo "Reasons: $REASONS"
exit 1
