#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage12_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
ABS_TOL=${STAGE12_7_INVARIANCE_ABS_TOL:-1.0e-12}
REL_TOL=${STAGE12_7_INVARIANCE_REL_TOL:-1.0e-14}
BASELINE_DAT="$OUTPUT_DIR/stage12_7_baseline_stage11_np1.dat"
CANDIDATE_DAT="$OUTPUT_DIR/stage12_7_candidate_np1.dat"
STAGE12_DIAG="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate.dat"
GATE_DAT="$OUTPUT_DIR/stage12_7_force_candidate_invariance_np1.dat"
BASELINE_LOG="$OUTPUT_DIR/stage12_7_baseline_stage11_np1.log"
CANDIDATE_LOG="$OUTPUT_DIR/stage12_7_candidate_np1.log"
BUILD_LOG="$OUTPUT_DIR/stage12_7_build.log"
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
        echo "stage12_7_requested_flag $requested_flag"
        echo "stage12_7_build_status $build_status"
        echo "stage12_7_baseline_run_status $baseline_run_status"
        echo "stage12_7_candidate_run_status $candidate_run_status"
        echo "stage12_7_hook_active_status $hook_active_status"
        echo "stage12_7_force_candidate_computed_status $force_candidate_computed_status"
        echo "stage12_7_force_candidate_finite_status $force_candidate_finite_status"
        echo "stage12_7_force_norm_finite_status $force_norm_finite_status"
        echo "stage12_7_power_diagnostics_finite_status $power_diagnostics_finite_status"
        echo "stage12_7_no_field_modification_status $no_field_modification_status"
        echo "stage12_7_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage12_7_no_eulerian_force_density_status $no_eulerian_force_density_status"
        echo "stage12_7_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage12_7_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage12_7_no_feedback_application_status $no_feedback_application_status"
        echo "stage12_7_no_twoway_force_status $no_twoway_force_status"
        echo "stage12_7_no_structure_advance_status $no_structure_advance_status"
        echo "stage12_7_np1_signature_invariance_status $np1_signature_invariance_status"
        echo "stage12_7_force_candidate_invariance_np1_status $force_candidate_invariance_np1_status"
    } > "$GATE_DAT"
}

compare_metric() {
    metric=$1
    reference=$(get_value "$metric" "$BASELINE_DAT" 2>/dev/null)
    candidate=$(get_value "$metric" "$CANDIDATE_DAT" 2>/dev/null)
    if [ -z "$reference" ] || [ -z "$candidate" ]; then
        echo "metric=$metric baseline=${reference:-missing} candidate=${candidate:-missing} delta=missing abs_tol=$ABS_TOL rel_tol=$REL_TOL effective_tol=missing FAIL"
        add_reason "$metric"
        return 1
    fi
    awk -v metric="$metric" -v ref="$reference" -v cand="$candidate" -v abs_tol="$ABS_TOL" -v rel_tol="$REL_TOL" '
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
        printf("metric=%s baseline=%.16e candidate=%.16e delta=%.16e abs_tol=%.16e rel_tol=%.16e effective_tol=%.16e %s\n", metric, ref, cand, delta, abs_tol, rel_tol, eff, status)
        exit(ok ? 0 : 1)
      }'
    rc=$?
    if [ $rc -ne 0 ]; then
        add_reason "$metric"
        return 1
    fi
    return 0
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
force_candidate_invariance_np1_status=0

if [ "${STAGE12_7_RUN_STAGE12_6:-0}" = "1" ]; then
    if ! bash stage12_checks/run_stage12_6_production_feedback_candidate.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage12_6_prerequisite_failed"
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
    rm -f "$BASELINE_DAT"
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
    if [ $? -eq 0 ] && [ -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat ]; then
        cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$BASELINE_DAT"
        baseline_run_status=1
    else
        add_reason "baseline_stage11_np1_run_failed"
    fi
else
    add_reason "baseline_stage11_np1_run_skipped"
fi

if [ "$build_ok" -eq 1 ]; then
    rm -f "$CANDIDATE_DAT" "$STAGE12_DIAG"
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
    if [ $? -eq 0 ] && [ -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat ]; then
        cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$CANDIDATE_DAT"
        candidate_run_status=1
    else
        add_reason "candidate_stage12_np1_run_failed"
    fi
else
    add_reason "candidate_stage12_np1_run_skipped"
fi

if [ -f "$BASELINE_DAT" ]; then
    check_key stage9_9_parallel_consistency_local_status 1 "$BASELINE_DAT" >/dev/null
    check_key stage9_9_decomposition_invariant_initial_state_status 1 "$BASELINE_DAT" >/dev/null
else
    add_reason "baseline_stage11_np1_dat_missing"
fi

if [ -f "$CANDIDATE_DAT" ]; then
    check_key stage9_9_parallel_consistency_local_status 1 "$CANDIDATE_DAT" >/dev/null
    check_key stage9_9_decomposition_invariant_initial_state_status 1 "$CANDIDATE_DAT" >/dev/null
else
    add_reason "candidate_stage12_np1_dat_missing"
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

if [ -f "$BASELINE_DAT" ] && [ -f "$CANDIDATE_DAT" ]; then
    invariance_ok=1
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
        if ! compare_metric "$metric"; then
            invariance_ok=0
        fi
    done
    if [ "$invariance_ok" -eq 1 ]; then
        np1_signature_invariance_status=1
    fi
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
   [ "$np1_signature_invariance_status" -eq 1 ]; then
    force_candidate_invariance_np1_status=1
fi

write_gate

if [ "$force_candidate_invariance_np1_status" -eq 1 ]; then
    echo 'STAGE 12.7 FORCE CANDIDATE INVARIANCE NP1 VERDICT: PASS'
    echo 'STAGE 12.7 FINAL VERDICT: PASS'
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage12_7_failure"
fi

echo 'STAGE 12.7 FORCE CANDIDATE INVARIANCE NP1 VERDICT: FAIL'
echo 'STAGE 12.7 FINAL VERDICT: FAIL'
echo "Reasons: $REASONS"
exit 1
