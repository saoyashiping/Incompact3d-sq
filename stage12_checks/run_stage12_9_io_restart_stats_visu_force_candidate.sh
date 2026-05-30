#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage12_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
STAGE12_DIAG="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate.dat"
GATE_DAT="$OUTPUT_DIR/stage12_9_io_restart_stats_visu_force_candidate.dat"
BUILD_LOG="$OUTPUT_DIR/stage12_9_build.log"
STAGE9_7_LOG="$OUTPUT_DIR/stage12_9_stage9_7_stats_visu_io_force_candidate.log"
STAGE9_8_LOG="$OUTPUT_DIR/stage12_9_stage9_8_restart_io_force_candidate.log"
STAGE9_7_DAT_COPY="$OUTPUT_DIR/stage12_9_stage9_7_stats_visu_io_force_candidate.dat"
STAGE9_8_DAT_COPY="$OUTPUT_DIR/stage12_9_stage9_8_restart_io_force_candidate.dat"
STAGE9_7_DIAG_COPY="$OUTPUT_DIR/stage12_9_stage9_7_feedback_candidate_diagnostics.dat"
STAGE9_8_DIAG_COPY="$OUTPUT_DIR/stage12_9_stage9_8_feedback_candidate_diagnostics.dat"
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

check_diag_key() {
    phase=$1
    key=$2
    expected=$3
    file=$4
    value=$(get_value "$key" "$file" 2>/dev/null)
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "${phase}_${key}"
    return 1
}

write_gate() {
    {
        echo "stage12_9_requested_flag $requested_flag"
        echo "stage12_9_build_status $build_status"
        echo "stage12_9_stage9_7_stats_visu_status $stage9_7_stats_visu_status"
        echo "stage12_9_stage9_8_restart_status $stage9_8_restart_status"
        echo "stage12_9_stage9_7_hook_active_status $stage9_7_hook_active_status"
        echo "stage12_9_stage9_8_hook_active_status $stage9_8_hook_active_status"
        echo "stage12_9_force_candidate_computed_status $force_candidate_computed_status"
        echo "stage12_9_force_candidate_finite_status $force_candidate_finite_status"
        echo "stage12_9_force_norm_finite_status $force_norm_finite_status"
        echo "stage12_9_power_diagnostics_finite_status $power_diagnostics_finite_status"
        echo "stage12_9_action_reaction_status $action_reaction_status"
        echo "stage12_9_pair_power_consistency_status $pair_power_consistency_status"
        echo "stage12_9_no_field_modification_status $no_field_modification_status"
        echo "stage12_9_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage12_9_no_eulerian_force_density_status $no_eulerian_force_density_status"
        echo "stage12_9_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage12_9_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage12_9_no_feedback_application_status $no_feedback_application_status"
        echo "stage12_9_no_twoway_force_status $no_twoway_force_status"
        echo "stage12_9_no_structure_advance_status $no_structure_advance_status"
        echo "stage12_9_no_restart_contamination_status $no_restart_contamination_status"
        echo "stage12_9_no_stats_visu_contamination_status $no_stats_visu_contamination_status"
        echo "stage12_9_io_restart_stats_visu_force_candidate_status $io_restart_stats_visu_force_candidate_status"
    } > "$GATE_DAT"
}

copy_optional_matches() {
    phase=$1
    first_copy=$2
    copied_first=0
    for pattern in "$@"; do
        if [ "$pattern" = "$phase" ] || [ "$pattern" = "$first_copy" ]; then
            continue
        fi
        for file in $pattern; do
            if [ -f "$file" ]; then
                if [ "$copied_first" -eq 0 ]; then
                    cp "$file" "$first_copy"
                    copied_first=1
                fi
                base=$(basename "$file")
                cp "$file" "$OUTPUT_DIR/stage12_9_${phase}_${base}"
            fi
        done
    done
    return 0
}

verify_stage12_diag() {
    phase=$1
    missing_reason=$2
    copy_file=$3
    diag_ok=1
    requested_ok=0
    readonly_ok=0
    initialized_ok=0
    sample_called_ok=0

    if [ ! -f "$STAGE12_DIAG" ]; then
        add_reason "$missing_reason"
        return 1
    fi

    if check_diag_key "$phase" stage12_6_requested_flag 1 "$STAGE12_DIAG"; then requested_ok=1; else diag_ok=0; fi
    if check_diag_key "$phase" stage12_6_readonly_mode_status 1 "$STAGE12_DIAG"; then readonly_ok=1; else diag_ok=0; fi
    if check_diag_key "$phase" stage12_6_hook_initialized_status 1 "$STAGE12_DIAG"; then initialized_ok=1; else diag_ok=0; fi
    if check_diag_key "$phase" stage12_6_hook_sample_called_status 1 "$STAGE12_DIAG"; then sample_called_ok=1; else diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_sampled_velocity_available_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_force_candidate_computed_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_force_candidate_finite_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_force_norm_finite_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_power_diagnostics_finite_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_action_reaction_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_pair_power_consistency_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_field_modified_status 0 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_rhs_modified_status 0 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_no_eulerian_force_density_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_no_rhs_injection_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_no_ibm_spreading_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_no_feedback_application_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_no_twoway_force_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_no_structure_advance_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi
    if ! check_diag_key "$phase" stage12_6_production_feedback_candidate_status 1 "$STAGE12_DIAG"; then diag_ok=0; fi

    cp "$STAGE12_DIAG" "$copy_file"

    if [ "$requested_ok" -eq 1 ] && [ "$readonly_ok" -eq 1 ] && [ "$initialized_ok" -eq 1 ] && [ "$sample_called_ok" -eq 1 ]; then
        if [ "$phase" = "stage9_7" ]; then
            stage9_7_hook_active_status=1
        else
            stage9_8_hook_active_status=1
        fi
    fi

    if [ "$diag_ok" -eq 0 ]; then
        return 1
    fi
    return 0
}

mkdir -p "$OUTPUT_DIR" "$STAGE11_OUTPUT_DIR" "$STAGE9_OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$STAGE9_7_LOG"
: > "$STAGE9_8_LOG"

requested_flag=1
build_status=0
stage9_7_stats_visu_status=0
stage9_8_restart_status=0
stage9_7_hook_active_status=0
stage9_8_hook_active_status=0
force_candidate_computed_status=0
force_candidate_finite_status=0
force_norm_finite_status=0
power_diagnostics_finite_status=0
action_reaction_status=0
pair_power_consistency_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_eulerian_force_density_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_restart_contamination_status=0
no_stats_visu_contamination_status=0
io_restart_stats_visu_force_candidate_status=0
stage9_7_diag_status=0
stage9_8_diag_status=0

if [ "${STAGE12_9_RUN_STAGE12_8:-0}" = "1" ]; then
    if ! bash stage12_checks/run_stage12_8_force_candidate_parallel_consistency.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage12_8_prerequisite_failed"
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
    rm -f "$STAGE12_DIAG"
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
    bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh >> "$STAGE9_7_LOG" 2>&1
    if [ $? -eq 0 ] && grep -q 'STAGE 9.7 FINAL VERDICT: PASS' "$STAGE9_7_LOG"; then
        stage9_7_stats_visu_status=1
    else
        add_reason "stage9_7_stats_visu_final_verdict"
    fi
    copy_optional_matches stage9_7 "$STAGE9_7_DAT_COPY" "$STAGE9_OUTPUT_DIR"/*stage9_7*.dat "$STAGE9_OUTPUT_DIR"/*stats*.dat "$STAGE9_OUTPUT_DIR"/*visu*.dat "$STAGE9_OUTPUT_DIR"/*coarse*.dat
    if verify_stage12_diag stage9_7 missing_stage12_6_feedback_candidate_diagnostics_after_stage9_7 "$STAGE9_7_DIAG_COPY"; then
        stage9_7_diag_status=1
    fi
else
    add_reason "stage9_7_stats_visu_run_skipped"
fi

if [ "$build_ok" -eq 1 ]; then
    rm -f "$STAGE12_DIAG"
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
    bash stage9_checks/run_stage9_8_restart_io_regression.sh >> "$STAGE9_8_LOG" 2>&1
    if [ $? -eq 0 ] && grep -q 'STAGE 9.8 FINAL VERDICT: PASS' "$STAGE9_8_LOG"; then
        stage9_8_restart_status=1
    else
        add_reason "stage9_8_restart_final_verdict"
    fi
    copy_optional_matches stage9_8 "$STAGE9_8_DAT_COPY" "$STAGE9_OUTPUT_DIR"/*stage9_8*.dat "$STAGE9_OUTPUT_DIR"/*restart*.dat
    if verify_stage12_diag stage9_8 missing_stage12_6_feedback_candidate_diagnostics_after_stage9_8 "$STAGE9_8_DIAG_COPY"; then
        stage9_8_diag_status=1
    fi
else
    add_reason "stage9_8_restart_run_skipped"
fi

if [ "$stage9_7_diag_status" -eq 1 ] && [ "$stage9_8_diag_status" -eq 1 ]; then
    force_candidate_computed_status=1
    force_candidate_finite_status=1
    force_norm_finite_status=1
    power_diagnostics_finite_status=1
    action_reaction_status=1
    pair_power_consistency_status=1
    no_field_modification_status=1
    no_rhs_modification_status=1
    no_eulerian_force_density_status=1
    no_rhs_injection_status=1
    no_ibm_spreading_status=1
    no_feedback_application_status=1
    no_twoway_force_status=1
    no_structure_advance_status=1
fi

if [ "$stage9_8_restart_status" -eq 1 ] && [ "$stage9_8_diag_status" -eq 1 ]; then
    no_restart_contamination_status=1
fi
if [ "$stage9_7_stats_visu_status" -eq 1 ] && [ "$stage9_7_diag_status" -eq 1 ]; then
    no_stats_visu_contamination_status=1
fi

if [ "$requested_flag" -eq 1 ] && \
   [ "$build_status" -eq 1 ] && \
   [ "$stage9_7_stats_visu_status" -eq 1 ] && \
   [ "$stage9_8_restart_status" -eq 1 ] && \
   [ "$stage9_7_hook_active_status" -eq 1 ] && \
   [ "$stage9_8_hook_active_status" -eq 1 ] && \
   [ "$force_candidate_computed_status" -eq 1 ] && \
   [ "$force_candidate_finite_status" -eq 1 ] && \
   [ "$force_norm_finite_status" -eq 1 ] && \
   [ "$power_diagnostics_finite_status" -eq 1 ] && \
   [ "$action_reaction_status" -eq 1 ] && \
   [ "$pair_power_consistency_status" -eq 1 ] && \
   [ "$no_field_modification_status" -eq 1 ] && \
   [ "$no_rhs_modification_status" -eq 1 ] && \
   [ "$no_eulerian_force_density_status" -eq 1 ] && \
   [ "$no_rhs_injection_status" -eq 1 ] && \
   [ "$no_ibm_spreading_status" -eq 1 ] && \
   [ "$no_feedback_application_status" -eq 1 ] && \
   [ "$no_twoway_force_status" -eq 1 ] && \
   [ "$no_structure_advance_status" -eq 1 ] && \
   [ "$no_restart_contamination_status" -eq 1 ] && \
   [ "$no_stats_visu_contamination_status" -eq 1 ]; then
    io_restart_stats_visu_force_candidate_status=1
fi

write_gate

if [ "$io_restart_stats_visu_force_candidate_status" -eq 1 ]; then
    echo 'STAGE 12.9 IO RESTART STATS VISU FORCE CANDIDATE VERDICT: PASS'
    echo 'STAGE 12.9 FINAL VERDICT: PASS'
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage12_9_failure"
fi

echo 'STAGE 12.9 IO RESTART STATS VISU FORCE CANDIDATE VERDICT: FAIL'
echo 'STAGE 12.9 FINAL VERDICT: FAIL'
echo "Reasons: $REASONS"
exit 1
