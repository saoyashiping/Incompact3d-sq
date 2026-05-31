#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
OUT_DAT="$OUTPUT_DIR/stage13_11_total_closure.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_11_total_closure_reasons.tmp"
STAGE13_DIAG="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate.dat"
STAGE13_CLOSED_FILE=stage13_checks/STAGE13_CLOSED.md

STAGE13_6_LOG="$OUTPUT_DIR/stage13_11_stage13_6_production_force_density_candidate.log"
STAGE13_7_LOG="$OUTPUT_DIR/stage13_11_stage13_7_force_density_invariance_np1.log"
STAGE13_8_LOG="$OUTPUT_DIR/stage13_11_stage13_8_force_density_parallel_consistency.log"
STAGE13_9_LOG="$OUTPUT_DIR/stage13_11_stage13_9_io_restart_stats_visu_force_density.log"
STAGE13_10_LOG="$OUTPUT_DIR/stage13_11_stage13_10_rhs_injection_contamination_audit.log"

mkdir -p stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$OUT_DAT" "$STAGE13_6_LOG" "$STAGE13_7_LOG" "$STAGE13_8_LOG" "$STAGE13_9_LOG" "$STAGE13_10_LOG"

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        if [ -n "${DECOMP2D_ROOT:-}" ]; then
            cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT"
        else
            cmake -S . -B "$BUILD_DIR"
        fi
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

write_output_dat() {
    final_status=$1
    {
        echo "stage13_11_requested_flag 1"
        echo "stage13_11_build_configure_status $configure_status"
        echo "stage13_11_build_status $build_status"
        echo "stage13_11_xcompact3d_build_status $xcompact3d_build_status"
        echo "stage13_11_stage11_check_targets_build_status $stage11_build_status"
        echo "stage13_11_stage12_check_targets_build_status $stage12_build_status"
        echo "stage13_11_stage13_1_6_check_targets_build_status $stage13_build_status"
        echo "stage13_11_stage13_6_hook_smoke_status $stage13_6_hook_smoke_status"
        echo "stage13_11_stage13_7_np1_invariance_status $stage13_7_np1_invariance_status"
        echo "stage13_11_stage13_8_parallel_consistency_status $stage13_8_parallel_consistency_status"
        echo "stage13_11_stage13_9_io_restart_status $stage13_9_io_restart_status"
        echo "stage13_11_stage13_10_contamination_audit_status $stage13_10_contamination_status"
        echo "stage13_11_stage13_6_diagnostic_dat_status $stage13_6_diagnostic_dat_status"
        echo "stage13_11_force_density_candidate_computed_status $force_density_candidate_computed_status"
        echo "stage13_11_force_density_candidate_finite_status $force_density_candidate_finite_status"
        echo "stage13_11_force_density_norm_finite_status $force_density_norm_finite_status"
        echo "stage13_11_integrated_force_finite_status $integrated_force_finite_status"
        echo "stage13_11_integrated_force_conservation_status $integrated_force_conservation_status"
        echo "stage13_11_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_11_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_11_field_modified_status $field_modified_status"
        echo "stage13_11_rhs_modified_status $rhs_modified_status"
        echo "stage13_11_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_11_no_production_ibm_forcing_status $no_production_ibm_forcing_status"
        echo "stage13_11_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_11_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_11_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_11_stage13_closed_file_written_status $closed_file_written_status"
        echo "stage13_11_total_closure_status $final_status"
    } > "$OUT_DAT"
}

build_target() {
    target=$1
    group=$2
    if cmake --build "$BUILD_DIR" --target "$target" -j; then
        return 0
    fi
    build_status=0
    add_reason "build failed for $target"
    case "$group" in
        xcompact3d) xcompact3d_build_status=0 ;;
        stage11) stage11_build_status=0 ;;
        stage12) stage12_build_status=0 ;;
        stage13) stage13_build_status=0 ;;
    esac
    return 1
}

run_stage_gate() {
    script_path=$1
    log_path=$2
    pass_line=$3
    reason_text=$4
    if sh "$script_path" > "$log_path" 2>&1 && grep "$pass_line" "$log_path" >/dev/null 2>&1; then
        return 0
    fi
    add_reason "$reason_text"
    return 1
}

validate_stage13_6_diagnostics() {
    if [ ! -f "$STAGE13_DIAG" ]; then
        add_reason "missing_stage13_6_force_density_candidate_diagnostics"
        return 1
    fi
    stage13_6_diagnostic_dat_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_force_density_candidate_computed_status 1 >/dev/null && force_density_candidate_computed_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_force_density_candidate_finite_status 1 >/dev/null && force_density_candidate_finite_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_force_density_norm_finite_status 1 >/dev/null && force_density_norm_finite_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_integrated_force_finite_status 1 >/dev/null && integrated_force_finite_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_integrated_force_conservation_status 1 >/dev/null && integrated_force_conservation_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_spreading_input_sign_status 1 >/dev/null && spreading_input_sign_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_wrong_sign_rejection_status 1 >/dev/null && wrong_sign_rejection_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_field_modified_status 0 >/dev/null && field_modified_status=0
    check_key_equals "$STAGE13_DIAG" stage13_6_rhs_modified_status 0 >/dev/null && rhs_modified_status=0
    check_key_equals "$STAGE13_DIAG" stage13_6_no_rhs_injection_status 1 >/dev/null && no_rhs_injection_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_no_production_ibm_forcing_status 1 >/dev/null && no_production_ibm_forcing_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_no_feedback_application_status 1 >/dev/null && no_feedback_application_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_no_twoway_force_status 1 >/dev/null && no_twoway_force_status=1
    check_key_equals "$STAGE13_DIAG" stage13_6_no_structure_advance_status 1 >/dev/null && no_structure_advance_status=1
}

write_stage13_closed_file() {
    cat > "$STAGE13_CLOSED_FILE" <<'CLOSED_EOF'
# Stage 13 CLOSED

Stage 13 is CLOSED after Stage 13.11 total closure validation.

The closure gate requires:

- `xcompact3d` build PASS;
- Stage 11 check target builds PASS;
- Stage 12 check target builds PASS;
- Stage 13.1–13.6 check target builds PASS;
- Stage 13.6 production force-density candidate hook smoke PASS;
- Stage 13.7 np=1 no-contamination invariance PASS;
- Stage 13.8 np=1/2/4 parallel consistency PASS;
- Stage 13.9 restart / stats / visu / coarse I/O compatibility PASS;
- Stage 13.10 RHS injection / production IBM forcing / structure contamination audit PASS;
- Stage 13.6 force-density diagnostic `.dat` exists and reports all required no-contamination statuses.

Do not modify closed Stage 13.0–13.10 files without opening a new stage.
CLOSED_EOF
}

configure_status=1
build_status=1
xcompact3d_build_status=1
stage11_build_status=1
stage12_build_status=1
stage13_build_status=1
stage13_6_hook_smoke_status=0
stage13_7_np1_invariance_status=0
stage13_8_parallel_consistency_status=0
stage13_9_io_restart_status=0
stage13_10_contamination_status=0
stage13_6_diagnostic_dat_status=0
force_density_candidate_computed_status=0
force_density_candidate_finite_status=0
force_density_norm_finite_status=0
integrated_force_finite_status=0
integrated_force_conservation_status=0
spreading_input_sign_status=0
wrong_sign_rejection_status=0
field_modified_status=1
rhs_modified_status=1
no_rhs_injection_status=0
no_production_ibm_forcing_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
closed_file_written_status=0

if ! ensure_build_dir; then
    configure_status=0
    build_status=0
    add_reason "build_stage9 configure failed"
fi

build_target xcompact3d xcompact3d || true

for target in \
    fibre_stage11_config_check \
    fibre_stage11_lagrangian_state_check \
    fibre_stage11_grid_metadata_check \
    fibre_stage11_oneway_interpolation_check \
    fibre_stage11_controlled_interpolation_check \
    fibre_stage11_production_oneway_hook_check; do
    build_target "$target" stage11 || true
done

for target in \
    fibre_stage12_config_check \
    fibre_stage12_force_buffer_check \
    fibre_stage12_prescribed_velocity_check \
    fibre_stage12_feedback_formula_check \
    fibre_stage12_sign_convention_audit_check \
    fibre_stage12_power_diagnostics_check \
    fibre_stage12_production_feedback_candidate_check; do
    build_target "$target" stage12 || true
done

for target in \
    fibre_stage13_config_check \
    fibre_stage13_force_density_buffer_check \
    fibre_stage13_spreading_kernel_check \
    fibre_stage13_volume_normalization_audit_check \
    fibre_stage13_conservation_sign_audit_check \
    fibre_stage13_production_force_density_candidate_check; do
    build_target "$target" stage13 || true
done

if [ "$build_status" = "1" ]; then
    STAGE13_6_RUN_STAGE13_5=0 \
        run_stage_gate stage13_checks/run_stage13_6_production_force_density_candidate.sh \
        "$STAGE13_6_LOG" 'STAGE 13.6 FINAL VERDICT: PASS' \
        'Stage 13.6 production force-density candidate hook smoke failed' && stage13_6_hook_smoke_status=1

    STAGE13_7_RUN_STAGE13_6=0 \
        run_stage_gate stage13_checks/run_stage13_7_force_density_invariance_np1.sh \
        "$STAGE13_7_LOG" 'STAGE 13.7 FINAL VERDICT: PASS' \
        'Stage 13.7 np=1 force-density no-contamination invariance failed' && stage13_7_np1_invariance_status=1

    STAGE13_8_RUN_STAGE13_7=0 \
        run_stage_gate stage13_checks/run_stage13_8_force_density_parallel_consistency.sh \
        "$STAGE13_8_LOG" 'STAGE 13.8 FINAL VERDICT: PASS' \
        'Stage 13.8 np=1/2/4 parallel force-density consistency failed' && stage13_8_parallel_consistency_status=1

    STAGE13_9_RUN_STAGE13_8=0 \
        run_stage_gate stage13_checks/run_stage13_9_io_restart_stats_visu_force_density.sh \
        "$STAGE13_9_LOG" 'STAGE 13.9 FINAL VERDICT: PASS' \
        'Stage 13.9 restart / stats / visu / coarse I/O compatibility failed' && stage13_9_io_restart_status=1

    STAGE13_10_RUN_STAGE13_9=0 \
        run_stage_gate stage13_checks/run_stage13_10_rhs_injection_contamination_audit.sh \
        "$STAGE13_10_LOG" 'STAGE 13.10 FINAL VERDICT: PASS' \
        'Stage 13.10 RHS injection / IBM / structure contamination audit failed' && stage13_10_contamination_status=1

    validate_stage13_6_diagnostics
else
    add_reason "Stage 13.11 gate run phase skipped because required build failed"
fi

if [ ! -s "$REASONS_FILE" ] && [ "$configure_status" = "1" ] && [ "$build_status" = "1" ] && \
   [ "$xcompact3d_build_status" = "1" ] && [ "$stage11_build_status" = "1" ] && \
   [ "$stage12_build_status" = "1" ] && [ "$stage13_build_status" = "1" ] && \
   [ "$stage13_6_hook_smoke_status" = "1" ] && [ "$stage13_7_np1_invariance_status" = "1" ] && \
   [ "$stage13_8_parallel_consistency_status" = "1" ] && [ "$stage13_9_io_restart_status" = "1" ] && \
   [ "$stage13_10_contamination_status" = "1" ] && [ "$stage13_6_diagnostic_dat_status" = "1" ] && \
   [ "$force_density_candidate_computed_status" = "1" ] && [ "$force_density_candidate_finite_status" = "1" ] && \
   [ "$force_density_norm_finite_status" = "1" ] && [ "$integrated_force_finite_status" = "1" ] && \
   [ "$integrated_force_conservation_status" = "1" ] && [ "$spreading_input_sign_status" = "1" ] && \
   [ "$wrong_sign_rejection_status" = "1" ] && [ "$field_modified_status" = "0" ] && \
   [ "$rhs_modified_status" = "0" ] && [ "$no_rhs_injection_status" = "1" ] && \
   [ "$no_production_ibm_forcing_status" = "1" ] && [ "$no_feedback_application_status" = "1" ] && \
   [ "$no_twoway_force_status" = "1" ] && [ "$no_structure_advance_status" = "1" ]; then
    write_stage13_closed_file
    closed_file_written_status=1
    write_output_dat 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.11 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.11 closure failed without a captured reason"
fi
write_output_dat 0
echo "STAGE 13.11 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
