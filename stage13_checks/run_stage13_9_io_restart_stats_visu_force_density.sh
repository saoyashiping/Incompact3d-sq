#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
STAGE13_DAT="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate.dat"
OUT_DAT="$OUTPUT_DIR/stage13_9_io_restart_stats_visu_force_density.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_9_io_restart_stats_visu_force_density_reasons.tmp"
STAGE9_7_LOG="$OUTPUT_DIR/stage13_9_stage9_7_stats_visu_io_force_density.log"
STAGE9_8_LOG="$OUTPUT_DIR/stage13_9_stage9_8_restart_io_force_density.log"
STAGE9_7_COPY="$OUTPUT_DIR/stage13_9_stage9_7_stats_visu_io_force_density.dat"
STAGE9_8_COPY="$OUTPUT_DIR/stage13_9_stage9_8_restart_io_force_density.dat"
STAGE9_7_DIAG_COPY="$OUTPUT_DIR/stage13_9_stage9_7_force_density_diagnostics.dat"
STAGE9_8_DIAG_COPY="$OUTPUT_DIR/stage13_9_stage9_8_force_density_diagnostics.dat"

mkdir -p stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$OUT_DAT" "$STAGE13_DAT" "$STAGE9_7_LOG" "$STAGE9_8_LOG" \
      "$STAGE9_7_COPY" "$STAGE9_8_COPY" "$STAGE9_7_DIAG_COPY" "$STAGE9_8_DIAG_COPY"

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

check_key_equals_phase() {
    file=$1
    key=$2
    expected=$3
    phase=$4
    value=$(get_dat_value "$file" "$key")
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$phase $key expected $expected but found ${value:-missing}"
    return 1
}

copy_matching_stage9_files() {
    destination=$1
    shift
    copied=0
    rm -f "$destination"
    for pattern in "$@"; do
        for file in $pattern; do
            if [ -f "$file" ]; then
                {
                    echo "# copied_from $file"
                    cat "$file"
                    echo ""
                } >> "$destination"
                copied=1
            fi
        done
    done
    if [ "$copied" = "0" ]; then
        rm -f "$destination"
    fi
    return 0
}

reset_diag_values() {
    diag_hook_active_status=0
    diag_candidate_computed_status=0
    diag_candidate_finite_status=0
    diag_norm_finite_status=0
    diag_integrated_finite_status=0
    diag_integrated_conservation_status=0
    diag_spreading_input_sign_status=0
    diag_wrong_sign_rejection_status=0
    diag_no_field_modification_status=0
    diag_no_rhs_modification_status=0
    diag_no_rhs_injection_status=0
    diag_no_production_ibm_forcing_status=0
    diag_no_feedback_application_status=0
    diag_no_twoway_force_status=0
    diag_no_structure_advance_status=0
    diag_production_force_density_status=0
}

validate_stage13_diagnostics() {
    phase=$1
    missing_reason=$2
    copy_path=$3
    status=1
    reset_diag_values
    if [ ! -f "$STAGE13_DAT" ]; then
        add_reason "$missing_reason"
        return 1
    fi
    cp "$STAGE13_DAT" "$copy_path"
    check_key_equals_phase "$STAGE13_DAT" stage13_6_requested_flag 1 "$phase" >/dev/null || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_readonly_mode_status 1 "$phase" >/dev/null || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_spreading_readonly_status 1 "$phase" >/dev/null || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_hook_initialized_status 1 "$phase" >/dev/null && diag_hook_active_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_hook_sample_called_status 1 "$phase" >/dev/null || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_sampled_velocity_available_status 1 "$phase" >/dev/null || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_force_density_candidate_computed_status 1 "$phase" >/dev/null && diag_candidate_computed_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_force_density_candidate_finite_status 1 "$phase" >/dev/null && diag_candidate_finite_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_force_density_norm_finite_status 1 "$phase" >/dev/null && diag_norm_finite_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_integrated_force_finite_status 1 "$phase" >/dev/null && diag_integrated_finite_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_integrated_force_conservation_status 1 "$phase" >/dev/null && diag_integrated_conservation_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_spreading_input_sign_status 1 "$phase" >/dev/null && diag_spreading_input_sign_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_wrong_sign_rejection_status 1 "$phase" >/dev/null && diag_wrong_sign_rejection_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_field_modified_status 0 "$phase" >/dev/null && diag_no_field_modification_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_rhs_modified_status 0 "$phase" >/dev/null && diag_no_rhs_modification_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_no_rhs_injection_status 1 "$phase" >/dev/null && diag_no_rhs_injection_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_no_production_ibm_forcing_status 1 "$phase" >/dev/null && diag_no_production_ibm_forcing_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_no_feedback_application_status 1 "$phase" >/dev/null && diag_no_feedback_application_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_no_twoway_force_status 1 "$phase" >/dev/null && diag_no_twoway_force_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_no_structure_advance_status 1 "$phase" >/dev/null && diag_no_structure_advance_status=1 || status=0
    check_key_equals_phase "$STAGE13_DAT" stage13_6_production_force_density_candidate_status 1 "$phase" >/dev/null && diag_production_force_density_status=1 || status=0
    return "$status"
}

write_output_dat() {
    final_status=$1
    {
        echo "stage13_9_requested_flag 1"
        echo "stage13_9_build_status $build_status"
        echo "stage13_9_stage9_7_stats_visu_status $stage9_7_stats_visu_status"
        echo "stage13_9_stage9_8_restart_status $stage9_8_restart_status"
        echo "stage13_9_stage9_7_hook_active_status $stage9_7_hook_active_status"
        echo "stage13_9_stage9_8_hook_active_status $stage9_8_hook_active_status"
        echo "stage13_9_force_density_candidate_computed_status $force_density_candidate_computed_status"
        echo "stage13_9_force_density_candidate_finite_status $force_density_candidate_finite_status"
        echo "stage13_9_force_density_norm_finite_status $force_density_norm_finite_status"
        echo "stage13_9_integrated_force_finite_status $integrated_force_finite_status"
        echo "stage13_9_integrated_force_conservation_status $integrated_force_conservation_status"
        echo "stage13_9_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_9_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_9_no_field_modification_status $no_field_modification_status"
        echo "stage13_9_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage13_9_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_9_no_production_ibm_forcing_status $no_production_ibm_forcing_status"
        echo "stage13_9_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_9_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_9_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_9_no_restart_contamination_status $no_restart_contamination_status"
        echo "stage13_9_no_stats_visu_contamination_status $no_stats_visu_contamination_status"
        echo "stage13_9_io_restart_stats_visu_force_density_status $final_status"
    } > "$OUT_DAT"
}

set_common_diagnostic_statuses() {
    if [ "$stage9_7_candidate_computed_status" = "1" ] && [ "$stage9_8_candidate_computed_status" = "1" ]; then
        force_density_candidate_computed_status=1
    fi
    if [ "$stage9_7_candidate_finite_status" = "1" ] && [ "$stage9_8_candidate_finite_status" = "1" ]; then
        force_density_candidate_finite_status=1
    fi
    if [ "$stage9_7_norm_finite_status" = "1" ] && [ "$stage9_8_norm_finite_status" = "1" ]; then
        force_density_norm_finite_status=1
    fi
    if [ "$stage9_7_integrated_finite_status" = "1" ] && [ "$stage9_8_integrated_finite_status" = "1" ]; then
        integrated_force_finite_status=1
    fi
    if [ "$stage9_7_integrated_conservation_status" = "1" ] && [ "$stage9_8_integrated_conservation_status" = "1" ]; then
        integrated_force_conservation_status=1
    fi
    if [ "$stage9_7_spreading_input_sign_status" = "1" ] && [ "$stage9_8_spreading_input_sign_status" = "1" ]; then
        spreading_input_sign_status=1
    fi
    if [ "$stage9_7_wrong_sign_rejection_status" = "1" ] && [ "$stage9_8_wrong_sign_rejection_status" = "1" ]; then
        wrong_sign_rejection_status=1
    fi
    if [ "$stage9_7_no_field_modification_status" = "1" ] && [ "$stage9_8_no_field_modification_status" = "1" ]; then
        no_field_modification_status=1
    fi
    if [ "$stage9_7_no_rhs_modification_status" = "1" ] && [ "$stage9_8_no_rhs_modification_status" = "1" ]; then
        no_rhs_modification_status=1
    fi
    if [ "$stage9_7_no_rhs_injection_status" = "1" ] && [ "$stage9_8_no_rhs_injection_status" = "1" ]; then
        no_rhs_injection_status=1
    fi
    if [ "$stage9_7_no_production_ibm_forcing_status" = "1" ] && [ "$stage9_8_no_production_ibm_forcing_status" = "1" ]; then
        no_production_ibm_forcing_status=1
    fi
    if [ "$stage9_7_no_feedback_application_status" = "1" ] && [ "$stage9_8_no_feedback_application_status" = "1" ]; then
        no_feedback_application_status=1
    fi
    if [ "$stage9_7_no_twoway_force_status" = "1" ] && [ "$stage9_8_no_twoway_force_status" = "1" ]; then
        no_twoway_force_status=1
    fi
    if [ "$stage9_7_no_structure_advance_status" = "1" ] && [ "$stage9_8_no_structure_advance_status" = "1" ]; then
        no_structure_advance_status=1
    fi
    if [ "$stage9_8_restart_status" = "1" ] && [ "$stage9_8_no_field_modification_status" = "1" ] && [ "$stage9_8_no_rhs_modification_status" = "1" ]; then
        no_restart_contamination_status=1
    fi
    if [ "$stage9_7_stats_visu_status" = "1" ] && [ "$stage9_7_no_field_modification_status" = "1" ] && [ "$stage9_7_no_rhs_modification_status" = "1" ]; then
        no_stats_visu_contamination_status=1
    fi
}

build_status=1
stage9_7_stats_visu_status=0
stage9_8_restart_status=0
stage9_7_hook_active_status=0
stage9_8_hook_active_status=0
force_density_candidate_computed_status=0
force_density_candidate_finite_status=0
force_density_norm_finite_status=0
integrated_force_finite_status=0
integrated_force_conservation_status=0
spreading_input_sign_status=0
wrong_sign_rejection_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_rhs_injection_status=0
no_production_ibm_forcing_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_restart_contamination_status=0
no_stats_visu_contamination_status=0

stage9_7_candidate_computed_status=0
stage9_7_candidate_finite_status=0
stage9_7_norm_finite_status=0
stage9_7_integrated_finite_status=0
stage9_7_integrated_conservation_status=0
stage9_7_spreading_input_sign_status=0
stage9_7_wrong_sign_rejection_status=0
stage9_7_no_field_modification_status=0
stage9_7_no_rhs_modification_status=0
stage9_7_no_rhs_injection_status=0
stage9_7_no_production_ibm_forcing_status=0
stage9_7_no_feedback_application_status=0
stage9_7_no_twoway_force_status=0
stage9_7_no_structure_advance_status=0

stage9_8_candidate_computed_status=0
stage9_8_candidate_finite_status=0
stage9_8_norm_finite_status=0
stage9_8_integrated_finite_status=0
stage9_8_integrated_conservation_status=0
stage9_8_spreading_input_sign_status=0
stage9_8_wrong_sign_rejection_status=0
stage9_8_no_field_modification_status=0
stage9_8_no_rhs_modification_status=0
stage9_8_no_rhs_injection_status=0
stage9_8_no_production_ibm_forcing_status=0
stage9_8_no_feedback_application_status=0
stage9_8_no_twoway_force_status=0
stage9_8_no_structure_advance_status=0

STAGE13_9_RUN_STAGE13_8=${STAGE13_9_RUN_STAGE13_8:-0}
if [ "$STAGE13_9_RUN_STAGE13_8" = "1" ]; then
    if ! sh stage13_checks/run_stage13_8_force_density_parallel_consistency.sh; then
        add_reason "optional Stage 13.8 prerequisite failed"
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
        bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh > "$STAGE9_7_LOG" 2>&1
    if [ $? -eq 0 ] && grep 'STAGE 9.7 FINAL VERDICT: PASS' "$STAGE9_7_LOG" >/dev/null 2>&1; then
        stage9_7_stats_visu_status=1
    else
        add_reason "Stage 9.7 stats/visu/coarse I/O did not report PASS under Stage 13 hook environment"
    fi
    copy_matching_stage9_files "$STAGE9_7_COPY" 'stage9_outputs/*stage9_7*' 'stage9_outputs/*stats*' 'stage9_outputs/*visu*' 'stage9_outputs/*coarse*'
    if validate_stage13_diagnostics stage9_7 missing_stage13_6_force_density_candidate_diagnostics_after_stage9_7 "$STAGE9_7_DIAG_COPY"; then
        :
    fi
    stage9_7_hook_active_status=$diag_hook_active_status
    stage9_7_candidate_computed_status=$diag_candidate_computed_status
    stage9_7_candidate_finite_status=$diag_candidate_finite_status
    stage9_7_norm_finite_status=$diag_norm_finite_status
    stage9_7_integrated_finite_status=$diag_integrated_finite_status
    stage9_7_integrated_conservation_status=$diag_integrated_conservation_status
    stage9_7_spreading_input_sign_status=$diag_spreading_input_sign_status
    stage9_7_wrong_sign_rejection_status=$diag_wrong_sign_rejection_status
    stage9_7_no_field_modification_status=$diag_no_field_modification_status
    stage9_7_no_rhs_modification_status=$diag_no_rhs_modification_status
    stage9_7_no_rhs_injection_status=$diag_no_rhs_injection_status
    stage9_7_no_production_ibm_forcing_status=$diag_no_production_ibm_forcing_status
    stage9_7_no_feedback_application_status=$diag_no_feedback_application_status
    stage9_7_no_twoway_force_status=$diag_no_twoway_force_status
    stage9_7_no_structure_advance_status=$diag_no_structure_advance_status

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
        bash stage9_checks/run_stage9_8_restart_io_regression.sh > "$STAGE9_8_LOG" 2>&1
    if [ $? -eq 0 ] && grep 'STAGE 9.8 FINAL VERDICT: PASS' "$STAGE9_8_LOG" >/dev/null 2>&1; then
        stage9_8_restart_status=1
    else
        add_reason "Stage 9.8 restart I/O regression did not report PASS under Stage 13 hook environment"
    fi
    copy_matching_stage9_files "$STAGE9_8_COPY" 'stage9_outputs/*stage9_8*' 'stage9_outputs/*restart*'
    if validate_stage13_diagnostics stage9_8 missing_stage13_6_force_density_candidate_diagnostics_after_stage9_8 "$STAGE9_8_DIAG_COPY"; then
        :
    fi
    stage9_8_hook_active_status=$diag_hook_active_status
    stage9_8_candidate_computed_status=$diag_candidate_computed_status
    stage9_8_candidate_finite_status=$diag_candidate_finite_status
    stage9_8_norm_finite_status=$diag_norm_finite_status
    stage9_8_integrated_finite_status=$diag_integrated_finite_status
    stage9_8_integrated_conservation_status=$diag_integrated_conservation_status
    stage9_8_spreading_input_sign_status=$diag_spreading_input_sign_status
    stage9_8_wrong_sign_rejection_status=$diag_wrong_sign_rejection_status
    stage9_8_no_field_modification_status=$diag_no_field_modification_status
    stage9_8_no_rhs_modification_status=$diag_no_rhs_modification_status
    stage9_8_no_rhs_injection_status=$diag_no_rhs_injection_status
    stage9_8_no_production_ibm_forcing_status=$diag_no_production_ibm_forcing_status
    stage9_8_no_feedback_application_status=$diag_no_feedback_application_status
    stage9_8_no_twoway_force_status=$diag_no_twoway_force_status
    stage9_8_no_structure_advance_status=$diag_no_structure_advance_status
else
    add_reason "run phase skipped because required build failed"
fi

set_common_diagnostic_statuses

if [ ! -s "$REASONS_FILE" ] && [ "$build_status" = "1" ] && [ "$stage9_7_stats_visu_status" = "1" ] && \
   [ "$stage9_8_restart_status" = "1" ] && [ "$stage9_7_hook_active_status" = "1" ] && \
   [ "$stage9_8_hook_active_status" = "1" ] && [ "$force_density_candidate_computed_status" = "1" ] && \
   [ "$force_density_candidate_finite_status" = "1" ] && [ "$force_density_norm_finite_status" = "1" ] && \
   [ "$integrated_force_finite_status" = "1" ] && [ "$integrated_force_conservation_status" = "1" ] && \
   [ "$spreading_input_sign_status" = "1" ] && [ "$wrong_sign_rejection_status" = "1" ] && \
   [ "$no_field_modification_status" = "1" ] && [ "$no_rhs_modification_status" = "1" ] && \
   [ "$no_rhs_injection_status" = "1" ] && [ "$no_production_ibm_forcing_status" = "1" ] && \
   [ "$no_feedback_application_status" = "1" ] && [ "$no_twoway_force_status" = "1" ] && \
   [ "$no_structure_advance_status" = "1" ] && [ "$no_restart_contamination_status" = "1" ] && \
   [ "$no_stats_visu_contamination_status" = "1" ]; then
    write_output_dat 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.9 IO RESTART STATS VISU FORCE DENSITY VERDICT: PASS"
    echo "STAGE 13.9 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.9 gate failed without a captured reason"
fi
write_output_dat 0
echo "STAGE 13.9 IO RESTART STATS VISU FORCE DENSITY VERDICT: FAIL"
echo "STAGE 13.9 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
