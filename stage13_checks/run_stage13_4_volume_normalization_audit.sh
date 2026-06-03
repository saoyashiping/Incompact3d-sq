#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
DAT_FILE="$OUTPUT_DIR/fibre_stage13_4_volume_normalization_audit.dat"
GATE_DAT="$OUTPUT_DIR/stage13_4_volume_normalization_audit_gate.dat"
BUILD_LOG="$OUTPUT_DIR/stage13_4_build.log"
CHECK_LOG="$OUTPUT_DIR/stage13_4_volume_normalization_audit_check.log"
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

check_dat_key() {
    key=$1
    expected=$2
    value=$(get_value "$key" "$DAT_FILE" 2>/dev/null)
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$key"
    return 1
}

check_scalar_le() {
    key=$1
    limit=$2
    value=$(get_value "$key" "$DAT_FILE" 2>/dev/null)
    if [ -n "$value" ] && awk -v value="$value" -v limit="$limit" 'BEGIN { exit !((value + 0.0) <= (limit + 0.0)) }'; then
        return 0
    fi
    add_reason "$key"
    return 1
}

write_gate() {
    {
        echo "stage13_4_gate_requested_flag $requested_flag"
        echo "stage13_4_gate_build_status $build_status"
        echo "stage13_4_gate_volume_normalization_audit_check_status $volume_normalization_audit_check_status"
        echo "stage13_4_gate_positive_volume_status $positive_volume_status"
        echo "stage13_4_gate_finite_volume_status $finite_volume_status"
        echo "stage13_4_gate_invalid_zero_volume_rejection_status $invalid_zero_volume_rejection_status"
        echo "stage13_4_gate_invalid_negative_volume_rejection_status $invalid_negative_volume_rejection_status"
        echo "stage13_4_gate_uniform_volume_conservation_status $uniform_volume_conservation_status"
        echo "stage13_4_gate_nonuniform_volume_conservation_status $nonuniform_volume_conservation_status"
        echo "stage13_4_gate_density_times_volume_integral_status $density_times_volume_integral_status"
        echo "stage13_4_gate_component_volume_normalization_status $component_volume_normalization_status"
        echo "stage13_4_gate_boundary_volume_normalization_status $boundary_volume_normalization_status"
        echo "stage13_4_gate_finite_force_density_status $finite_force_density_status"
        echo "stage13_4_gate_clear_status $clear_status"
        echo "stage13_4_gate_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_4_gate_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage13_4_gate_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_4_gate_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_4_gate_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_4_gate_no_fluid_field_access_status $no_fluid_field_access_status"
        echo "stage13_4_gate_no_fluid_field_modification_status $no_fluid_field_modification_status"
        echo "stage13_4_gate_status $gate_status"
    } > "$GATE_DAT"
}

mkdir -p "$OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$CHECK_LOG"

requested_flag=1
build_status=0
volume_normalization_audit_check_status=0
initialized_status=0
positive_volume_status=0
finite_volume_status=0
invalid_zero_volume_rejection_status=0
invalid_negative_volume_rejection_status=0
uniform_volume_conservation_status=0
nonuniform_volume_conservation_status=0
density_times_volume_integral_status=0
component_volume_normalization_status=0
boundary_volume_normalization_status=0
finite_force_density_status=0
clear_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_fluid_field_access_status=0
no_fluid_field_modification_status=0
scalar_error_status=0
gate_status=0

if [ "${STAGE13_4_RUN_STAGE13_3:-0}" = "1" ]; then
    if ! bash stage13_checks/run_stage13_3_spreading_kernel.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage13_3_prerequisite_failed"
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
        fibre_stage12_production_feedback_candidate_check \
        fibre_stage13_config_check \
        fibre_stage13_force_density_buffer_check \
        fibre_stage13_spreading_kernel_check \
        fibre_stage13_volume_normalization_audit_check
    do
        if ! cmake --build "$BUILD_DIR" --target "$target" -j >> "$BUILD_LOG" 2>&1; then
            build_ok=0
            add_reason "build_failed_$target"
        fi
    done
fi
build_status=$build_ok

if [ "$build_status" -eq 1 ]; then
    rm -f "$DAT_FILE"
    exe="$BUILD_DIR/bin/fibre_stage13_volume_normalization_audit_check"
    if [ ! -x "$exe" ]; then
        exe="$BUILD_DIR/src/fibre_stage13_volume_normalization_audit_check"
    fi
    if [ ! -x "$exe" ]; then
        add_reason "missing_fibre_stage13_volume_normalization_audit_check_executable"
    elif X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1 \
         X3D_STAGE13_FORCE_READONLY=1 \
         X3D_STAGE13_SPREADING_READONLY=1 \
         X3D_STAGE13_MAX_POINTS=8 \
         X3D_STAGE13_MAX_EULERIAN_POINTS=64 \
         X3D_STAGE13_SPREADING_NORMALIZATION=conservative \
         "$exe" > "$CHECK_LOG" 2>&1; then
        if grep 'STAGE 13.4 VOLUME NORMALIZATION AUDIT VERDICT: PASS' "$CHECK_LOG" >/dev/null 2>&1; then
            volume_normalization_audit_check_status=1
        else
            add_reason "stage13_4_volume_normalization_audit_pass_line_missing"
        fi
    else
        add_reason "fibre_stage13_volume_normalization_audit_check_failed"
    fi

    if [ ! -f "$DAT_FILE" ]; then
        add_reason "missing_fibre_stage13_4_volume_normalization_audit_dat"
    else
        if check_dat_key stage13_4_requested_flag 1; then :; fi
        if check_dat_key stage13_4_readonly_mode_status 1; then :; fi
        if check_dat_key stage13_4_spreading_readonly_status 1; then :; fi
        if check_dat_key stage13_4_initialized_status 1; then initialized_status=1; fi
        if check_dat_key stage13_4_positive_volume_status 1; then positive_volume_status=1; fi
        if check_dat_key stage13_4_finite_volume_status 1; then finite_volume_status=1; fi
        if check_dat_key stage13_4_invalid_zero_volume_rejection_status 1; then invalid_zero_volume_rejection_status=1; fi
        if check_dat_key stage13_4_invalid_negative_volume_rejection_status 1; then invalid_negative_volume_rejection_status=1; fi
        if check_dat_key stage13_4_uniform_volume_conservation_status 1; then uniform_volume_conservation_status=1; fi
        if check_dat_key stage13_4_nonuniform_volume_conservation_status 1; then nonuniform_volume_conservation_status=1; fi
        if check_dat_key stage13_4_density_times_volume_integral_status 1; then density_times_volume_integral_status=1; fi
        if check_dat_key stage13_4_component_volume_normalization_status 1; then component_volume_normalization_status=1; fi
        if check_dat_key stage13_4_boundary_volume_normalization_status 1; then boundary_volume_normalization_status=1; fi
        if check_dat_key stage13_4_finite_force_density_status 1; then finite_force_density_status=1; fi
        if check_dat_key stage13_4_clear_status 1; then clear_status=1; fi
        if check_dat_key stage13_4_no_rhs_injection_status 1; then no_rhs_injection_status=1; fi
        if check_dat_key stage13_4_no_ibm_spreading_status 1; then no_ibm_spreading_status=1; fi
        if check_dat_key stage13_4_no_feedback_application_status 1; then no_feedback_application_status=1; fi
        if check_dat_key stage13_4_no_twoway_force_status 1; then no_twoway_force_status=1; fi
        if check_dat_key stage13_4_no_structure_advance_status 1; then no_structure_advance_status=1; fi
        if check_dat_key stage13_4_no_fluid_field_access_status 1; then no_fluid_field_access_status=1; fi
        if check_dat_key stage13_4_no_fluid_field_modification_status 1; then no_fluid_field_modification_status=1; fi
        if check_dat_key stage13_4_volume_normalization_audit_status 1; then :; fi

        scalar_error_status=1
        if ! check_scalar_le stage13_4_uniform_integrated_force_error_l2 1.0e-12; then scalar_error_status=0; fi
        if ! check_scalar_le stage13_4_nonuniform_integrated_force_error_l2 1.0e-12; then scalar_error_status=0; fi
        if ! check_scalar_le stage13_4_boundary_integrated_force_error_l2 1.0e-12; then scalar_error_status=0; fi
        if ! check_scalar_le stage13_4_max_abs_force_density_norm_after_clear 1.0e-12; then scalar_error_status=0; fi
    fi
fi

if [ -z "$REASONS" ] && \
   [ "$build_status" -eq 1 ] && \
   [ "$volume_normalization_audit_check_status" -eq 1 ] && \
   [ "$initialized_status" -eq 1 ] && \
   [ "$positive_volume_status" -eq 1 ] && \
   [ "$finite_volume_status" -eq 1 ] && \
   [ "$invalid_zero_volume_rejection_status" -eq 1 ] && \
   [ "$invalid_negative_volume_rejection_status" -eq 1 ] && \
   [ "$uniform_volume_conservation_status" -eq 1 ] && \
   [ "$nonuniform_volume_conservation_status" -eq 1 ] && \
   [ "$density_times_volume_integral_status" -eq 1 ] && \
   [ "$component_volume_normalization_status" -eq 1 ] && \
   [ "$boundary_volume_normalization_status" -eq 1 ] && \
   [ "$finite_force_density_status" -eq 1 ] && \
   [ "$clear_status" -eq 1 ] && \
   [ "$no_rhs_injection_status" -eq 1 ] && \
   [ "$no_ibm_spreading_status" -eq 1 ] && \
   [ "$no_feedback_application_status" -eq 1 ] && \
   [ "$no_twoway_force_status" -eq 1 ] && \
   [ "$no_structure_advance_status" -eq 1 ] && \
   [ "$no_fluid_field_access_status" -eq 1 ] && \
   [ "$no_fluid_field_modification_status" -eq 1 ] && \
   [ "$scalar_error_status" -eq 1 ]; then
    gate_status=1
fi

write_gate

if [ "$gate_status" -eq 1 ]; then
    echo "STAGE 13.4 FINAL VERDICT: PASS"
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage13_4_failure"
fi

echo "STAGE 13.4 FINAL VERDICT: FAIL"
echo "Reasons: $REASONS"
exit 1
