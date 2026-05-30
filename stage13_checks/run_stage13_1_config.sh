#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
DAT_FILE="$OUTPUT_DIR/fibre_stage13_1_config.dat"
GATE_DAT="$OUTPUT_DIR/stage13_1_config_gate.dat"
BUILD_LOG="$OUTPUT_DIR/stage13_1_build.log"
CHECK_LOG="$OUTPUT_DIR/stage13_1_config_check.log"
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

write_gate() {
    {
        echo "stage13_1_gate_requested_flag $requested_flag"
        echo "stage13_1_gate_build_status $build_status"
        echo "stage13_1_gate_config_check_status $config_check_status"
        echo "stage13_1_gate_readonly_mode_status $readonly_mode_status"
        echo "stage13_1_gate_spreading_readonly_status $spreading_readonly_status"
        echo "stage13_1_gate_max_points_status $max_points_status"
        echo "stage13_1_gate_max_eulerian_points_status $max_eulerian_points_status"
        echo "stage13_1_gate_normalization_mode_status $normalization_mode_status"
        echo "stage13_1_gate_no_force_density_allocation_status $no_force_density_allocation_status"
        echo "stage13_1_gate_no_spreading_status $no_spreading_status"
        echo "stage13_1_gate_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_1_gate_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage13_1_gate_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_1_gate_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_1_gate_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_1_gate_no_fluid_field_access_status $no_fluid_field_access_status"
        echo "stage13_1_gate_no_fluid_field_modification_status $no_fluid_field_modification_status"
        echo "stage13_1_gate_status $gate_status"
    } > "$GATE_DAT"
}

mkdir -p "$OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$CHECK_LOG"

requested_flag=1
build_status=0
config_check_status=0
readonly_mode_status=0
spreading_readonly_status=0
max_points_status=0
max_eulerian_points_status=0
normalization_mode_status=0
no_force_density_allocation_status=0
no_spreading_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_fluid_field_access_status=0
no_fluid_field_modification_status=0
gate_status=0

if [ "${STAGE13_1_RUN_STAGE13_0:-0}" = "1" ]; then
    if ! bash stage13_checks/run_stage13_0_preflight_closure_integrity.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage13_0_prerequisite_failed"
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
        fibre_stage13_config_check
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
    exe="$BUILD_DIR/bin/fibre_stage13_config_check"
    if [ ! -x "$exe" ]; then
        exe="$BUILD_DIR/src/fibre_stage13_config_check"
    fi
    if [ ! -x "$exe" ]; then
        add_reason "missing_fibre_stage13_config_check_executable"
    elif X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1 \
         X3D_STAGE13_FORCE_READONLY=1 \
         X3D_STAGE13_SPREADING_READONLY=1 \
         X3D_STAGE13_MAX_POINTS=8 \
         X3D_STAGE13_MAX_EULERIAN_POINTS=64 \
         X3D_STAGE13_SPREADING_NORMALIZATION=conservative \
         "$exe" > "$CHECK_LOG" 2>&1; then
        if grep 'STAGE 13.1 CONFIG VERDICT: PASS' "$CHECK_LOG" >/dev/null 2>&1; then
            config_check_status=1
        else
            add_reason "stage13_1_config_pass_line_missing"
        fi
    else
        add_reason "fibre_stage13_config_check_failed"
    fi

    if [ ! -f "$DAT_FILE" ]; then
        add_reason "missing_fibre_stage13_1_config_dat"
    else
        if check_dat_key stage13_1_requested_flag 1; then :; fi
        if check_dat_key stage13_1_readonly_mode_status 1; then readonly_mode_status=1; fi
        if check_dat_key stage13_1_spreading_readonly_status 1; then spreading_readonly_status=1; fi
        if check_dat_key stage13_1_max_points_status 1; then max_points_status=1; fi
        if check_dat_key stage13_1_max_eulerian_points_status 1; then max_eulerian_points_status=1; fi
        if check_dat_key stage13_1_normalization_mode_status 1; then normalization_mode_status=1; fi
        if check_dat_key stage13_1_no_force_density_allocation_status 1; then no_force_density_allocation_status=1; fi
        if check_dat_key stage13_1_no_spreading_status 1; then no_spreading_status=1; fi
        if check_dat_key stage13_1_no_rhs_injection_status 1; then no_rhs_injection_status=1; fi
        if check_dat_key stage13_1_no_ibm_spreading_status 1; then no_ibm_spreading_status=1; fi
        if check_dat_key stage13_1_no_feedback_application_status 1; then no_feedback_application_status=1; fi
        if check_dat_key stage13_1_no_twoway_force_status 1; then no_twoway_force_status=1; fi
        if check_dat_key stage13_1_no_structure_advance_status 1; then no_structure_advance_status=1; fi
        if check_dat_key stage13_1_no_fluid_field_access_status 1; then no_fluid_field_access_status=1; fi
        if check_dat_key stage13_1_no_fluid_field_modification_status 1; then no_fluid_field_modification_status=1; fi
        if check_dat_key stage13_1_config_status 1; then :; fi
    fi
fi

if [ "$build_status" -eq 1 ] && \
   [ "$config_check_status" -eq 1 ] && \
   [ "$readonly_mode_status" -eq 1 ] && \
   [ "$spreading_readonly_status" -eq 1 ] && \
   [ "$max_points_status" -eq 1 ] && \
   [ "$max_eulerian_points_status" -eq 1 ] && \
   [ "$normalization_mode_status" -eq 1 ] && \
   [ "$no_force_density_allocation_status" -eq 1 ] && \
   [ "$no_spreading_status" -eq 1 ] && \
   [ "$no_rhs_injection_status" -eq 1 ] && \
   [ "$no_ibm_spreading_status" -eq 1 ] && \
   [ "$no_feedback_application_status" -eq 1 ] && \
   [ "$no_twoway_force_status" -eq 1 ] && \
   [ "$no_structure_advance_status" -eq 1 ] && \
   [ "$no_fluid_field_access_status" -eq 1 ] && \
   [ "$no_fluid_field_modification_status" -eq 1 ]; then
    gate_status=1
fi

write_gate

if [ "$gate_status" -eq 1 ]; then
    echo "STAGE 13.1 FINAL VERDICT: PASS"
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage13_1_failure"
fi

echo "STAGE 13.1 FINAL VERDICT: FAIL"
echo "Reasons: $REASONS"
exit 1
