#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
LOG_FILE="$OUTPUT_DIR/stage13_6_production_force_density_candidate_check.log"
SMOKE_LOG="$OUTPUT_DIR/stage13_6_production_force_density_candidate_smoke.log"
CHECK_DAT="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate_check.dat"
PROD_DAT="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate.dat"
GATE_FILE="$OUTPUT_DIR/stage13_6_production_force_density_candidate_gate.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_6_production_force_density_candidate_reasons.tmp"

mkdir -p stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"

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

write_gate_file() {
    gate_status=$1
    {
        echo "stage13_6_gate_requested_flag $requested_flag"
        echo "stage13_6_gate_build_status $build_status"
        echo "stage13_6_gate_standalone_check_status $standalone_check_status"
        echo "stage13_6_gate_production_smoke_status $production_smoke_status"
        echo "stage13_6_gate_hook_active_status $hook_active_status"
        echo "stage13_6_gate_sampled_velocity_available_status $sampled_velocity_available_status"
        echo "stage13_6_gate_force_density_candidate_computed_status $candidate_computed_status"
        echo "stage13_6_gate_force_density_candidate_finite_status $candidate_finite_status"
        echo "stage13_6_gate_force_density_norm_finite_status $norm_finite_status"
        echo "stage13_6_gate_integrated_force_finite_status $integrated_finite_status"
        echo "stage13_6_gate_integrated_force_conservation_status $integrated_conservation_status"
        echo "stage13_6_gate_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_6_gate_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_6_gate_field_modified_status $field_modified_status"
        echo "stage13_6_gate_rhs_modified_status $rhs_modified_status"
        echo "stage13_6_gate_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_6_gate_no_production_ibm_forcing_status $no_production_ibm_forcing_status"
        echo "stage13_6_gate_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_6_gate_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_6_gate_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_6_gate_status $gate_status"
    } > "$GATE_FILE"
}

requested_flag=1
build_status=1
standalone_check_status=0
production_smoke_status=0
hook_active_status=0
sampled_velocity_available_status=0
candidate_computed_status=0
candidate_finite_status=0
norm_finite_status=0
integrated_finite_status=0
integrated_conservation_status=0
spreading_input_sign_status=0
wrong_sign_rejection_status=0
field_modified_status=1
rhs_modified_status=1
no_rhs_injection_status=0
no_production_ibm_forcing_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0

STAGE13_6_RUN_STAGE13_5=${STAGE13_6_RUN_STAGE13_5:-0}
if [ "$STAGE13_6_RUN_STAGE13_5" = "1" ]; then
    if ! sh stage13_checks/run_stage13_5_conservation_sign_audit.sh; then
        add_reason "optional Stage 13.5 prerequisite failed"
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

if [ -x "$BUILD_DIR/bin/fibre_stage13_production_force_density_candidate_check" ]; then
    EXE="$BUILD_DIR/bin/fibre_stage13_production_force_density_candidate_check"
elif [ -x "$BUILD_DIR/src/fibre_stage13_production_force_density_candidate_check" ]; then
    EXE="$BUILD_DIR/src/fibre_stage13_production_force_density_candidate_check"
else
    EXE=""
    add_reason "fibre_stage13_production_force_density_candidate_check executable not found"
fi

if [ -n "$EXE" ]; then
    X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1 \
    X3D_STAGE13_FORCE_READONLY=1 \
    X3D_STAGE13_SPREADING_READONLY=1 \
    X3D_STAGE13_MAX_POINTS=8 \
    X3D_STAGE13_MAX_EULERIAN_POINTS=64 \
    X3D_STAGE13_SPREADING_NORMALIZATION=conservative \
        "$EXE" > "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        add_reason "standalone Stage 13.6 check execution failed"
    fi
else
    : > "$LOG_FILE"
fi

if grep -q "STAGE 13.6 PRODUCTION FORCE DENSITY CANDIDATE CHECK VERDICT: PASS" "$LOG_FILE"; then
    standalone_check_status=1
else
    add_reason "PASS verdict missing from Stage 13.6 standalone check log"
fi

if [ ! -f "$CHECK_DAT" ]; then
    add_reason "$CHECK_DAT missing"
else
    check_key_equals "$CHECK_DAT" stage13_6_check_requested_flag 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_readonly_mode_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_spreading_readonly_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_hook_initialized_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_hook_sample_called_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_sampled_velocity_available_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_force_density_candidate_computed_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_force_density_candidate_finite_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_force_density_norm_finite_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_integrated_force_finite_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_integrated_force_conservation_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_spreading_input_sign_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_wrong_sign_rejection_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_field_unchanged_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_rhs_modified_status 0 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_no_rhs_injection_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_no_production_ibm_forcing_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_no_feedback_application_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_no_twoway_force_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_no_structure_advance_status 1 >/dev/null || true
    check_key_equals "$CHECK_DAT" stage13_6_check_production_force_density_candidate_status 1 >/dev/null || true
fi

rm -f "$PROD_DAT"

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
X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$SMOKE_LOG" 2>&1
if [ $? -eq 0 ]; then
    production_smoke_status=1
else
    add_reason "Stage 9.9 deterministic production smoke failed under Stage 11/12/13 hook environment"
fi

if [ ! -f "$PROD_DAT" ]; then
    add_reason "$PROD_DAT missing"
else
    check_key_equals "$PROD_DAT" stage13_6_requested_flag 1 >/dev/null && requested_flag=1
    check_key_equals "$PROD_DAT" stage13_6_readonly_mode_status 1 >/dev/null || true
    check_key_equals "$PROD_DAT" stage13_6_spreading_readonly_status 1 >/dev/null || true
    check_key_equals "$PROD_DAT" stage13_6_hook_initialized_status 1 >/dev/null && hook_active_status=1
    check_key_equals "$PROD_DAT" stage13_6_hook_sample_called_status 1 >/dev/null || true
    check_key_equals "$PROD_DAT" stage13_6_sampled_velocity_available_status 1 >/dev/null && sampled_velocity_available_status=1
    check_key_equals "$PROD_DAT" stage13_6_force_density_candidate_computed_status 1 >/dev/null && candidate_computed_status=1
    check_key_equals "$PROD_DAT" stage13_6_force_density_candidate_finite_status 1 >/dev/null && candidate_finite_status=1
    check_key_equals "$PROD_DAT" stage13_6_force_density_norm_finite_status 1 >/dev/null && norm_finite_status=1
    check_key_equals "$PROD_DAT" stage13_6_integrated_force_finite_status 1 >/dev/null && integrated_finite_status=1
    check_key_equals "$PROD_DAT" stage13_6_integrated_force_conservation_status 1 >/dev/null && integrated_conservation_status=1
    check_key_equals "$PROD_DAT" stage13_6_spreading_input_sign_status 1 >/dev/null && spreading_input_sign_status=1
    check_key_equals "$PROD_DAT" stage13_6_wrong_sign_rejection_status 1 >/dev/null && wrong_sign_rejection_status=1
    check_key_equals "$PROD_DAT" stage13_6_field_modified_status 0 >/dev/null && field_modified_status=0
    check_key_equals "$PROD_DAT" stage13_6_rhs_modified_status 0 >/dev/null && rhs_modified_status=0
    check_key_equals "$PROD_DAT" stage13_6_no_rhs_injection_status 1 >/dev/null && no_rhs_injection_status=1
    check_key_equals "$PROD_DAT" stage13_6_no_production_ibm_forcing_status 1 >/dev/null && no_production_ibm_forcing_status=1
    check_key_equals "$PROD_DAT" stage13_6_no_feedback_application_status 1 >/dev/null && no_feedback_application_status=1
    check_key_equals "$PROD_DAT" stage13_6_no_twoway_force_status 1 >/dev/null && no_twoway_force_status=1
    check_key_equals "$PROD_DAT" stage13_6_no_structure_advance_status 1 >/dev/null && no_structure_advance_status=1
    check_key_equals "$PROD_DAT" stage13_6_production_force_density_candidate_status 1 >/dev/null || true
fi

if [ ! -s "$REASONS_FILE" ] && [ "$build_status" = "1" ] && [ "$standalone_check_status" = "1" ] && \
   [ "$production_smoke_status" = "1" ]; then
    write_gate_file 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.6 PRODUCTION FORCE DENSITY CANDIDATE VERDICT: PASS"
    echo "STAGE 13.6 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.6 gate failed without a captured reason"
fi
write_gate_file 0
echo "STAGE 13.6 PRODUCTION FORCE DENSITY CANDIDATE VERDICT: FAIL"
echo "STAGE 13.6 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
