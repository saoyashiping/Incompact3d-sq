#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage12_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
LOG_FILE="$OUTPUT_DIR/stage12_6_production_feedback_candidate_check.log"
SMOKE_LOG="$OUTPUT_DIR/stage12_6_stage9_9_production_smoke.log"
DAT_FILE="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate.dat"
CHECK_DAT_FILE="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate_check.dat"
GATE_FILE="$OUTPUT_DIR/stage12_6_production_feedback_candidate_gate.dat"
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
    add_reason "$key expected $expected but found ${value:-missing}"
    return 1
}

check_scalar_le() {
    key=$1
    limit=$2
    file=$3
    value=$(get_value "$key" "$file" 2>/dev/null)
    if [ -z "$value" ]; then
        add_reason "$key missing"
        return 1
    fi
    awk -v value="$value" -v limit="$limit" 'BEGIN { exit !(value + 0.0 <= limit + 0.0) }'
    if [ $? -eq 0 ]; then
        return 0
    fi
    add_reason "$key expected <= $limit but found $value"
    return 1
}

write_gate() {
    {
        echo "stage12_6_gate_requested_flag $gate_requested_flag"
        echo "stage12_6_gate_build_status $gate_build_status"
        echo "stage12_6_gate_standalone_check_status $gate_standalone_check_status"
        echo "stage12_6_gate_production_smoke_status $gate_production_smoke_status"
        echo "stage12_6_gate_hook_active_status $gate_hook_active_status"
        echo "stage12_6_gate_sampled_velocity_available_status $gate_sampled_velocity_available_status"
        echo "stage12_6_gate_force_candidate_computed_status $gate_force_candidate_computed_status"
        echo "stage12_6_gate_force_candidate_finite_status $gate_force_candidate_finite_status"
        echo "stage12_6_gate_force_norm_finite_status $gate_force_norm_finite_status"
        echo "stage12_6_gate_power_diagnostics_finite_status $gate_power_diagnostics_finite_status"
        echo "stage12_6_gate_action_reaction_status $gate_action_reaction_status"
        echo "stage12_6_gate_pair_power_consistency_status $gate_pair_power_consistency_status"
        echo "stage12_6_gate_field_modified_status $gate_field_modified_status"
        echo "stage12_6_gate_rhs_modified_status $gate_rhs_modified_status"
        echo "stage12_6_gate_no_eulerian_force_density_status $gate_no_eulerian_force_density_status"
        echo "stage12_6_gate_no_rhs_injection_status $gate_no_rhs_injection_status"
        echo "stage12_6_gate_no_ibm_spreading_status $gate_no_ibm_spreading_status"
        echo "stage12_6_gate_no_feedback_application_status $gate_no_feedback_application_status"
        echo "stage12_6_gate_no_twoway_force_status $gate_no_twoway_force_status"
        echo "stage12_6_gate_no_structure_advance_status $gate_no_structure_advance_status"
        echo "stage12_6_gate_status $gate_status"
    } > "$GATE_FILE"
}

mkdir -p "$OUTPUT_DIR" "$STAGE11_OUTPUT_DIR" "$STAGE9_OUTPUT_DIR"
: > "$LOG_FILE"
: > "$SMOKE_LOG"

gate_requested_flag=1
gate_build_status=0
gate_standalone_check_status=0
gate_production_smoke_status=0
gate_hook_active_status=0
gate_sampled_velocity_available_status=0
gate_force_candidate_computed_status=0
gate_force_candidate_finite_status=0
gate_force_norm_finite_status=0
gate_power_diagnostics_finite_status=0
gate_action_reaction_status=0
gate_pair_power_consistency_status=0
gate_field_modified_status=1
gate_rhs_modified_status=1
gate_no_eulerian_force_density_status=0
gate_no_rhs_injection_status=0
gate_no_ibm_spreading_status=0
gate_no_feedback_application_status=0
gate_no_twoway_force_status=0
gate_no_structure_advance_status=0
gate_status=0

if [ "${STAGE12_6_RUN_STAGE12_5:-0}" = "1" ]; then
    if ! bash stage12_checks/run_stage12_5_power_diagnostics.sh >> "$LOG_FILE" 2>&1; then
        add_reason "optional Stage 12.5 prerequisite failed"
    fi
fi

if ensure_build_dir >> "$LOG_FILE" 2>&1; then
    build_ok=1
else
    build_ok=0
    add_reason "cmake configure failed"
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
        if ! cmake --build "$BUILD_DIR" --target "$target" -j >> "$LOG_FILE" 2>&1; then
            build_ok=0
            add_reason "build failed for $target"
        fi
    done
fi

gate_build_status=$build_ok

exe="$BUILD_DIR/bin/fibre_stage12_production_feedback_candidate_check"
if [ ! -x "$exe" ]; then
    exe="$BUILD_DIR/src/fibre_stage12_production_feedback_candidate_check"
fi

if [ "$build_ok" -eq 1 ] && [ -x "$exe" ]; then
    if X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
       X3D_STAGE12_FORCE_READONLY=1 \
       X3D_STAGE12_FEEDBACK_GAIN=1.0 \
       X3D_STAGE12_FORCE_SIGN=1 \
       X3D_STAGE12_MAX_POINTS=8 \
       "$exe" >> "$LOG_FILE" 2>&1; then
        if grep -q 'STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE CHECK VERDICT: PASS' "$LOG_FILE"; then
            gate_standalone_check_status=1
        else
            add_reason "standalone check PASS verdict missing"
        fi
    else
        add_reason "standalone Stage 12.6 check executable failed"
    fi
else
    add_reason "standalone Stage 12.6 check executable not available"
fi

if [ "$build_ok" -eq 1 ]; then
    if X3D_STAGE11_ONEWAY_HOOK=1 \
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
       bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh >> "$SMOKE_LOG" 2>&1; then
        gate_production_smoke_status=1
    else
        add_reason "Stage 9.9 deterministic no-fibre production smoke failed"
    fi
else
    add_reason "production smoke skipped because build failed"
fi

if [ -f "$CHECK_DAT_FILE" ]; then
    check_key stage12_6_check_requested_flag 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_readonly_mode_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_hook_initialized_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_hook_sample_called_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_sampled_velocity_available_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_force_candidate_computed_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_force_candidate_finite_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_force_norm_finite_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_power_diagnostics_finite_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_action_reaction_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_pair_power_consistency_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_field_unchanged_status 1 "$CHECK_DAT_FILE" >/dev/null
    check_key stage12_6_check_rhs_modified_status 0 "$CHECK_DAT_FILE" >/dev/null
else
    add_reason "standalone check dat file missing"
fi

if [ -f "$DAT_FILE" ]; then
    hook_initialized_ok=0
    hook_sample_called_ok=0
    check_key stage12_6_requested_flag 1 "$DAT_FILE" >/dev/null
    check_key stage12_6_readonly_mode_status 1 "$DAT_FILE" >/dev/null
    check_key stage12_6_hook_initialized_status 1 "$DAT_FILE" && hook_initialized_ok=1
    check_key stage12_6_hook_sample_called_status 1 "$DAT_FILE" && hook_sample_called_ok=1
    if [ "$hook_initialized_ok" -eq 1 ] && [ "$hook_sample_called_ok" -eq 1 ]; then
        gate_hook_active_status=1
    fi
    check_key stage12_6_sampled_velocity_available_status 1 "$DAT_FILE" && gate_sampled_velocity_available_status=1
    check_key stage12_6_force_candidate_computed_status 1 "$DAT_FILE" && gate_force_candidate_computed_status=1
    check_key stage12_6_force_candidate_finite_status 1 "$DAT_FILE" && gate_force_candidate_finite_status=1
    check_key stage12_6_force_norm_finite_status 1 "$DAT_FILE" && gate_force_norm_finite_status=1
    check_key stage12_6_power_diagnostics_finite_status 1 "$DAT_FILE" && gate_power_diagnostics_finite_status=1
    check_key stage12_6_action_reaction_status 1 "$DAT_FILE" && gate_action_reaction_status=1
    check_key stage12_6_pair_power_consistency_status 1 "$DAT_FILE" && gate_pair_power_consistency_status=1
    check_key stage12_6_field_modified_status 0 "$DAT_FILE" && gate_field_modified_status=0
    check_key stage12_6_rhs_modified_status 0 "$DAT_FILE" && gate_rhs_modified_status=0
    check_key stage12_6_no_eulerian_force_density_status 1 "$DAT_FILE" && gate_no_eulerian_force_density_status=1
    check_key stage12_6_no_rhs_injection_status 1 "$DAT_FILE" && gate_no_rhs_injection_status=1
    check_key stage12_6_no_ibm_spreading_status 1 "$DAT_FILE" && gate_no_ibm_spreading_status=1
    check_key stage12_6_no_feedback_application_status 1 "$DAT_FILE" && gate_no_feedback_application_status=1
    check_key stage12_6_no_twoway_force_status 1 "$DAT_FILE" && gate_no_twoway_force_status=1
    check_key stage12_6_no_structure_advance_status 1 "$DAT_FILE" && gate_no_structure_advance_status=1
    check_key stage12_6_production_feedback_candidate_status 1 "$DAT_FILE" >/dev/null
    check_scalar_le stage12_6_pair_power_consistency_error 1.0e-12 "$DAT_FILE" >/dev/null
    check_scalar_le stage12_6_action_reaction_error 1.0e-12 "$DAT_FILE" >/dev/null
else
    add_reason "production Stage 12.6 diagnostics dat file missing"
fi

if [ "$gate_build_status" -eq 1 ] && \
   [ "$gate_standalone_check_status" -eq 1 ] && \
   [ "$gate_production_smoke_status" -eq 1 ] && \
   [ "$gate_hook_active_status" -eq 1 ] && \
   [ "$gate_sampled_velocity_available_status" -eq 1 ] && \
   [ "$gate_force_candidate_computed_status" -eq 1 ] && \
   [ "$gate_force_candidate_finite_status" -eq 1 ] && \
   [ "$gate_force_norm_finite_status" -eq 1 ] && \
   [ "$gate_power_diagnostics_finite_status" -eq 1 ] && \
   [ "$gate_action_reaction_status" -eq 1 ] && \
   [ "$gate_pair_power_consistency_status" -eq 1 ] && \
   [ "$gate_field_modified_status" -eq 0 ] && \
   [ "$gate_rhs_modified_status" -eq 0 ] && \
   [ "$gate_no_eulerian_force_density_status" -eq 1 ] && \
   [ "$gate_no_rhs_injection_status" -eq 1 ] && \
   [ "$gate_no_ibm_spreading_status" -eq 1 ] && \
   [ "$gate_no_feedback_application_status" -eq 1 ] && \
   [ "$gate_no_twoway_force_status" -eq 1 ] && \
   [ "$gate_no_structure_advance_status" -eq 1 ]; then
    gate_status=1
fi

write_gate

if [ "$gate_status" -eq 1 ]; then
    echo 'STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE VERDICT: PASS'
    echo 'STAGE 12.6 FINAL VERDICT: PASS'
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown Stage 12.6 gate failure"
fi

echo 'STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE VERDICT: FAIL'
echo 'STAGE 12.6 FINAL VERDICT: FAIL'
echo "Reasons: $REASONS"
exit 1
