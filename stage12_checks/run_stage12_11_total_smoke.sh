#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage12_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
STAGE12_DIAG="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate.dat"
STAGE12_10_DAT="$OUTPUT_DIR/stage12_10_rhs_spreading_structure_contamination_audit.dat"
GATE_DAT="$OUTPUT_DIR/stage12_11_total_smoke.dat"
BUILD_LOG="$OUTPUT_DIR/stage12_11_build.log"
STAGE12_6_LOG="$OUTPUT_DIR/stage12_11_stage12_6.log"
STAGE12_7_LOG="$OUTPUT_DIR/stage12_11_stage12_7.log"
STAGE12_8_LOG="$OUTPUT_DIR/stage12_11_stage12_8.log"
STAGE12_9_LOG="$OUTPUT_DIR/stage12_11_stage12_9.log"
STAGE12_10_LOG="$OUTPUT_DIR/stage12_11_stage12_10.log"
CLOSED_FILE=stage12_checks/STAGE12_CLOSED.md
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
    reason_prefix=$4
    value=$(get_value "$key" "$file" 2>/dev/null)
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "${reason_prefix}${key}"
    return 1
}

write_gate() {
    {
        echo "stage12_11_requested_flag $requested_flag"
        echo "stage12_11_build_status $build_status"
        echo "stage12_11_stage12_6_status $stage12_6_status"
        echo "stage12_11_stage12_7_status $stage12_7_status"
        echo "stage12_11_stage12_8_status $stage12_8_status"
        echo "stage12_11_stage12_9_status $stage12_9_status"
        echo "stage12_11_stage12_10_status $stage12_10_status"
        echo "stage12_11_hook_active_status $hook_active_status"
        echo "stage12_11_force_candidate_computed_status $force_candidate_computed_status"
        echo "stage12_11_force_candidate_finite_status $force_candidate_finite_status"
        echo "stage12_11_force_norm_finite_status $force_norm_finite_status"
        echo "stage12_11_power_diagnostics_finite_status $power_diagnostics_finite_status"
        echo "stage12_11_action_reaction_status $action_reaction_status"
        echo "stage12_11_pair_power_consistency_status $pair_power_consistency_status"
        echo "stage12_11_no_field_modification_status $no_field_modification_status"
        echo "stage12_11_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage12_11_no_eulerian_force_density_status $no_eulerian_force_density_status"
        echo "stage12_11_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage12_11_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage12_11_no_feedback_application_status $no_feedback_application_status"
        echo "stage12_11_no_twoway_force_status $no_twoway_force_status"
        echo "stage12_11_no_structure_advance_status $no_structure_advance_status"
        echo "stage12_11_total_smoke_status $total_smoke_status"
    } > "$GATE_DAT"
}


verify_stage12_diag() {
    diag_ok=1
    requested_ok=0
    readonly_ok=0
    initialized_ok=0
    sample_called_ok=0

    if [ ! -f "$STAGE12_DIAG" ]; then
        add_reason "missing_stage12_6_feedback_candidate_diagnostics"
        return 1
    fi

    if check_key stage12_6_requested_flag 1 "$STAGE12_DIAG" stage12_diag_; then requested_ok=1; else diag_ok=0; fi
    if check_key stage12_6_readonly_mode_status 1 "$STAGE12_DIAG" stage12_diag_; then readonly_ok=1; else diag_ok=0; fi
    if check_key stage12_6_hook_initialized_status 1 "$STAGE12_DIAG" stage12_diag_; then initialized_ok=1; else diag_ok=0; fi
    if check_key stage12_6_hook_sample_called_status 1 "$STAGE12_DIAG" stage12_diag_; then sample_called_ok=1; else diag_ok=0; fi
    if check_key stage12_6_sampled_velocity_available_status 1 "$STAGE12_DIAG" stage12_diag_; then :; else diag_ok=0; fi
    if check_key stage12_6_force_candidate_computed_status 1 "$STAGE12_DIAG" stage12_diag_; then force_candidate_computed_status=1; else diag_ok=0; fi
    if check_key stage12_6_force_candidate_finite_status 1 "$STAGE12_DIAG" stage12_diag_; then force_candidate_finite_status=1; else diag_ok=0; fi
    if check_key stage12_6_force_norm_finite_status 1 "$STAGE12_DIAG" stage12_diag_; then force_norm_finite_status=1; else diag_ok=0; fi
    if check_key stage12_6_power_diagnostics_finite_status 1 "$STAGE12_DIAG" stage12_diag_; then power_diagnostics_finite_status=1; else diag_ok=0; fi
    if check_key stage12_6_action_reaction_status 1 "$STAGE12_DIAG" stage12_diag_; then action_reaction_status=1; else diag_ok=0; fi
    if check_key stage12_6_pair_power_consistency_status 1 "$STAGE12_DIAG" stage12_diag_; then pair_power_consistency_status=1; else diag_ok=0; fi
    if check_key stage12_6_field_modified_status 0 "$STAGE12_DIAG" stage12_diag_; then no_field_modification_status=1; else diag_ok=0; fi
    if check_key stage12_6_rhs_modified_status 0 "$STAGE12_DIAG" stage12_diag_; then no_rhs_modification_status=1; else diag_ok=0; fi
    if check_key stage12_6_no_eulerian_force_density_status 1 "$STAGE12_DIAG" stage12_diag_; then no_eulerian_force_density_status=1; else diag_ok=0; fi
    if check_key stage12_6_no_rhs_injection_status 1 "$STAGE12_DIAG" stage12_diag_; then no_rhs_injection_status=1; else diag_ok=0; fi
    if check_key stage12_6_no_ibm_spreading_status 1 "$STAGE12_DIAG" stage12_diag_; then no_ibm_spreading_status=1; else diag_ok=0; fi
    if check_key stage12_6_no_feedback_application_status 1 "$STAGE12_DIAG" stage12_diag_; then no_feedback_application_status=1; else diag_ok=0; fi
    if check_key stage12_6_no_twoway_force_status 1 "$STAGE12_DIAG" stage12_diag_; then no_twoway_force_status=1; else diag_ok=0; fi
    if check_key stage12_6_no_structure_advance_status 1 "$STAGE12_DIAG" stage12_diag_; then no_structure_advance_status=1; else diag_ok=0; fi
    if check_key stage12_6_production_feedback_candidate_status 1 "$STAGE12_DIAG" stage12_diag_; then :; else diag_ok=0; fi

    if [ "$requested_ok" -eq 1 ] && [ "$readonly_ok" -eq 1 ] && \
       [ "$initialized_ok" -eq 1 ] && [ "$sample_called_ok" -eq 1 ]; then
        hook_active_status=1
    fi

    if [ "$diag_ok" -eq 1 ]; then
        return 0
    fi
    return 1
}

verify_stage12_10_dat() {
    audit_ok=1

    if [ ! -f "$STAGE12_10_DAT" ]; then
        add_reason "missing_stage12_10_contamination_audit_dat"
        return 1
    fi

    for key in \
        stage12_10_static_audit_status \
        stage12_10_velocity_intent_readonly_status \
        stage12_10_no_velocity_write_static_status \
        stage12_10_no_pressure_write_static_status \
        stage12_10_no_rhs_write_static_status \
        stage12_10_no_eulerian_force_density_static_status \
        stage12_10_no_rhs_injection_static_status \
        stage12_10_no_ibm_spreading_static_status \
        stage12_10_no_feedback_application_static_status \
        stage12_10_no_twoway_force_static_status \
        stage12_10_no_structure_advance_static_status \
        stage12_10_runtime_smoke_status \
        stage12_10_hook_active_status \
        stage12_10_force_candidate_computed_status \
        stage12_10_force_candidate_finite_status \
        stage12_10_force_norm_finite_status \
        stage12_10_power_diagnostics_finite_status \
        stage12_10_no_field_modification_status \
        stage12_10_no_rhs_modification_status \
        stage12_10_no_eulerian_force_density_status \
        stage12_10_no_rhs_injection_status \
        stage12_10_no_ibm_spreading_status \
        stage12_10_no_feedback_application_status \
        stage12_10_no_twoway_force_status \
        stage12_10_no_structure_advance_status \
        stage12_10_rhs_spreading_structure_contamination_audit_status
    do
        if ! check_key "$key" 1 "$STAGE12_10_DAT" stage12_10_audit_; then
            audit_ok=0
        fi
    done

    if [ "$audit_ok" -eq 1 ]; then
        return 0
    fi
    return 1
}

generate_stage12_closed() {
    cat > "$CLOSED_FILE" <<'CLOSED_EOF'
# Stage 12 Closed

Stage 12 is closed after the Stage 12.11 total smoke gate passed.

## Purpose

Stage 12 established the production Lagrangian feedback force candidate diagnostic path. It computes and audits the Lagrangian force candidate without enabling Eulerian force density, IBM spreading, RHS injection, two-way force density, or fibre structure advance.

## Closed sub-stages

- Stage 12.0 config and global switches.
- Stage 12.1 Lagrangian force candidate buffer skeleton.
- Stage 12.2 prescribed fibre/control-point velocity skeleton.
- Stage 12.3 controlled feedback force formula checks.
- Stage 12.4 feedback sign convention audit.
- Stage 12.5 force work / power diagnostic candidate.
- Stage 12.6 production feedback candidate hook skeleton.
- Stage 12.7 np=1 force-candidate no-contamination invariance.
- Stage 12.8 np=1/2/4 parallel force-candidate consistency.
- Stage 12.9 restart / stats / visu / coarse I/O compatibility.
- Stage 12.10 RHS / spreading / structure contamination audit.
- Stage 12.11 total smoke closure.

## Governing diagnostic model

```text
U_f = I_h[u](X_f)
slip = U_f - V_f
F_fluid_to_fibre_cand = alpha * slip
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
P_pair + P_slip = 0
f_fsi = 0
RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Closure statement

Stage 12 closes the Lagrangian feedback force candidate diagnostic path only. No Eulerian force density, no IBM spreading, no RHS injection, no two-way force density, and no fibre structure advance are activated.

Real Eulerian force-density construction begins in Stage 13.

## Next recommended stage

Stage 13 production two-way force-density RHS candidate.
CLOSED_EOF
}

mkdir -p "$OUTPUT_DIR" "$STAGE11_OUTPUT_DIR" "$STAGE9_OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$STAGE12_6_LOG"
: > "$STAGE12_7_LOG"
: > "$STAGE12_8_LOG"
: > "$STAGE12_9_LOG"
: > "$STAGE12_10_LOG"

requested_flag=1
build_status=0
stage12_6_status=0
stage12_7_status=0
stage12_8_status=0
stage12_9_status=0
stage12_10_status=0
hook_active_status=0
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
total_smoke_status=0
final_diag_status=0
stage12_10_dat_status=0

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
    if STAGE12_6_RUN_STAGE12_5=0 bash stage12_checks/run_stage12_6_production_feedback_candidate.sh > "$STAGE12_6_LOG" 2>&1 && \
       grep 'STAGE 12.6 FINAL VERDICT: PASS' "$STAGE12_6_LOG" >/dev/null 2>&1; then
        stage12_6_status=1
    else
        add_reason "stage12_6_total_smoke_failed"
    fi

    if STAGE12_7_RUN_STAGE12_6=0 bash stage12_checks/run_stage12_7_force_candidate_invariance_np1.sh > "$STAGE12_7_LOG" 2>&1 && \
       grep 'STAGE 12.7 FINAL VERDICT: PASS' "$STAGE12_7_LOG" >/dev/null 2>&1; then
        stage12_7_status=1
    else
        add_reason "stage12_7_total_smoke_failed"
    fi

    if STAGE12_8_RUN_STAGE12_7=0 bash stage12_checks/run_stage12_8_force_candidate_parallel_consistency.sh > "$STAGE12_8_LOG" 2>&1 && \
       grep 'STAGE 12.8 FINAL VERDICT: PASS' "$STAGE12_8_LOG" >/dev/null 2>&1; then
        stage12_8_status=1
    else
        add_reason "stage12_8_total_smoke_failed"
    fi

    if STAGE12_9_RUN_STAGE12_8=0 bash stage12_checks/run_stage12_9_io_restart_stats_visu_force_candidate.sh > "$STAGE12_9_LOG" 2>&1 && \
       grep 'STAGE 12.9 FINAL VERDICT: PASS' "$STAGE12_9_LOG" >/dev/null 2>&1; then
        stage12_9_status=1
    else
        add_reason "stage12_9_total_smoke_failed"
    fi

    if STAGE12_10_RUN_STAGE12_9=0 bash stage12_checks/run_stage12_10_rhs_spreading_structure_contamination_audit.sh > "$STAGE12_10_LOG" 2>&1 && \
       grep 'STAGE 12.10 FINAL VERDICT: PASS' "$STAGE12_10_LOG" >/dev/null 2>&1; then
        stage12_10_status=1
    else
        add_reason "stage12_10_total_smoke_failed"
    fi

    if verify_stage12_diag; then
        final_diag_status=1
    fi
    if verify_stage12_10_dat; then
        stage12_10_dat_status=1
    fi
fi

if [ "$build_status" -eq 1 ] && \
   [ "$stage12_6_status" -eq 1 ] && \
   [ "$stage12_7_status" -eq 1 ] && \
   [ "$stage12_8_status" -eq 1 ] && \
   [ "$stage12_9_status" -eq 1 ] && \
   [ "$stage12_10_status" -eq 1 ] && \
   [ "$hook_active_status" -eq 1 ] && \
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
   [ "$final_diag_status" -eq 1 ] && \
   [ "$stage12_10_dat_status" -eq 1 ]; then
    total_smoke_status=1
fi

write_gate

if [ "$total_smoke_status" -eq 1 ]; then
    generate_stage12_closed
    echo "STAGE 12.11 TOTAL SMOKE VERDICT: PASS"
    echo "STAGE 12.11 FINAL VERDICT: PASS"
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage12_11_failure"
fi

echo "STAGE 12.11 TOTAL SMOKE VERDICT: FAIL"
echo "STAGE 12.11 FINAL VERDICT: FAIL"
echo "Reasons: $REASONS"
exit 1
