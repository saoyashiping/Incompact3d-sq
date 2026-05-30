#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}

mkdir -p stage11_outputs stage9_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

build_status=1
stage11_5_status=1
stage11_6_status=1
stage11_7_status=1
stage11_8_status=1
stage11_9_status=1
hook_active_status=0
sample_performed_status=0
sampled_velocity_finite_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_force_status=0
no_twoway_force_status=0
no_structure_advance_status=0
requested_flag=0
final_status=1
reasons="init"

for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check; do
    cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done

run_gate() {
    gate_name=$1
    log_file=$2
    shift 2
    if "$@" > "$log_file" 2>&1; then
        grep -q "STAGE ${gate_name} FINAL VERDICT: PASS" "$log_file" || return 1
        return 0
    fi
    return 1
}

if [ "$build_status" -eq 1 ]; then
    run_gate "11.5" stage11_outputs/stage11_10_stage11_5.log env STAGE11_5_RUN_STAGE11_4=0 bash stage11_checks/run_stage11_5_production_oneway_hook.sh || stage11_5_status=0
    run_gate "11.6" stage11_outputs/stage11_10_stage11_6.log env STAGE11_6_RUN_STAGE11_5=0 bash stage11_checks/run_stage11_6_oneway_sampling_invariance_np1.sh || stage11_6_status=0
    run_gate "11.7" stage11_outputs/stage11_10_stage11_7.log env STAGE11_7_RUN_STAGE11_6=0 bash stage11_checks/run_stage11_7_oneway_sampling_parallel_consistency.sh || stage11_7_status=0
    run_gate "11.8" stage11_outputs/stage11_10_stage11_8.log env STAGE11_8_RUN_STAGE11_7=0 bash stage11_checks/run_stage11_8_io_restart_stats_visu_oneway.sh || stage11_8_status=0
    run_gate "11.9" stage11_outputs/stage11_10_stage11_9.log env STAGE11_9_RUN_STAGE11_8=0 bash stage11_checks/run_stage11_9_rhs_coupling_contamination_audit.sh || stage11_9_status=0
else
    stage11_5_status=0
    stage11_6_status=0
    stage11_7_status=0
    stage11_8_status=0
    stage11_9_status=0
fi

HOOK=stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
AUDIT=stage11_outputs/stage11_9_rhs_coupling_contamination_audit.dat
get_val(){ awk -v k="$1" '$1==k{print $2}' "$2"; }

if [ -f "$HOOK" ]; then
    requested_flag=$(get_val stage11_5_requested_flag "$HOOK")
    hook_active_status=$(get_val stage11_5_hook_sample_called_status "$HOOK")
    sample_performed_status=$(get_val stage11_5_sample_performed_status "$HOOK")
    sampled_velocity_finite_status=$(get_val stage11_5_sampled_velocity_finite_status "$HOOK")
    [ "$(get_val stage11_5_field_modified_status "$HOOK")" = "0" ] && no_field_modification_status=1 || no_field_modification_status=0
    [ "$(get_val stage11_5_rhs_modified_status "$HOOK")" = "0" ] && no_rhs_modification_status=1 || no_rhs_modification_status=0
    no_rhs_injection_status=$(get_val stage11_5_no_rhs_injection_status "$HOOK")
    no_ibm_spreading_status=$(get_val stage11_5_no_ibm_spreading_status "$HOOK")
    no_feedback_force_status=$(get_val stage11_5_no_feedback_force_status "$HOOK")
    no_twoway_force_status=$(get_val stage11_5_no_twoway_force_status "$HOOK")
    no_structure_advance_status=$(get_val stage11_5_no_structure_advance_status "$HOOK")
else
    final_status=0
    reasons="missing_stage11_5_hook_diagnostics"
fi

if [ -f "$AUDIT" ]; then
    for k in \
      stage11_9_static_audit_status \
      stage11_9_velocity_intent_readonly_status \
      stage11_9_no_velocity_write_static_status \
      stage11_9_no_rhs_write_static_status \
      stage11_9_no_rhs_injection_static_status \
      stage11_9_no_ibm_spreading_static_status \
      stage11_9_no_feedback_force_static_status \
      stage11_9_no_twoway_force_static_status \
      stage11_9_no_structure_advance_static_status \
      stage11_9_runtime_smoke_status \
      stage11_9_hook_active_status \
      stage11_9_sample_performed_status \
      stage11_9_sampled_velocity_finite_status \
      stage11_9_no_field_modification_status \
      stage11_9_no_rhs_modification_status \
      stage11_9_no_rhs_injection_status \
      stage11_9_no_ibm_spreading_status \
      stage11_9_no_feedback_force_status \
      stage11_9_no_twoway_force_status \
      stage11_9_no_structure_advance_status \
      stage11_9_rhs_coupling_contamination_audit_status; do
      [ "$(get_val "$k" "$AUDIT")" = "1" ] || { final_status=0; [ "$reasons" = "init" ] && reasons="$k"; }
    done
else
    final_status=0
    [ "$reasons" = "init" ] && reasons="missing_stage11_9_audit_dat"
fi

if [ "$final_status" -eq 1 ]; then
  if [ "$build_status" -ne 1 ] || [ "$stage11_5_status" -ne 1 ] || [ "$stage11_6_status" -ne 1 ] || [ "$stage11_7_status" -ne 1 ] || [ "$stage11_8_status" -ne 1 ] || [ "$stage11_9_status" -ne 1 ] || [ "$hook_active_status" -ne 1 ] || [ "$sample_performed_status" -ne 1 ] || [ "$sampled_velocity_finite_status" -ne 1 ] || [ "$no_field_modification_status" -ne 1 ] || [ "$no_rhs_modification_status" -ne 1 ] || [ "$no_rhs_injection_status" -ne 1 ] || [ "$no_ibm_spreading_status" -ne 1 ] || [ "$no_feedback_force_status" -ne 1 ] || [ "$no_twoway_force_status" -ne 1 ] || [ "$no_structure_advance_status" -ne 1 ]; then
    final_status=0
    [ "$reasons" = "init" ] && reasons="stage11_10_total_smoke_checks_failed"
  fi
fi

cat > stage11_outputs/stage11_10_total_smoke.dat <<EOD
stage11_10_requested_flag $requested_flag
stage11_10_build_status $build_status
stage11_10_stage11_5_status $stage11_5_status
stage11_10_stage11_6_status $stage11_6_status
stage11_10_stage11_7_status $stage11_7_status
stage11_10_stage11_8_status $stage11_8_status
stage11_10_stage11_9_status $stage11_9_status
stage11_10_hook_active_status $hook_active_status
stage11_10_sample_performed_status $sample_performed_status
stage11_10_sampled_velocity_finite_status $sampled_velocity_finite_status
stage11_10_no_field_modification_status $no_field_modification_status
stage11_10_no_rhs_modification_status $no_rhs_modification_status
stage11_10_no_rhs_injection_status $no_rhs_injection_status
stage11_10_no_ibm_spreading_status $no_ibm_spreading_status
stage11_10_no_feedback_force_status $no_feedback_force_status
stage11_10_no_twoway_force_status $no_twoway_force_status
stage11_10_no_structure_advance_status $no_structure_advance_status
stage11_10_total_smoke_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
    cat > stage11_checks/STAGE11_CLOSED.md <<'EOC'
# STAGE11_CLOSED

## Stage 11 purpose
Production one-way fluid-to-fibre sampling path.

## Closed sub-stages
- Stage 11.0 config and global switches.
- Stage 11.1 Lagrangian/control-point state skeleton.
- Stage 11.2 grid metadata bridge skeleton.
- Stage 11.3 interpolation API compile-only skeleton.
- Stage 11.4 controlled analytic interpolation checks.
- Stage 11.5 production one-way sampling hook.
- Stage 11.6 np=1 one-way sampling invariance.
- Stage 11.7 np=1/2/4 parallel consistency.
- Stage 11.8 restart / stats / visu / coarse I/O compatibility.
- Stage 11.9 RHS / coupling contamination audit.
- Stage 11.10 total smoke closure.

## Governing read-only model
- div(u)=0
- du/dt + u·grad(u) = -grad(p) + nu laplacian(u) + f_drive + f_fsi
- f_fsi=0
- RHS_stage11=RHS_stage10=RHS_stage9
- U_f = I_h[u](X_f)

## Explicit statement
Stage 11 closes the production one-way fluid-to-fibre sampled-velocity data path only.
No real feedback force, no IBM spreading, no RHS injection, no two-way force density, and no fibre structure advance are activated.
Real feedback force candidate begins in Stage 12.

## Next recommended stage
Stage 12 production feedback force candidate.
EOC
    echo "STAGE 11.10 TOTAL SMOKE VERDICT: PASS"
    echo "STAGE 11.10 FINAL VERDICT: PASS"
else
    [ "$reasons" = "init" ] && reasons="unknown_failure"
    echo "STAGE 11.10 TOTAL SMOKE VERDICT: FAIL"
    echo "STAGE 11.10 FINAL VERDICT: FAIL"
    echo "Reasons:$reasons"
    exit 1
fi
