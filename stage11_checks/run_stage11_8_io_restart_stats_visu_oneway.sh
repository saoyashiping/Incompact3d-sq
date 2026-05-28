#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_8_RUN_STAGE11_7=${STAGE11_8_RUN_STAGE11_7:-0}

mkdir -p stage11_outputs stage9_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_8_RUN_STAGE11_7" = "1" ]; then
    bash stage11_checks/run_stage11_7_oneway_sampling_parallel_consistency.sh
fi

build_status=1
stage9_7_status=1
stage9_8_status=1
stage9_7_hook_active_status=0
stage9_8_hook_active_status=0
sample_performed_status=0
sampled_velocity_finite_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_force_status=0
no_twoway_force_status=0
no_structure_advance_status=0
no_restart_contamination_status=0
no_stats_visu_contamination_status=0
final_status=1
requested_flag=0
reasons="init"

for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check; do
    cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done

rm -f stage11_outputs/fibre_stage11_5_production_oneway_hook.dat

LOG97=stage11_outputs/stage11_8_stage9_7_stats_visu_io_oneway.log
if [ "$build_status" -eq 1 ]; then
    X3D_STAGE11_ONEWAY_HOOK=1 \
    X3D_STAGE11_FORCE_READONLY=1 \
    X3D_STAGE11_MAX_POINTS=8 \
    X3D_STAGE11_MAX_STEPS=3 \
    STAGE9_SKIP_PREREQS=1 \
    bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh > "$LOG97" 2>&1 || stage9_7_status=0
else
    stage9_7_status=0
fi

if [ "$stage9_7_status" -eq 1 ]; then
    grep -q "STAGE 9.7 FINAL VERDICT: PASS" "$LOG97" || stage9_7_status=0
fi

for f in stage9_outputs/*stage9_7* stage9_outputs/*stats* stage9_outputs/*visu* stage9_outputs/*coarse*; do
    [ -e "$f" ] && cp "$f" stage11_outputs/ || true
done
[ -f stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat ] && cp stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat stage11_outputs/stage11_8_stage9_7_stats_visu_io_oneway.dat || true

HOOK=stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
get_val(){ awk -v k="$1" '$1==k{print $2}' "$2"; }

check_hook_diag() {
    phase=$1
    if [ ! -f "$HOOK" ]; then
        reasons="missing_stage11_5_hook_diagnostics_after_${phase}"
        return 1
    fi
    for k in stage11_5_requested_flag stage11_5_readonly_mode_status stage11_5_hook_initialized_status stage11_5_hook_sample_called_status stage11_5_sample_performed_status stage11_5_sample_count_status stage11_5_sampled_velocity_finite_status stage11_5_sampled_velocity_signature_status stage11_5_field_modified_status stage11_5_rhs_modified_status stage11_5_no_rhs_injection_status stage11_5_no_ibm_spreading_status stage11_5_no_feedback_force_status stage11_5_no_twoway_force_status stage11_5_no_structure_advance_status stage11_5_production_oneway_hook_status; do
        v=$(get_val "$k" "$HOOK")
        case "$k" in
            stage11_5_field_modified_status|stage11_5_rhs_modified_status) exp=0 ;;
            *) exp=1 ;;
        esac
        if [ "$v" != "$exp" ]; then
            reasons="${k}_${phase}"
            return 1
        fi
    done
    return 0
}

if check_hook_diag stage9_7; then
    stage9_7_hook_active_status=$(get_val stage11_5_hook_sample_called_status "$HOOK")
    sample_performed_status=$(get_val stage11_5_sample_performed_status "$HOOK")
    sampled_velocity_finite_status=$(get_val stage11_5_sampled_velocity_finite_status "$HOOK")
    [ "$(get_val stage11_5_field_modified_status "$HOOK")" = "0" ] && no_field_modification_status=1 || no_field_modification_status=0
    [ "$(get_val stage11_5_rhs_modified_status "$HOOK")" = "0" ] && no_rhs_modification_status=1 || no_rhs_modification_status=0
    no_rhs_injection_status=$(get_val stage11_5_no_rhs_injection_status "$HOOK")
    no_ibm_spreading_status=$(get_val stage11_5_no_ibm_spreading_status "$HOOK")
    no_feedback_force_status=$(get_val stage11_5_no_feedback_force_status "$HOOK")
    no_twoway_force_status=$(get_val stage11_5_no_twoway_force_status "$HOOK")
    no_structure_advance_status=$(get_val stage11_5_no_structure_advance_status "$HOOK")
    cp "$HOOK" stage11_outputs/stage11_8_stage9_7_hook_diagnostics.dat
else
    stage9_7_hook_active_status=0
    final_status=0
fi

rm -f stage11_outputs/fibre_stage11_5_production_oneway_hook.dat

LOG98=stage11_outputs/stage11_8_stage9_8_restart_io_oneway.log
if [ "$build_status" -eq 1 ]; then
    X3D_STAGE11_ONEWAY_HOOK=1 \
    X3D_STAGE11_FORCE_READONLY=1 \
    X3D_STAGE11_MAX_POINTS=8 \
    X3D_STAGE11_MAX_STEPS=3 \
    STAGE9_SKIP_PREREQS=1 \
    bash stage9_checks/run_stage9_8_restart_io_regression.sh > "$LOG98" 2>&1 || stage9_8_status=0
else
    stage9_8_status=0
fi

if [ "$stage9_8_status" -eq 1 ]; then
    grep -q "STAGE 9.8 FINAL VERDICT: PASS" "$LOG98" || stage9_8_status=0
fi

for f in stage9_outputs/*stage9_8* stage9_outputs/*restart*; do
    [ -e "$f" ] && cp "$f" stage11_outputs/ || true
done
[ -f stage9_outputs/fibre_stage9_8_restart_io_regression.dat ] && cp stage9_outputs/fibre_stage9_8_restart_io_regression.dat stage11_outputs/stage11_8_stage9_8_restart_io_oneway.dat || true

if check_hook_diag stage9_8; then
    stage9_8_hook_active_status=$(get_val stage11_5_hook_sample_called_status "$HOOK")
    cp "$HOOK" stage11_outputs/stage11_8_stage9_8_hook_diagnostics.dat
else
    stage9_8_hook_active_status=0
    final_status=0
fi

[ "$stage9_8_status" -eq 1 ] && no_restart_contamination_status=1 || no_restart_contamination_status=0
[ "$stage9_7_status" -eq 1 ] && no_stats_visu_contamination_status=1 || no_stats_visu_contamination_status=0

if [ "$stage9_7_hook_active_status" = "1" ] || [ "$stage9_8_hook_active_status" = "1" ]; then
    requested_flag=1
fi

if [ "$final_status" -eq 1 ]; then
  if [ "$build_status" -ne 1 ] || [ "$stage9_7_status" -ne 1 ] || [ "$stage9_8_status" -ne 1 ] || [ "$stage9_7_hook_active_status" -ne 1 ] || [ "$stage9_8_hook_active_status" -ne 1 ] || [ "$sample_performed_status" -ne 1 ] || [ "$sampled_velocity_finite_status" -ne 1 ] || [ "$no_field_modification_status" -ne 1 ] || [ "$no_rhs_modification_status" -ne 1 ] || [ "$no_rhs_injection_status" -ne 1 ] || [ "$no_ibm_spreading_status" -ne 1 ] || [ "$no_feedback_force_status" -ne 1 ] || [ "$no_twoway_force_status" -ne 1 ] || [ "$no_structure_advance_status" -ne 1 ] || [ "$no_restart_contamination_status" -ne 1 ] || [ "$no_stats_visu_contamination_status" -ne 1 ]; then
    final_status=0
    [ "$reasons" = "init" ] && reasons="io_restart_or_hook_checks_failed"
  fi
fi

cat > stage11_outputs/stage11_8_io_restart_stats_visu_oneway.dat <<EOD
stage11_8_requested_flag $requested_flag
stage11_8_build_status $build_status
stage11_8_stage9_7_stats_visu_status $stage9_7_status
stage11_8_stage9_8_restart_status $stage9_8_status
stage11_8_stage9_7_hook_active_status $stage9_7_hook_active_status
stage11_8_stage9_8_hook_active_status $stage9_8_hook_active_status
stage11_8_sample_performed_status $sample_performed_status
stage11_8_sampled_velocity_finite_status $sampled_velocity_finite_status
stage11_8_no_field_modification_status $no_field_modification_status
stage11_8_no_rhs_modification_status $no_rhs_modification_status
stage11_8_no_rhs_injection_status $no_rhs_injection_status
stage11_8_no_ibm_spreading_status $no_ibm_spreading_status
stage11_8_no_feedback_force_status $no_feedback_force_status
stage11_8_no_twoway_force_status $no_twoway_force_status
stage11_8_no_structure_advance_status $no_structure_advance_status
stage11_8_no_restart_contamination_status $no_restart_contamination_status
stage11_8_no_stats_visu_contamination_status $no_stats_visu_contamination_status
stage11_8_io_restart_stats_visu_oneway_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 11.8 IO RESTART STATS VISU ONEWAY VERDICT: PASS"
  echo "STAGE 11.8 FINAL VERDICT: PASS"
else
  [ "$reasons" = "init" ] && reasons="unknown_failure"
  echo "STAGE 11.8 IO RESTART STATS VISU ONEWAY VERDICT: FAIL"
  echo "STAGE 11.8 FINAL VERDICT: FAIL"
  echo "Reasons:$reasons"
  exit 1
fi
