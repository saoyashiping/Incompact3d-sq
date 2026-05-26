#!/bin/sh
set -eu
BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_5_RUN_STAGE11_4=${STAGE11_5_RUN_STAGE11_4:-0}
mkdir -p stage11_outputs stage9_outputs
ensure_build_dir(){ if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then cmake -S . -B "$BUILD_DIR"; fi; }
ensure_build_dir
[ "$STAGE11_5_RUN_STAGE11_4" = "1" ] && sh stage11_checks/run_stage11_4_controlled_interpolation.sh
build_status=1; standalone_check_status=1; production_smoke_status=1; hook_diagnostics_status=1; gate_status=1; reasons=""
for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check; do
  cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done
LOG=stage11_outputs/stage11_5_production_oneway_hook_check.log
EXE="$BUILD_DIR/bin/fibre_stage11_production_oneway_hook_check"; [ -x "$EXE" ] || EXE="$BUILD_DIR/src/fibre_stage11_production_oneway_hook_check"
if [ "$build_status" -eq 1 ]; then
  X3D_STAGE11_ONEWAY_HOOK=1 X3D_STAGE11_FORCE_READONLY=1 X3D_STAGE11_MAX_POINTS=8 X3D_STAGE11_MAX_STEPS=3 "$EXE" > "$LOG" 2>&1 || standalone_check_status=0
else standalone_check_status=0; fi
[ "$standalone_check_status" -eq 1 ] && grep -q "STAGE 11.5 PRODUCTION ONEWAY HOOK CHECK VERDICT: PASS" "$LOG" || standalone_check_status=0
if [ "$build_status" -eq 1 ]; then
  X3D_STAGE11_ONEWAY_HOOK=1 X3D_STAGE11_FORCE_READONLY=1 X3D_STAGE11_MAX_POINTS=8 X3D_STAGE11_MAX_STEPS=3 STAGE9_SKIP_PREREQS=1 X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 X3D_STAGE9_9_DETERMINISTIC_INIT=1 X3D_STAGE9_9_MAX_STEPS=3 bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh || production_smoke_status=0
else production_smoke_status=0; fi
DAT=stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
getv(){ awk -v k="$1" '$1==k{print $2}' "$DAT"; }
hook_active_status=0; sample_performed_status=0; sampled_velocity_finite_status=0; field_modified_status=1; rhs_modified_status=1; no_rhs_injection_status=0; no_ibm_spreading_status=0; no_feedback_force_status=0; no_twoway_force_status=0; no_structure_advance_status=0; requested_flag=0
if [ -f "$DAT" ]; then
 requested_flag=$(getv stage11_5_requested_flag); hook_active_status=$(getv stage11_5_hook_sample_called_status); sample_performed_status=$(getv stage11_5_sample_performed_status); sampled_velocity_finite_status=$(getv stage11_5_sampled_velocity_finite_status); field_modified_status=$(getv stage11_5_field_modified_status); rhs_modified_status=$(getv stage11_5_rhs_modified_status); no_rhs_injection_status=$(getv stage11_5_no_rhs_injection_status); no_ibm_spreading_status=$(getv stage11_5_no_ibm_spreading_status); no_feedback_force_status=$(getv stage11_5_no_feedback_force_status); no_twoway_force_status=$(getv stage11_5_no_twoway_force_status); no_structure_advance_status=$(getv stage11_5_no_structure_advance_status)
fi
if [ ! -f "$DAT" ]; then
  hook_diagnostics_status=0
  reasons="$reasons missing_production_hook_dat;"
elif [ "$requested_flag" != "1" ] || [ "$hook_active_status" != "1" ] || [ "$sample_performed_status" != "1" ] || [ "$sampled_velocity_finite_status" != "1" ] || [ "$field_modified_status" != "0" ] || [ "$rhs_modified_status" != "0" ] || [ "$no_rhs_injection_status" != "1" ] || [ "$no_ibm_spreading_status" != "1" ] || [ "$no_feedback_force_status" != "1" ] || [ "$no_twoway_force_status" != "1" ] || [ "$no_structure_advance_status" != "1" ]; then
  hook_diagnostics_status=0
  reasons="$reasons hook_diagnostics_failed;"
fi
[ "$build_status" -eq 1 ] || { reasons="$reasons build_failed;"; gate_status=0; }
[ "$standalone_check_status" -eq 1 ] || { reasons="$reasons standalone_check_failed;"; gate_status=0; }
[ "$production_smoke_status" -eq 1 ] || { reasons="$reasons production_smoke_failed;"; gate_status=0; }
[ "$hook_diagnostics_status" -eq 1 ] || gate_status=0
cat > stage11_outputs/stage11_5_production_oneway_hook_gate.dat <<EOD
stage11_5_gate_requested_flag $requested_flag
stage11_5_gate_build_status $build_status
stage11_5_gate_standalone_check_status $standalone_check_status
stage11_5_gate_production_smoke_status $production_smoke_status
stage11_5_gate_hook_active_status $hook_active_status
stage11_5_gate_sample_performed_status $sample_performed_status
stage11_5_gate_sampled_velocity_finite_status $sampled_velocity_finite_status
stage11_5_gate_field_modified_status $field_modified_status
stage11_5_gate_rhs_modified_status $rhs_modified_status
stage11_5_gate_no_rhs_injection_status $no_rhs_injection_status
stage11_5_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage11_5_gate_no_feedback_force_status $no_feedback_force_status
stage11_5_gate_no_twoway_force_status $no_twoway_force_status
stage11_5_gate_no_structure_advance_status $no_structure_advance_status
stage11_5_gate_status $gate_status
EOD
if [ "$gate_status" -eq 1 ]; then
 echo "STAGE 11.5 PRODUCTION ONEWAY HOOK VERDICT: PASS"
 echo "STAGE 11.5 FINAL VERDICT: PASS"
else
 echo "STAGE 11.5 FINAL VERDICT: FAIL"
 echo "Reasons:$reasons"
 exit 1
fi
