#!/bin/sh
set -eu
BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_3_RUN_STAGE11_2=${STAGE11_3_RUN_STAGE11_2:-0}
mkdir -p stage11_outputs
ensure_build_dir(){ if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then cmake -S . -B "$BUILD_DIR"; fi; }
ensure_build_dir
if [ "$STAGE11_3_RUN_STAGE11_2" = "1" ]; then sh stage11_checks/run_stage11_2_grid_metadata.sh; fi
build_status=1
api_check_status=1
gate_status=1
reasons=""
for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check; do
  if ! cmake --build "$BUILD_DIR" --target "$tgt" -j; then build_status=0; fi
done
LOG_FILE=stage11_outputs/stage11_3_oneway_interpolation_api.log
EXE_PATH="$BUILD_DIR/bin/fibre_stage11_oneway_interpolation_check"
if [ ! -x "$EXE_PATH" ]; then EXE_PATH="$BUILD_DIR/src/fibre_stage11_oneway_interpolation_check"; fi
if [ "$build_status" -eq 1 ]; then
  if ! X3D_STAGE11_ONEWAY_HOOK=1 X3D_STAGE11_FORCE_READONLY=1 X3D_STAGE11_MAX_POINTS=8 X3D_STAGE11_MAX_STEPS=3 "$EXE_PATH" > "$LOG_FILE" 2>&1; then api_check_status=0; fi
else
  api_check_status=0
fi
if [ "$api_check_status" -eq 1 ]; then if ! grep -q "STAGE 11.3 ONEWAY INTERPOLATION API VERDICT: PASS" "$LOG_FILE"; then api_check_status=0; fi; fi
DAT=stage11_outputs/fibre_stage11_3_oneway_interpolation_api.dat
requested_flag=0; lagrangian_state_available_status=0; grid_metadata_available_status=0; interpolation_initialized_status=0; interface_available_status=0; prepare_called_status=0; sample_interface_called_status=0; sample_not_performed_status=0; velocity_placeholder_unchanged_status=0; no_fluid_field_access_status=0; no_velocity_sampling_status=0; no_fluid_field_modification_status=0; no_rhs_injection_status=0; no_ibm_spreading_status=0; no_feedback_force_status=0; no_twoway_force_status=0; no_structure_advance_status=0; oneway_interpolation_api_status=0
if [ "$api_check_status" -eq 1 ] && [ -f "$DAT" ]; then
 get_val(){ awk -v k="$1" '$1==k{print $2}' "$DAT"; }
 requested_flag=$(get_val stage11_3_requested_flag)
 readonly_mode_status=$(get_val stage11_3_readonly_mode_status)
 lagrangian_state_available_status=$(get_val stage11_3_lagrangian_state_available_status)
 grid_metadata_available_status=$(get_val stage11_3_grid_metadata_available_status)
 interpolation_initialized_status=$(get_val stage11_3_interpolation_initialized_status)
 interface_available_status=$(get_val stage11_3_interface_available_status)
 prepare_called_status=$(get_val stage11_3_prepare_called_status)
 sample_interface_called_status=$(get_val stage11_3_sample_interface_called_status)
 sample_not_performed_status=$(get_val stage11_3_sample_not_performed_status)
 lagrangian_state_input_status=$(get_val stage11_3_lagrangian_state_input_status)
 grid_metadata_input_status=$(get_val stage11_3_grid_metadata_input_status)
 velocity_placeholder_unchanged_status=$(get_val stage11_3_velocity_placeholder_unchanged_status)
 no_fluid_field_access_status=$(get_val stage11_3_no_fluid_field_access_status)
 no_velocity_sampling_status=$(get_val stage11_3_no_velocity_sampling_status)
 no_fluid_field_modification_status=$(get_val stage11_3_no_fluid_field_modification_status)
 no_rhs_injection_status=$(get_val stage11_3_no_rhs_injection_status)
 no_ibm_spreading_status=$(get_val stage11_3_no_ibm_spreading_status)
 no_feedback_force_status=$(get_val stage11_3_no_feedback_force_status)
 no_twoway_force_status=$(get_val stage11_3_no_twoway_force_status)
 no_structure_advance_status=$(get_val stage11_3_no_structure_advance_status)
 oneway_interpolation_api_status=$(get_val stage11_3_oneway_interpolation_api_status)
 if [ "$requested_flag" != "1" ] || [ "$readonly_mode_status" != "1" ] || [ "$lagrangian_state_available_status" != "1" ] || [ "$grid_metadata_available_status" != "1" ] || [ "$interpolation_initialized_status" != "1" ] || [ "$interface_available_status" != "1" ] || [ "$prepare_called_status" != "1" ] || [ "$sample_interface_called_status" != "1" ] || [ "$sample_not_performed_status" != "1" ] || [ "$lagrangian_state_input_status" != "1" ] || [ "$grid_metadata_input_status" != "1" ] || [ "$velocity_placeholder_unchanged_status" != "1" ] || [ "$no_fluid_field_access_status" != "1" ] || [ "$no_velocity_sampling_status" != "1" ] || [ "$no_fluid_field_modification_status" != "1" ] || [ "$no_rhs_injection_status" != "1" ] || [ "$no_ibm_spreading_status" != "1" ] || [ "$no_feedback_force_status" != "1" ] || [ "$no_twoway_force_status" != "1" ] || [ "$no_structure_advance_status" != "1" ] || [ "$oneway_interpolation_api_status" != "1" ]; then gate_status=0; fi
else gate_status=0; fi
if [ "$build_status" -ne 1 ]; then reasons="$reasons build_failed;"; gate_status=0; fi
if [ "$api_check_status" -ne 1 ]; then reasons="$reasons interpolation_api_check_failed;"; gate_status=0; fi
cat > stage11_outputs/stage11_3_oneway_interpolation_api_gate.dat <<EOD
stage11_3_gate_requested_flag $requested_flag
stage11_3_gate_build_status $build_status
stage11_3_gate_interpolation_api_check_status $api_check_status
stage11_3_gate_lagrangian_state_available_status $lagrangian_state_available_status
stage11_3_gate_grid_metadata_available_status $grid_metadata_available_status
stage11_3_gate_interpolation_initialized_status $interpolation_initialized_status
stage11_3_gate_interface_available_status $interface_available_status
stage11_3_gate_prepare_called_status $prepare_called_status
stage11_3_gate_sample_interface_called_status $sample_interface_called_status
stage11_3_gate_sample_not_performed_status $sample_not_performed_status
stage11_3_gate_velocity_placeholder_unchanged_status $velocity_placeholder_unchanged_status
stage11_3_gate_no_fluid_field_access_status $no_fluid_field_access_status
stage11_3_gate_no_velocity_sampling_status $no_velocity_sampling_status
stage11_3_gate_no_fluid_field_modification_status $no_fluid_field_modification_status
stage11_3_gate_no_rhs_injection_status $no_rhs_injection_status
stage11_3_gate_no_ibm_spreading_status $no_ibm_spreading_status
stage11_3_gate_no_feedback_force_status $no_feedback_force_status
stage11_3_gate_no_twoway_force_status $no_twoway_force_status
stage11_3_gate_no_structure_advance_status $no_structure_advance_status
stage11_3_gate_status $gate_status
EOD
if [ "$gate_status" -eq 1 ]; then echo "STAGE 11.3 FINAL VERDICT: PASS"; else echo "STAGE 11.3 FINAL VERDICT: FAIL"; echo "Reasons:$reasons"; exit 1; fi
