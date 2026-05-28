#!/bin/sh
set -eu
BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_4_RUN_STAGE11_3=${STAGE11_4_RUN_STAGE11_3:-0}
mkdir -p stage11_outputs
ensure_build_dir(){ if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then cmake -S . -B "$BUILD_DIR"; fi; }
ensure_build_dir
[ "$STAGE11_4_RUN_STAGE11_3" = "1" ] && sh stage11_checks/run_stage11_3_oneway_interpolation_api.sh
build_status=1; check_status=1; gate_status=1; reasons=""
for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check; do
  cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done
LOG=stage11_outputs/stage11_4_controlled_interpolation.log
EXE="$BUILD_DIR/bin/fibre_stage11_controlled_interpolation_check"; [ -x "$EXE" ] || EXE="$BUILD_DIR/src/fibre_stage11_controlled_interpolation_check"
if [ "$build_status" -eq 1 ]; then X3D_STAGE11_ONEWAY_HOOK=1 X3D_STAGE11_FORCE_READONLY=1 X3D_STAGE11_MAX_POINTS=8 X3D_STAGE11_MAX_STEPS=3 "$EXE" > "$LOG" 2>&1 || check_status=0; else check_status=0; fi
[ "$check_status" -eq 1 ] && grep -q "STAGE 11.4 CONTROLLED INTERPOLATION VERDICT: PASS" "$LOG" || check_status=0
DAT=stage11_outputs/fibre_stage11_4_controlled_interpolation.dat
getv(){ awk -v k="$1" '$1==k{print $2}' "$DAT"; }
if [ "$check_status" -eq 1 ] && [ -f "$DAT" ]; then
  requested=$(getv stage11_4_requested_flag); readonly=$(getv stage11_4_readonly_mode_status); cstat=$(getv stage11_4_controlled_interpolation_status)
  conste=$(getv stage11_4_constant_max_error); lineare=$(getv stage11_4_linear_max_error); sheare=$(getv stage11_4_shear_max_error); poise=$(getv stage11_4_poiseuille_max_error); wsum=$(getv stage11_4_weight_sum_max_error)
  awk -v a="$conste" -v b="$lineare" -v c="$sheare" -v d="$wsum" -v e="$poise" 'BEGIN{ok=(a<=1e-12 && b<=1e-12 && c<=1e-12 && d<=1e-12 && e<=5e-3); exit(ok?0:1)}' || gate_status=0
  [ "$requested" = "1" ] && [ "$readonly" = "1" ] && [ "$cstat" = "1" ] || gate_status=0
else gate_status=0; fi
[ "$build_status" -eq 1 ] || { reasons="$reasons build_failed;"; gate_status=0; }
[ "$check_status" -eq 1 ] || { reasons="$reasons controlled_interpolation_check_failed;"; gate_status=0; }
cat > stage11_outputs/stage11_4_controlled_interpolation_gate.dat <<EOD
stage11_4_gate_requested_flag ${requested:-0}
stage11_4_gate_build_status $build_status
stage11_4_gate_controlled_interpolation_check_status $check_status
stage11_4_gate_constant_field_status $(getv stage11_4_constant_field_status 2>/dev/null || echo 0)
stage11_4_gate_linear_field_status $(getv stage11_4_linear_field_status 2>/dev/null || echo 0)
stage11_4_gate_shear_field_status $(getv stage11_4_shear_field_status 2>/dev/null || echo 0)
stage11_4_gate_poiseuille_field_status $(getv stage11_4_poiseuille_field_status 2>/dev/null || echo 0)
stage11_4_gate_weight_sum_status $(getv stage11_4_weight_sum_status 2>/dev/null || echo 0)
stage11_4_gate_periodic_safety_status $(getv stage11_4_periodic_safety_status 2>/dev/null || echo 0)
stage11_4_gate_near_wall_safety_status $(getv stage11_4_near_wall_safety_status 2>/dev/null || echo 0)
stage11_4_gate_out_of_domain_rejection_status $(getv stage11_4_out_of_domain_rejection_status 2>/dev/null || echo 0)
stage11_4_gate_no_production_fluid_access_status $(getv stage11_4_no_production_fluid_access_status 2>/dev/null || echo 0)
stage11_4_gate_no_fluid_field_modification_status $(getv stage11_4_no_fluid_field_modification_status 2>/dev/null || echo 0)
stage11_4_gate_no_rhs_injection_status $(getv stage11_4_no_rhs_injection_status 2>/dev/null || echo 0)
stage11_4_gate_no_ibm_spreading_status $(getv stage11_4_no_ibm_spreading_status 2>/dev/null || echo 0)
stage11_4_gate_no_feedback_force_status $(getv stage11_4_no_feedback_force_status 2>/dev/null || echo 0)
stage11_4_gate_no_twoway_force_status $(getv stage11_4_no_twoway_force_status 2>/dev/null || echo 0)
stage11_4_gate_no_structure_advance_status $(getv stage11_4_no_structure_advance_status 2>/dev/null || echo 0)
stage11_4_gate_status $gate_status
EOD
[ "$gate_status" -eq 1 ] && echo "STAGE 11.4 FINAL VERDICT: PASS" || { echo "STAGE 11.4 FINAL VERDICT: FAIL"; echo "Reasons:$reasons"; exit 1; }
