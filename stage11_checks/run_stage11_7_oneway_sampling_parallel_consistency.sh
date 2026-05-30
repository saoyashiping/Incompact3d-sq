#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_7_RUN_STAGE11_6=${STAGE11_7_RUN_STAGE11_6:-0}
STAGE11_7_INVARIANCE_ABS_TOL=${STAGE11_7_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE11_7_INVARIANCE_REL_TOL=${STAGE11_7_INVARIANCE_REL_TOL:-1.0e-14}

mkdir -p stage11_outputs stage9_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_7_RUN_STAGE11_6" = "1" ]; then
    bash stage11_checks/run_stage11_6_oneway_sampling_invariance_np1.sh
fi

build_status=1
baseline_run_status=1
hook_run_status=1
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
np1_signature_invariance_status=1
np2_signature_invariance_status=1
np4_signature_invariance_status=1
parallel_consistency_status=0
final_status=1
reasons="init"

for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check; do
    cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done

BASE_LOG=stage11_outputs/stage11_7_baseline.log
HOOK_LOG=stage11_outputs/stage11_7_hook.log

if [ "$build_status" -eq 1 ]; then
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$BASE_LOG" 2>&1 || baseline_run_status=0
else
    baseline_run_status=0
fi

if [ "$baseline_run_status" -eq 1 ]; then
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat stage11_outputs/stage11_7_baseline_np1.dat || baseline_run_status=0
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np2.dat stage11_outputs/stage11_7_baseline_np2.dat || baseline_run_status=0
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np4.dat stage11_outputs/stage11_7_baseline_np4.dat || baseline_run_status=0
fi

rm -f stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
if [ "$build_status" -eq 1 ]; then
    X3D_STAGE11_ONEWAY_HOOK=1 \
    X3D_STAGE11_FORCE_READONLY=1 \
    X3D_STAGE11_MAX_POINTS=8 \
    X3D_STAGE11_MAX_STEPS=3 \
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$HOOK_LOG" 2>&1 || hook_run_status=0
else
    hook_run_status=0
fi

if [ "$hook_run_status" -eq 1 ]; then
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat stage11_outputs/stage11_7_hook_np1.dat || hook_run_status=0
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np2.dat stage11_outputs/stage11_7_hook_np2.dat || hook_run_status=0
  cp stage9_outputs/fibre_stage9_9_parallel_consistency_np4.dat stage11_outputs/stage11_7_hook_np4.dat || hook_run_status=0
fi

if [ "$hook_run_status" -eq 1 ]; then
  grep -q "STAGE 9.9 FINAL VERDICT: PASS" "$HOOK_LOG" || parallel_consistency_status=0
  if grep -q "STAGE 9.9 FINAL VERDICT: PASS" "$HOOK_LOG"; then
    parallel_consistency_status=1
  fi
fi

H11=stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
get_val() { awk -v k="$1" '$1==k{print $2}' "$2"; }

if [ ! -f "$H11" ]; then
  reasons="missing_stage11_5_hook_diagnostics"
  final_status=0
else
  for k in stage11_5_requested_flag stage11_5_readonly_mode_status stage11_5_hook_initialized_status stage11_5_hook_sample_called_status stage11_5_sample_performed_status stage11_5_sample_count_status stage11_5_sampled_velocity_finite_status stage11_5_sampled_velocity_signature_status stage11_5_field_modified_status stage11_5_rhs_modified_status stage11_5_no_rhs_injection_status stage11_5_no_ibm_spreading_status stage11_5_no_feedback_force_status stage11_5_no_twoway_force_status stage11_5_no_structure_advance_status stage11_5_production_oneway_hook_status; do
    v=$(get_val "$k" "$H11")
    case "$k" in
      stage11_5_field_modified_status|stage11_5_rhs_modified_status) exp=0 ;;
      *) exp=1 ;;
    esac
    if [ "$v" != "$exp" ]; then
      reasons="$k"
      final_status=0
      break
    fi
  done
  hook_active_status=$(get_val stage11_5_hook_sample_called_status "$H11")
  sample_performed_status=$(get_val stage11_5_sample_performed_status "$H11")
  sampled_velocity_finite_status=$(get_val stage11_5_sampled_velocity_finite_status "$H11")
  [ "$(get_val stage11_5_field_modified_status "$H11")" = "0" ] && no_field_modification_status=1
  [ "$(get_val stage11_5_rhs_modified_status "$H11")" = "0" ] && no_rhs_modification_status=1
  no_rhs_injection_status=$(get_val stage11_5_no_rhs_injection_status "$H11")
  no_ibm_spreading_status=$(get_val stage11_5_no_ibm_spreading_status "$H11")
  no_feedback_force_status=$(get_val stage11_5_no_feedback_force_status "$H11")
  no_twoway_force_status=$(get_val stage11_5_no_twoway_force_status "$H11")
  no_structure_advance_status=$(get_val stage11_5_no_structure_advance_status "$H11")
fi

compare_np_metric() {
  np=$1; metric=$2
  base="stage11_outputs/stage11_7_baseline_np${np}.dat"
  hook="stage11_outputs/stage11_7_hook_np${np}.dat"
  ref=$(get_val "$metric" "$base")
  cur=$(get_val "$metric" "$hook")
  res=$(awk -v r="$ref" -v c="$cur" -v at="$STAGE11_7_INVARIANCE_ABS_TOL" -v rt="$STAGE11_7_INVARIANCE_REL_TOL" 'BEGIN{d=c-r; ad=(d<0?-d:d); ar=(r<0?-r:r); eff=at; tmp=rt*((ar>1)?ar:1); if(tmp>eff)eff=tmp; pass=(ad<=eff)?1:0; printf "%.17e %.17e %.17e %.17e %d", d,at,rt,eff,pass}')
  d=$(echo "$res"|awk '{print $1}'); at=$(echo "$res"|awk '{print $2}'); rt=$(echo "$res"|awk '{print $3}'); et=$(echo "$res"|awk '{print $4}'); pass=$(echo "$res"|awk '{print $5}')
  echo "np=$np metric=$metric baseline=$ref hook=$cur delta=$d abs_tol=$at rel_tol=$rt effective_tolerance=$et result=$( [ "$pass" = "1" ] && echo PASS || echo FAIL )"
  [ "$pass" = "1" ]
}

run_np_compare() {
  np=$1
  ok=1
  for m in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
    compare_np_metric "$np" "$m" || ok=0
  done
  [ "$ok" = "1" ]
}

if [ "$baseline_run_status" -eq 1 ] && [ "$hook_run_status" -eq 1 ]; then
  run_np_compare 1 || np1_signature_invariance_status=0
  run_np_compare 2 || np2_signature_invariance_status=0
  run_np_compare 4 || np4_signature_invariance_status=0
else
  np1_signature_invariance_status=0
  np2_signature_invariance_status=0
  np4_signature_invariance_status=0
fi

requested_flag=0
[ "$hook_active_status" = "1" ] && requested_flag=1

if [ "$final_status" -eq 1 ]; then
  if [ "$build_status" -ne 1 ] || [ "$baseline_run_status" -ne 1 ] || [ "$hook_run_status" -ne 1 ] || [ "$hook_active_status" -ne 1 ] || [ "$sample_performed_status" -ne 1 ] || [ "$sampled_velocity_finite_status" -ne 1 ] || [ "$no_field_modification_status" -ne 1 ] || [ "$no_rhs_modification_status" -ne 1 ] || [ "$no_rhs_injection_status" -ne 1 ] || [ "$no_ibm_spreading_status" -ne 1 ] || [ "$no_feedback_force_status" -ne 1 ] || [ "$no_twoway_force_status" -ne 1 ] || [ "$no_structure_advance_status" -ne 1 ] || [ "$np1_signature_invariance_status" -ne 1 ] || [ "$np2_signature_invariance_status" -ne 1 ] || [ "$np4_signature_invariance_status" -ne 1 ] || [ "$parallel_consistency_status" -ne 1 ]; then
    final_status=0
    [ "$reasons" = "init" ] && reasons="invariance_or_hook_checks_failed"
  fi
fi

cat > stage11_outputs/stage11_7_oneway_sampling_parallel_consistency.dat <<EOD
stage11_7_requested_flag $requested_flag
stage11_7_build_status $build_status
stage11_7_baseline_run_status $baseline_run_status
stage11_7_hook_run_status $hook_run_status
stage11_7_hook_active_status $hook_active_status
stage11_7_sample_performed_status $sample_performed_status
stage11_7_sampled_velocity_finite_status $sampled_velocity_finite_status
stage11_7_no_field_modification_status $no_field_modification_status
stage11_7_no_rhs_modification_status $no_rhs_modification_status
stage11_7_no_rhs_injection_status $no_rhs_injection_status
stage11_7_no_ibm_spreading_status $no_ibm_spreading_status
stage11_7_no_feedback_force_status $no_feedback_force_status
stage11_7_no_twoway_force_status $no_twoway_force_status
stage11_7_no_structure_advance_status $no_structure_advance_status
stage11_7_np1_signature_invariance_status $np1_signature_invariance_status
stage11_7_np2_signature_invariance_status $np2_signature_invariance_status
stage11_7_np4_signature_invariance_status $np4_signature_invariance_status
stage11_7_parallel_consistency_status $parallel_consistency_status
stage11_7_oneway_sampling_parallel_consistency_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 11.7 ONEWAY SAMPLING PARALLEL CONSISTENCY VERDICT: PASS"
  echo "STAGE 11.7 FINAL VERDICT: PASS"
else
  [ "$reasons" = "init" ] && reasons="unknown_failure"
  echo "STAGE 11.7 ONEWAY SAMPLING PARALLEL CONSISTENCY VERDICT: FAIL"
  echo "STAGE 11.7 FINAL VERDICT: FAIL"
  echo "Reasons:$reasons"
  exit 1
fi
