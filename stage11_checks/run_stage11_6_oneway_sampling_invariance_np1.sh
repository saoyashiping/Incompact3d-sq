#!/bin/sh
set -eu
BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_6_RUN_STAGE11_5=${STAGE11_6_RUN_STAGE11_5:-0}
STAGE11_6_INVARIANCE_ABS_TOL=${STAGE11_6_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE11_6_INVARIANCE_REL_TOL=${STAGE11_6_INVARIANCE_REL_TOL:-1.0e-14}
mkdir -p stage11_outputs stage9_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_6_RUN_STAGE11_5" = "1" ]; then
    bash stage11_checks/run_stage11_5_production_oneway_hook.sh
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
final_status=1
reasons=""

for tgt in xcompact3d fibre_stage10_config_check fibre_stage10_noop_hook_check fibre_stage11_config_check fibre_stage11_lagrangian_state_check fibre_stage11_grid_metadata_check fibre_stage11_oneway_interpolation_check fibre_stage11_controlled_interpolation_check fibre_stage11_production_oneway_hook_check; do
    cmake --build "$BUILD_DIR" --target "$tgt" -j || build_status=0
done

if [ "$build_status" -eq 1 ]; then
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh || baseline_run_status=0
else
    baseline_run_status=0
fi

if [ "$baseline_run_status" -eq 1 ]; then
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat stage11_outputs/stage11_6_baseline_np1.dat || baseline_run_status=0
fi

if [ "$build_status" -eq 1 ]; then
    X3D_STAGE11_ONEWAY_HOOK=1 \
    X3D_STAGE11_FORCE_READONLY=1 \
    X3D_STAGE11_MAX_POINTS=8 \
    X3D_STAGE11_MAX_STEPS=3 \
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh || hook_run_status=0
else
    hook_run_status=0
fi

if [ "$hook_run_status" -eq 1 ]; then
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat stage11_outputs/stage11_6_hook_np1.dat || hook_run_status=0
fi

BASE=stage11_outputs/stage11_6_baseline_np1.dat
HOOK=stage11_outputs/stage11_6_hook_np1.dat
H11=stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
get_val() { awk -v k="$1" '$1==k{print $2}' "$2"; }

if [ -f "$H11" ]; then
    hook_active_status=$(get_val stage11_5_hook_sample_called_status "$H11")
    sample_performed_status=$(get_val stage11_5_sample_performed_status "$H11")
    sampled_velocity_finite_status=$(get_val stage11_5_sampled_velocity_finite_status "$H11")
    fmod=$(get_val stage11_5_field_modified_status "$H11")
    rmod=$(get_val stage11_5_rhs_modified_status "$H11")
    [ "$fmod" = "0" ] && no_field_modification_status=1 || no_field_modification_status=0
    [ "$rmod" = "0" ] && no_rhs_modification_status=1 || no_rhs_modification_status=0
    no_rhs_injection_status=$(get_val stage11_5_no_rhs_injection_status "$H11")
    no_ibm_spreading_status=$(get_val stage11_5_no_ibm_spreading_status "$H11")
    no_feedback_force_status=$(get_val stage11_5_no_feedback_force_status "$H11")
    no_twoway_force_status=$(get_val stage11_5_no_twoway_force_status "$H11")
    no_structure_advance_status=$(get_val stage11_5_no_structure_advance_status "$H11")
else
    np1_signature_invariance_status=0
fi

compare_metric() {
    metric=$1
    ref=$(get_val "$metric" "$BASE")
    cur=$(get_val "$metric" "$HOOK")
    res=$(awk -v r="$ref" -v c="$cur" -v at="$STAGE11_6_INVARIANCE_ABS_TOL" -v rt="$STAGE11_6_INVARIANCE_REL_TOL" 'BEGIN{d=c-r; ad=(d<0?-d:d); ar=(r<0?-r:r); eff=at; tmp=rt*((ar>1)?ar:1); if(tmp>eff)eff=tmp; pass=(ad<=eff)?1:0; printf "%.17e %.17e %.17e %.17e %d", d,at,tmp,eff,pass}')
    delta=$(echo "$res" | awk '{print $1}')
    abs_tol=$(echo "$res" | awk '{print $2}')
    rel_term=$(echo "$res" | awk '{print $3}')
    eff_tol=$(echo "$res" | awk '{print $4}')
    pass=$(echo "$res" | awk '{print $5}')
    echo "metric=$metric baseline=$ref hook=$cur delta=$delta abs_tol=$abs_tol rel_tol=$STAGE11_6_INVARIANCE_REL_TOL effective_tol=$eff_tol result=$( [ "$pass" = "1" ] && echo PASS || echo FAIL )"
    [ "$pass" = "1" ] || np1_signature_invariance_status=0
}

if [ "$baseline_run_status" -eq 1 ] && [ "$hook_run_status" -eq 1 ] && [ -f "$BASE" ] && [ -f "$HOOK" ]; then
    for m in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        compare_metric "$m"
    done
else
    np1_signature_invariance_status=0
fi

requested_flag=0
[ "$hook_active_status" = "1" ] && requested_flag=1

if [ "$build_status" -ne 1 ] || [ "$baseline_run_status" -ne 1 ] || [ "$hook_run_status" -ne 1 ] || [ "$hook_active_status" -ne 1 ] || [ "$sample_performed_status" -ne 1 ] || [ "$sampled_velocity_finite_status" -ne 1 ] || [ "$no_field_modification_status" -ne 1 ] || [ "$no_rhs_modification_status" -ne 1 ] || [ "$no_rhs_injection_status" -ne 1 ] || [ "$no_ibm_spreading_status" -ne 1 ] || [ "$no_feedback_force_status" -ne 1 ] || [ "$no_twoway_force_status" -ne 1 ] || [ "$no_structure_advance_status" -ne 1 ] || [ "$np1_signature_invariance_status" -ne 1 ]; then
    final_status=0
fi

cat > stage11_outputs/stage11_6_oneway_sampling_invariance_np1.dat <<EOD
stage11_6_requested_flag $requested_flag
stage11_6_build_status $build_status
stage11_6_baseline_run_status $baseline_run_status
stage11_6_hook_run_status $hook_run_status
stage11_6_hook_active_status $hook_active_status
stage11_6_sample_performed_status $sample_performed_status
stage11_6_sampled_velocity_finite_status $sampled_velocity_finite_status
stage11_6_no_field_modification_status $no_field_modification_status
stage11_6_no_rhs_modification_status $no_rhs_modification_status
stage11_6_no_rhs_injection_status $no_rhs_injection_status
stage11_6_no_ibm_spreading_status $no_ibm_spreading_status
stage11_6_no_feedback_force_status $no_feedback_force_status
stage11_6_no_twoway_force_status $no_twoway_force_status
stage11_6_no_structure_advance_status $no_structure_advance_status
stage11_6_np1_signature_invariance_status $np1_signature_invariance_status
stage11_6_oneway_sampling_invariance_np1_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
    echo "STAGE 11.6 ONEWAY SAMPLING INVARIANCE NP1 VERDICT: PASS"
    echo "STAGE 11.6 FINAL VERDICT: PASS"
else
    echo "STAGE 11.6 ONEWAY SAMPLING INVARIANCE NP1 VERDICT: FAIL"
    echo "STAGE 11.6 FINAL VERDICT: FAIL"
    exit 1
fi
