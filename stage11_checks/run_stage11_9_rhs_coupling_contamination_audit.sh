#!/bin/sh
set -eu

BUILD_DIR=${BUILD_DIR:-build_stage9}
STAGE11_9_RUN_STAGE11_8=${STAGE11_9_RUN_STAGE11_8:-0}

mkdir -p stage11_outputs stage9_outputs

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

ensure_build_dir

if [ "$STAGE11_9_RUN_STAGE11_8" = "1" ]; then
    bash stage11_checks/run_stage11_8_io_restart_stats_visu_oneway.sh
fi

build_status=1
static_audit_status=1
velocity_intent_readonly_status=1
no_velocity_write_static_status=1
no_rhs_write_static_status=1
no_rhs_injection_static_status=1
no_ibm_spreading_static_status=1
no_feedback_force_static_status=1
no_twoway_force_static_status=1
no_structure_advance_static_status=1
runtime_smoke_status=1
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

active_tmp=$(mktemp)
viol_tmp=$(mktemp)
trap 'rm -f "$active_tmp" "$viol_tmp"' EXIT

extract_active() {
    in_file=$1
    out_file=$2
    awk '{
      line=$0
      if (line ~ /^[[:space:]]*!/) next
      sub(/!.*/, "", line)
      if (line ~ /^[[:space:]]*$/) next
      print tolower(line)
    }' "$in_file" > "$out_file"
}

contains_word() {
    word=$1
    text=$2
    echo "$text" | awk -v w="$word" '{ if ($0 ~ ("(^|[^a-z0-9_])" w "([^a-z0-9_]|$)")) exit 0; else exit 1 }'
}

check_forbidden_calls_uses() {
    in_file=$1
    file_tag=$2
    while IFS= read -r line; do
        case "$line" in
            *no_rhs_injection_status*|*no_ibm_spreading_status*|*no_feedback_force_status*|*no_twoway_force_status*|*no_structure_advance_status*)
                continue
                ;;
        esac

        case "$line" in
            *stage11_*rhs*|*rhs_injection*|*spread*|*spreading*|*feedback_force*|*twoway*|*two_way*|*structure_advance*)
                # only count executable/import path lines
                case "$line" in
                    *" use "*|use\ *|*" call "*|call\ *)
                        printf '%s|%s|%s\n' "$file_tag" "unknown" "$line" >> "$viol_tmp"
                        ;;
                esac
                ;;
        esac

        for token in fibre_ibm_spreading fibre_ibm_feedback fibre_ibm_force_buffer fibre_ibm_structure_coupling fibre_stage8_twoway_force_density fibre_stage8_feedback_candidate fibre_stage8_oneway_forcing fibre_stage7_force_spreading fibre_stage7_rhs_candidate fibre_stage5_rhs_wrapper fibre_stage5_spreading_candidate fibre_stage5_rhs_injection fibre_structure_integrator fibre_structure_advance fibre_tension_solver; do
            if contains_word "$token" "$line"; then
                case "$line" in
                    *" use "*|use\ *|*" call "*|call\ *)
                        printf '%s|%s|%s\n' "$file_tag" "unknown" "$line" >> "$viol_tmp"
                        ;;
                esac
            fi
        done
    done < "$in_file"
}

# preprocess active statements (read-only audit)
extract_active src/xcompact3d.f90 "$active_tmp"
check_forbidden_calls_uses "$active_tmp" "src/xcompact3d.f90"

extract_active src/fibre_stage11_production_oneway_hook.f90 "$active_tmp"
check_forbidden_calls_uses "$active_tmp" "src/fibre_stage11_production_oneway_hook.f90"

# intent check for production velocity arrays
awk '
{
  line=$0
  if (line ~ /^[[:space:]]*!/) next
  sub(/!.*/, "", line)
  line=tolower(line)
  if (line ~ /intent[[:space:]]*\([[:space:]]*inout[[:space:]]*\)/ || line ~ /intent[[:space:]]*\([[:space:]]*out[[:space:]]*\)/) {
    if (line ~ /\<ux\>/ || line ~ /\<uy\>/ || line ~ /\<uz\>/) {
      printf "src/fibre_stage11_production_oneway_hook.f90|%d|%s|velocity_intent_not_readonly\n", NR, $0
    }
  }
}' src/fibre_stage11_production_oneway_hook.f90 >> "$viol_tmp"

# write check: assignment to restricted LHS in production hook
awk '
{
  raw=$0
  if (raw ~ /^[[:space:]]*!/) next
  sub(/!.*/, "", raw)
  line=tolower(raw)
  if (line ~ /^[[:space:]]*$/) next
  if (line ~ /::/) next
  if (line ~ /=>/) next
  split(line, arr, "=")
  lhs=arr[1]
  if (length(arr) > 1) {
    gsub(/[[:space:]]+/, "", lhs)
    sub(/\(.*/, "", lhs)
    if (lhs ~ /^(ux|uy|uz|pressure|pp3|phi|dpdyx1|rhs|gx|gy|gz|force|fpx|fpy|fpz)$/) {
      printf "src/fibre_stage11_production_oneway_hook.f90|%d|%s|forbidden_write_lhs\n", NR, $0
    }
  }
}' src/fibre_stage11_production_oneway_hook.f90 >> "$viol_tmp"

if [ -s "$viol_tmp" ]; then
    static_audit_status=0
    velocity_intent_readonly_status=0
    no_velocity_write_static_status=0
    no_rhs_write_static_status=0
    no_rhs_injection_static_status=0
    no_ibm_spreading_static_status=0
    no_feedback_force_static_status=0
    no_twoway_force_static_status=0
    no_structure_advance_static_status=0
    final_status=0
    reasons="static_audit_violation"
    echo "Static audit violations detected:"
    while IFS='|' read -r f l txt why; do
        echo "file=$f line=$l offending_line=$txt reason=$why"
    done < "$viol_tmp"
else
    velocity_intent_readonly_status=1
    no_velocity_write_static_status=1
    no_rhs_write_static_status=1
    no_rhs_injection_static_status=1
    no_ibm_spreading_static_status=1
    no_feedback_force_static_status=1
    no_twoway_force_static_status=1
    no_structure_advance_static_status=1
fi

rm -f stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
LOG99=stage11_outputs/stage11_9_stage9_9_runtime_hook.log
if [ "$build_status" -eq 1 ]; then
    X3D_STAGE11_ONEWAY_HOOK=1 \
    X3D_STAGE11_FORCE_READONLY=1 \
    X3D_STAGE11_MAX_POINTS=8 \
    X3D_STAGE11_MAX_STEPS=3 \
    STAGE9_SKIP_PREREQS=1 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
    bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$LOG99" 2>&1 || runtime_smoke_status=0
else
    runtime_smoke_status=0
fi

if [ "$runtime_smoke_status" -eq 1 ]; then
    grep -q "STAGE 9.9 FINAL VERDICT: PASS" "$LOG99" || runtime_smoke_status=0
fi

HOOK=stage11_outputs/fibre_stage11_5_production_oneway_hook.dat
get_val(){ awk -v k="$1" '$1==k{print $2}' "$2"; }

if [ ! -f "$HOOK" ]; then
    final_status=0
    [ "$reasons" = "init" ] && reasons="missing_stage11_5_hook_diagnostics"
else
    for k in stage11_5_requested_flag stage11_5_readonly_mode_status stage11_5_hook_initialized_status stage11_5_hook_sample_called_status stage11_5_sample_performed_status stage11_5_sample_count_status stage11_5_sampled_velocity_finite_status stage11_5_sampled_velocity_signature_status stage11_5_field_modified_status stage11_5_rhs_modified_status stage11_5_no_rhs_injection_status stage11_5_no_ibm_spreading_status stage11_5_no_feedback_force_status stage11_5_no_twoway_force_status stage11_5_no_structure_advance_status stage11_5_production_oneway_hook_status; do
        v=$(get_val "$k" "$HOOK")
        case "$k" in
            stage11_5_field_modified_status|stage11_5_rhs_modified_status) exp=0 ;;
            *) exp=1 ;;
        esac
        if [ "$v" != "$exp" ]; then
            final_status=0
            [ "$reasons" = "init" ] && reasons="$k"
        fi
    done

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
fi

if [ "$final_status" -eq 1 ]; then
  if [ "$build_status" -ne 1 ] || [ "$static_audit_status" -ne 1 ] || [ "$runtime_smoke_status" -ne 1 ] || [ "$velocity_intent_readonly_status" -ne 1 ] || [ "$no_velocity_write_static_status" -ne 1 ] || [ "$no_rhs_write_static_status" -ne 1 ] || [ "$no_rhs_injection_static_status" -ne 1 ] || [ "$no_ibm_spreading_static_status" -ne 1 ] || [ "$no_feedback_force_static_status" -ne 1 ] || [ "$no_twoway_force_static_status" -ne 1 ] || [ "$no_structure_advance_static_status" -ne 1 ] || [ "$hook_active_status" -ne 1 ] || [ "$sample_performed_status" -ne 1 ] || [ "$sampled_velocity_finite_status" -ne 1 ] || [ "$no_field_modification_status" -ne 1 ] || [ "$no_rhs_modification_status" -ne 1 ] || [ "$no_rhs_injection_status" -ne 1 ] || [ "$no_ibm_spreading_status" -ne 1 ] || [ "$no_feedback_force_status" -ne 1 ] || [ "$no_twoway_force_status" -ne 1 ] || [ "$no_structure_advance_status" -ne 1 ]; then
    final_status=0
    [ "$reasons" = "init" ] && reasons="rhs_coupling_contamination_audit_failed"
  fi
fi

cat > stage11_outputs/stage11_9_rhs_coupling_contamination_audit.dat <<EOD
stage11_9_requested_flag $requested_flag
stage11_9_build_status $build_status
stage11_9_static_audit_status $static_audit_status
stage11_9_velocity_intent_readonly_status $velocity_intent_readonly_status
stage11_9_no_velocity_write_static_status $no_velocity_write_static_status
stage11_9_no_rhs_write_static_status $no_rhs_write_static_status
stage11_9_no_rhs_injection_static_status $no_rhs_injection_static_status
stage11_9_no_ibm_spreading_static_status $no_ibm_spreading_static_status
stage11_9_no_feedback_force_static_status $no_feedback_force_static_status
stage11_9_no_twoway_force_static_status $no_twoway_force_static_status
stage11_9_no_structure_advance_static_status $no_structure_advance_static_status
stage11_9_runtime_smoke_status $runtime_smoke_status
stage11_9_hook_active_status $hook_active_status
stage11_9_sample_performed_status $sample_performed_status
stage11_9_sampled_velocity_finite_status $sampled_velocity_finite_status
stage11_9_no_field_modification_status $no_field_modification_status
stage11_9_no_rhs_modification_status $no_rhs_modification_status
stage11_9_no_rhs_injection_status $no_rhs_injection_status
stage11_9_no_ibm_spreading_status $no_ibm_spreading_status
stage11_9_no_feedback_force_status $no_feedback_force_status
stage11_9_no_twoway_force_status $no_twoway_force_status
stage11_9_no_structure_advance_status $no_structure_advance_status
stage11_9_rhs_coupling_contamination_audit_status $final_status
EOD

if [ "$final_status" -eq 1 ]; then
  echo "STAGE 11.9 RHS COUPLING CONTAMINATION AUDIT VERDICT: PASS"
  echo "STAGE 11.9 FINAL VERDICT: PASS"
else
  [ "$reasons" = "init" ] && reasons="unknown_failure"
  echo "STAGE 11.9 RHS COUPLING CONTAMINATION AUDIT VERDICT: FAIL"
  echo "STAGE 11.9 FINAL VERDICT: FAIL"
  echo "Reasons:$reasons"
  exit 1
fi
