#!/bin/sh
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
STAGE13_DAT="$OUTPUT_DIR/fibre_stage13_6_production_force_density_candidate.dat"
OUT_DAT="$OUTPUT_DIR/stage13_10_rhs_injection_contamination_audit.dat"
REASONS_FILE="$OUTPUT_DIR/stage13_10_rhs_injection_contamination_audit_reasons.tmp"
RUNTIME_LOG="$OUTPUT_DIR/stage13_10_stage9_9_rhs_injection_contamination_audit.log"
STATIC_FILES="src/xcompact3d.f90 src/fibre_stage13_production_force_density_candidate.f90"
PROD_FILE=src/fibre_stage13_production_force_density_candidate.f90
XCOMPACT_FILE=src/xcompact3d.f90

mkdir -p stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$OUT_DAT" "$STAGE13_DAT" "$RUNTIME_LOG"

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

add_reason() {
    echo "$1" >> "$REASONS_FILE"
}

get_dat_value() {
    file=$1
    key=$2
    awk -v key="$key" '$1 == key { value=$2 } END { if (value != "") print value }' "$file"
}

check_key_equals() {
    file=$1
    key=$2
    expected=$3
    value=$(get_dat_value "$file" "$key")
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$key expected $expected but found ${value:-missing}"
    return 1
}

write_output_dat() {
    final_status=$1
    {
        echo "stage13_10_requested_flag 1"
        echo "stage13_10_build_status $build_status"
        echo "stage13_10_static_audit_status $static_audit_status"
        echo "stage13_10_velocity_intent_readonly_status $velocity_intent_readonly_status"
        echo "stage13_10_no_velocity_write_static_status $no_velocity_write_static_status"
        echo "stage13_10_no_pressure_write_static_status $no_pressure_write_static_status"
        echo "stage13_10_no_rhs_write_static_status $no_rhs_write_static_status"
        echo "stage13_10_no_rhs_injection_static_status $no_rhs_injection_static_status"
        echo "stage13_10_no_production_ibm_forcing_static_status $no_production_ibm_forcing_static_status"
        echo "stage13_10_no_feedback_application_static_status $no_feedback_application_static_status"
        echo "stage13_10_no_twoway_force_static_status $no_twoway_force_static_status"
        echo "stage13_10_no_structure_advance_static_status $no_structure_advance_static_status"
        echo "stage13_10_hook_placement_status $hook_placement_status"
        echo "stage13_10_runtime_smoke_status $runtime_smoke_status"
        echo "stage13_10_hook_active_status $hook_active_status"
        echo "stage13_10_force_density_candidate_computed_status $force_density_candidate_computed_status"
        echo "stage13_10_force_density_candidate_finite_status $force_density_candidate_finite_status"
        echo "stage13_10_force_density_norm_finite_status $force_density_norm_finite_status"
        echo "stage13_10_integrated_force_finite_status $integrated_force_finite_status"
        echo "stage13_10_integrated_force_conservation_status $integrated_force_conservation_status"
        echo "stage13_10_spreading_input_sign_status $spreading_input_sign_status"
        echo "stage13_10_wrong_sign_rejection_status $wrong_sign_rejection_status"
        echo "stage13_10_no_field_modification_status $no_field_modification_status"
        echo "stage13_10_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage13_10_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_10_no_production_ibm_forcing_status $no_production_ibm_forcing_status"
        echo "stage13_10_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_10_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_10_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_10_rhs_injection_contamination_audit_status $final_status"
    } > "$OUT_DAT"
}

run_forbidden_call_import_audit() {
    for file in $STATIC_FILES; do
        awk -v file="$file" '
          function trim(s) { sub(/^[[:space:]]+/, "", s); sub(/[[:space:]]+$/, "", s); return s }
          function has_forbidden(s,    i) {
            for (i = 1; i <= n; i++) if (index(s, forbidden[i]) > 0) return forbidden[i]
            return ""
          }
          BEGIN {
            n = split("rhs_injection inject_rhs add_to_rhs apply_rhs production_ibm_forcing ibm_forcing_application apply_ibm_force apply_force_to_fluid feedback_application apply_feedback twoway two_way two_way_force fibre_stage8_twoway_force_density fibre_stage8_feedback_candidate fibre_stage8_oneway_forcing fibre_stage7_rhs_candidate fibre_stage7_force_spreading fibre_stage5_rhs_wrapper fibre_stage5_rhs_injection fibre_structure_integrator fibre_structure_advance advance_structure structure_advance fibre_tension_solver", forbidden, " ")
            found = 0
          }
          {
            raw = $0
            sub(/!.*/, "", raw)
            stmt = trim(tolower(raw))
            if (stmt == "") next
            if (stmt ~ /^(use|call)[[:space:]]/) {
              hit = has_forbidden(stmt)
              if (hit != "") {
                printf("file: %s line number: %d offending line: %s reason: forbidden active use/call token %s\n", file, FNR, stmt, hit)
                found = 1
              }
            }
          }
          END { exit(found ? 1 : 0) }
        ' "$file" >> "$REASONS_FILE" || static_audit_status=0
    done
}

run_velocity_intent_audit() {
    awk -v file="$PROD_FILE" '
      function trim(s) { sub(/^[[:space:]]+/, "", s); sub(/[[:space:]]+$/, "", s); return s }
      function has_token(s, tok) { return s ~ "(^|[^a-z0-9_])" tok "([^a-z0-9_]|$)" }
      BEGIN { found = 0; ux_in = 0; uy_in = 0; uz_in = 0 }
      {
        raw = $0
        sub(/!.*/, "", raw)
        stmt = trim(tolower(raw))
        if (stmt == "") next
        if (stmt ~ /intent[[:space:]]*\([[:space:]]*(inout|out)[[:space:]]*\)/ && \
            (has_token(stmt, "ux") || has_token(stmt, "uy") || has_token(stmt, "uz"))) {
          printf("file: %s line number: %d offending line: %s reason: production velocity dummy has writable intent\n", file, FNR, stmt)
          found = 1
        }
        if (stmt ~ /intent[[:space:]]*\([[:space:]]*in[[:space:]]*\)/ && has_token(stmt, "ux")) ux_in = 1
        if (stmt ~ /intent[[:space:]]*\([[:space:]]*in[[:space:]]*\)/ && has_token(stmt, "uy")) uy_in = 1
        if (stmt ~ /intent[[:space:]]*\([[:space:]]*in[[:space:]]*\)/ && has_token(stmt, "uz")) uz_in = 1
      }
      END {
        if (ux_in != 1) { printf("file: %s line number: 0 offending line: ux reason: missing intent(in) production velocity declaration\n", file); found = 1 }
        if (uy_in != 1) { printf("file: %s line number: 0 offending line: uy reason: missing intent(in) production velocity declaration\n", file); found = 1 }
        if (uz_in != 1) { printf("file: %s line number: 0 offending line: uz reason: missing intent(in) production velocity declaration\n", file); found = 1 }
        exit(found ? 1 : 0)
      }
    ' "$PROD_FILE" >> "$REASONS_FILE" || velocity_intent_readonly_status=0
}

run_production_field_write_audit() {
    awk -v file="$PROD_FILE" '
      function trim(s) { sub(/^[[:space:]]+/, "", s); sub(/[[:space:]]+$/, "", s); return s }
      function has_token(s, tok) { return s ~ "(^|[^a-z0-9_])" tok "([^a-z0-9_]|$)" }
      BEGIN {
        found = 0
        n = split("ux uy uz pressure pp3 phi dpdyx1 rhs gx gy gz fpx fpy fpz eulerian_force production_force ibm_force", fields, " ")
      }
      {
        raw = $0
        sub(/!.*/, "", raw)
        stmt = trim(tolower(raw))
        if (stmt == "") next
        if (index(stmt, "=") == 0) next
        lhs = stmt
        sub(/=.*/, "", lhs)
        lhs = trim(lhs)
        if (lhs ~ /::/) next
        if (lhs ~ /^(if|do|where|else|select)[[:space:]\(]/) next
        if (lhs ~ /[<>\/]/) next
        for (i = 1; i <= n; i++) {
          if (has_token(lhs, fields[i])) {
            printf("file: %s line number: %d offending line: %s reason: active assignment to production field token %s\n", file, FNR, stmt, fields[i])
            found = 1
            if (fields[i] == "ux" || fields[i] == "uy" || fields[i] == "uz") velocity = 1
            if (fields[i] == "pressure" || fields[i] == "pp3" || fields[i] == "phi" || fields[i] == "dpdyx1") pressure = 1
            if (fields[i] == "rhs" || fields[i] == "gx" || fields[i] == "gy" || fields[i] == "gz" || fields[i] == "fpx" || fields[i] == "fpy" || fields[i] == "fpz" || fields[i] == "eulerian_force" || fields[i] == "production_force" || fields[i] == "ibm_force") rhs = 1
          }
        }
      }
      END {
        if (velocity) print "STATIC_FLAG velocity_write"
        if (pressure) print "STATIC_FLAG pressure_write"
        if (rhs) print "STATIC_FLAG rhs_write"
        exit(found ? 1 : 0)
      }
    ' "$PROD_FILE" > "$OUTPUT_DIR/stage13_10_static_write_audit.tmp"
    rc=$?
    if [ $rc -ne 0 ]; then
        sed '/^STATIC_FLAG /d' "$OUTPUT_DIR/stage13_10_static_write_audit.tmp" >> "$REASONS_FILE"
        if grep '^STATIC_FLAG velocity_write' "$OUTPUT_DIR/stage13_10_static_write_audit.tmp" >/dev/null 2>&1; then
            no_velocity_write_static_status=0
        fi
        if grep '^STATIC_FLAG pressure_write' "$OUTPUT_DIR/stage13_10_static_write_audit.tmp" >/dev/null 2>&1; then
            no_pressure_write_static_status=0
        fi
        if grep '^STATIC_FLAG rhs_write' "$OUTPUT_DIR/stage13_10_static_write_audit.tmp" >/dev/null 2>&1; then
            no_rhs_write_static_status=0
        fi
        static_audit_status=0
    fi
    rm -f "$OUTPUT_DIR/stage13_10_static_write_audit.tmp"
}

run_hook_placement_audit() {
    awk -v file="$XCOMPACT_FILE" '
      function trim(s) { sub(/^[[:space:]]+/, "", s); sub(/[[:space:]]+$/, "", s); return s }
      BEGIN { found = 0; stage12_sample = 0; stage13_sample = 0; in_rhs = 0; guard_seen = 0 }
      {
        raw = $0
        sub(/!.*/, "", raw)
        stmt = trim(tolower(raw))
        if (stmt == "") next
        active[FNR] = stmt
        text = text " " stmt
        if (index(stmt, "call stage10_hook_pre_rhs") > 0) in_rhs = 1
        if (index(stmt, "call stage10_hook_post_projection") > 0) in_rhs = 0
        if (index(stmt, "call stage12_production_feedback_candidate_sample") > 0 && stage12_sample == 0) stage12_sample = FNR
        if (index(stmt, "call stage13_production_force_density_candidate_sample") > 0 && stage13_sample == 0) stage13_sample = FNR
        if (index(stmt, "stage13_force_density_reg = stage13_requested()") > 0) guard_seen = 1
        if (index(stmt, "call stage13_production_force_density_candidate_") > 0) {
          if (stmt !~ /^if[[:space:]]*\([[:space:]]*stage13_force_density_reg[[:space:]]*\)[[:space:]]*call[[:space:]]+stage13_production_force_density_candidate_/) {
            printf("file: %s line number: %d offending line: %s reason: Stage 13 hook call is not guarded by stage13_force_density_reg\n", file, FNR, stmt)
            found = 1
          }
          if (in_rhs) {
            printf("file: %s line number: %d offending line: %s reason: Stage 13 hook call appears inside RHS/projection region\n", file, FNR, stmt)
            found = 1
          }
        }
      }
      END {
        if (text !~ /stage13_requested\(\)/ || text !~ /stage13_readonly_mode\(\)/ || text !~ /stage13_spreading_readonly_mode\(\)/ || guard_seen != 1) {
          printf("file: %s line number: 0 offending line: stage13_force_density_reg reason: missing required Stage 13 read-only guard predicates\n", file)
          found = 1
        }
        if (stage12_sample == 0 || stage13_sample == 0 || stage13_sample <= stage12_sample) {
          printf("file: %s line number: %d offending line: stage13_production_force_density_candidate_sample reason: Stage 13 sample is not after Stage 12 sample\n", file, stage13_sample)
          found = 1
        }
        exit(found ? 1 : 0)
      }
    ' "$XCOMPACT_FILE" >> "$REASONS_FILE" || hook_placement_status=0
}

run_static_audit() {
    run_forbidden_call_import_audit
    run_velocity_intent_audit
    run_production_field_write_audit
    run_hook_placement_audit
    if [ "$velocity_intent_readonly_status" = "0" ] || [ "$no_velocity_write_static_status" = "0" ] || \
       [ "$no_pressure_write_static_status" = "0" ] || [ "$no_rhs_write_static_status" = "0" ] || \
       [ "$hook_placement_status" = "0" ]; then
        static_audit_status=0
    fi
    if [ "$static_audit_status" = "1" ]; then
        no_rhs_injection_static_status=1
        no_production_ibm_forcing_static_status=1
        no_feedback_application_static_status=1
        no_twoway_force_static_status=1
        no_structure_advance_static_status=1
    else
        no_rhs_injection_static_status=0
        no_production_ibm_forcing_static_status=0
        no_feedback_application_static_status=0
        no_twoway_force_static_status=0
        no_structure_advance_static_status=0
    fi
}

validate_runtime_diagnostics() {
    if [ ! -f "$STAGE13_DAT" ]; then
        add_reason "missing_stage13_6_force_density_candidate_diagnostics"
        return 1
    fi
    check_key_equals "$STAGE13_DAT" stage13_6_requested_flag 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_readonly_mode_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_spreading_readonly_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_hook_initialized_status 1 >/dev/null && hook_active_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_hook_sample_called_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_sampled_velocity_available_status 1 >/dev/null || true
    check_key_equals "$STAGE13_DAT" stage13_6_force_density_candidate_computed_status 1 >/dev/null && force_density_candidate_computed_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_force_density_candidate_finite_status 1 >/dev/null && force_density_candidate_finite_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_force_density_norm_finite_status 1 >/dev/null && force_density_norm_finite_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_integrated_force_finite_status 1 >/dev/null && integrated_force_finite_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_integrated_force_conservation_status 1 >/dev/null && integrated_force_conservation_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_spreading_input_sign_status 1 >/dev/null && spreading_input_sign_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_wrong_sign_rejection_status 1 >/dev/null && wrong_sign_rejection_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_field_modified_status 0 >/dev/null && no_field_modification_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_rhs_modified_status 0 >/dev/null && no_rhs_modification_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_rhs_injection_status 1 >/dev/null && no_rhs_injection_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_production_ibm_forcing_status 1 >/dev/null && no_production_ibm_forcing_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_feedback_application_status 1 >/dev/null && no_feedback_application_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_twoway_force_status 1 >/dev/null && no_twoway_force_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_no_structure_advance_status 1 >/dev/null && no_structure_advance_status=1
    check_key_equals "$STAGE13_DAT" stage13_6_production_force_density_candidate_status 1 >/dev/null || true
}

build_status=1
static_audit_status=1
velocity_intent_readonly_status=1
no_velocity_write_static_status=1
no_pressure_write_static_status=1
no_rhs_write_static_status=1
no_rhs_injection_static_status=1
no_production_ibm_forcing_static_status=1
no_feedback_application_static_status=1
no_twoway_force_static_status=1
no_structure_advance_static_status=1
hook_placement_status=1
runtime_smoke_status=0
hook_active_status=0
force_density_candidate_computed_status=0
force_density_candidate_finite_status=0
force_density_norm_finite_status=0
integrated_force_finite_status=0
integrated_force_conservation_status=0
spreading_input_sign_status=0
wrong_sign_rejection_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_rhs_injection_status=0
no_production_ibm_forcing_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0

STAGE13_10_RUN_STAGE13_9=${STAGE13_10_RUN_STAGE13_9:-0}
if [ "$STAGE13_10_RUN_STAGE13_9" = "1" ]; then
    if ! sh stage13_checks/run_stage13_9_io_restart_stats_visu_force_density.sh; then
        add_reason "optional Stage 13.9 prerequisite failed"
    fi
fi

run_static_audit
ensure_build_dir

targets="xcompact3d \
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
fibre_stage12_production_feedback_candidate_check \
fibre_stage13_config_check \
fibre_stage13_force_density_buffer_check \
fibre_stage13_spreading_kernel_check \
fibre_stage13_volume_normalization_audit_check \
fibre_stage13_conservation_sign_audit_check \
fibre_stage13_production_force_density_candidate_check"

for target in $targets; do
    if ! cmake --build "$BUILD_DIR" --target "$target" -j; then
        build_status=0
        add_reason "build failed for $target"
    fi
done

if [ "$build_status" = "1" ]; then
    rm -f "$STAGE13_DAT"
    X3D_STAGE11_ONEWAY_HOOK=1 \
    X3D_STAGE11_FORCE_READONLY=1 \
    X3D_STAGE11_MAX_POINTS=8 \
    X3D_STAGE11_MAX_STEPS=3 \
    X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
    X3D_STAGE12_FORCE_READONLY=1 \
    X3D_STAGE12_FEEDBACK_GAIN=1.0 \
    X3D_STAGE12_FORCE_SIGN=1 \
    X3D_STAGE12_MAX_POINTS=8 \
    X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1 \
    X3D_STAGE13_FORCE_READONLY=1 \
    X3D_STAGE13_SPREADING_READONLY=1 \
    X3D_STAGE13_MAX_POINTS=8 \
    X3D_STAGE13_MAX_EULERIAN_POINTS=64 \
    X3D_STAGE13_SPREADING_NORMALIZATION=conservative \
    STAGE9_SKIP_PREREQS=1 \
    STAGE9_9_MAX_STEPS=3 \
    X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
    X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
    X3D_STAGE9_9_MAX_STEPS=3 \
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$RUNTIME_LOG" 2>&1
    if [ $? -eq 0 ] && grep 'STAGE 9.9 FINAL VERDICT: PASS' "$RUNTIME_LOG" >/dev/null 2>&1; then
        runtime_smoke_status=1
    else
        add_reason "Stage 9.9 deterministic runtime smoke did not report PASS under Stage 13 hook environment"
    fi
    validate_runtime_diagnostics
else
    add_reason "runtime smoke skipped because required build failed"
fi

if [ ! -s "$REASONS_FILE" ] && [ "$build_status" = "1" ] && [ "$static_audit_status" = "1" ] && \
   [ "$velocity_intent_readonly_status" = "1" ] && [ "$no_velocity_write_static_status" = "1" ] && \
   [ "$no_pressure_write_static_status" = "1" ] && [ "$no_rhs_write_static_status" = "1" ] && \
   [ "$no_rhs_injection_static_status" = "1" ] && [ "$no_production_ibm_forcing_static_status" = "1" ] && \
   [ "$no_feedback_application_static_status" = "1" ] && [ "$no_twoway_force_static_status" = "1" ] && \
   [ "$no_structure_advance_static_status" = "1" ] && [ "$hook_placement_status" = "1" ] && \
   [ "$runtime_smoke_status" = "1" ] && [ "$hook_active_status" = "1" ] && \
   [ "$force_density_candidate_computed_status" = "1" ] && [ "$force_density_candidate_finite_status" = "1" ] && \
   [ "$force_density_norm_finite_status" = "1" ] && [ "$integrated_force_finite_status" = "1" ] && \
   [ "$integrated_force_conservation_status" = "1" ] && [ "$spreading_input_sign_status" = "1" ] && \
   [ "$wrong_sign_rejection_status" = "1" ] && [ "$no_field_modification_status" = "1" ] && \
   [ "$no_rhs_modification_status" = "1" ] && [ "$no_rhs_injection_status" = "1" ] && \
   [ "$no_production_ibm_forcing_status" = "1" ] && [ "$no_feedback_application_status" = "1" ] && \
   [ "$no_twoway_force_status" = "1" ] && [ "$no_structure_advance_status" = "1" ]; then
    write_output_dat 1
    rm -f "$REASONS_FILE"
    echo "STAGE 13.10 RHS INJECTION CONTAMINATION AUDIT VERDICT: PASS"
    echo "STAGE 13.10 FINAL VERDICT: PASS"
    exit 0
fi

if [ ! -s "$REASONS_FILE" ]; then
    add_reason "Stage 13.10 gate failed without a captured reason"
fi
write_output_dat 0
echo "STAGE 13.10 RHS INJECTION CONTAMINATION AUDIT VERDICT: FAIL"
echo "STAGE 13.10 FINAL VERDICT: FAIL"
echo "Reasons:"
sed 's/^/ - /' "$REASONS_FILE"
rm -f "$REASONS_FILE"
exit 1
