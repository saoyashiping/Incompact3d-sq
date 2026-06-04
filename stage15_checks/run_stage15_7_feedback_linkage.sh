#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_7_RUN_STAGE15_6=${STAGE15_7_RUN_STAGE15_6:-0}
STAGE15_7_REQUIRE_STAGE14_CLOSED=${STAGE15_7_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_7_REQUIRE_STAGE15_6=${STAGE15_7_REQUIRE_STAGE15_6:-1}
STAGE15_7_ENABLE=${STAGE15_7_ENABLE:-1}
STAGE15_7_CONTROLLED_STEP_ENABLE=${STAGE15_7_CONTROLLED_STEP_ENABLE:-1}
STAGE15_7_STRUCTURE_ADVANCE_ENABLE=${STAGE15_7_STRUCTURE_ADVANCE_ENABLE:-1}
STAGE15_7_DIAGNOSTIC_ONLY=${STAGE15_7_DIAGNOSTIC_ONLY:-1}
STAGE15_7_NP=${STAGE15_7_NP:-1}
STAGE15_7_NPTS=${STAGE15_7_NPTS:-8}
STAGE15_7_DT=${STAGE15_7_DT:-1.0e-4}
STAGE15_7_RHO_TILDE=${STAGE15_7_RHO_TILDE:-1.0}
STAGE15_7_TEST_FORCE=${STAGE15_7_TEST_FORCE:-1.0e-6}
STAGE15_7_FEEDBACK_ALPHA=${STAGE15_7_FEEDBACK_ALPHA:-1.0}
STAGE15_7_LAMBDA=${STAGE15_7_LAMBDA:-1.0e-8}
STAGE15_7_MAX_VELOCITY_UPDATE=${STAGE15_7_MAX_VELOCITY_UPDATE:-1.0e-9}
STAGE15_7_MAX_SLIP_ERROR=${STAGE15_7_MAX_SLIP_ERROR:-1.0e-14}
STAGE15_7_MAX_FORCE_ERROR=${STAGE15_7_MAX_FORCE_ERROR:-1.0e-14}
STAGE15_7_MAX_FORCE_RESPONSE=${STAGE15_7_MAX_FORCE_RESPONSE:-1.0e-8}
STAGE15_7_MAX_RHS_RESPONSE=${STAGE15_7_MAX_RHS_RESPONSE:-1.0e-12}
STAGE15_7_RUN_PRODUCTION_SMOKE=${STAGE15_7_RUN_PRODUCTION_SMOKE:-0}

OUTPUT_DIR=stage15_outputs
OUT_DAT=$OUTPUT_DIR/stage15_7_feedback_linkage.dat
CHECK_DAT=$OUTPUT_DIR/fibre_stage15_7_feedback_linkage.dat
CHECK_LOG=$OUTPUT_DIR/stage15_7_feedback_linkage.log
REASONS_FILE=$OUTPUT_DIR/stage15_7_feedback_linkage_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_7_feedback_linkage_static_audit_report.txt
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_6_SCRIPT=stage15_checks/run_stage15_6_controlled_structure_step_np1.sh
STAGE15_6_DOC=stage15_checks/stage15_6_controlled_structure_step_np1.md

BUILD_STATUS=1
RUN_STATUS=0
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
XCOMPACT_HOOK_STATUS=0
CONTROL_GUARD_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_6_STATUS=0
STAGE14_LAMBDA_GATE_ABSENT_STATUS=0
STAGE11_DIAGNOSTIC_STATUS=0
STAGE13_DIAGNOSTIC_STATUS=0
STAGE14_DIAGNOSTIC_STATUS=0
STAGE15_1_DIAGNOSTIC_STATUS=0
STAGE15_2_DIAGNOSTIC_STATUS=0
STAGE15_3_DIAGNOSTIC_STATUS=0
STAGE15_4_DIAGNOSTIC_STATUS=0
STAGE15_5_DIAGNOSTIC_STATUS=0
STAGE15_6_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOCS_AND_TARGET_STATUS=0
DIAGNOSTIC_STATUS=0
PRODUCTION_SMOKE_STATUS=0
PRODUCTION_SMOKE_DEFERRED_STATUS=0
FINAL_STATUS=0

mkdir -p "$OUTPUT_DIR" stage14_outputs stage13_outputs stage11_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$CHECK_DAT" "$CHECK_LOG"

add_reason() { echo "$1" >> "$REASONS_FILE"; }
static_note() { echo "$1" >> "$STATIC_REPORT"; }
search_silent() { pattern=$1; shift; grep -Eq "$pattern" "$@"; }
search_report() { pattern=$1; shift; grep -En "$pattern" "$@"; }

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        if [ -n "$DECOMP2D_ROOT" ]; then
            cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT"
        else
            cmake -S . -B "$BUILD_DIR"
        fi
    fi
}

build_target() {
    target=$1
    ensure_build_dir || return 1
    cmake --build "$BUILD_DIR" --target "$target" -j
}

stage15_7_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage15_feedback_linkage_check" \
               "$BUILD_DIR/src/fibre_stage15_feedback_linkage_check" \
               "$BUILD_DIR/fibre_stage15_feedback_linkage_check"; do
        if [ -x "$exe" ]; then
            echo "$exe"
            return 0
        fi
    done
    return 1
}

get_value() {
    file=$1
    key=$2
    awk -v key="$key" '$1 == key { print $2; found=1; exit } END { if (!found) exit 1 }' "$file"
}

require_key_value() {
    file=$1
    key=$2
    expected=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || { add_reason "missing_${key}_in_${file}"; return 1; }
    if [ "$value" != "$expected" ]; then
        add_reason "${key}_expected_${expected}_got_${value}_in_${file}"
        return 1
    fi
    return 0
}

require_finite_key() {
    file=$1
    key=$2
    value=$(get_value "$file" "$key" 2>/dev/null) || { add_reason "missing_${key}_in_${file}"; return 1; }
    awk -v v="$value" 'BEGIN { if (v ~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([EeDd][+-]?[0-9]+)?$/) exit 0; exit 1 }' || {
        add_reason "nonfinite_or_non_numeric_${key}_${value}_in_${file}"
        return 1
    }
    return 0
}

require_real_le_key() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || { add_reason "missing_${key}_in_${file}"; return 1; }
    awk -v value="$value" -v limit="$limit" 'BEGIN { d=value+0.0; if (d < 0.0) d=-d; if (d <= limit+0.0) exit 0; exit 1 }' || {
        add_reason "${key}_expected_abs_le_${limit}_got_${value}_in_${file}"
        return 1
    }
    return 0
}

scan_stage15_7_source_guards() {
    status=0
    files="src/fibre_stage15_feedback_linkage_check.f90"
    static_note "Stage 15.7 source guard audit"

    if awk '
      /^[[:space:]]*!/ { next }
      {
        raw=$0
        line=tolower($0)
        sub(/!.*/, "", line)
        if (line ~ /^[[:space:]]*call[[:space:]]+/) {
          routine=line
          sub(/^[[:space:]]*call[[:space:]]+/, "", routine)
          sub(/[[:space:]]*\(.*/, "", routine)
          gsub(/[[:space:]]/, "", routine)
          if (routine ~ /(bending|bend|tension|wall|contact|multifibre|multi_fibre)/ || \
              routine ~ /(poisson|projection|pressure|rk3|channel_forcing|xcompact3d)/) {
            print FILENAME ":" FNR ":" raw
            bad=1
          }
        }
      }
      END { exit(bad ? 1 : 0) }
    ' $files >> "$STATIC_REPORT"; then
        :
    else
        add_reason "stage15_7_forbidden_production_call_found"
        status=1
    fi

    if search_report '^[[:space:]]*use[[:space:]]+(fibre_stage14|xcompact3d|navier|poisson|schemes|derive)' \
       $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_7_forbidden_source_module_connection_found"
        status=1
    fi
    if search_report 'rg[[:space:]]|grep[[:space:]]+-R' stage15_checks/run_stage15_7_feedback_linkage.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_7_wrapper_must_not_require_ripgrep_or_recursive_grep"
        status=1
    fi

    if [ "$status" = "0" ]; then
        SOURCE_GUARD_STATUS=1
    fi
    return $status
}

check_xcompact_hook_insertion() {
    status=0
    use_count=$(grep -Ec 'use[[:space:]]+fibre_stage15_production_structure_hook' src/xcompact3d.f90 || true)
    reset_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_reset' src/xcompact3d.f90 || true)
    register_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_register\(1\)' src/xcompact3d.f90 || true)
    apply_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_apply\(itime,[[:space:]]*0\)' src/xcompact3d.f90 || true)
    linkage_count=$(grep -Ec 'stage15_feedback_linkage|fibre_stage15_feedback_linkage' src/xcompact3d.f90 || true)
    forbidden_count=$(grep -Ec 'call[[:space:]]+stage15_.*full.*advance|call[[:space:]]+.*bending|call[[:space:]]+.*tension' src/xcompact3d.f90 || true)
    printf 'xcompact3d hook counts: use=%s reset=%s register=%s apply=%s linkage=%s forbidden=%s\n' \
        "$use_count" "$reset_count" "$register_count" "$apply_count" "$linkage_count" "$forbidden_count" >> "$STATIC_REPORT"
    [ "$use_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_use_count_${use_count}"; status=1; }
    [ "$reset_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_reset_count_${reset_count}"; status=1; }
    [ "$register_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_register_count_${register_count}"; status=1; }
    [ "$apply_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_apply_count_${apply_count}"; status=1; }
    [ "$linkage_count" = "0" ] || { add_reason "xcompact3d_stage15_7_unapproved_linkage_connection_found"; status=1; }
    [ "$forbidden_count" = "0" ] || { add_reason "xcompact3d_forbidden_structure_solve_call_found"; status=1; }
    if [ "$status" = "0" ]; then
        XCOMPACT_HOOK_STATUS=1
        CONTROL_GUARD_STATUS=1
    fi
    return $status
}

run_static_audit() {
    status=0
    [ -f stage15_checks/run_stage15_7_feedback_linkage.sh ] || { add_reason "missing_stage15_7_wrapper"; status=1; }
    [ -f stage15_checks/stage15_7_feedback_linkage.md ] || { add_reason "missing_stage15_7_doc"; status=1; }
    if grep -q 'fibre_stage15_feedback_linkage_check' src/CMakeLists.txt && \
       grep -q 'fibre_stage15_feedback_linkage_check.f90' src/CMakeLists.txt; then
        DOCS_AND_TARGET_STATUS=1
    else
        add_reason "stage15_7_check_target_missing"
        status=1
    fi

    scan_stage15_7_source_guards || status=1
    check_xcompact_hook_insertion || status=1

    if grep -Eq 'stage14_get_injection_gain\(\)[[:space:]]*==[[:space:]]*0\.0|stage14_get_injection_gain\(\)[[:space:]]*\.eq\.[[:space:]]*0\.0' src/xcompact3d.f90 src/fibre_stage14_production_rhs_injection.f90; then
        add_reason "stage14_forbidden_lambda_zero_hook_registration_gate_found"
        status=1
    else
        STAGE14_LAMBDA_GATE_ABSENT_STATUS=1
    fi

    if search_silent 'fibre_stage11_5_production_oneway_hook.dat' src/fibre_stage11_production_oneway_hook.f90 && \
       search_silent 'stage11_5_production_oneway_hook_status' src/fibre_stage11_production_oneway_hook.f90; then
        STAGE11_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage11_5_production_oneway_hook_diagnostics_missing"
        status=1
    fi

    if search_silent 'fibre_stage13_6_production_force_density_candidate.dat' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_silent 'stage13_6_production_force_density_candidate_status' src/fibre_stage13_production_force_density_candidate.f90; then
        STAGE13_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage13_6_production_force_density_candidate_diagnostics_missing"
        status=1
    fi

    if search_silent 'fibre_stage14_5_production_rhs_hook.dat' src/fibre_stage14_production_rhs_injection.f90 && \
       search_silent 'stage14_5_nonzero_lambda_blocked_status' src/fibre_stage14_production_rhs_injection.f90 && \
       search_silent 'stage14_5_production_rhs_hook_status' src/fibre_stage14_production_rhs_injection.f90; then
        STAGE14_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage14_small_lambda_production_rhs_hook_diagnostics_missing"
        status=1
    fi

    if search_silent 'stage15_structure_state_write_diagnostics' src/fibre_stage15_structure_state.f90 && \
       search_silent 'fibre_stage15_1_structure_state_buffer.dat' src/fibre_stage15_structure_state_check.f90; then
        STAGE15_1_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_1_structure_state_diagnostics_missing"
        status=1
    fi

    if search_silent 'stage15_velocity_source_write_diagnostics' src/fibre_stage15_velocity_source_adapter.f90 && \
       search_silent 'fibre_stage15_2_velocity_source_adapter.dat' src/fibre_stage15_velocity_source_adapter_check.f90; then
        STAGE15_2_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_2_velocity_source_adapter_diagnostics_missing"
        status=1
    fi

    if search_silent 'stage15_structure_advance_formula_write_diagnostics' src/fibre_stage15_structure_advance_formula.f90 && \
       search_silent 'fibre_stage15_3_structure_advance_formula.dat' src/fibre_stage15_structure_advance_formula_check.f90; then
        STAGE15_3_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_3_structure_advance_formula_diagnostics_missing"
        status=1
    fi

    if search_silent 'stage15_production_structure_hook_write_diagnostics' src/fibre_stage15_production_structure_hook.f90 && \
       search_silent 'fibre_stage15_4_production_structure_hook.dat' src/fibre_stage15_production_structure_hook.f90; then
        STAGE15_4_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_4_production_structure_hook_diagnostics_missing"
        status=1
    fi

    if [ -f stage15_checks/run_stage15_5_structure_noop_invariance.sh ] && \
       search_silent 'fibre_stage15_5_structure_noop_invariance.dat' stage15_checks/run_stage15_5_structure_noop_invariance.sh; then
        STAGE15_5_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_5_structure_noop_invariance_diagnostics_missing"
        status=1
    fi

    if [ -f "$STAGE15_6_SCRIPT" ] && [ -s "$STAGE15_6_DOC" ] && \
       search_silent 'fibre_stage15_6_controlled_structure_step_np1.dat' "$STAGE15_6_SCRIPT"; then
        STAGE15_6_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_6_controlled_structure_step_diagnostics_missing"
        status=1
    fi

    if search_silent 'rank0_write_allowed' src/fibre_stage11_production_oneway_hook.f90 && \
       search_silent 'rank0_write_allowed' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_silent 'rank0_write_allowed' src/fibre_stage14_production_rhs_injection.f90 && \
       search_silent 'rank0_write_allowed' src/fibre_stage15_production_structure_hook.f90; then
        RANK0_DIAGNOSTIC_STATUS=1
    else
        add_reason "rank0_safe_diagnostic_writing_regressed"
        status=1
    fi

    if search_silent 'lbound\(ux, 1\) \+ 2' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_silent 'np=1/2/4' src/fibre_stage13_production_force_density_candidate.f90; then
        STAGE13_SAMPLING_REPAIR_STATUS=1
    else
        add_reason "stage13_force_density_sampling_repair_missing"
        status=1
    fi

    if [ "$STAGE15_7_REQUIRE_STAGE14_CLOSED" = "0" ]; then
        STAGE14_CLOSED_STATUS=1
    elif [ -s "$STAGE14_CLOSED_FILE" ]; then
        STAGE14_CLOSED_STATUS=1
    else
        add_reason "missing_stage14_checks_STAGE14_CLOSED_md"
        status=1
    fi

    if [ "$STAGE15_7_RUN_STAGE15_6" = "1" ]; then
        bash "$STAGE15_6_SCRIPT" || { add_reason "stage15_6_optional_prerequisite_run_failed"; status=1; }
    fi
    if [ "$STAGE15_7_REQUIRE_STAGE15_6" = "0" ]; then
        STAGE15_6_STATUS=1
    elif [ -s "$STAGE15_6_DOC" ] && [ -f "$STAGE15_6_SCRIPT" ] && [ -f src/fibre_stage15_controlled_structure_step_check.f90 ]; then
        STAGE15_6_STATUS=1
    else
        add_reason "stage15_6_closed_pass_evidence_missing"
        status=1
    fi

    if [ "$STAGE15_7_NP" = "1" ]; then
        :
    else
        add_reason "stage15_7_requires_np1_got_${STAGE15_7_NP}"
        status=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

run_check() {
    exe=$(stage15_7_exe) || { add_reason "missing_fibre_stage15_feedback_linkage_check_executable"; return 1; }
    export STAGE15_7_ENABLE STAGE15_7_CONTROLLED_STEP_ENABLE STAGE15_7_STRUCTURE_ADVANCE_ENABLE
    export STAGE15_7_DIAGNOSTIC_ONLY STAGE15_7_NP STAGE15_7_NPTS STAGE15_7_DT STAGE15_7_RHO_TILDE
    export STAGE15_7_TEST_FORCE STAGE15_7_FEEDBACK_ALPHA STAGE15_7_LAMBDA STAGE15_7_MAX_VELOCITY_UPDATE
    export STAGE15_7_MAX_SLIP_ERROR STAGE15_7_MAX_FORCE_ERROR STAGE15_7_MAX_FORCE_RESPONSE STAGE15_7_MAX_RHS_RESPONSE
    # shellcheck disable=SC2086
    "$MPIEXEC" $MPIEXEC_FLAGS -n 1 "$exe" > "$CHECK_LOG" 2>&1
}

handle_production_smoke() {
    if [ "$STAGE15_7_RUN_PRODUCTION_SMOKE" = "0" ]; then
        PRODUCTION_SMOKE_DEFERRED_STATUS=1
        PRODUCTION_SMOKE_STATUS=1
        echo "production_smoke_deferred: standalone feedback-linkage validation is authoritative for Stage 15.7" >> "$STATIC_REPORT"
    else
        PRODUCTION_SMOKE_STATUS=0
        add_reason "stage15_7_production_smoke_requested_but_not_enabled_for_feedback_linkage"
    fi
}

verify_diagnostics() {
    status=0
    [ -s "$CHECK_DAT" ] || { add_reason "missing_stage15_7_feedback_linkage_dat"; return 1; }
    require_key_value "$CHECK_DAT" stage15_7_requested_status 1 || status=1
    require_key_value "$CHECK_DAT" controlled_step_enabled_status 1 || status=1
    require_key_value "$CHECK_DAT" structure_advance_enable_status 1 || status=1
    require_key_value "$CHECK_DAT" diagnostic_only_status 1 || status=1
    require_key_value "$CHECK_DAT" np 1 || status=1
    require_key_value "$CHECK_DAT" npts "$STAGE15_7_NPTS" || status=1
    require_key_value "$CHECK_DAT" velocity_update_nonzero_status 1 || status=1
    require_key_value "$CHECK_DAT" slip_change_nonzero_status 1 || status=1
    require_key_value "$CHECK_DAT" slip_consistency_status 1 || status=1
    require_key_value "$CHECK_DAT" feedback_force_change_nonzero_status 1 || status=1
    require_key_value "$CHECK_DAT" feedback_force_consistency_status 1 || status=1
    require_key_value "$CHECK_DAT" force_response_bounded_status 1 || status=1
    require_key_value "$CHECK_DAT" rhs_response_bounded_status 1 || status=1
    require_key_value "$CHECK_DAT" controlled_update_count 1 || status=1
    require_key_value "$CHECK_DAT" production_full_structure_advance_count 0 || status=1
    require_key_value "$CHECK_DAT" bending_solve_count 0 || status=1
    require_key_value "$CHECK_DAT" tension_solve_count 0 || status=1
    require_key_value "$CHECK_DAT" wall_contact_count 0 || status=1
    require_key_value "$CHECK_DAT" multifibre_count 0 || status=1
    require_key_value "$CHECK_DAT" rhs_injection_connection_count 0 || status=1
    require_key_value "$CHECK_DAT" approved_stage12_13_14_chain_status 1 || status=1
    require_key_value "$CHECK_DAT" no_fluid_rhs_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_pressure_projection_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_poisson_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_rk3_channel_forcing_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" final_status 1 || status=1

    for key in dt rho_tilde test_force_magnitude feedback_alpha lambda_value max_velocity_update max_slip_change \
               slip_error max_feedback_force_change feedback_force_error; do
        require_finite_key "$CHECK_DAT" "$key" || status=1
    done
    require_real_le_key "$CHECK_DAT" max_velocity_update "$STAGE15_7_MAX_VELOCITY_UPDATE" || status=1
    require_real_le_key "$CHECK_DAT" slip_error "$STAGE15_7_MAX_SLIP_ERROR" || status=1
    require_real_le_key "$CHECK_DAT" feedback_force_error "$STAGE15_7_MAX_FORCE_ERROR" || status=1
    require_real_le_key "$CHECK_DAT" max_feedback_force_change "$STAGE15_7_MAX_FORCE_RESPONSE" || status=1

    if [ "$status" = "0" ]; then
        DIAGNOSTIC_STATUS=1
    fi
    return $status
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage15_7_requested_status 1
stage15_7_build_status $BUILD_STATUS
stage15_7_run_status $RUN_STATUS
stage15_7_static_audit_status $STATIC_AUDIT_STATUS
stage15_7_source_guard_status $SOURCE_GUARD_STATUS
stage15_7_xcompact_hook_status $XCOMPACT_HOOK_STATUS
stage15_7_control_guard_status $CONTROL_GUARD_STATUS
stage15_7_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_7_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_7_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_7_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_7_stage15_1_diagnostic_status $STAGE15_1_DIAGNOSTIC_STATUS
stage15_7_stage15_2_diagnostic_status $STAGE15_2_DIAGNOSTIC_STATUS
stage15_7_stage15_3_diagnostic_status $STAGE15_3_DIAGNOSTIC_STATUS
stage15_7_stage15_4_diagnostic_status $STAGE15_4_DIAGNOSTIC_STATUS
stage15_7_stage15_5_diagnostic_status $STAGE15_5_DIAGNOSTIC_STATUS
stage15_7_stage15_6_diagnostic_status $STAGE15_6_DIAGNOSTIC_STATUS
stage15_7_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_7_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_7_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_7_stage15_6_status $STAGE15_6_STATUS
stage15_7_docs_and_target_status $DOCS_AND_TARGET_STATUS
stage15_7_diagnostic_status $DIAGNOSTIC_STATUS
stage15_7_production_smoke_status $PRODUCTION_SMOKE_STATUS
stage15_7_production_smoke_deferred_status $PRODUCTION_SMOKE_DEFERRED_STATUS
stage15_7_np_value $STAGE15_7_NP
stage15_7_npts_value $STAGE15_7_NPTS
stage15_7_dt_value $STAGE15_7_DT
stage15_7_rho_tilde_value $STAGE15_7_RHO_TILDE
stage15_7_test_force_value $STAGE15_7_TEST_FORCE
stage15_7_feedback_alpha_value $STAGE15_7_FEEDBACK_ALPHA
stage15_7_lambda_value $STAGE15_7_LAMBDA
stage15_7_final_status $final_status
EOF_DAT
}

run_static_audit || true
handle_production_smoke

if ! build_target fibre_stage15_feedback_linkage_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage15_feedback_linkage_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_check; then
        RUN_STATUS=1
        verify_diagnostics || true
    else
        add_reason "run_failed_fibre_stage15_feedback_linkage_check"
    fi
else
    add_reason "stage15_7_feedback_linkage_check_skipped_due_to_build_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$RUN_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && \
   [ "$SOURCE_GUARD_STATUS" = "1" ] && [ "$XCOMPACT_HOOK_STATUS" = "1" ] && \
   [ "$CONTROL_GUARD_STATUS" = "1" ] && [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && \
   [ "$STAGE11_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE13_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE14_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_1_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_2_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_3_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_4_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_5_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_6_DIAGNOSTIC_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$STAGE15_6_STATUS" = "1" ] && [ "$DOCS_AND_TARGET_STATUS" = "1" ] && \
   [ "$DIAGNOSTIC_STATUS" = "1" ] && [ "$PRODUCTION_SMOKE_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_output_dat 1
    echo 'STAGE 15.7 FEEDBACK LINKAGE VERDICT: PASS'
    echo 'STAGE 15.7 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 15.7 FEEDBACK LINKAGE VERDICT: FAIL'
echo 'STAGE 15.7 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_7_failure'
fi
exit 1
