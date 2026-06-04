#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_8_RUN_STAGE15_7=${STAGE15_8_RUN_STAGE15_7:-0}
STAGE15_8_REQUIRE_STAGE14_CLOSED=${STAGE15_8_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_8_REQUIRE_STAGE15_7=${STAGE15_8_REQUIRE_STAGE15_7:-1}
STAGE15_8_ENABLE=${STAGE15_8_ENABLE:-1}
STAGE15_8_CONTROLLED_STEP_ENABLE=${STAGE15_8_CONTROLLED_STEP_ENABLE:-1}
STAGE15_8_STRUCTURE_ADVANCE_ENABLE=${STAGE15_8_STRUCTURE_ADVANCE_ENABLE:-1}
STAGE15_8_DIAGNOSTIC_ONLY=${STAGE15_8_DIAGNOSTIC_ONLY:-1}
STAGE15_8_NP_LIST=${STAGE15_8_NP_LIST:-"1 2 4"}
STAGE15_8_NPTS=${STAGE15_8_NPTS:-8}
STAGE15_8_DT=${STAGE15_8_DT:-1.0e-4}
STAGE15_8_RHO_TILDE=${STAGE15_8_RHO_TILDE:-1.0}
STAGE15_8_TEST_FORCE=${STAGE15_8_TEST_FORCE:-1.0e-6}
STAGE15_8_FEEDBACK_ALPHA=${STAGE15_8_FEEDBACK_ALPHA:-1.0}
STAGE15_8_LAMBDA=${STAGE15_8_LAMBDA:-1.0e-8}
STAGE15_8_MAX_VELOCITY_UPDATE=${STAGE15_8_MAX_VELOCITY_UPDATE:-1.0e-9}
STAGE15_8_MAX_SLIP_ERROR=${STAGE15_8_MAX_SLIP_ERROR:-1.0e-14}
STAGE15_8_MAX_FORCE_ERROR=${STAGE15_8_MAX_FORCE_ERROR:-1.0e-14}
STAGE15_8_MAX_FORCE_RESPONSE=${STAGE15_8_MAX_FORCE_RESPONSE:-1.0e-8}
STAGE15_8_MAX_RHS_RESPONSE=${STAGE15_8_MAX_RHS_RESPONSE:-1.0e-12}
STAGE15_8_MAX_PARALLEL_VELOCITY_DIFF=${STAGE15_8_MAX_PARALLEL_VELOCITY_DIFF:-1.0e-14}
STAGE15_8_MAX_PARALLEL_SLIP_DIFF=${STAGE15_8_MAX_PARALLEL_SLIP_DIFF:-1.0e-14}
STAGE15_8_MAX_PARALLEL_FORCE_DIFF=${STAGE15_8_MAX_PARALLEL_FORCE_DIFF:-1.0e-14}
STAGE15_8_MAX_PARALLEL_RHS_DIFF=${STAGE15_8_MAX_PARALLEL_RHS_DIFF:-1.0e-14}
STAGE15_8_RUN_PRODUCTION_SMOKE=${STAGE15_8_RUN_PRODUCTION_SMOKE:-0}

OUTPUT_DIR=stage15_outputs
SUMMARY_DAT=$OUTPUT_DIR/fibre_stage15_8_controlled_structure_parallel_consistency.dat
OUT_DAT=$OUTPUT_DIR/stage15_8_controlled_structure_parallel_consistency.dat
REASONS_FILE=$OUTPUT_DIR/stage15_8_controlled_structure_parallel_consistency_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_8_controlled_structure_parallel_static_audit_report.txt
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_7_SCRIPT=stage15_checks/run_stage15_7_feedback_linkage.sh
STAGE15_7_DOC=stage15_checks/stage15_7_feedback_linkage.md
CHECK_DAT=$OUTPUT_DIR/fibre_stage15_7_feedback_linkage.dat

BUILD_STATUS=1
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
XCOMPACT_HOOK_STATUS=0
CONTROL_GUARD_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_7_STATUS=0
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
STAGE15_7_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOCS_AND_TARGET_STATUS=0
PRODUCTION_SMOKE_STATUS=0
PRODUCTION_SMOKE_DEFERRED_STATUS=0
NP1_RUN_STATUS=0
NP2_RUN_STATUS=0
NP4_RUN_STATUS=0
NP1_FINAL_STATUS=0
NP2_FINAL_STATUS=0
NP4_FINAL_STATUS=0
PARALLEL_VELOCITY_STATUS=0
PARALLEL_SLIP_STATUS=0
PARALLEL_FORCE_STATUS=0
PARALLEL_RHS_STATUS=0
CONTROLLED_UPDATE_COUNT_CONSISTENCY_STATUS=0
NO_RANK_CORRUPTION_STATUS=0
DIAGNOSTIC_STATUS=0
FINAL_STATUS=0

NP1_MAX_VELOCITY_UPDATE=0.0
NP2_MAX_VELOCITY_UPDATE=0.0
NP4_MAX_VELOCITY_UPDATE=0.0
NP1_MAX_SLIP_CHANGE=0.0
NP2_MAX_SLIP_CHANGE=0.0
NP4_MAX_SLIP_CHANGE=0.0
NP1_MAX_FEEDBACK_FORCE_CHANGE=0.0
NP2_MAX_FEEDBACK_FORCE_CHANGE=0.0
NP4_MAX_FEEDBACK_FORCE_CHANGE=0.0
NP1_RHS_RESPONSE=0.0
NP2_RHS_RESPONSE=0.0
NP4_RHS_RESPONSE=0.0
PARALLEL_VELOCITY_DIFF_NP2=0.0
PARALLEL_VELOCITY_DIFF_NP4=0.0
PARALLEL_SLIP_DIFF_NP2=0.0
PARALLEL_SLIP_DIFF_NP4=0.0
PARALLEL_FORCE_DIFF_NP2=0.0
PARALLEL_FORCE_DIFF_NP4=0.0
PARALLEL_RHS_DIFF_NP2=0.0
PARALLEL_RHS_DIFF_NP4=0.0
BENDING_SOLVE_COUNT=0
TENSION_SOLVE_COUNT=0
WALL_CONTACT_COUNT=0
MULTIFIBRE_COUNT=0
RHS_INJECTION_CONNECTION_COUNT=0
APPROVED_STAGE12_13_14_CHAIN_STATUS=0
NO_FLUID_RHS_MODIFICATION_STATUS=0
NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=0

mkdir -p "$OUTPUT_DIR" stage14_outputs stage13_outputs stage11_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$SUMMARY_DAT" "$OUT_DAT" "$CHECK_DAT" \
      "$OUTPUT_DIR"/stage15_8_np1_feedback_linkage.dat \
      "$OUTPUT_DIR"/stage15_8_np2_feedback_linkage.dat \
      "$OUTPUT_DIR"/stage15_8_np4_feedback_linkage.dat \
      "$OUTPUT_DIR"/stage15_8_np1_feedback_linkage.log \
      "$OUTPUT_DIR"/stage15_8_np2_feedback_linkage.log \
      "$OUTPUT_DIR"/stage15_8_np4_feedback_linkage.log

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

abs_diff() {
    awk -v a="$1" -v b="$2" 'BEGIN { d=(a+0.0)-(b+0.0); if (d<0.0) d=-d; printf "%.17e", d }'
}

le_value() {
    value=$1
    limit=$2
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value+0.0) <= (limit+0.0)) exit 0; exit 1 }'
}

compute_rhs_response() {
    file=$1
    force=$(get_value "$file" max_feedback_force_change 2>/dev/null || echo 0.0)
    lambda=$(get_value "$file" lambda_value 2>/dev/null || echo 0.0)
    awk -v f="$force" -v l="$lambda" 'BEGIN { if (l<0.0) l=-l; if (f<0.0) f=-f; printf "%.17e", l*f }'
}

scan_stage15_8_source_guards() {
    status=0
    files="stage15_checks/run_stage15_8_controlled_structure_parallel_consistency.sh stage15_checks/stage15_8_controlled_structure_parallel_consistency.md"
    static_note "Stage 15.8 source guard audit"

    if search_report 'rg[[:space:]]|grep[[:space:]]+-R' stage15_checks/run_stage15_8_controlled_structure_parallel_consistency.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_8_wrapper_must_not_require_ripgrep_or_recursive_grep"
        status=1
    fi
    if search_report '^[[:space:]]*call[[:space:]]+.*(bending|tension|wall|contact|multifibre|multi_fibre|poisson|projection|pressure|rk3|channel_forcing)' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_8_forbidden_production_call_marker_found"
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
    linkage_count=$(grep -Ec 'stage15_feedback_linkage|stage15_8_controlled_structure_parallel' src/xcompact3d.f90 || true)
    forbidden_count=$(grep -Ec 'call[[:space:]]+stage15_.*full.*advance|call[[:space:]]+.*bending|call[[:space:]]+.*tension' src/xcompact3d.f90 || true)
    printf 'xcompact3d hook counts: use=%s reset=%s register=%s apply=%s linkage=%s forbidden=%s\n' \
        "$use_count" "$reset_count" "$register_count" "$apply_count" "$linkage_count" "$forbidden_count" >> "$STATIC_REPORT"
    [ "$use_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_use_count_${use_count}"; status=1; }
    [ "$reset_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_reset_count_${reset_count}"; status=1; }
    [ "$register_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_register_count_${register_count}"; status=1; }
    [ "$apply_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_apply_count_${apply_count}"; status=1; }
    [ "$linkage_count" = "0" ] || { add_reason "xcompact3d_stage15_8_unapproved_linkage_connection_found"; status=1; }
    [ "$forbidden_count" = "0" ] || { add_reason "xcompact3d_forbidden_structure_solve_call_found"; status=1; }
    if [ "$status" = "0" ]; then
        XCOMPACT_HOOK_STATUS=1
        CONTROL_GUARD_STATUS=1
    fi
    return $status
}

run_static_audit() {
    status=0
    [ -f stage15_checks/run_stage15_8_controlled_structure_parallel_consistency.sh ] || { add_reason "missing_stage15_8_wrapper"; status=1; }
    [ -f stage15_checks/stage15_8_controlled_structure_parallel_consistency.md ] || { add_reason "missing_stage15_8_doc"; status=1; }
    if grep -q 'fibre_stage15_feedback_linkage_check' src/CMakeLists.txt && [ -f src/fibre_stage15_feedback_linkage_check.f90 ]; then
        DOCS_AND_TARGET_STATUS=1
    else
        add_reason "stage15_8_required_stage15_7_check_target_missing"
        status=1
    fi

    scan_stage15_8_source_guards || status=1
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

    if [ -f stage15_checks/run_stage15_6_controlled_structure_step_np1.sh ] && \
       search_silent 'fibre_stage15_6_controlled_structure_step_np1.dat' stage15_checks/run_stage15_6_controlled_structure_step_np1.sh; then
        STAGE15_6_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_6_controlled_structure_step_diagnostics_missing"
        status=1
    fi

    if [ -f "$STAGE15_7_SCRIPT" ] && [ -s "$STAGE15_7_DOC" ] && \
       search_silent 'fibre_stage15_7_feedback_linkage.dat' "$STAGE15_7_SCRIPT"; then
        STAGE15_7_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage15_7_feedback_linkage_diagnostics_missing"
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

    if [ "$STAGE15_8_REQUIRE_STAGE14_CLOSED" = "0" ]; then
        STAGE14_CLOSED_STATUS=1
    elif [ -s "$STAGE14_CLOSED_FILE" ]; then
        STAGE14_CLOSED_STATUS=1
    else
        add_reason "missing_stage14_checks_STAGE14_CLOSED_md"
        status=1
    fi

    if [ "$STAGE15_8_RUN_STAGE15_7" = "1" ]; then
        bash "$STAGE15_7_SCRIPT" || { add_reason "stage15_7_optional_prerequisite_run_failed"; status=1; }
    fi
    if [ "$STAGE15_8_REQUIRE_STAGE15_7" = "0" ]; then
        STAGE15_7_STATUS=1
    elif [ -s "$STAGE15_7_DOC" ] && [ -f "$STAGE15_7_SCRIPT" ] && [ -f src/fibre_stage15_feedback_linkage_check.f90 ]; then
        STAGE15_7_STATUS=1
    else
        add_reason "stage15_7_closed_pass_evidence_missing"
        status=1
    fi

    if [ "$STAGE15_8_NP_LIST" = "1 2 4" ]; then
        :
    else
        add_reason "stage15_8_np_list_must_be_1_2_4_got_${STAGE15_8_NP_LIST}"
        status=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

run_np_case() {
    np=$1
    exe=$(stage15_7_exe) || { add_reason "missing_fibre_stage15_feedback_linkage_check_executable"; return 1; }
    out_file=$OUTPUT_DIR/stage15_8_np${np}_feedback_linkage.dat
    log_file=$OUTPUT_DIR/stage15_8_np${np}_feedback_linkage.log
    rm -f "$CHECK_DAT" "$out_file" "$log_file"
    export STAGE15_7_ENABLE="$STAGE15_8_ENABLE"
    export STAGE15_7_CONTROLLED_STEP_ENABLE="$STAGE15_8_CONTROLLED_STEP_ENABLE"
    export STAGE15_7_STRUCTURE_ADVANCE_ENABLE="$STAGE15_8_STRUCTURE_ADVANCE_ENABLE"
    export STAGE15_7_DIAGNOSTIC_ONLY="$STAGE15_8_DIAGNOSTIC_ONLY"
    export STAGE15_7_NP=1
    export STAGE15_7_NPTS="$STAGE15_8_NPTS"
    export STAGE15_7_DT="$STAGE15_8_DT"
    export STAGE15_7_RHO_TILDE="$STAGE15_8_RHO_TILDE"
    export STAGE15_7_TEST_FORCE="$STAGE15_8_TEST_FORCE"
    export STAGE15_7_FEEDBACK_ALPHA="$STAGE15_8_FEEDBACK_ALPHA"
    export STAGE15_7_LAMBDA="$STAGE15_8_LAMBDA"
    export STAGE15_7_MAX_VELOCITY_UPDATE="$STAGE15_8_MAX_VELOCITY_UPDATE"
    export STAGE15_7_MAX_SLIP_ERROR="$STAGE15_8_MAX_SLIP_ERROR"
    export STAGE15_7_MAX_FORCE_ERROR="$STAGE15_8_MAX_FORCE_ERROR"
    export STAGE15_7_MAX_FORCE_RESPONSE="$STAGE15_8_MAX_FORCE_RESPONSE"
    export STAGE15_7_MAX_RHS_RESPONSE="$STAGE15_8_MAX_RHS_RESPONSE"
    # shellcheck disable=SC2086
    "$MPIEXEC" $MPIEXEC_FLAGS -n "$np" "$exe" > "$log_file" 2>&1 || {
        add_reason "stage15_8_np${np}_feedback_linkage_run_failed"
        return 1
    }
    [ -s "$CHECK_DAT" ] || { add_reason "stage15_8_np${np}_feedback_linkage_dat_missing"; return 1; }
    cp "$CHECK_DAT" "$out_file"
    verify_np_dat "$np" "$out_file" || return 1
    return 0
}

verify_np_dat() {
    np=$1
    file=$2
    status=0
    require_key_value "$file" stage15_7_requested_status 1 || status=1
    require_key_value "$file" controlled_step_enabled_status 1 || status=1
    require_key_value "$file" structure_advance_enable_status 1 || status=1
    require_key_value "$file" diagnostic_only_status 1 || status=1
    require_key_value "$file" velocity_update_nonzero_status 1 || status=1
    require_key_value "$file" slip_change_nonzero_status 1 || status=1
    require_key_value "$file" slip_consistency_status 1 || status=1
    require_key_value "$file" feedback_force_change_nonzero_status 1 || status=1
    require_key_value "$file" feedback_force_consistency_status 1 || status=1
    require_key_value "$file" force_response_bounded_status 1 || status=1
    require_key_value "$file" rhs_response_bounded_status 1 || status=1
    require_key_value "$file" controlled_update_count 1 || status=1
    require_key_value "$file" production_full_structure_advance_count 0 || status=1
    require_key_value "$file" bending_solve_count 0 || status=1
    require_key_value "$file" tension_solve_count 0 || status=1
    require_key_value "$file" wall_contact_count 0 || status=1
    require_key_value "$file" multifibre_count 0 || status=1
    require_key_value "$file" rhs_injection_connection_count 0 || status=1
    require_key_value "$file" approved_stage12_13_14_chain_status 1 || status=1
    require_key_value "$file" no_fluid_rhs_modification_status 1 || status=1
    require_key_value "$file" no_pressure_projection_modification_status 1 || status=1
    require_key_value "$file" no_poisson_modification_status 1 || status=1
    require_key_value "$file" no_rk3_channel_forcing_modification_status 1 || status=1
    require_key_value "$file" final_status 1 || status=1
    for key in max_velocity_update max_slip_change slip_error max_feedback_force_change feedback_force_error lambda_value; do
        require_finite_key "$file" "$key" || status=1
    done
    final_count=$(awk '$1 == "final_status" { count++ } END { print count+0 }' "$file")
    [ "$final_count" = "1" ] || { add_reason "stage15_8_np${np}_rank_corrupted_final_status_count_${final_count}"; status=1; }
    return $status
}

collect_np_values() {
    file1=$OUTPUT_DIR/stage15_8_np1_feedback_linkage.dat
    file2=$OUTPUT_DIR/stage15_8_np2_feedback_linkage.dat
    file4=$OUTPUT_DIR/stage15_8_np4_feedback_linkage.dat
    NP1_FINAL_STATUS=$(get_value "$file1" final_status 2>/dev/null || echo 0)
    NP2_FINAL_STATUS=$(get_value "$file2" final_status 2>/dev/null || echo 0)
    NP4_FINAL_STATUS=$(get_value "$file4" final_status 2>/dev/null || echo 0)
    NP1_MAX_VELOCITY_UPDATE=$(get_value "$file1" max_velocity_update 2>/dev/null || echo 0.0)
    NP2_MAX_VELOCITY_UPDATE=$(get_value "$file2" max_velocity_update 2>/dev/null || echo 0.0)
    NP4_MAX_VELOCITY_UPDATE=$(get_value "$file4" max_velocity_update 2>/dev/null || echo 0.0)
    NP1_MAX_SLIP_CHANGE=$(get_value "$file1" max_slip_change 2>/dev/null || echo 0.0)
    NP2_MAX_SLIP_CHANGE=$(get_value "$file2" max_slip_change 2>/dev/null || echo 0.0)
    NP4_MAX_SLIP_CHANGE=$(get_value "$file4" max_slip_change 2>/dev/null || echo 0.0)
    NP1_MAX_FEEDBACK_FORCE_CHANGE=$(get_value "$file1" max_feedback_force_change 2>/dev/null || echo 0.0)
    NP2_MAX_FEEDBACK_FORCE_CHANGE=$(get_value "$file2" max_feedback_force_change 2>/dev/null || echo 0.0)
    NP4_MAX_FEEDBACK_FORCE_CHANGE=$(get_value "$file4" max_feedback_force_change 2>/dev/null || echo 0.0)
    NP1_RHS_RESPONSE=$(compute_rhs_response "$file1")
    NP2_RHS_RESPONSE=$(compute_rhs_response "$file2")
    NP4_RHS_RESPONSE=$(compute_rhs_response "$file4")
    PARALLEL_VELOCITY_DIFF_NP2=$(abs_diff "$NP2_MAX_VELOCITY_UPDATE" "$NP1_MAX_VELOCITY_UPDATE")
    PARALLEL_VELOCITY_DIFF_NP4=$(abs_diff "$NP4_MAX_VELOCITY_UPDATE" "$NP1_MAX_VELOCITY_UPDATE")
    PARALLEL_SLIP_DIFF_NP2=$(abs_diff "$NP2_MAX_SLIP_CHANGE" "$NP1_MAX_SLIP_CHANGE")
    PARALLEL_SLIP_DIFF_NP4=$(abs_diff "$NP4_MAX_SLIP_CHANGE" "$NP1_MAX_SLIP_CHANGE")
    PARALLEL_FORCE_DIFF_NP2=$(abs_diff "$NP2_MAX_FEEDBACK_FORCE_CHANGE" "$NP1_MAX_FEEDBACK_FORCE_CHANGE")
    PARALLEL_FORCE_DIFF_NP4=$(abs_diff "$NP4_MAX_FEEDBACK_FORCE_CHANGE" "$NP1_MAX_FEEDBACK_FORCE_CHANGE")
    PARALLEL_RHS_DIFF_NP2=$(abs_diff "$NP2_RHS_RESPONSE" "$NP1_RHS_RESPONSE")
    PARALLEL_RHS_DIFF_NP4=$(abs_diff "$NP4_RHS_RESPONSE" "$NP1_RHS_RESPONSE")
    BENDING_SOLVE_COUNT=$(get_value "$file1" bending_solve_count 2>/dev/null || echo 1)
    TENSION_SOLVE_COUNT=$(get_value "$file1" tension_solve_count 2>/dev/null || echo 1)
    WALL_CONTACT_COUNT=$(get_value "$file1" wall_contact_count 2>/dev/null || echo 1)
    MULTIFIBRE_COUNT=$(get_value "$file1" multifibre_count 2>/dev/null || echo 1)
    RHS_INJECTION_CONNECTION_COUNT=$(get_value "$file1" rhs_injection_connection_count 2>/dev/null || echo 1)
    APPROVED_STAGE12_13_14_CHAIN_STATUS=$(get_value "$file1" approved_stage12_13_14_chain_status 2>/dev/null || echo 0)
    NO_FLUID_RHS_MODIFICATION_STATUS=$(get_value "$file1" no_fluid_rhs_modification_status 2>/dev/null || echo 0)
    NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=$(get_value "$file1" no_pressure_projection_modification_status 2>/dev/null || echo 0)
    NO_POISSON_MODIFICATION_STATUS=$(get_value "$file1" no_poisson_modification_status 2>/dev/null || echo 0)
    NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value "$file1" no_rk3_channel_forcing_modification_status 2>/dev/null || echo 0)
}

compare_parallel_values() {
    status=0
    collect_np_values
    le_value "$PARALLEL_VELOCITY_DIFF_NP2" "$STAGE15_8_MAX_PARALLEL_VELOCITY_DIFF" && \
    le_value "$PARALLEL_VELOCITY_DIFF_NP4" "$STAGE15_8_MAX_PARALLEL_VELOCITY_DIFF" || { add_reason "stage15_8_parallel_velocity_diff_exceeded"; status=1; }
    le_value "$PARALLEL_SLIP_DIFF_NP2" "$STAGE15_8_MAX_PARALLEL_SLIP_DIFF" && \
    le_value "$PARALLEL_SLIP_DIFF_NP4" "$STAGE15_8_MAX_PARALLEL_SLIP_DIFF" || { add_reason "stage15_8_parallel_slip_diff_exceeded"; status=1; }
    le_value "$PARALLEL_FORCE_DIFF_NP2" "$STAGE15_8_MAX_PARALLEL_FORCE_DIFF" && \
    le_value "$PARALLEL_FORCE_DIFF_NP4" "$STAGE15_8_MAX_PARALLEL_FORCE_DIFF" || { add_reason "stage15_8_parallel_force_diff_exceeded"; status=1; }
    le_value "$PARALLEL_RHS_DIFF_NP2" "$STAGE15_8_MAX_PARALLEL_RHS_DIFF" && \
    le_value "$PARALLEL_RHS_DIFF_NP4" "$STAGE15_8_MAX_PARALLEL_RHS_DIFF" || { add_reason "stage15_8_parallel_rhs_diff_exceeded"; status=1; }

    if [ "$status" = "0" ]; then
        PARALLEL_VELOCITY_STATUS=1
        PARALLEL_SLIP_STATUS=1
        PARALLEL_FORCE_STATUS=1
        PARALLEL_RHS_STATUS=1
    fi

    c1=$(get_value "$OUTPUT_DIR/stage15_8_np1_feedback_linkage.dat" controlled_update_count 2>/dev/null || echo -1)
    c2=$(get_value "$OUTPUT_DIR/stage15_8_np2_feedback_linkage.dat" controlled_update_count 2>/dev/null || echo -2)
    c4=$(get_value "$OUTPUT_DIR/stage15_8_np4_feedback_linkage.dat" controlled_update_count 2>/dev/null || echo -4)
    if [ "$c1" = "1" ] && [ "$c2" = "1" ] && [ "$c4" = "1" ]; then
        CONTROLLED_UPDATE_COUNT_CONSISTENCY_STATUS=1
    else
        add_reason "stage15_8_controlled_update_counts_not_consistent_${c1}_${c2}_${c4}"
        status=1
    fi

    if [ "$NP1_FINAL_STATUS" = "1" ] && [ "$NP2_FINAL_STATUS" = "1" ] && [ "$NP4_FINAL_STATUS" = "1" ]; then
        NO_RANK_CORRUPTION_STATUS=1
    else
        add_reason "stage15_8_np_final_status_not_all_pass"
        status=1
    fi

    if [ "$status" = "0" ]; then
        DIAGNOSTIC_STATUS=1
    fi
    return $status
}

handle_production_smoke() {
    if [ "$STAGE15_8_RUN_PRODUCTION_SMOKE" = "0" ]; then
        PRODUCTION_SMOKE_DEFERRED_STATUS=1
        PRODUCTION_SMOKE_STATUS=1
        echo "production_smoke_deferred: MPI-launched standalone feedback-linkage validation is authoritative for Stage 15.8" >> "$STATIC_REPORT"
    else
        PRODUCTION_SMOKE_STATUS=0
        add_reason "stage15_8_production_smoke_requested_but_not_enabled_for_parallel_linkage"
    fi
}

write_summary_dat() {
    final_status=$1
    cat > "$SUMMARY_DAT" <<EOF_DAT
stage15_8_requested_status 1
np_list 1_2_4
np1_run_status $NP1_RUN_STATUS
np2_run_status $NP2_RUN_STATUS
np4_run_status $NP4_RUN_STATUS
np1_final_status $NP1_FINAL_STATUS
np2_final_status $NP2_FINAL_STATUS
np4_final_status $NP4_FINAL_STATUS
np1_max_velocity_update $NP1_MAX_VELOCITY_UPDATE
np2_max_velocity_update $NP2_MAX_VELOCITY_UPDATE
np4_max_velocity_update $NP4_MAX_VELOCITY_UPDATE
np1_max_slip_change $NP1_MAX_SLIP_CHANGE
np2_max_slip_change $NP2_MAX_SLIP_CHANGE
np4_max_slip_change $NP4_MAX_SLIP_CHANGE
np1_max_feedback_force_change $NP1_MAX_FEEDBACK_FORCE_CHANGE
np2_max_feedback_force_change $NP2_MAX_FEEDBACK_FORCE_CHANGE
np4_max_feedback_force_change $NP4_MAX_FEEDBACK_FORCE_CHANGE
parallel_velocity_diff_np2 $PARALLEL_VELOCITY_DIFF_NP2
parallel_velocity_diff_np4 $PARALLEL_VELOCITY_DIFF_NP4
parallel_slip_diff_np2 $PARALLEL_SLIP_DIFF_NP2
parallel_slip_diff_np4 $PARALLEL_SLIP_DIFF_NP4
parallel_force_diff_np2 $PARALLEL_FORCE_DIFF_NP2
parallel_force_diff_np4 $PARALLEL_FORCE_DIFF_NP4
parallel_rhs_diff_np2 $PARALLEL_RHS_DIFF_NP2
parallel_rhs_diff_np4 $PARALLEL_RHS_DIFF_NP4
parallel_velocity_status $PARALLEL_VELOCITY_STATUS
parallel_slip_status $PARALLEL_SLIP_STATUS
parallel_force_status $PARALLEL_FORCE_STATUS
parallel_rhs_status $PARALLEL_RHS_STATUS
controlled_update_count_consistency_status $CONTROLLED_UPDATE_COUNT_CONSISTENCY_STATUS
no_rank_corruption_status $NO_RANK_CORRUPTION_STATUS
bending_solve_count $BENDING_SOLVE_COUNT
tension_solve_count $TENSION_SOLVE_COUNT
wall_contact_count $WALL_CONTACT_COUNT
multifibre_count $MULTIFIBRE_COUNT
rhs_injection_connection_count $RHS_INJECTION_CONNECTION_COUNT
approved_stage12_13_14_chain_status $APPROVED_STAGE12_13_14_CHAIN_STATUS
no_fluid_rhs_modification_status $NO_FLUID_RHS_MODIFICATION_STATUS
no_pressure_projection_modification_status $NO_PRESSURE_PROJECTION_MODIFICATION_STATUS
no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
no_rk3_channel_forcing_modification_status $NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS
final_status $final_status
EOF_DAT
}

write_gate_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage15_8_requested_status 1
stage15_8_build_status $BUILD_STATUS
stage15_8_static_audit_status $STATIC_AUDIT_STATUS
stage15_8_source_guard_status $SOURCE_GUARD_STATUS
stage15_8_xcompact_hook_status $XCOMPACT_HOOK_STATUS
stage15_8_control_guard_status $CONTROL_GUARD_STATUS
stage15_8_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_8_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_8_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_8_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_8_stage15_1_diagnostic_status $STAGE15_1_DIAGNOSTIC_STATUS
stage15_8_stage15_2_diagnostic_status $STAGE15_2_DIAGNOSTIC_STATUS
stage15_8_stage15_3_diagnostic_status $STAGE15_3_DIAGNOSTIC_STATUS
stage15_8_stage15_4_diagnostic_status $STAGE15_4_DIAGNOSTIC_STATUS
stage15_8_stage15_5_diagnostic_status $STAGE15_5_DIAGNOSTIC_STATUS
stage15_8_stage15_6_diagnostic_status $STAGE15_6_DIAGNOSTIC_STATUS
stage15_8_stage15_7_diagnostic_status $STAGE15_7_DIAGNOSTIC_STATUS
stage15_8_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_8_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_8_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_8_stage15_7_status $STAGE15_7_STATUS
stage15_8_docs_and_target_status $DOCS_AND_TARGET_STATUS
stage15_8_diagnostic_status $DIAGNOSTIC_STATUS
stage15_8_production_smoke_status $PRODUCTION_SMOKE_STATUS
stage15_8_production_smoke_deferred_status $PRODUCTION_SMOKE_DEFERRED_STATUS
stage15_8_final_status $final_status
EOF_DAT
}

run_static_audit || true
handle_production_smoke

if ! build_target fibre_stage15_feedback_linkage_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage15_feedback_linkage_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_np_case 1; then NP1_RUN_STATUS=1; fi
    if run_np_case 2; then NP2_RUN_STATUS=1; fi
    if run_np_case 4; then NP4_RUN_STATUS=1; fi
else
    add_reason "stage15_8_parallel_runs_skipped_due_to_build_failure"
fi

if [ "$NP1_RUN_STATUS" = "1" ] && [ "$NP2_RUN_STATUS" = "1" ] && [ "$NP4_RUN_STATUS" = "1" ]; then
    compare_parallel_values || true
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$SOURCE_GUARD_STATUS" = "1" ] && \
   [ "$XCOMPACT_HOOK_STATUS" = "1" ] && [ "$CONTROL_GUARD_STATUS" = "1" ] && \
   [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && [ "$STAGE11_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE14_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_1_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_2_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_3_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_4_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_5_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_6_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_7_DIAGNOSTIC_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$STAGE15_7_STATUS" = "1" ] && [ "$DOCS_AND_TARGET_STATUS" = "1" ] && \
   [ "$NP1_RUN_STATUS" = "1" ] && [ "$NP2_RUN_STATUS" = "1" ] && [ "$NP4_RUN_STATUS" = "1" ] && \
   [ "$DIAGNOSTIC_STATUS" = "1" ] && [ "$PRODUCTION_SMOKE_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_summary_dat 1
    write_gate_dat 1
    echo 'STAGE 15.8 CONTROLLED STRUCTURE PARALLEL CONSISTENCY VERDICT: PASS'
    echo 'STAGE 15.8 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_summary_dat 0
write_gate_dat 0
echo 'STAGE 15.8 CONTROLLED STRUCTURE PARALLEL CONSISTENCY VERDICT: FAIL'
echo 'STAGE 15.8 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_8_failure'
fi
exit 1
