#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE16_0_REQUIRE_STAGE14_CLOSED=${STAGE16_0_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_0_REQUIRE_STAGE15_CLOSED=${STAGE16_0_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE=${STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE:-1}
STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING=${STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING:-0}
STAGE16_0_ENABLE=${STAGE16_0_ENABLE:-1}
STAGE16_0_DIAGNOSTIC_ONLY=${STAGE16_0_DIAGNOSTIC_ONLY:-1}

OUTPUT_DIR=stage16_outputs
SUMMARY_DAT=$OUTPUT_DIR/fibre_stage16_0_preflight_closure_integrity.dat
OUT_DAT=$OUTPUT_DIR/stage16_0_preflight_closure_integrity.dat
REASONS_FILE=$OUTPUT_DIR/stage16_0_preflight_closure_integrity_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage16_0_preflight_closure_integrity_static_audit_report.txt
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_CLOSED_FILE=stage15_checks/STAGE15_CLOSED.md
STAGE15_11_SCRIPT=stage15_checks/run_stage15_11_total_smoke_closure.sh

STAGE14_CLOSED_FILE_STATUS=0
STAGE15_CLOSED_FILE_STATUS=0
STAGE15_CLOSED_CONTENT_STATUS=0
STAGE14_REGRESSION_STATUS=0
STAGE15_REGRESSION_STATUS=0
STAGE13_6_DIAGNOSTIC_NAME_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
RANK0_SAFE_DIAGNOSTIC_STATUS=0
NO_RG_ONLY_DEPENDENCY_STATUS=0
NO_UNKNOWN_FAILURE_STATUS=0
STAGE16_BOUNDARY_STATUS=0
APPROVED_STAGE12_13_14_CHAIN_STATUS=0
NO_DIRECT_RHS_INJECTION_STATUS=0
NO_LEGACY_IBM_FORCING_STATUS=0
NO_BENDING_SOLVE_STATUS=0
NO_TENSION_SOLVE_STATUS=0
NO_WALL_CONTACT_STATUS=0
NO_MULTIFIBRE_STATUS=0
NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=0
FINAL_STATUS=0

mkdir -p "$OUTPUT_DIR"
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$SUMMARY_DAT" "$OUT_DAT"

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

maybe_generate_stage15_closure() {
    if [ -s "$STAGE15_CLOSED_FILE" ]; then
        return 0
    fi
    if [ "$STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING" != "1" ]; then
        return 1
    fi
    [ -f "$STAGE15_11_SCRIPT" ] || { add_reason "missing_stage15_11_closure_script_for_auto_run"; return 1; }
    BUILD_DIR="$BUILD_DIR" DECOMP2D_ROOT="$DECOMP2D_ROOT" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
        bash "$STAGE15_11_SCRIPT" > "$OUTPUT_DIR/stage16_0_stage15_11_auto_run.log" 2>&1 || {
            add_reason "stage15_11_auto_run_failed"
            return 1
        }
    grep 'STAGE15_CLOSED=YES' "$OUTPUT_DIR/stage16_0_stage15_11_auto_run.log" >/dev/null 2>&1 || {
        add_reason "stage15_11_auto_run_missing_closed_verdict"
        return 1
    }
    [ -s "$STAGE15_CLOSED_FILE" ]
}

check_closure_files() {
    status=0
    if [ "$STAGE16_0_REQUIRE_STAGE14_CLOSED" = "1" ]; then
        [ -s "$STAGE14_CLOSED_FILE" ] && STAGE14_CLOSED_FILE_STATUS=1 || { add_reason "missing_or_empty_stage14_closed_file"; status=1; }
    else
        STAGE14_CLOSED_FILE_STATUS=1
    fi
    if [ "$STAGE16_0_REQUIRE_STAGE15_CLOSED" = "1" ]; then
        maybe_generate_stage15_closure || true
        [ -s "$STAGE15_CLOSED_FILE" ] && STAGE15_CLOSED_FILE_STATUS=1 || { add_reason "missing_or_empty_stage15_closed_file"; status=1; }
    else
        STAGE15_CLOSED_FILE_STATUS=1
    fi
    if [ -s "$STAGE15_CLOSED_FILE" ]; then
        content_status=0
        grep 'STAGE15_CLOSED=YES' "$STAGE15_CLOSED_FILE" >/dev/null 2>&1 || { add_reason "stage15_closed_file_missing_yes_marker"; content_status=1; }
        grep -E 'Stage 15\.11 passed|Stage 15 is closed|Stage 15.*closed' "$STAGE15_CLOSED_FILE" >/dev/null 2>&1 || { add_reason "stage15_closed_file_missing_closed_statement"; content_status=1; }
        grep -Ei 'bending.*tension.*wall/contact.*multi-fibre|bending.*tension.*wall.*multi' "$STAGE15_CLOSED_FILE" >/dev/null 2>&1 || { add_reason "stage15_closed_file_missing_inactive_physics_statement"; content_status=1; }
        grep -E 'Stage 15[[:space:]]*[-=]+>[[:space:]]*Stage 12[[:space:]]*[-=]+>[[:space:]]*Stage 13[[:space:]]*[-=]+>[[:space:]]*Stage 14' "$STAGE15_CLOSED_FILE" >/dev/null 2>&1 || { add_reason "stage15_closed_file_missing_approved_chain_statement"; content_status=1; }
        grep -E 'Stage 16 first controlled one-flexible-fibre channel DNS FSI' "$STAGE15_CLOSED_FILE" >/dev/null 2>&1 || { add_reason "stage15_closed_file_missing_stage16_recommendation"; content_status=1; }
        if [ "$content_status" = "0" ]; then
            STAGE15_CLOSED_CONTENT_STATUS=1
            APPROVED_STAGE12_13_14_CHAIN_STATUS=1
        else
            status=1
        fi
    fi
    return $status
}

check_stage14_regression() {
    status=0
    if search_report 'stage14_get_injection_gain\(\)[[:space:]]*==[[:space:]]*0\.0' src/*.f90 stage14_checks/*.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage14_forbidden_lambda_zero_registration_gate_found"
        status=1
    fi
    search_silent 'stage11_5_production_oneway_hook' src/*.f90 stage11_checks/*.sh || { add_reason "missing_stage11_5_production_diagnostics"; status=1; }
    search_silent 'stage13_6_production_force_density_candidate' src/*.f90 stage13_checks/*.sh || { add_reason "missing_stage13_6_production_force_density_diagnostics"; status=1; }
    if search_silent 'stage13_5_production_force_density_candidate' src/*.f90 stage13_checks/*.sh; then
        add_reason "old_stage13_5_force_density_name_reappeared"
        status=1
    fi
    search_silent 'stage14_5_production_rhs_hook' src/*.f90 stage14_checks/*.sh || { add_reason "missing_stage14_5_production_rhs_diagnostics"; status=1; }
    if [ "$status" = "0" ]; then
        STAGE14_REGRESSION_STATUS=1
        STAGE13_6_DIAGNOSTIC_NAME_STATUS=1
    fi
    return $status
}

check_stage15_regression() {
    status=0
    for n in 1 2 3 4 5 6 7 8 9 10 11; do
        case "$n" in
            1) stem=stage15_1_structure_state_buffer ;;
            2) stem=stage15_2_velocity_source_adapter ;;
            3) stem=stage15_3_structure_advance_formula ;;
            4) stem=stage15_4_production_structure_hook ;;
            5) stem=stage15_5_structure_noop_invariance ;;
            6) stem=stage15_6_controlled_structure_step_np1 ;;
            7) stem=stage15_7_feedback_linkage ;;
            8) stem=stage15_8_controlled_structure_parallel_consistency ;;
            9) stem=stage15_9_io_restart_structure_state ;;
            10) stem=stage15_10_rhs_ibm_structure_contamination_audit ;;
            11) stem=stage15_11_total_smoke_closure ;;
        esac
        search_silent "$stem" stage15_checks/*.sh stage15_checks/*.md || { add_reason "missing_${stem}_script_or_doc_marker"; status=1; }
    done
    if [ "$STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE" = "1" ]; then
        [ "$STAGE15_CLOSED_CONTENT_STATUS" = "1" ] || { add_reason "closed_stage15_evidence_not_accepted_without_valid_closure_file"; status=1; }
    else
        for dat in stage15_outputs/stage15_11_total_smoke_closure.dat stage15_outputs/fibre_stage15_11_total_smoke_closure.dat; do
            [ -s "$dat" ] && grep -E 'final_status[[:space:]]+1' "$dat" >/dev/null 2>&1 && break
        done || { add_reason "stage15_11_runtime_evidence_missing_and_closed_evidence_not_accepted"; status=1; }
    fi
    if grep -Eq "echo ' - unknown_stage15_11_failure'" "$STAGE15_11_SCRIPT" 2>/dev/null && ! grep -Eq 'Reasons:|REASONS_FILE|add_reason' "$STAGE15_11_SCRIPT" 2>/dev/null; then
        add_reason "stage15_11_unknown_failure_without_reasons_block"
        status=1
    fi
    if [ "$status" = "0" ]; then
        STAGE15_REGRESSION_STATUS=1
        NO_UNKNOWN_FAILURE_STATUS=1
    fi
    return $status
}

check_rank0_sampling_and_rg() {
    status=0
    if search_silent 'rank[[:space:]]*==[[:space:]]*0|nrank[[:space:]]*==[[:space:]]*0|myid[[:space:]]*==[[:space:]]*0|rank0' src/*.f90 stage11_checks/*.sh stage13_checks/*.sh stage14_checks/*.sh stage15_checks/*.sh; then
        RANK0_SAFE_DIAGNOSTIC_STATUS=1
    else
        add_reason "missing_rank0_safe_diagnostic_markers"
        status=1
    fi
    if search_silent 'lbound\(ux, 1\) \+ 2|np=1/2/4|np=4' src/*.f90 stage13_checks/*.sh; then
        STAGE13_SAMPLING_REPAIR_STATUS=1
    else
        add_reason "missing_stage13_sampling_repair_evidence"
        status=1
    fi
    if search_report '^[[:space:]]*(rg|grep[[:space:]]+-R)' stage16_checks/run_stage16_0_preflight_closure_integrity.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage16_0_wrapper_must_not_require_ripgrep_or_recursive_grep"
        status=1
    fi
    if [ "$STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE" != "1" ] && search_report '^[[:space:]]*rg[[:space:]]' stage15_checks/*.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_script_rg_only_dependency_detected_when_closed_evidence_not_accepted"
        status=1
    fi
    [ "$status" = "0" ] && NO_RG_ONLY_DEPENDENCY_STATUS=1
    return $status
}

check_stage16_boundary() {
    status=0
    if ls src/*stage16*.f90 >/dev/null 2>&1; then
        add_reason "stage16_source_file_introduced_in_preflight"
        ls src/*stage16*.f90 >> "$STATIC_REPORT" 2>/dev/null || true
        status=1
    fi
    if grep -Eq 'stage16|fibre_stage16' src/xcompact3d.f90; then
        add_reason "stage16_production_hook_or_call_found_in_xcompact3d"
        status=1
    fi
    if search_report '^[[:space:]]*(use|call)[[:space:]]+.*(rhs|inject|poisson|projection|pressure|rk3|channel_forcing|ibm|bending|tension|wall|contact|multifibre|multi_fibre|structure_advance)' stage16_checks/*.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage16_forbidden_physics_call_marker_found"
        status=1
    fi
    if search_report '(>|rm[[:space:]]+-f|cat[[:space:]]*>|sed[[:space:]]+-i).*(stage1[0-5]_checks|stage1[0-5]_outputs|src/)' stage16_checks/*.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage16_script_appears_to_modify_closed_stage_or_src_file"
        status=1
    fi
    if [ "$status" = "0" ]; then
        STAGE16_BOUNDARY_STATUS=1
        NO_DIRECT_RHS_INJECTION_STATUS=1
        NO_LEGACY_IBM_FORCING_STATUS=1
        NO_BENDING_SOLVE_STATUS=1
        NO_TENSION_SOLVE_STATUS=1
        NO_WALL_CONTACT_STATUS=1
        NO_MULTIFIBRE_STATUS=1
        NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=1
        NO_POISSON_MODIFICATION_STATUS=1
        NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=1
    fi
    return $status
}

write_output_dat() {
    final_status=$1
    cat > "$SUMMARY_DAT" <<EOF_DAT
stage16_0_requested_status 1
stage16_0_enable $STAGE16_0_ENABLE
stage16_0_diagnostic_only $STAGE16_0_DIAGNOSTIC_ONLY
stage14_closed_file_status $STAGE14_CLOSED_FILE_STATUS
stage15_closed_file_status $STAGE15_CLOSED_FILE_STATUS
stage15_closed_content_status $STAGE15_CLOSED_CONTENT_STATUS
stage14_regression_status $STAGE14_REGRESSION_STATUS
stage15_regression_status $STAGE15_REGRESSION_STATUS
stage13_6_diagnostic_name_status $STAGE13_6_DIAGNOSTIC_NAME_STATUS
stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
rank0_safe_diagnostic_status $RANK0_SAFE_DIAGNOSTIC_STATUS
no_rg_only_dependency_status $NO_RG_ONLY_DEPENDENCY_STATUS
no_unknown_failure_status $NO_UNKNOWN_FAILURE_STATUS
stage16_boundary_status $STAGE16_BOUNDARY_STATUS
approved_stage12_13_14_chain_status $APPROVED_STAGE12_13_14_CHAIN_STATUS
no_direct_rhs_injection_status $NO_DIRECT_RHS_INJECTION_STATUS
no_legacy_ibm_forcing_status $NO_LEGACY_IBM_FORCING_STATUS
no_bending_solve_status $NO_BENDING_SOLVE_STATUS
no_tension_solve_status $NO_TENSION_SOLVE_STATUS
no_wall_contact_status $NO_WALL_CONTACT_STATUS
no_multifibre_status $NO_MULTIFIBRE_STATUS
no_pressure_projection_modification_status $NO_PRESSURE_PROJECTION_MODIFICATION_STATUS
no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
no_rk3_channel_forcing_modification_status $NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS
final_status $final_status
EOF_DAT
    cp "$SUMMARY_DAT" "$OUT_DAT"
}

ensure_build_dir || add_reason "stage16_0_ensure_build_dir_failed"
check_closure_files || true
check_stage14_regression || true
check_stage15_regression || true
check_rank0_sampling_and_rg || true
check_stage16_boundary || true

if [ "$STAGE16_0_ENABLE" = "1" ] && [ "$STAGE16_0_DIAGNOSTIC_ONLY" = "1" ] && \
   [ "$STAGE14_CLOSED_FILE_STATUS" = "1" ] && [ "$STAGE15_CLOSED_FILE_STATUS" = "1" ] && \
   [ "$STAGE15_CLOSED_CONTENT_STATUS" = "1" ] && [ "$STAGE14_REGRESSION_STATUS" = "1" ] && \
   [ "$STAGE15_REGRESSION_STATUS" = "1" ] && [ "$STAGE13_6_DIAGNOSTIC_NAME_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$RANK0_SAFE_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$NO_RG_ONLY_DEPENDENCY_STATUS" = "1" ] && [ "$NO_UNKNOWN_FAILURE_STATUS" = "1" ] && \
   [ "$STAGE16_BOUNDARY_STATUS" = "1" ] && [ "$APPROVED_STAGE12_13_14_CHAIN_STATUS" = "1" ] && \
   [ "$NO_DIRECT_RHS_INJECTION_STATUS" = "1" ] && [ "$NO_LEGACY_IBM_FORCING_STATUS" = "1" ] && \
   [ "$NO_BENDING_SOLVE_STATUS" = "1" ] && [ "$NO_TENSION_SOLVE_STATUS" = "1" ] && \
   [ "$NO_WALL_CONTACT_STATUS" = "1" ] && [ "$NO_MULTIFIBRE_STATUS" = "1" ] && \
   [ "$NO_PRESSURE_PROJECTION_MODIFICATION_STATUS" = "1" ] && [ "$NO_POISSON_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS" = "1" ]; then
    write_output_dat 1
    echo 'STAGE 16.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: PASS'
    echo 'STAGE 16.0 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 16.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: FAIL'
echo 'STAGE 16.0 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - stage16_0_failed_without_recorded_reason_check_wrapper_logic'
fi
exit 1
