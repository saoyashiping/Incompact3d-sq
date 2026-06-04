#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_10_RUN_STAGE15_9=${STAGE15_10_RUN_STAGE15_9:-0}
STAGE15_10_REQUIRE_STAGE14_CLOSED=${STAGE15_10_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_10_REQUIRE_STAGE15_9=${STAGE15_10_REQUIRE_STAGE15_9:-1}
STAGE15_10_ENABLE=${STAGE15_10_ENABLE:-1}
STAGE15_10_CONTROLLED_STEP_ENABLE=${STAGE15_10_CONTROLLED_STEP_ENABLE:-1}
STAGE15_10_STRUCTURE_ADVANCE_ENABLE=${STAGE15_10_STRUCTURE_ADVANCE_ENABLE:-1}
STAGE15_10_DIAGNOSTIC_ONLY=${STAGE15_10_DIAGNOSTIC_ONLY:-1}
STAGE15_10_NP=${STAGE15_10_NP:-2}
STAGE15_10_NPTS=${STAGE15_10_NPTS:-8}
STAGE15_10_DT=${STAGE15_10_DT:-1.0e-4}
STAGE15_10_RHO_TILDE=${STAGE15_10_RHO_TILDE:-1.0}
STAGE15_10_TEST_FORCE=${STAGE15_10_TEST_FORCE:-1.0e-6}
STAGE15_10_FEEDBACK_ALPHA=${STAGE15_10_FEEDBACK_ALPHA:-1.0}
STAGE15_10_LAMBDA=${STAGE15_10_LAMBDA:-1.0e-8}
STAGE15_10_MAX_FORCE_RESPONSE=${STAGE15_10_MAX_FORCE_RESPONSE:-1.0e-8}
STAGE15_10_MAX_RHS_RESPONSE=${STAGE15_10_MAX_RHS_RESPONSE:-1.0e-12}
STAGE15_10_MAX_STAGE14_RHS_INCREMENT=${STAGE15_10_MAX_STAGE14_RHS_INCREMENT:-1.0e-4}
STAGE15_10_MAX_FLUID_DELTA=${STAGE15_10_MAX_FLUID_DELTA:-1.0e-8}
STAGE15_10_RUN_PRODUCTION_SMOKE=${STAGE15_10_RUN_PRODUCTION_SMOKE:-1}

OUTPUT_DIR=stage15_outputs
SUMMARY_DAT=$OUTPUT_DIR/fibre_stage15_10_rhs_ibm_structure_contamination_audit.dat
OUT_DAT=$OUTPUT_DIR/stage15_10_rhs_ibm_structure_contamination_audit.dat
REASONS_FILE=$OUTPUT_DIR/stage15_10_rhs_ibm_structure_contamination_audit_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_10_rhs_ibm_structure_contamination_static_audit_report.txt
STAGE15_7_DAT=$OUTPUT_DIR/fibre_stage15_7_feedback_linkage.dat
STAGE15_10_LINKAGE_DAT=$OUTPUT_DIR/stage15_10_feedback_linkage.dat
STAGE14_10_DAT=stage14_outputs/stage14_10_rhs_ibm_structure_contamination_audit.dat
STAGE15_10_STAGE14_DAT=$OUTPUT_DIR/stage15_10_stage14_rhs_ibm_structure_contamination_audit.dat
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md

BUILD_STATUS=1
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
XCOMPACT_HOOK_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_9_STATUS=0
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
STAGE15_8_DIAGNOSTIC_STATUS=0
STAGE15_9_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOCS_AND_TARGET_STATUS=0
PRODUCTION_SMOKE_STATUS=0
PRODUCTION_SMOKE_DEFERRED_STATUS=0

CONTROLLED_UPDATE_STATUS=0
FEEDBACK_LINKAGE_STATUS=0
STAGE13_FORCE_DENSITY_STATUS=0
STAGE14_RHS_STATUS=0
FORCE_RESPONSE_BOUNDED_STATUS=0
RHS_RESPONSE_BOUNDED_STATUS=0
STAGE14_RHS_INCREMENT_BOUNDED_STATUS=0
FLUID_SIGNATURE_STATUS=0
APPROVED_STAGE12_13_14_CHAIN_STATUS=0
NO_FLUID_RHS_MODIFICATION_STATUS=0
NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=0
NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
NO_NAN_INF_STATUS=0
FINAL_STATUS=0

FLUID_SIGNATURE_DELTA=0.0
DIRECT_RHS_INJECTION_CONNECTION_COUNT=0
LEGACY_IBM_FORCING_COUNT=0
PRODUCTION_IBM_FORCING_OUTSIDE_APPROVED_CHAIN_COUNT=0
PRODUCTION_FULL_STRUCTURE_ADVANCE_COUNT=0
CONTROLLED_STRUCTURE_UPDATE_COUNT=0
BENDING_SOLVE_COUNT=0
TENSION_SOLVE_COUNT=0
IMPLICIT_STRUCTURE_SOLVE_COUNT=0
WALL_CONTACT_COUNT=0
MULTIFIBRE_COUNT=0

mkdir -p "$OUTPUT_DIR" stage14_outputs stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$SUMMARY_DAT" "$OUT_DAT" "$STAGE15_10_LINKAGE_DAT" "$STAGE15_10_STAGE14_DAT"

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

is_finite_number() {
    value=$1
    awk -v value="$value" 'BEGIN { v=value+0.0; if (v==v && v<1.0e300 && v>-1.0e300) exit 0; exit 1 }'
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
    is_finite_number "$value" || { add_reason "nonfinite_${key}_${value}_in_${file}"; return 1; }
}

require_real_le_key() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || { add_reason "missing_${key}_in_${file}"; return 1; }
    is_finite_number "$value" || { add_reason "nonfinite_${key}_${value}_in_${file}"; return 1; }
    awk -v value="$value" -v limit="$limit" 'BEGIN { exit !((value+0.0) <= (limit+0.0)) }' || {
        add_reason "${key}_expected_le_${limit}_got_${value}_in_${file}"
        return 1
    }
}

validate_controls() {
    status=0
    for value in "$STAGE15_10_DT" "$STAGE15_10_RHO_TILDE" "$STAGE15_10_TEST_FORCE" \
                 "$STAGE15_10_FEEDBACK_ALPHA" "$STAGE15_10_LAMBDA" \
                 "$STAGE15_10_MAX_FORCE_RESPONSE" "$STAGE15_10_MAX_RHS_RESPONSE" \
                 "$STAGE15_10_MAX_STAGE14_RHS_INCREMENT" "$STAGE15_10_MAX_FLUID_DELTA"; do
        is_finite_number "$value" || status=1
    done
    [ "$status" = "0" ] || add_reason "stage15_10_numeric_control_not_finite"
    case "$STAGE15_10_NP" in 1|2|4) ;; *) add_reason "stage15_10_np_must_be_1_2_or_4"; status=1 ;; esac
    return $status
}

scan_active_stage15_forbidden() {
    status=0
    files="src/fibre_stage15_structure_state.f90 src/fibre_stage15_velocity_source_adapter.f90 src/fibre_stage15_structure_advance_formula.f90 src/fibre_stage15_controlled_structure_step.f90 src/fibre_stage15_feedback_linkage_check.f90 src/fibre_stage15_production_structure_hook.f90"
    if search_report '^[[:space:]]*(rg|grep[[:space:]]+-R)' stage15_checks/run_stage15_10_rhs_ibm_structure_contamination_audit.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_10_wrapper_must_not_require_ripgrep_or_recursive_grep"
        status=1
    fi
    if search_report '^[[:space:]]*(use|call)[[:space:]]+.*(poisson|projection|pressure|rk3|channel_forcing|production_ibm|legacy_ibm|ibm_forcing|apply_ibm|bending|tension|wall_contact|contact|multifibre|multi_fibre|implicit_structure|full_structure|stage14.*rhs.*inject)' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_active_forbidden_solver_or_structure_call_found"
        status=1
    fi
    if search_report '^[[:space:]]*(rhs|ux|uy|uz)[[:space:]]*=' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_direct_fluid_field_assignment_found"
        status=1
    fi
    [ "$status" = "0" ] && SOURCE_GUARD_STATUS=1
    return $status
}

check_xcompact_hook() {
    status=0
    search_silent 'stage15_production_structure_hook_register' src/xcompact3d.f90 || { add_reason "xcompact3d_missing_stage15_register_hook"; status=1; }
    search_silent 'stage15_production_structure_hook_apply' src/xcompact3d.f90 || { add_reason "xcompact3d_missing_stage15_apply_hook"; status=1; }
    search_silent 'stage15_production_structure_hook_finalize' src/xcompact3d.f90 || { add_reason "xcompact3d_missing_stage15_finalize_hook"; status=1; }
    if search_report 'stage15_10.*call|call.*stage15_10|stage15.*production.*advance|legacy_ibm|production_ibm' src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "xcompact3d_unapproved_stage15_10_full_advance_or_ibm_call"
        status=1
    fi
    [ "$status" = "0" ] && XCOMPACT_HOOK_STATUS=1
    return $status
}

check_regression_markers() {
    status=0
    if [ "$STAGE15_10_REQUIRE_STAGE14_CLOSED" = "1" ]; then
        [ -f "$STAGE14_CLOSED_FILE" ] && STAGE14_CLOSED_STATUS=1 || { add_reason "missing_stage14_closed_file"; status=1; }
    else
        STAGE14_CLOSED_STATUS=1
    fi
    if search_report 'stage14_get_injection_gain\(\)[[:space:]]*==[[:space:]]*0\.0' src/*.f90 stage14_checks/*.sh >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage14_forbidden_lambda_zero_registration_gate_found"
        status=1
    else
        STAGE14_LAMBDA_GATE_ABSENT_STATUS=1
    fi
    search_silent 'stage11_5_production_oneway_hook' src/*.f90 stage11_checks/*.sh && STAGE11_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage11_5_production_diagnostics"; status=1; }
    search_silent 'stage13_6_production_force_density_candidate' src/*.f90 stage13_checks/*.sh && STAGE13_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage13_6_production_force_density_diagnostics"; status=1; }
    if search_silent 'stage13_5_production_force_density_candidate' src/*.f90 stage13_checks/*.sh; then
        add_reason "old_stage13_5_force_density_name_reappeared"
        status=1
    fi
    search_silent 'stage14_5_production_rhs_hook' src/*.f90 stage14_checks/*.sh && STAGE14_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage14_5_production_rhs_diagnostics"; status=1; }
    search_silent 'stage15_1_structure_state_buffer' stage15_checks/*.sh stage15_checks/*.md && STAGE15_1_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_1_diagnostics"; status=1; }
    search_silent 'stage15_2_velocity_source_adapter' stage15_checks/*.sh stage15_checks/*.md && STAGE15_2_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_2_diagnostics"; status=1; }
    search_silent 'stage15_3_structure_advance_formula' stage15_checks/*.sh stage15_checks/*.md && STAGE15_3_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_3_diagnostics"; status=1; }
    search_silent 'stage15_4_production_structure_hook' stage15_checks/*.sh stage15_checks/*.md && STAGE15_4_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_4_diagnostics"; status=1; }
    search_silent 'stage15_5_structure_noop_invariance' stage15_checks/*.sh stage15_checks/*.md && STAGE15_5_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_5_diagnostics"; status=1; }
    search_silent 'stage15_6_controlled_structure_step_np1' stage15_checks/*.sh stage15_checks/*.md && STAGE15_6_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_6_diagnostics"; status=1; }
    search_silent 'stage15_7_feedback_linkage' stage15_checks/*.sh stage15_checks/*.md src/*.f90 && STAGE15_7_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_7_diagnostics"; status=1; }
    search_silent 'stage15_8_controlled_structure_parallel_consistency' stage15_checks/*.sh stage15_checks/*.md && STAGE15_8_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_8_diagnostics"; status=1; }
    search_silent 'stage15_9_io_restart_structure_state' stage15_checks/*.sh stage15_checks/*.md && STAGE15_9_DIAGNOSTIC_STATUS=1 || { add_reason "missing_stage15_9_diagnostics"; status=1; }
    if search_silent 'rank[[:space:]]*==[[:space:]]*0|nrank[[:space:]]*==[[:space:]]*0|myid[[:space:]]*==[[:space:]]*0|rank0' src/*.f90 stage11_checks/*.sh stage13_checks/*.sh stage14_checks/*.sh stage15_checks/*.sh; then
        RANK0_DIAGNOSTIC_STATUS=1
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
    if [ "$STAGE15_10_REQUIRE_STAGE15_9" = "1" ]; then
        [ -f stage15_checks/run_stage15_9_io_restart_structure_state.sh ] || { add_reason "missing_stage15_9_wrapper"; status=1; }
        [ -f stage15_checks/stage15_9_io_restart_structure_state.md ] || { add_reason "missing_stage15_9_doc"; status=1; }
    fi
    [ -f stage15_checks/run_stage15_10_rhs_ibm_structure_contamination_audit.sh ] || { add_reason "missing_stage15_10_wrapper"; status=1; }
    [ -f stage15_checks/stage15_10_rhs_ibm_structure_contamination_audit.md ] || { add_reason "missing_stage15_10_doc"; status=1; }
    [ -f stage14_checks/run_stage14_10_rhs_ibm_structure_contamination_audit.sh ] || { add_reason "missing_stage14_10_contamination_wrapper"; status=1; }
    search_silent 'fibre_stage15_feedback_linkage_check' src/CMakeLists.txt || { add_reason "missing_stage15_feedback_linkage_check_target"; status=1; }
    [ "$status" = "0" ] && DOCS_AND_TARGET_STATUS=1
    return $status
}

run_static_audit() {
    status=0
    scan_active_stage15_forbidden || status=1
    check_xcompact_hook || status=1
    check_regression_markers || status=1
    [ "$status" = "0" ] && STATIC_AUDIT_STATUS=1
    return $status
}

run_stage15_linkage_audit() {
    status=0
    exe=$(stage15_7_exe) || { add_reason "missing_fibre_stage15_feedback_linkage_check_executable"; return 1; }
    rm -f "$STAGE15_7_DAT" "$STAGE15_10_LINKAGE_DAT"
    export STAGE15_7_ENABLE="$STAGE15_10_ENABLE"
    export STAGE15_7_CONTROLLED_STEP_ENABLE="$STAGE15_10_CONTROLLED_STEP_ENABLE"
    export STAGE15_7_STRUCTURE_ADVANCE_ENABLE="$STAGE15_10_STRUCTURE_ADVANCE_ENABLE"
    export STAGE15_7_DIAGNOSTIC_ONLY="$STAGE15_10_DIAGNOSTIC_ONLY"
    export STAGE15_7_NP=1
    export STAGE15_7_NPTS="$STAGE15_10_NPTS"
    export STAGE15_7_DT="$STAGE15_10_DT"
    export STAGE15_7_RHO_TILDE="$STAGE15_10_RHO_TILDE"
    export STAGE15_7_TEST_FORCE="$STAGE15_10_TEST_FORCE"
    export STAGE15_7_FEEDBACK_ALPHA="$STAGE15_10_FEEDBACK_ALPHA"
    export STAGE15_7_LAMBDA="$STAGE15_10_LAMBDA"
    export STAGE15_7_MAX_VELOCITY_UPDATE=1.0e-9
    export STAGE15_7_MAX_SLIP_ERROR=1.0e-14
    export STAGE15_7_MAX_FORCE_ERROR=1.0e-14
    export STAGE15_7_MAX_FORCE_RESPONSE="$STAGE15_10_MAX_FORCE_RESPONSE"
    export STAGE15_7_MAX_RHS_RESPONSE="$STAGE15_10_MAX_RHS_RESPONSE"
    # shellcheck disable=SC2086
    "$MPIEXEC" $MPIEXEC_FLAGS -n "$STAGE15_10_NP" "$exe" > "$OUTPUT_DIR/stage15_10_feedback_linkage.log" 2>&1 || { add_reason "stage15_10_feedback_linkage_run_failed"; return 1; }
    [ -s "$STAGE15_7_DAT" ] || { add_reason "missing_stage15_10_feedback_linkage_dat"; return 1; }
    cp "$STAGE15_7_DAT" "$STAGE15_10_LINKAGE_DAT"
    require_key_value "$STAGE15_10_LINKAGE_DAT" final_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" velocity_update_nonzero_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" feedback_force_consistency_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" force_response_bounded_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" rhs_response_bounded_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" controlled_update_count 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" production_full_structure_advance_count 0 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" bending_solve_count 0 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" tension_solve_count 0 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" wall_contact_count 0 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" multifibre_count 0 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" rhs_injection_connection_count 0 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" approved_stage12_13_14_chain_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" no_fluid_rhs_modification_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" no_pressure_projection_modification_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" no_poisson_modification_status 1 || status=1
    require_key_value "$STAGE15_10_LINKAGE_DAT" no_rk3_channel_forcing_modification_status 1 || status=1
    require_real_le_key "$STAGE15_10_LINKAGE_DAT" max_feedback_force_change "$STAGE15_10_MAX_FORCE_RESPONSE" || status=1
    for key in max_velocity_update max_slip_change slip_error max_feedback_force_change feedback_force_error; do
        require_finite_key "$STAGE15_10_LINKAGE_DAT" "$key" || status=1
    done
    if [ "$status" = "0" ]; then
        CONTROLLED_UPDATE_STATUS=1
        FEEDBACK_LINKAGE_STATUS=1
        FORCE_RESPONSE_BOUNDED_STATUS=1
        RHS_RESPONSE_BOUNDED_STATUS=1
        CONTROLLED_STRUCTURE_UPDATE_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" controlled_update_count)
        PRODUCTION_FULL_STRUCTURE_ADVANCE_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" production_full_structure_advance_count)
        BENDING_SOLVE_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" bending_solve_count)
        TENSION_SOLVE_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" tension_solve_count)
        WALL_CONTACT_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" wall_contact_count)
        MULTIFIBRE_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" multifibre_count)
        DIRECT_RHS_INJECTION_CONNECTION_COUNT=$(get_value "$STAGE15_10_LINKAGE_DAT" rhs_injection_connection_count)
        APPROVED_STAGE12_13_14_CHAIN_STATUS=$(get_value "$STAGE15_10_LINKAGE_DAT" approved_stage12_13_14_chain_status)
        NO_FLUID_RHS_MODIFICATION_STATUS=$(get_value "$STAGE15_10_LINKAGE_DAT" no_fluid_rhs_modification_status)
        NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=$(get_value "$STAGE15_10_LINKAGE_DAT" no_pressure_projection_modification_status)
        NO_POISSON_MODIFICATION_STATUS=$(get_value "$STAGE15_10_LINKAGE_DAT" no_poisson_modification_status)
        NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value "$STAGE15_10_LINKAGE_DAT" no_rk3_channel_forcing_modification_status)
        NO_CHANNEL_FORCING_MODIFICATION_STATUS=1
    fi
    return $status
}

run_stage14_10_contamination_audit() {
    status=0
    if [ "$STAGE15_10_RUN_PRODUCTION_SMOKE" != "1" ]; then
        PRODUCTION_SMOKE_DEFERRED_STATUS=1
        add_reason "stage15_10_production_smoke_disabled_static_and_standalone_only"
        return 1
    fi
    STAGE14_10_NP="$STAGE15_10_NP" \
    STAGE14_10_SMALL_LAMBDA="$STAGE15_10_LAMBDA" \
    STAGE14_10_MAX_RHS_INCREMENT="$STAGE15_10_MAX_STAGE14_RHS_INCREMENT" \
    STAGE14_10_MAX_FLUID_DELTA="$STAGE15_10_MAX_FLUID_DELTA" \
    BUILD_DIR="$BUILD_DIR" DECOMP2D_ROOT="$DECOMP2D_ROOT" MPIEXEC="$MPIEXEC" MPIEXEC_FLAGS="$MPIEXEC_FLAGS" \
        bash stage14_checks/run_stage14_10_rhs_ibm_structure_contamination_audit.sh > "$OUTPUT_DIR/stage15_10_stage14_10_contamination.log" 2>&1 || {
            add_reason "stage15_10_stage14_10_contamination_wrapper_failed"
            status=1
        }
    [ -s "$STAGE14_10_DAT" ] || { add_reason "missing_stage14_10_contamination_dat"; return 1; }
    cp "$STAGE14_10_DAT" "$STAGE15_10_STAGE14_DAT"
    require_key_value "$STAGE14_10_DAT" stage14_10_rhs_ibm_structure_contamination_audit_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_stage13_force_density_active_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_stage14_hook_active_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_rhs_increment_bounded_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_fluid_response_bounded_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_nan_inf_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_pressure_projection_contamination_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_poisson_modification_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_rk3_channel_forcing_contamination_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_bending_solve_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_tension_solve_status 1 || status=1
    require_key_value "$STAGE14_10_DAT" stage14_10_no_wall_contact_status 1 || status=1
    require_real_le_key "$STAGE14_10_DAT" stage14_10_max_rhs_increment "$STAGE15_10_MAX_STAGE14_RHS_INCREMENT" || status=1
    require_real_le_key "$STAGE14_10_DAT" stage14_10_max_fluid_signature_delta "$STAGE15_10_MAX_FLUID_DELTA" || status=1
    FLUID_SIGNATURE_DELTA=$(get_value "$STAGE14_10_DAT" stage14_10_max_fluid_signature_delta 2>/dev/null || echo 0.0)
    if [ "$status" = "0" ]; then
        PRODUCTION_SMOKE_STATUS=1
        STAGE13_FORCE_DENSITY_STATUS=1
        STAGE14_RHS_STATUS=1
        STAGE14_RHS_INCREMENT_BOUNDED_STATUS=1
        FLUID_SIGNATURE_STATUS=1
        NO_NAN_INF_STATUS=1
        LEGACY_IBM_FORCING_COUNT=0
        PRODUCTION_IBM_FORCING_OUTSIDE_APPROVED_CHAIN_COUNT=0
        IMPLICIT_STRUCTURE_SOLVE_COUNT=0
    fi
    return $status
}

write_output_dat() {
    final_status=$1
    cat > "$SUMMARY_DAT" <<EOF_DAT
stage15_10_requested_status 1
np $STAGE15_10_NP
npts $STAGE15_10_NPTS
dt $STAGE15_10_DT
rho_tilde $STAGE15_10_RHO_TILDE
test_force_magnitude $STAGE15_10_TEST_FORCE
feedback_alpha $STAGE15_10_FEEDBACK_ALPHA
lambda_value $STAGE15_10_LAMBDA
stage15_10_build_status $BUILD_STATUS
stage15_10_static_audit_status $STATIC_AUDIT_STATUS
stage15_10_source_guard_status $SOURCE_GUARD_STATUS
stage15_10_xcompact_hook_status $XCOMPACT_HOOK_STATUS
stage15_10_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_10_stage15_9_status $STAGE15_9_STATUS
stage15_10_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_10_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_10_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_10_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_10_stage15_1_diagnostic_status $STAGE15_1_DIAGNOSTIC_STATUS
stage15_10_stage15_2_diagnostic_status $STAGE15_2_DIAGNOSTIC_STATUS
stage15_10_stage15_3_diagnostic_status $STAGE15_3_DIAGNOSTIC_STATUS
stage15_10_stage15_4_diagnostic_status $STAGE15_4_DIAGNOSTIC_STATUS
stage15_10_stage15_5_diagnostic_status $STAGE15_5_DIAGNOSTIC_STATUS
stage15_10_stage15_6_diagnostic_status $STAGE15_6_DIAGNOSTIC_STATUS
stage15_10_stage15_7_diagnostic_status $STAGE15_7_DIAGNOSTIC_STATUS
stage15_10_stage15_8_diagnostic_status $STAGE15_8_DIAGNOSTIC_STATUS
stage15_10_stage15_9_diagnostic_status $STAGE15_9_DIAGNOSTIC_STATUS
stage15_10_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_10_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_10_docs_and_target_status $DOCS_AND_TARGET_STATUS
production_smoke_status $PRODUCTION_SMOKE_STATUS
production_smoke_deferred_status $PRODUCTION_SMOKE_DEFERRED_STATUS
controlled_update_status $CONTROLLED_UPDATE_STATUS
feedback_linkage_status $FEEDBACK_LINKAGE_STATUS
stage13_force_density_status $STAGE13_FORCE_DENSITY_STATUS
stage14_rhs_status $STAGE14_RHS_STATUS
force_response_bounded_status $FORCE_RESPONSE_BOUNDED_STATUS
rhs_response_bounded_status $RHS_RESPONSE_BOUNDED_STATUS
stage14_rhs_increment_bounded_status $STAGE14_RHS_INCREMENT_BOUNDED_STATUS
fluid_signature_delta $FLUID_SIGNATURE_DELTA
fluid_signature_status $FLUID_SIGNATURE_STATUS
approved_stage12_13_14_chain_status $APPROVED_STAGE12_13_14_CHAIN_STATUS
direct_rhs_injection_connection_count $DIRECT_RHS_INJECTION_CONNECTION_COUNT
legacy_ibm_forcing_count $LEGACY_IBM_FORCING_COUNT
production_ibm_forcing_outside_approved_chain_count $PRODUCTION_IBM_FORCING_OUTSIDE_APPROVED_CHAIN_COUNT
production_full_structure_advance_count $PRODUCTION_FULL_STRUCTURE_ADVANCE_COUNT
controlled_structure_update_count $CONTROLLED_STRUCTURE_UPDATE_COUNT
bending_solve_count $BENDING_SOLVE_COUNT
tension_solve_count $TENSION_SOLVE_COUNT
implicit_structure_solve_count $IMPLICIT_STRUCTURE_SOLVE_COUNT
wall_contact_count $WALL_CONTACT_COUNT
multifibre_count $MULTIFIBRE_COUNT
no_fluid_rhs_modification_status $NO_FLUID_RHS_MODIFICATION_STATUS
no_pressure_projection_modification_status $NO_PRESSURE_PROJECTION_MODIFICATION_STATUS
no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
no_rk3_channel_forcing_modification_status $NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS
no_channel_forcing_modification_status $NO_CHANNEL_FORCING_MODIFICATION_STATUS
no_nan_inf_status $NO_NAN_INF_STATUS
final_status $final_status
EOF_DAT
    cp "$SUMMARY_DAT" "$OUT_DAT"
}

if ! validate_controls; then
    BUILD_STATUS=0
fi

if [ "$BUILD_STATUS" = "1" ]; then
    build_target fibre_stage15_feedback_linkage_check || { BUILD_STATUS=0; add_reason "build_failed_fibre_stage15_feedback_linkage_check"; }
fi

run_static_audit || true

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$STAGE15_10_RUN_STAGE15_9" = "1" ]; then
    bash stage15_checks/run_stage15_9_io_restart_structure_state.sh > "$OUTPUT_DIR/stage15_10_stage15_9_prereq.log" 2>&1
    if [ $? -eq 0 ] && grep 'STAGE 15.9 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage15_10_stage15_9_prereq.log" >/dev/null 2>&1; then
        STAGE15_9_STATUS=1
    else
        add_reason "stage15_10_stage15_9_prerequisite_failed"
    fi
elif [ "$STAGE15_10_REQUIRE_STAGE15_9" = "1" ] && [ -f stage15_checks/run_stage15_9_io_restart_structure_state.sh ] && [ -f stage15_checks/stage15_9_io_restart_structure_state.md ]; then
    STAGE15_9_STATUS=1
elif [ "$STAGE15_10_REQUIRE_STAGE15_9" != "1" ]; then
    STAGE15_9_STATUS=1
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$STAGE15_9_STATUS" = "1" ]; then
    run_stage15_linkage_audit || true
else
    add_reason "stage15_10_linkage_audit_skipped_due_to_build_static_or_prereq_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$STAGE15_9_STATUS" = "1" ] && [ "$CONTROLLED_UPDATE_STATUS" = "1" ]; then
    run_stage14_10_contamination_audit || true
else
    add_reason "stage15_10_stage14_contamination_audit_skipped_due_to_missing_linkage_evidence"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$STAGE15_9_STATUS" = "1" ] && \
   [ "$CONTROLLED_UPDATE_STATUS" = "1" ] && [ "$FEEDBACK_LINKAGE_STATUS" = "1" ] && \
   [ "$STAGE13_FORCE_DENSITY_STATUS" = "1" ] && [ "$STAGE14_RHS_STATUS" = "1" ] && \
   [ "$FORCE_RESPONSE_BOUNDED_STATUS" = "1" ] && [ "$RHS_RESPONSE_BOUNDED_STATUS" = "1" ] && \
   [ "$STAGE14_RHS_INCREMENT_BOUNDED_STATUS" = "1" ] && [ "$FLUID_SIGNATURE_STATUS" = "1" ] && \
   [ "$APPROVED_STAGE12_13_14_CHAIN_STATUS" = "1" ] && [ "$DIRECT_RHS_INJECTION_CONNECTION_COUNT" = "0" ] && \
   [ "$LEGACY_IBM_FORCING_COUNT" = "0" ] && [ "$PRODUCTION_IBM_FORCING_OUTSIDE_APPROVED_CHAIN_COUNT" = "0" ] && \
   [ "$PRODUCTION_FULL_STRUCTURE_ADVANCE_COUNT" = "0" ] && [ "$CONTROLLED_STRUCTURE_UPDATE_COUNT" = "1" ] && \
   [ "$BENDING_SOLVE_COUNT" = "0" ] && [ "$TENSION_SOLVE_COUNT" = "0" ] && [ "$IMPLICIT_STRUCTURE_SOLVE_COUNT" = "0" ] && \
   [ "$WALL_CONTACT_COUNT" = "0" ] && [ "$MULTIFIBRE_COUNT" = "0" ] && \
   [ "$NO_FLUID_RHS_MODIFICATION_STATUS" = "1" ] && [ "$NO_PRESSURE_PROJECTION_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_POISSON_MODIFICATION_STATUS" = "1" ] && [ "$NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_CHANNEL_FORCING_MODIFICATION_STATUS" = "1" ] && [ "$NO_NAN_INF_STATUS" = "1" ]; then
    write_output_dat 1
    echo 'STAGE 15.10 RHS/IBM/STRUCTURE CONTAMINATION AUDIT VERDICT: PASS'
    echo 'STAGE 15.10 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 15.10 RHS/IBM/STRUCTURE CONTAMINATION AUDIT VERDICT: FAIL'
echo 'STAGE 15.10 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_10_failure'
fi
exit 1
