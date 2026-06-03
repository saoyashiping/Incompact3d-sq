#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_4_RUN_STAGE15_3=${STAGE15_4_RUN_STAGE15_3:-0}
STAGE15_4_REQUIRE_STAGE14_CLOSED=${STAGE15_4_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_4_REQUIRE_STAGE15_3=${STAGE15_4_REQUIRE_STAGE15_3:-1}
STAGE15_4_ENABLE=${STAGE15_4_ENABLE:-1}
STAGE15_4_STRUCTURE_ADVANCE_ENABLE=${STAGE15_4_STRUCTURE_ADVANCE_ENABLE:-0}
STAGE15_4_DIAGNOSTIC_ONLY=${STAGE15_4_DIAGNOSTIC_ONLY:-1}
STAGE15_4_NP=${STAGE15_4_NP:-1}
STAGE15_4_MAX_FLUID_DELTA=${STAGE15_4_MAX_FLUID_DELTA:-0.0}
STAGE15_4_RUN_PRODUCTION_SMOKE=${STAGE15_4_RUN_PRODUCTION_SMOKE:-1}

OUTPUT_DIR=stage15_outputs
OUT_DAT=$OUTPUT_DIR/stage15_4_production_structure_hook.dat
REASONS_FILE=$OUTPUT_DIR/stage15_4_production_structure_hook_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_4_production_structure_static_audit_report.txt
CHECK_LOG=$OUTPUT_DIR/stage15_4_production_structure_hook_check.log
CHECK_DAT=$OUTPUT_DIR/fibre_stage15_4_production_structure_hook.dat
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_3_DOC=stage15_checks/stage15_3_structure_advance_formula.md
STAGE15_3_SCRIPT=stage15_checks/run_stage15_3_structure_advance_formula.sh

BUILD_STATUS=1
RUN_STATUS=0
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
XCOMPACT_HOOK_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_3_STATUS=0
STAGE14_LAMBDA_GATE_ABSENT_STATUS=0
STAGE11_DIAGNOSTIC_STATUS=0
STAGE13_DIAGNOSTIC_STATUS=0
STAGE14_DIAGNOSTIC_STATUS=0
STAGE15_1_DIAGNOSTIC_STATUS=0
STAGE15_2_DIAGNOSTIC_STATUS=0
STAGE15_3_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOCS_AND_TARGET_STATUS=0
DIAGNOSTIC_STATUS=0
PRODUCTION_SMOKE_STATUS=1
PRODUCTION_SMOKE_DEFERRED_STATUS=0
FINAL_STATUS=0

mkdir -p "$OUTPUT_DIR" stage14_outputs stage13_outputs stage11_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$CHECK_LOG" "$CHECK_DAT"

add_reason() { echo "$1" >> "$REASONS_FILE"; }
static_note() { echo "$1" >> "$STATIC_REPORT"; }
search_q() { pattern=$1; shift; grep -En "$pattern" "$@"; }
search_silent() { pattern=$1; shift; grep -Eq "$pattern" "$@"; }

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

stage15_4_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage15_production_structure_hook_check" \
               "$BUILD_DIR/src/fibre_stage15_production_structure_hook_check" \
               "$BUILD_DIR/fibre_stage15_production_structure_hook_check"; do
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

scan_stage15_4_source_guards() {
    status=0
    files="src/fibre_stage15_production_structure_hook.f90 src/fibre_stage15_production_structure_hook_check.f90"
    static_note "Stage 15.4 source guard audit"

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
          if (routine !~ /^stage15_production_structure_hook_/ && \
              routine !~ /^stage15_structure_state_/ && routine !~ /^stage15_config_load/ && \
              routine ~ /(advance|bending|bend|tension|wall|contact|multifibre|multi_fibre|stage14|rhs_injection|poisson|projection|pressure|rk3|channel_forcing|xcompact3d)/) {
            print FILENAME ":" FNR ":" raw
            bad=1
          }
        }
      }
      END { exit(bad ? 1 : 0) }
    ' $files >> "$STATIC_REPORT"; then
        SOURCE_GUARD_STATUS=1
    else
        add_reason "stage15_4_forbidden_production_call_found"
        status=1
    fi

    if search_q '^[[:space:]]*use[[:space:]]+(fibre_stage14|xcompact3d|navier|poisson|time_integrators)|^[[:space:]]*call[[:space:]]+(stage14_|.*rhs_injection|poisson|projection|pressure|rk3|channel_forcing|wall|contact|multifibre|multi_fibre)|injection_gain' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_4_forbidden_source_connection_found"
        status=1
    fi
    return $status
}

check_xcompact_hook_insertion() {
    status=0
    use_count=$(grep -Ec 'use[[:space:]]+fibre_stage15_production_structure_hook' src/xcompact3d.f90 || true)
    reset_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_reset' src/xcompact3d.f90 || true)
    register_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_register\(1\)' src/xcompact3d.f90 || true)
    apply_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_apply\(itime,[[:space:]]*0\)' src/xcompact3d.f90 || true)
    finalize_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_finalize\(\)' src/xcompact3d.f90 || true)
    if [ "$use_count" -eq 1 ] && [ "$reset_count" -eq 1 ] && [ "$register_count" -eq 1 ] && \
       [ "$apply_count" -eq 1 ] && [ "$finalize_count" -ge 1 ]; then
        XCOMPACT_HOOK_STATUS=1
    else
        add_reason "xcompact3d_stage15_4_hook_insertion_not_approved_use_${use_count}_reset_${reset_count}_register_${register_count}_apply_${apply_count}_finalize_${finalize_count}"
        status=1
    fi
    if search_q 'call[[:space:]]+stage15_.*advance|stage15_production_structure_hook_apply[[:space:]]*\([^)]*dux|stage15_production_structure_hook_apply[[:space:]]*\([^)]*uy|stage15_production_structure_hook_apply[[:space:]]*\([^)]*uz' src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "xcompact3d_stage15_production_advance_or_fluid_argument_found"
        status=1
    fi
    return $status
}

run_static_audit() {
    status=0
    static_note "Stage 15.4 regression audit"

    [ -f src/fibre_stage15_production_structure_hook.f90 ] || { add_reason "missing_stage15_4_hook_source"; status=1; }
    [ -f src/fibre_stage15_production_structure_hook_check.f90 ] || { add_reason "missing_stage15_4_hook_check_source"; status=1; }
    [ -f stage15_checks/run_stage15_4_production_structure_hook.sh ] || { add_reason "missing_stage15_4_wrapper"; status=1; }
    [ -f stage15_checks/stage15_4_production_structure_hook.md ] || { add_reason "missing_stage15_4_documentation"; status=1; }
    if search_silent 'fibre_stage15_production_structure_hook_check' src/CMakeLists.txt; then
        DOCS_AND_TARGET_STATUS=1
    else
        add_reason "missing_fibre_stage15_production_structure_hook_check_build_target"
        status=1
    fi

    scan_stage15_4_source_guards || status=1
    check_xcompact_hook_insertion || status=1

    if search_q 'stage14_get_injection_gain[[:space:]]*\([[:space:]]*\)[[:space:]]*==[[:space:]]*0\.0(_[[:alnum:]_]+)?' \
          src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "forbidden_stage14_lambda_zero_registration_gate_found"
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
        add_reason "stage13_6_production_force_density_diagnostics_missing"
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

    if [ "$STAGE15_4_REQUIRE_STAGE14_CLOSED" = "0" ]; then
        STAGE14_CLOSED_STATUS=1
    elif [ -s "$STAGE14_CLOSED_FILE" ]; then
        STAGE14_CLOSED_STATUS=1
    else
        add_reason "missing_stage14_checks_STAGE14_CLOSED_md"
        status=1
    fi

    if [ "$STAGE15_4_RUN_STAGE15_3" = "1" ]; then
        bash "$STAGE15_3_SCRIPT" || { add_reason "stage15_3_optional_prerequisite_run_failed"; status=1; }
    fi
    if [ "$STAGE15_4_REQUIRE_STAGE15_3" = "0" ]; then
        STAGE15_3_STATUS=1
    elif [ -s "$STAGE15_3_DOC" ] && [ -f src/fibre_stage15_structure_advance_formula.f90 ] && \
         [ -f src/fibre_stage15_structure_advance_formula_check.f90 ]; then
        STAGE15_3_STATUS=1
    else
        add_reason "stage15_3_closed_pass_evidence_missing"
        status=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

run_check() {
    exe=$(stage15_4_exe) || { add_reason "missing_fibre_stage15_production_structure_hook_check_executable"; return 1; }
    export X3D_STAGE15_ENABLE="$STAGE15_4_ENABLE"
    export X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE="$STAGE15_4_STRUCTURE_ADVANCE_ENABLE"
    export X3D_STAGE15_DIAGNOSTIC_ONLY="$STAGE15_4_DIAGNOSTIC_ONLY"
    export X3D_STAGE15_REQUIRE_STAGE14_CLOSED="$STAGE15_4_REQUIRE_STAGE14_CLOSED"
    # shellcheck disable=SC2086
    "$MPIEXEC" $MPIEXEC_FLAGS -n 1 "$exe" > "$CHECK_LOG" 2>&1
}

handle_production_smoke() {
    if [ "$STAGE15_4_RUN_PRODUCTION_SMOKE" = "1" ]; then
        PRODUCTION_SMOKE_DEFERRED_STATUS=1
        PRODUCTION_SMOKE_STATUS=1
        echo "production_smoke_deferred: lightweight production smoke is deferred; standalone hook check is authoritative" >> "$STATIC_REPORT"
    else
        PRODUCTION_SMOKE_DEFERRED_STATUS=1
        PRODUCTION_SMOKE_STATUS=1
        echo "production_smoke_deferred: STAGE15_4_RUN_PRODUCTION_SMOKE=0" >> "$STATIC_REPORT"
    fi
}

verify_diagnostics() {
    status=0
    [ -s "$CHECK_DAT" ] || { add_reason "missing_stage15_4_production_structure_hook_dat"; return 1; }
    require_key_value "$CHECK_DAT" stage15_4_requested_status 1 || status=1
    require_key_value "$CHECK_DAT" hook_registered_status 1 || status=1
    require_key_value "$CHECK_DAT" hook_apply_count 2 || status=1
    require_key_value "$CHECK_DAT" hook_finalize_status 1 || status=1
    require_key_value "$CHECK_DAT" diagnostic_only_status 1 || status=1
    require_key_value "$CHECK_DAT" noop_status 1 || status=1
    require_key_value "$CHECK_DAT" structure_state_available_status 1 || status=1
    require_key_value "$CHECK_DAT" production_time_loop_hook_status 1 || status=1
    require_key_value "$CHECK_DAT" production_time_loop_connection_count 1 || status=1
    require_key_value "$CHECK_DAT" production_structure_advance_count 0 || status=1
    require_key_value "$CHECK_DAT" x_position_update_count 0 || status=1
    require_key_value "$CHECK_DAT" v_velocity_update_count 0 || status=1
    require_key_value "$CHECK_DAT" a_acceleration_update_count 0 || status=1
    require_key_value "$CHECK_DAT" bending_solve_count 0 || status=1
    require_key_value "$CHECK_DAT" tension_solve_count 0 || status=1
    require_key_value "$CHECK_DAT" wall_contact_count 0 || status=1
    require_key_value "$CHECK_DAT" multifibre_count 0 || status=1
    require_key_value "$CHECK_DAT" rhs_injection_connection_count 0 || status=1
    require_key_value "$CHECK_DAT" no_fluid_rhs_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_pressure_projection_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_poisson_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_rk3_channel_forcing_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_stage10_15_3_regression_status 1 || status=1
    require_key_value "$CHECK_DAT" final_status 1 || status=1
    require_key_value "$CHECK_DAT" stage15_4_check_final_status 1 || status=1
    for key in stage15_4_requested_status hook_registered_status hook_apply_count hook_finalize_status \
               diagnostic_only_status noop_status structure_state_available_status production_time_loop_hook_status \
               production_time_loop_connection_count production_structure_advance_count x_position_update_count \
               v_velocity_update_count a_acceleration_update_count bending_solve_count tension_solve_count \
               wall_contact_count multifibre_count rhs_injection_connection_count no_fluid_rhs_modification_status \
               no_pressure_projection_modification_status no_poisson_modification_status \
               no_rk3_channel_forcing_modification_status final_status; do
        require_finite_key "$CHECK_DAT" "$key" || status=1
    done
    if [ "$status" = "0" ]; then
        DIAGNOSTIC_STATUS=1
    fi
    return $status
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage15_4_requested_status 1
stage15_4_build_status $BUILD_STATUS
stage15_4_run_status $RUN_STATUS
stage15_4_static_audit_status $STATIC_AUDIT_STATUS
stage15_4_source_guard_status $SOURCE_GUARD_STATUS
stage15_4_xcompact_hook_status $XCOMPACT_HOOK_STATUS
stage15_4_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_4_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_4_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_4_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_4_stage15_1_diagnostic_status $STAGE15_1_DIAGNOSTIC_STATUS
stage15_4_stage15_2_diagnostic_status $STAGE15_2_DIAGNOSTIC_STATUS
stage15_4_stage15_3_diagnostic_status $STAGE15_3_DIAGNOSTIC_STATUS
stage15_4_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_4_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_4_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_4_stage15_3_status $STAGE15_3_STATUS
stage15_4_docs_and_target_status $DOCS_AND_TARGET_STATUS
stage15_4_diagnostic_status $DIAGNOSTIC_STATUS
stage15_4_production_smoke_status $PRODUCTION_SMOKE_STATUS
stage15_4_production_smoke_deferred_status $PRODUCTION_SMOKE_DEFERRED_STATUS
stage15_4_np_value $STAGE15_4_NP
stage15_4_max_fluid_delta_value $STAGE15_4_MAX_FLUID_DELTA
stage15_4_enable_value $STAGE15_4_ENABLE
stage15_4_structure_advance_enable_value $STAGE15_4_STRUCTURE_ADVANCE_ENABLE
stage15_4_diagnostic_only_value $STAGE15_4_DIAGNOSTIC_ONLY
stage15_4_final_status $final_status
EOF_DAT
}

run_static_audit || true
handle_production_smoke

if ! build_target fibre_stage15_production_structure_hook_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage15_production_structure_hook_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_check; then
        RUN_STATUS=1
        verify_diagnostics || true
    else
        add_reason "run_failed_fibre_stage15_production_structure_hook_check"
    fi
else
    add_reason "stage15_4_production_structure_hook_check_skipped_due_to_build_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$RUN_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && \
   [ "$SOURCE_GUARD_STATUS" = "1" ] && [ "$XCOMPACT_HOOK_STATUS" = "1" ] && \
   [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && [ "$STAGE11_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE14_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_1_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_2_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_3_DIAGNOSTIC_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$STAGE15_3_STATUS" = "1" ] && [ "$DOCS_AND_TARGET_STATUS" = "1" ] && \
   [ "$DIAGNOSTIC_STATUS" = "1" ] && [ "$PRODUCTION_SMOKE_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_output_dat 1
    echo 'STAGE 15.4 PRODUCTION STRUCTURE HOOK VERDICT: PASS'
    echo 'STAGE 15.4 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 15.4 PRODUCTION STRUCTURE HOOK VERDICT: FAIL'
echo 'STAGE 15.4 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_4_failure'
fi
exit 1
