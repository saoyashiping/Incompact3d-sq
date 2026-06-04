#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_2_RUN_STAGE15_1=${STAGE15_2_RUN_STAGE15_1:-0}
STAGE15_2_REQUIRE_STAGE14_CLOSED=${STAGE15_2_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_2_REQUIRE_STAGE15_1=${STAGE15_2_REQUIRE_STAGE15_1:-1}
STAGE15_2_NPTS=${STAGE15_2_NPTS:-8}
STAGE15_2_ENABLE=${STAGE15_2_ENABLE:-1}
STAGE15_2_STRUCTURE_ADVANCE_ENABLE=${STAGE15_2_STRUCTURE_ADVANCE_ENABLE:-0}
STAGE15_2_DIAGNOSTIC_ONLY=${STAGE15_2_DIAGNOSTIC_ONLY:-1}
STAGE15_2_MAX_VELOCITY_DIFF=${STAGE15_2_MAX_VELOCITY_DIFF:-1.0e-14}
STAGE15_2_MAX_FORCE_DIFF=${STAGE15_2_MAX_FORCE_DIFF:-1.0e-14}

OUTPUT_DIR=stage15_outputs
OUT_DAT=$OUTPUT_DIR/stage15_2_velocity_source_adapter.dat
REASONS_FILE=$OUTPUT_DIR/stage15_2_velocity_source_adapter_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_2_velocity_source_static_audit_report.txt
CHECK_LOG=$OUTPUT_DIR/stage15_2_velocity_source_adapter_check.log
CHECK_DAT=$OUTPUT_DIR/fibre_stage15_2_velocity_source_adapter.dat
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_1_DOC=stage15_checks/stage15_1_structure_state_buffer.md
STAGE15_1_SCRIPT=stage15_checks/run_stage15_1_structure_state_buffer.sh

BUILD_STATUS=1
RUN_STATUS=0
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_1_STATUS=0
STAGE14_LAMBDA_GATE_ABSENT_STATUS=0
STAGE11_DIAGNOSTIC_STATUS=0
STAGE13_DIAGNOSTIC_STATUS=0
STAGE14_DIAGNOSTIC_STATUS=0
STAGE15_1_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOCS_AND_TARGET_STATUS=0
DIAGNOSTIC_STATUS=0
FINAL_STATUS=0

mkdir -p "$OUTPUT_DIR" stage14_outputs stage13_outputs stage11_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$CHECK_LOG" "$CHECK_DAT"

add_reason() { echo "$1" >> "$REASONS_FILE"; }
static_note() { echo "$1" >> "$STATIC_REPORT"; }
search_report() { pattern=$1; shift; grep -En "$pattern" "$@"; }
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

stage15_2_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage15_velocity_source_adapter_check" \
               "$BUILD_DIR/src/fibre_stage15_velocity_source_adapter_check" \
               "$BUILD_DIR/fibre_stage15_velocity_source_adapter_check"; do
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

require_numeric_le() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || { add_reason "missing_${key}_in_${file}"; return 1; }
    awk -v v="$value" -v l="$limit" 'BEGIN { gsub(/[dD]/,"E",v); gsub(/[dD]/,"E",l); if ((v+0) <= (l+0)) exit 0; exit 1 }' || {
        add_reason "${key}_${value}_exceeds_${limit}_in_${file}"
        return 1
    }
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

scan_stage15_2_source_guards() {
    status=0
    files="src/fibre_stage15_velocity_source_adapter.f90 src/fibre_stage15_velocity_source_adapter_check.f90"
    static_note "Stage 15.2 source guard audit"

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
          if (routine ~ /^(structure_advance|advance_structure|bending|bend_solve|solve_bending|tension_solve|solve_tension|update_position|update_velocity|wall|contact|multifibre|multi_fibre|poisson|projection|pressure|rk3|channel_forcing)/) {
            print FILENAME ":" FNR ":" raw
            bad=1
          }
        }
      }
      END { exit(bad ? 1 : 0) }
    ' $files >> "$STATIC_REPORT"; then
        SOURCE_GUARD_STATUS=1
    else
        add_reason "stage15_2_forbidden_structure_or_solver_call_found"
        status=1
    fi

    for file in src/poisson.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/Case-Channel.f90; do
        if [ -f "$file" ] && search_report 'fibre_stage15_velocity_source|stage15_2_' "$file" >> "$STATIC_REPORT" 2>/dev/null; then
            add_reason "stage15_2_static_match_in_unrelated_dns_file_${file}"
            status=1
        fi
    done
    return $status
}

run_static_audit() {
    status=0
    static_note "Stage 15.2 regression audit"

    [ -f src/fibre_stage15_velocity_source_adapter.f90 ] || { add_reason "missing_stage15_2_adapter_source"; status=1; }
    [ -f src/fibre_stage15_velocity_source_adapter_check.f90 ] || { add_reason "missing_stage15_2_adapter_check_source"; status=1; }
    [ -f stage15_checks/run_stage15_2_velocity_source_adapter.sh ] || { add_reason "missing_stage15_2_wrapper"; status=1; }
    [ -f stage15_checks/stage15_2_velocity_source_adapter.md ] || { add_reason "missing_stage15_2_documentation"; status=1; }
    if search_silent 'fibre_stage15_velocity_source_adapter_check' src/CMakeLists.txt; then
        DOCS_AND_TARGET_STATUS=1
    else
        add_reason "missing_fibre_stage15_velocity_source_adapter_check_build_target"
        status=1
    fi

    scan_stage15_2_source_guards || status=1

    if search_report 'stage14_get_injection_gain[[:space:]]*\([[:space:]]*\)[[:space:]]*==[[:space:]]*0\.0(_[[:alnum:]_]+)?' \
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
        add_reason "stage13_production_force_density_diagnostics_missing"
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

    if search_silent 'rank0_write_allowed' src/fibre_stage11_production_oneway_hook.f90 && \
       search_silent 'rank0_write_allowed' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_silent 'rank0_write_allowed' src/fibre_stage14_production_rhs_injection.f90; then
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

    if [ "$STAGE15_2_REQUIRE_STAGE14_CLOSED" = "0" ]; then
        STAGE14_CLOSED_STATUS=1
    elif [ -s "$STAGE14_CLOSED_FILE" ]; then
        STAGE14_CLOSED_STATUS=1
    else
        add_reason "missing_stage14_checks_STAGE14_CLOSED_md"
        status=1
    fi

    if [ "$STAGE15_2_RUN_STAGE15_1" = "1" ]; then
        bash "$STAGE15_1_SCRIPT" || { add_reason "stage15_1_optional_prerequisite_run_failed"; status=1; }
    fi
    if [ "$STAGE15_2_REQUIRE_STAGE15_1" = "0" ]; then
        STAGE15_1_STATUS=1
    elif [ -s "$STAGE15_1_DOC" ] && [ -f src/fibre_stage15_structure_state.f90 ] && \
         [ -f src/fibre_stage15_structure_state_check.f90 ]; then
        STAGE15_1_STATUS=1
    else
        add_reason "stage15_1_closed_pass_evidence_missing"
        status=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

run_check() {
    exe=$(stage15_2_exe) || { add_reason "missing_fibre_stage15_velocity_source_adapter_check_executable"; return 1; }
    export X3D_STAGE15_ENABLE="$STAGE15_2_ENABLE"
    export X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE="$STAGE15_2_STRUCTURE_ADVANCE_ENABLE"
    export X3D_STAGE15_DIAGNOSTIC_ONLY="$STAGE15_2_DIAGNOSTIC_ONLY"
    export X3D_STAGE15_REQUIRE_STAGE14_CLOSED="$STAGE15_2_REQUIRE_STAGE14_CLOSED"
    export STAGE15_2_NPTS STAGE15_2_MAX_VELOCITY_DIFF STAGE15_2_MAX_FORCE_DIFF
    # shellcheck disable=SC2086
    "$MPIEXEC" $MPIEXEC_FLAGS -n 1 "$exe" > "$CHECK_LOG" 2>&1
}

verify_diagnostics() {
    status=0
    [ -s "$CHECK_DAT" ] || { add_reason "missing_stage15_2_velocity_source_adapter_dat"; return 1; }
    require_key_value "$CHECK_DAT" stage15_2_requested_status 1 || status=1
    require_key_value "$CHECK_DAT" stage15_structure_owned_velocity_status 1 || status=1
    require_key_value "$CHECK_DAT" prescribed_velocity_reference_status 1 || status=1
    require_key_value "$CHECK_DAT" prescribed_stage12_reference_status 1 || status=1
    require_key_value "$CHECK_DAT" velocity_source_adapter_status 1 || status=1
    require_key_value "$CHECK_DAT" npts "$STAGE15_2_NPTS" || status=1
    require_numeric_le "$CHECK_DAT" max_velocity_source_diff "$STAGE15_2_MAX_VELOCITY_DIFF" || status=1
    require_numeric_le "$CHECK_DAT" max_feedback_force_diff "$STAGE15_2_MAX_FORCE_DIFF" || status=1
    require_key_value "$CHECK_DAT" velocity_equivalence_status 1 || status=1
    require_key_value "$CHECK_DAT" feedback_force_equivalence_status 1 || status=1
    require_key_value "$CHECK_DAT" zero_slip_status 1 || status=1
    require_key_value "$CHECK_DAT" finite_value_status 1 || status=1
    require_key_value "$CHECK_DAT" structure_advance_count 0 || status=1
    require_key_value "$CHECK_DAT" bending_solve_count 0 || status=1
    require_key_value "$CHECK_DAT" tension_solve_count 0 || status=1
    require_key_value "$CHECK_DAT" position_time_update_count 0 || status=1
    require_key_value "$CHECK_DAT" velocity_time_update_count 0 || status=1
    require_key_value "$CHECK_DAT" no_fluid_rhs_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_pressure_projection_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_poisson_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_rk3_channel_forcing_modification_status 1 || status=1
    require_key_value "$CHECK_DAT" no_stage11_14_regression_status 1 || status=1
    require_key_value "$CHECK_DAT" final_status 1 || status=1
    require_key_value "$CHECK_DAT" stage15_2_check_final_status 1 || status=1
    for key in stage15_2_requested_status stage15_structure_owned_velocity_status prescribed_velocity_reference_status \
               velocity_source_adapter_status npts max_velocity_source_diff max_feedback_force_diff \
               velocity_equivalence_status feedback_force_equivalence_status zero_slip_status finite_value_status \
               structure_advance_count bending_solve_count tension_solve_count position_time_update_count \
               velocity_time_update_count no_fluid_rhs_modification_status no_pressure_projection_modification_status \
               no_poisson_modification_status no_rk3_channel_forcing_modification_status final_status; do
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
stage15_2_requested_status 1
stage15_2_build_status $BUILD_STATUS
stage15_2_run_status $RUN_STATUS
stage15_2_static_audit_status $STATIC_AUDIT_STATUS
stage15_2_source_guard_status $SOURCE_GUARD_STATUS
stage15_2_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_2_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_2_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_2_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_2_stage15_1_diagnostic_status $STAGE15_1_DIAGNOSTIC_STATUS
stage15_2_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_2_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_2_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_2_stage15_1_status $STAGE15_1_STATUS
stage15_2_docs_and_target_status $DOCS_AND_TARGET_STATUS
stage15_2_diagnostic_status $DIAGNOSTIC_STATUS
stage15_2_npts_value $STAGE15_2_NPTS
stage15_2_max_velocity_diff_value $STAGE15_2_MAX_VELOCITY_DIFF
stage15_2_max_force_diff_value $STAGE15_2_MAX_FORCE_DIFF
stage15_2_enable_value $STAGE15_2_ENABLE
stage15_2_structure_advance_enable_value $STAGE15_2_STRUCTURE_ADVANCE_ENABLE
stage15_2_diagnostic_only_value $STAGE15_2_DIAGNOSTIC_ONLY
stage15_2_final_status $final_status
EOF_DAT
}

run_static_audit || true

if ! build_target fibre_stage15_velocity_source_adapter_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage15_velocity_source_adapter_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_check; then
        RUN_STATUS=1
        verify_diagnostics || true
    else
        add_reason "run_failed_fibre_stage15_velocity_source_adapter_check"
    fi
else
    add_reason "stage15_2_velocity_source_adapter_check_skipped_due_to_build_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$RUN_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && \
   [ "$SOURCE_GUARD_STATUS" = "1" ] && [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && \
   [ "$STAGE11_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE13_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE14_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_1_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && \
   [ "$STAGE14_CLOSED_STATUS" = "1" ] && [ "$STAGE15_1_STATUS" = "1" ] && \
   [ "$DOCS_AND_TARGET_STATUS" = "1" ] && [ "$DIAGNOSTIC_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_output_dat 1
    echo 'STAGE 15.2 VELOCITY SOURCE ADAPTER VERDICT: PASS'
    echo 'STAGE 15.2 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 15.2 VELOCITY SOURCE ADAPTER VERDICT: FAIL'
echo 'STAGE 15.2 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_2_failure'
fi
exit 1
