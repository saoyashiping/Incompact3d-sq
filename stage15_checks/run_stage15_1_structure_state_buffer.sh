#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_1_RUN_STAGE15_0=${STAGE15_1_RUN_STAGE15_0:-0}
STAGE15_1_REQUIRE_STAGE14_CLOSED=${STAGE15_1_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_1_REQUIRE_STAGE15_0=${STAGE15_1_REQUIRE_STAGE15_0:-1}
STAGE15_1_NPTS=${STAGE15_1_NPTS:-8}
STAGE15_1_ENABLE=${STAGE15_1_ENABLE:-0}
STAGE15_1_STRUCTURE_ADVANCE_ENABLE=${STAGE15_1_STRUCTURE_ADVANCE_ENABLE:-0}
STAGE15_1_DIAGNOSTIC_ONLY=${STAGE15_1_DIAGNOSTIC_ONLY:-1}

OUTPUT_DIR=stage15_outputs
OUT_DAT=$OUTPUT_DIR/stage15_1_structure_state_buffer.dat
REASONS_FILE=$OUTPUT_DIR/stage15_1_structure_state_buffer_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_1_structure_state_static_audit_report.txt
CHECK_LOG=$OUTPUT_DIR/stage15_1_structure_state_check.log
CHECK_DAT=$OUTPUT_DIR/fibre_stage15_1_structure_state_buffer.dat
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_0_DOC=stage15_checks/stage15_0_config.md

BUILD_STATUS=1
RUN_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_0_STATUS=0
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
NO_STRUCTURE_ADVANCE_STATIC_STATUS=0
NO_BENDING_SOLVE_STATIC_STATUS=0
NO_TENSION_SOLVE_STATIC_STATUS=0
NO_POSITION_UPDATE_STATIC_STATUS=0
NO_VELOCITY_UPDATE_STATIC_STATUS=0
NO_STAGE12_CONNECT_STATIC_STATUS=0
NO_STAGE14_CONNECT_STATIC_STATUS=0
NO_WALL_CONTACT_STATIC_STATUS=0
NO_MULTIFIBRE_STATIC_STATUS=0
NO_FLUID_RHS_STATIC_STATUS=0
NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=0
NO_RK3_CHANNEL_FORCING_STATIC_STATUS=0
STAGE14_LAMBDA_GATE_ABSENT_STATUS=0
STAGE11_DIAGNOSTIC_STATUS=0
STAGE13_DIAGNOSTIC_STATUS=0
STAGE14_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOC_STATUS=0
DIAGNOSTIC_STATUS=0
FINAL_STATUS=0

mkdir -p "$OUTPUT_DIR" stage14_outputs stage13_outputs stage11_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$CHECK_LOG" "$CHECK_DAT"

add_reason() { echo "$1" >> "$REASONS_FILE"; }
static_note() { echo "$1" >> "$STATIC_REPORT"; }

search_regex() {
    pattern=$1
    shift
    if command -v rg >/dev/null 2>&1; then
        rg -n "$pattern" "$@"
    else
        grep -En "$pattern" "$@"
    fi
}

search_regex_quiet() {
    pattern=$1
    shift
    if command -v rg >/dev/null 2>&1; then
        rg -q "$pattern" "$@"
    else
        grep -Eq "$pattern" "$@"
    fi
}

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

stage15_1_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage15_structure_state_check" \
               "$BUILD_DIR/src/fibre_stage15_structure_state_check" \
               "$BUILD_DIR/fibre_stage15_structure_state_check"; do
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
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    if [ "$value" != "$expected" ]; then
        add_reason "${key}_expected_${expected}_got_${value}_in_${file}"
        return 1
    fi
    return 0
}

require_finite_key() {
    file=$1
    key=$2
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    awk -v v="$value" 'BEGIN { if (v ~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([EeDd][+-]?[0-9]+)?$/) exit 0; exit 1 }' || {
        add_reason "nonfinite_or_non_numeric_${key}_${value}_in_${file}"
        return 1
    }
    return 0
}

scan_stage15_1_source_guards() {
    status=0
    files="src/fibre_stage15_structure_state.f90 src/fibre_stage15_structure_state_check.f90"
    static_note "Stage 15.1 source guard audit"

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
          if (routine ~ /(advance|bending|bend|tension_solve|solve_tension|position_update|velocity_update|wall|contact|multifibre|multi_fibre|poisson|projection|pressure|rk3|channel_forcing)/) {
            print FILENAME ":" FNR ":" raw
            bad=1
          }
        }
      }
      END { exit(bad ? 1 : 0) }
    ' $files >> "$STATIC_REPORT"; then
        NO_STRUCTURE_ADVANCE_STATIC_STATUS=1
        NO_BENDING_SOLVE_STATIC_STATUS=1
        NO_TENSION_SOLVE_STATIC_STATUS=1
        NO_POSITION_UPDATE_STATIC_STATUS=1
        NO_VELOCITY_UPDATE_STATIC_STATUS=1
        NO_WALL_CONTACT_STATIC_STATUS=1
        NO_MULTIFIBRE_STATIC_STATUS=1
        NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=1
        NO_RK3_CHANNEL_FORCING_STATIC_STATUS=1
    else
        add_reason "stage15_1_forbidden_structure_or_solver_call_found"
        status=1
    fi

    if search_regex "use[[:space:]]+fibre_stage12|stage12_|feedback_force|feedback_candidate" $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_1_stage12_feedback_connection_found"
        status=1
    else
        NO_STAGE12_CONNECT_STATIC_STATUS=1
    fi

    if search_regex "use[[:space:]]+fibre_stage14|stage14_|rhs_injection|injection_gain" $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_1_stage14_rhs_injection_connection_found"
        status=1
    else
        NO_STAGE14_CONNECT_STATIC_STATUS=1
    fi

    # Audit only real production-fluid/RHS access, not diagnostic status names such as
    # no_fluid_rhs_modification_status.  The previous broad grep matched those
    # negative diagnostic fields and caused a false Stage 15.1 failure even though
    # the structure-state skeleton does not access fluid arrays or RHS fields.
    if awk '
      /^[[:space:]]*!/ { next }
      {
        raw=$0
        line=tolower($0)
        sub(/!.*/, "", line)
        # calls to real fluid/RHS routines are forbidden in Stage 15.1
        if (line ~ /^[[:space:]]*call[[:space:]]+/) {
          routine=line
          sub(/^[[:space:]]*call[[:space:]]+/, "", routine)
          sub(/[[:space:]]*\(.*/, "", routine)
          gsub(/[[:space:]]/, "", routine)
          if (routine ~ /(fluid_rhs|rhs_injection|rhs_accumulator|navier|transeq|pressure|projection|poisson|rk3|channel_forcing)/) {
            print FILENAME ":" FNR ":" raw
            bad=1
          }
        }
        # direct access to production velocity/RHS arrays is also forbidden.
        # This deliberately does not match diagnostic names like no_fluid_rhs_*.
        if (line ~ /(^|[^[:alnum:]_])(ux|uy|uz|gx|gy|gz)[[:space:]]*\(/) {
          print FILENAME ":" FNR ":" raw
          bad=1
        }
        if (line ~ /^[[:space:]]*use[[:space:]]+(navier|transeq|poisson|variables|xcompact3d|fibre_stage14_production_rhs_injection)/) {
          print FILENAME ":" FNR ":" raw
          bad=1
        }
      }
      END { exit(bad ? 1 : 0) }
    ' $files >> "$STATIC_REPORT"; then
        NO_FLUID_RHS_STATIC_STATUS=1
    else
        add_reason "stage15_1_fluid_rhs_or_field_access_found"
        status=1
    fi

    if [ "$status" = "0" ]; then
        SOURCE_GUARD_STATUS=1
    fi
    return $status
}

run_static_audit() {
    status=0
    static_note "Stage 15.1 regression audit"

    [ -f src/fibre_stage15_structure_state.f90 ] || { add_reason "missing_src_fibre_stage15_structure_state_f90"; status=1; }
    [ -f src/fibre_stage15_structure_state_check.f90 ] || { add_reason "missing_src_fibre_stage15_structure_state_check_f90"; status=1; }
    [ -f stage15_checks/run_stage15_1_structure_state_buffer.sh ] || { add_reason "missing_stage15_1_wrapper"; status=1; }
    [ -f stage15_checks/stage15_1_structure_state_buffer.md ] || { add_reason "missing_stage15_1_documentation"; status=1; }
    if search_regex_quiet "fibre_stage15_structure_state_check" src/CMakeLists.txt; then
        DOC_STATUS=1
    else
        add_reason "missing_fibre_stage15_structure_state_check_build_target"
        status=1
    fi

    scan_stage15_1_source_guards || status=1

    if search_regex 'stage14_get_injection_gain[[:space:]]*\([[:space:]]*\)[[:space:]]*==[[:space:]]*0\.0(_[[:alnum:]_]+)?' src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "forbidden_stage14_lambda_zero_registration_gate_found"
        status=1
    else
        STAGE14_LAMBDA_GATE_ABSENT_STATUS=1
    fi

    if search_regex_quiet 'fibre_stage11_5_production_oneway_hook.dat' src/fibre_stage11_production_oneway_hook.f90 && \
       search_regex_quiet 'stage11_5_production_oneway_hook_status' src/fibre_stage11_production_oneway_hook.f90; then
        STAGE11_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage11_5_production_oneway_hook_diagnostics_missing"
        status=1
    fi

    if search_regex_quiet 'fibre_stage13_6_production_force_density_candidate.dat' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_regex_quiet 'stage13_6_production_force_density_candidate_status' src/fibre_stage13_production_force_density_candidate.f90; then
        STAGE13_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage13_production_force_density_diagnostics_missing"
        status=1
    fi

    if search_regex_quiet 'fibre_stage14_5_production_rhs_hook.dat' src/fibre_stage14_production_rhs_injection.f90 && \
       search_regex_quiet 'stage14_5_nonzero_lambda_blocked_status' src/fibre_stage14_production_rhs_injection.f90 && \
       search_regex_quiet 'stage14_5_production_rhs_hook_status' src/fibre_stage14_production_rhs_injection.f90; then
        STAGE14_DIAGNOSTIC_STATUS=1
    else
        add_reason "stage14_small_lambda_production_rhs_hook_diagnostics_missing"
        status=1
    fi

    if search_regex_quiet 'rank0_write_allowed' src/fibre_stage11_production_oneway_hook.f90 && \
       search_regex_quiet 'rank0_write_allowed' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_regex_quiet 'rank0_write_allowed' src/fibre_stage14_production_rhs_injection.f90; then
        RANK0_DIAGNOSTIC_STATUS=1
    else
        add_reason "rank0_safe_diagnostic_writing_regressed"
        status=1
    fi

    if search_regex_quiet 'lbound\(ux, 1\)[[:space:]]*\+[[:space:]]*2' src/fibre_stage13_production_force_density_candidate.f90 && \
       search_regex_quiet 'np=1/2/4' src/fibre_stage13_production_force_density_candidate.f90; then
        STAGE13_SAMPLING_REPAIR_STATUS=1
    else
        add_reason "stage13_force_density_sampling_repair_missing"
        status=1
    fi

    if [ "$STAGE15_1_REQUIRE_STAGE14_CLOSED" = "0" ]; then
        STAGE14_CLOSED_STATUS=1
    elif [ -s "$STAGE14_CLOSED_FILE" ]; then
        STAGE14_CLOSED_STATUS=1
    else
        add_reason "missing_stage14_checks_STAGE14_CLOSED_md"
        status=1
    fi

    if [ "$STAGE15_1_RUN_STAGE15_0" = "1" ]; then
        bash stage15_checks/run_stage15_0_config.sh || { add_reason "stage15_0_optional_prerequisite_run_failed"; status=1; }
    fi
    if [ "$STAGE15_1_REQUIRE_STAGE15_0" = "0" ]; then
        STAGE15_0_STATUS=1
    elif [ -s "$STAGE15_0_DOC" ] && [ -f src/fibre_stage15_config.f90 ] && [ -f src/fibre_stage15_config_check.f90 ]; then
        STAGE15_0_STATUS=1
    else
        add_reason "stage15_0_closed_pass_evidence_missing"
        status=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

run_check() {
    exe=$(stage15_1_exe) || {
        add_reason "missing_fibre_stage15_structure_state_check_executable"
        return 1
    }
    export X3D_STAGE15_ENABLE="$STAGE15_1_ENABLE"
    export X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE="$STAGE15_1_STRUCTURE_ADVANCE_ENABLE"
    export X3D_STAGE15_DIAGNOSTIC_ONLY="$STAGE15_1_DIAGNOSTIC_ONLY"
    export X3D_STAGE15_REQUIRE_STAGE14_CLOSED="$STAGE15_1_REQUIRE_STAGE14_CLOSED"
    export STAGE15_1_NPTS
    # shellcheck disable=SC2086
    "$MPIEXEC" $MPIEXEC_FLAGS -n 1 "$exe" > "$CHECK_LOG" 2>&1
}

verify_diagnostics() {
    status=0
    [ -s "$CHECK_DAT" ] || { add_reason "missing_stage15_1_structure_state_diagnostic_dat"; return 1; }
    require_key_value "$CHECK_DAT" stage15_1_requested_status 1 || status=1
    require_key_value "$CHECK_DAT" allocation_status 1 || status=1
    require_key_value "$CHECK_DAT" initialization_status 1 || status=1
    require_key_value "$CHECK_DAT" clear_status 1 || status=1
    require_key_value "$CHECK_DAT" validation_status 1 || status=1
    require_key_value "$CHECK_DAT" npts "$STAGE15_1_NPTS" || status=1
    require_key_value "$CHECK_DAT" x_finite_status 1 || status=1
    require_key_value "$CHECK_DAT" v_finite_status 1 || status=1
    require_key_value "$CHECK_DAT" optional_a_or_rhs_finite_status 1 || status=1
    require_key_value "$CHECK_DAT" optional_tension_finite_status 1 || status=1
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
    require_key_value "$CHECK_DAT" stage15_1_check_final_status 1 || status=1
    for key in stage15_1_requested_status allocation_status initialization_status clear_status validation_status npts \
               x_finite_status v_finite_status optional_a_or_rhs_finite_status optional_tension_finite_status \
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
stage15_1_requested_status 1
stage15_1_build_status $BUILD_STATUS
stage15_1_run_status $RUN_STATUS
stage15_1_static_audit_status $STATIC_AUDIT_STATUS
stage15_1_source_guard_status $SOURCE_GUARD_STATUS
stage15_1_no_structure_advance_static_status $NO_STRUCTURE_ADVANCE_STATIC_STATUS
stage15_1_no_bending_solve_static_status $NO_BENDING_SOLVE_STATIC_STATUS
stage15_1_no_tension_solve_static_status $NO_TENSION_SOLVE_STATIC_STATUS
stage15_1_no_position_update_static_status $NO_POSITION_UPDATE_STATIC_STATUS
stage15_1_no_velocity_update_static_status $NO_VELOCITY_UPDATE_STATIC_STATUS
stage15_1_no_stage12_connect_static_status $NO_STAGE12_CONNECT_STATIC_STATUS
stage15_1_no_stage14_connect_static_status $NO_STAGE14_CONNECT_STATIC_STATUS
stage15_1_no_wall_contact_static_status $NO_WALL_CONTACT_STATIC_STATUS
stage15_1_no_multifibre_static_status $NO_MULTIFIBRE_STATIC_STATUS
stage15_1_no_fluid_rhs_static_status $NO_FLUID_RHS_STATIC_STATUS
stage15_1_no_pressure_projection_poisson_static_status $NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS
stage15_1_no_rk3_channel_forcing_static_status $NO_RK3_CHANNEL_FORCING_STATIC_STATUS
stage15_1_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_1_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_1_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_1_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_1_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_1_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_1_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_1_stage15_0_status $STAGE15_0_STATUS
stage15_1_docs_and_target_status $DOC_STATUS
stage15_1_diagnostic_status $DIAGNOSTIC_STATUS
stage15_1_npts_value $STAGE15_1_NPTS
stage15_1_enable_value $STAGE15_1_ENABLE
stage15_1_structure_advance_enable_value $STAGE15_1_STRUCTURE_ADVANCE_ENABLE
stage15_1_diagnostic_only_value $STAGE15_1_DIAGNOSTIC_ONLY
stage15_1_final_status $final_status
EOF_DAT
}

run_static_audit || true

if ! build_target fibre_stage15_structure_state_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage15_structure_state_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_check; then
        RUN_STATUS=1
        verify_diagnostics || true
    else
        add_reason "run_failed_fibre_stage15_structure_state_check"
    fi
else
    add_reason "stage15_1_structure_state_check_skipped_due_to_build_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$RUN_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && \
   [ "$SOURCE_GUARD_STATUS" = "1" ] && [ "$NO_STRUCTURE_ADVANCE_STATIC_STATUS" = "1" ] && \
   [ "$NO_BENDING_SOLVE_STATIC_STATUS" = "1" ] && [ "$NO_TENSION_SOLVE_STATIC_STATUS" = "1" ] && \
   [ "$NO_POSITION_UPDATE_STATIC_STATUS" = "1" ] && [ "$NO_VELOCITY_UPDATE_STATIC_STATUS" = "1" ] && \
   [ "$NO_STAGE12_CONNECT_STATIC_STATUS" = "1" ] && [ "$NO_STAGE14_CONNECT_STATIC_STATUS" = "1" ] && \
   [ "$NO_WALL_CONTACT_STATIC_STATUS" = "1" ] && [ "$NO_MULTIFIBRE_STATIC_STATUS" = "1" ] && \
   [ "$NO_FLUID_RHS_STATIC_STATUS" = "1" ] && [ "$NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS" = "1" ] && \
   [ "$NO_RK3_CHANNEL_FORCING_STATIC_STATUS" = "1" ] && [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && \
   [ "$STAGE11_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE13_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE14_DIAGNOSTIC_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$STAGE15_0_STATUS" = "1" ] && [ "$DOC_STATUS" = "1" ] && [ "$DIAGNOSTIC_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_output_dat 1
    echo 'STAGE 15.1 STRUCTURE STATE BUFFER VERDICT: PASS'
    echo 'STAGE 15.1 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 15.1 STRUCTURE STATE BUFFER VERDICT: FAIL'
echo 'STAGE 15.1 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_1_failure'
fi
exit 1
