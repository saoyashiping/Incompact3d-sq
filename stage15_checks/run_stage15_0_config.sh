#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_0_RUN_STAGE14_11=${STAGE15_0_RUN_STAGE14_11:-0}
STAGE15_0_REQUIRE_STAGE14_CLOSED=${STAGE15_0_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_0_ENABLE=${STAGE15_0_ENABLE:-0}
STAGE15_0_STRUCTURE_ADVANCE_ENABLE=${STAGE15_0_STRUCTURE_ADVANCE_ENABLE:-0}
STAGE15_0_DIAGNOSTIC_ONLY=${STAGE15_0_DIAGNOSTIC_ONLY:-1}

OUTPUT_DIR=stage15_outputs
OUT_DAT=$OUTPUT_DIR/stage15_0_config.dat
REASONS_FILE=$OUTPUT_DIR/stage15_0_config_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_0_static_audit_report.txt
CONFIG_LOG=$OUTPUT_DIR/stage15_0_config_check.log
CONFIG_DAT=$OUTPUT_DIR/fibre_stage15_0_config.dat
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md

BUILD_STATUS=1
STAGE14_11_PREREQ_STATUS=1
STAGE14_CLOSED_STATUS=0
STATIC_AUDIT_STATUS=0
STAGE15_SOURCE_STATIC_STATUS=0
NO_STRUCTURE_ADVANCE_STATIC_STATUS=0
NO_BENDING_SOLVE_STATIC_STATUS=0
NO_TENSION_SOLVE_STATIC_STATUS=0
NO_FIBRE_POSITION_UPDATE_STATIC_STATUS=0
NO_FIBRE_VELOCITY_UPDATE_STATIC_STATUS=0
NO_WALL_CONTACT_STATIC_STATUS=0
NO_MULTIFIBRE_STATIC_STATUS=0
NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=0
NO_RK3_CHANNEL_FORCING_STATIC_STATUS=0
NO_PRODUCTION_IBM_STATIC_STATUS=0
STAGE14_LAMBDA_GATE_ABSENT_STATUS=0
STAGE11_DIAGNOSTIC_PATH_STATUS=0
STAGE13_DIAGNOSTIC_PATH_STATUS=0
STAGE14_DIAGNOSTIC_PATH_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
CONFIG_CHECK_RUN_STATUS=0
CONFIG_DEFAULT_SAFE_STATUS=0
CONFIG_REQUEST_PARSE_STATUS=0
CONFIG_STRUCTURE_ADVANCE_DISABLED_STATUS=0
CONFIG_DIAGNOSTIC_ONLY_STATUS=0
CONFIG_REQUIRE_STAGE14_CLOSED_STATUS=0
CONFIG_NO_PRODUCTION_BEHAVIOR_CHANGE_STATUS=0
CONFIG_NO_FLUID_RHS_CHANGE_STATUS=0
CONFIG_STATUS=0
FINAL_STATUS=0
STATIC_MATCH_COUNT=0

mkdir -p stage15_outputs stage14_outputs stage13_outputs stage11_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$CONFIG_LOG" "$CONFIG_DAT"

add_reason() {
    echo "$1" >> "$REASONS_FILE"
}

static_note() {
    echo "$1" >> "$STATIC_REPORT"
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

stage15_config_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage15_config_check" \
               "$BUILD_DIR/src/fibre_stage15_config_check" \
               "$BUILD_DIR/fibre_stage15_config_check"; do
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

require_file() {
    file=$1
    reason=$2
    if [ ! -s "$file" ]; then
        add_reason "$reason"
        return 1
    fi
    return 0
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

scan_active_use_call_forbidden() {
    label=$1
    shift
    awk -v label="$label" '
      BEGIN { bad=0 }
      /^[[:space:]]*!/ { next }
      {
        line=tolower($0)
        if (line ~ /^[[:space:]]*call[[:space:]]+/) {
          if (line ~ /(structure_advance|advance_structure|bending|bend_solve|tension|position_update|velocity_update|wall|contact|multifibre|multi_fibre|ibm|immersed|poisson|projection|pressure|rk3|channel_forcing)/) {
            print label ":" FILENAME ":" FNR ":" $0
            bad=1
          }
        }
      }
      END { exit(bad ? 1 : 0) }
    ' "$@" >> "$STATIC_REPORT"
}

run_static_audit() {
    status=0
    static_note "Stage 15.0 static audit"

    if scan_active_use_call_forbidden stage15_forbidden_use_call src/fibre_stage15_config.f90 src/fibre_stage15_config_check.f90; then
        STAGE15_SOURCE_STATIC_STATUS=1
        NO_STRUCTURE_ADVANCE_STATIC_STATUS=1
        NO_BENDING_SOLVE_STATIC_STATUS=1
        NO_TENSION_SOLVE_STATIC_STATUS=1
        NO_FIBRE_POSITION_UPDATE_STATIC_STATUS=1
        NO_FIBRE_VELOCITY_UPDATE_STATIC_STATUS=1
        NO_WALL_CONTACT_STATIC_STATUS=1
        NO_MULTIFIBRE_STATIC_STATUS=1
        NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=1
        NO_RK3_CHANNEL_FORCING_STATIC_STATUS=1
        NO_PRODUCTION_IBM_STATIC_STATUS=1
    else
        add_reason "stage15_forbidden_structure_fluid_or_ibm_use_call_found"
        status=1
    fi

    static_matches=0
    for file in src/poisson.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/Case-Channel.f90; do
        if [ -f "$file" ]; then
            count=$(grep -E 'fibre_stage15|stage15_' "$file" 2>/dev/null | wc -l | awk '{print $1}')
            static_matches=$((static_matches + count))
            if [ "$count" -gt 0 ]; then
                grep -nE 'fibre_stage15|stage15_' "$file" >> "$STATIC_REPORT" 2>/dev/null || true
            fi
        fi
    done
    STATIC_MATCH_COUNT=$static_matches
    if [ "$STATIC_MATCH_COUNT" -ne 0 ]; then
        NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=0
        NO_RK3_CHANNEL_FORCING_STATIC_STATUS=0
        add_reason "stage15_static_matches_in_pressure_projection_poisson_rk3_channel_files_${STATIC_MATCH_COUNT}"
        status=1
    fi

    if grep -nE 'stage14_get_injection_gain[[:space:]]*\([[:space:]]*\)[[:space:]]*==[[:space:]]*0\.0(_[[:alnum:]_]+)?' src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "forbidden_stage14_lambda_zero_registration_gate_found"
        status=1
    else
        STAGE14_LAMBDA_GATE_ABSENT_STATUS=1
    fi

    if grep -q 'fibre_stage11_5_production_oneway_hook.dat' src/fibre_stage11_production_oneway_hook.f90 && \
       grep -q 'stage11_5_production_oneway_hook_status' src/fibre_stage11_production_oneway_hook.f90; then
        STAGE11_DIAGNOSTIC_PATH_STATUS=1
    else
        add_reason "stage11_5_production_oneway_diagnostic_path_missing"
        status=1
    fi

    if grep -q 'fibre_stage13_6_production_force_density_candidate.dat' src/fibre_stage13_production_force_density_candidate.f90 && \
       grep -q 'stage13_6_production_force_density_candidate_status' src/fibre_stage13_production_force_density_candidate.f90; then
        STAGE13_DIAGNOSTIC_PATH_STATUS=1
    else
        add_reason "stage13_6_production_force_density_diagnostic_path_missing"
        status=1
    fi

    if grep -q 'fibre_stage14_5_production_rhs_hook.dat' src/fibre_stage14_production_rhs_injection.f90 && \
       grep -q 'stage14_5_production_rhs_hook_status' src/fibre_stage14_production_rhs_injection.f90; then
        STAGE14_DIAGNOSTIC_PATH_STATUS=1
    else
        add_reason "stage14_5_production_rhs_hook_diagnostic_path_missing"
        status=1
    fi

    if { grep -q 'rank0_write_allowed' src/fibre_stage11_production_oneway_hook.f90 || \
         grep -q 'MPI_COMM_RANK' src/fibre_stage11_production_oneway_hook.f90 || \
         grep -q 'nrank' src/fibre_stage11_production_oneway_hook.f90; } && \
       { grep -q 'rank0_write_allowed' src/fibre_stage13_production_force_density_candidate.f90 || \
         grep -q 'MPI_COMM_RANK' src/fibre_stage13_production_force_density_candidate.f90 || \
         grep -q 'nrank' src/fibre_stage13_production_force_density_candidate.f90; } && \
       { grep -q 'rank0_write_allowed' src/fibre_stage14_production_rhs_injection.f90 || \
         grep -q 'MPI_COMM_RANK' src/fibre_stage14_production_rhs_injection.f90 || \
         grep -q 'nrank' src/fibre_stage14_production_rhs_injection.f90; }; then
        RANK0_DIAGNOSTIC_STATUS=1
    else
        add_reason "rank0_safe_stage11_stage13_or_stage14_diagnostic_guard_missing"
        status=1
    fi

    if grep -nE 'i0[[:space:]]*=[[:space:]]*\(lbound\(ux,[[:space:]]*1\)[[:space:]]*\+[[:space:]]*ubound\(ux,[[:space:]]*1\)\)[[:space:]]*/[[:space:]]*2|j0[[:space:]]*=[[:space:]]*\(lbound\(ux,[[:space:]]*2\)[[:space:]]*\+[[:space:]]*ubound\(ux,[[:space:]]*2\)\)[[:space:]]*/[[:space:]]*2|k0[[:space:]]*=[[:space:]]*\(lbound\(ux,[[:space:]]*3\)[[:space:]]*\+[[:space:]]*ubound\(ux,[[:space:]]*3\)\)[[:space:]]*/[[:space:]]*2' src/fibre_stage13_production_force_density_candidate.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage13_force_density_local_subdomain_center_sampling_detected"
        status=1
    else
        STAGE13_SAMPLING_REPAIR_STATUS=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

verify_stage14_closed() {
    if [ "$STAGE15_0_REQUIRE_STAGE14_CLOSED" != "1" ]; then
        STAGE14_CLOSED_STATUS=1
        return 0
    fi
    if [ "$STAGE15_0_RUN_STAGE14_11" = "1" ]; then
        bash stage14_checks/run_stage14_11_total_smoke_closure.sh > "$OUTPUT_DIR/stage15_0_stage14_11_prereq.log" 2>&1
        if [ $? -ne 0 ] || ! grep 'STAGE 14.11 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage15_0_stage14_11_prereq.log" >/dev/null 2>&1; then
            STAGE14_11_PREREQ_STATUS=0
            add_reason "stage14_11_prerequisite_failed"
            return 1
        fi
    fi
    require_file "$STAGE14_CLOSED_FILE" "missing_stage14_closed_file" || return 1
    STAGE14_CLOSED_STATUS=1
    return 0
}

run_config_check() {
    exe=$(stage15_config_exe) || {
        add_reason "missing_fibre_stage15_config_check_executable"
        return 1
    }
    X3D_STAGE15_ENABLE="$STAGE15_0_ENABLE" \
    X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE="$STAGE15_0_STRUCTURE_ADVANCE_ENABLE" \
    X3D_STAGE15_DIAGNOSTIC_ONLY="$STAGE15_0_DIAGNOSTIC_ONLY" \
    X3D_STAGE15_REQUIRE_STAGE14_CLOSED="$STAGE15_0_REQUIRE_STAGE14_CLOSED" \
        "$exe" > "$CONFIG_LOG" 2>&1
    if [ $? -ne 0 ]; then
        add_reason "fibre_stage15_config_check_failed"
        return 1
    fi
    if ! grep -q 'STAGE 15.0 CONFIG VERDICT: PASS' "$CONFIG_LOG"; then
        add_reason "missing_stage15_0_config_pass_line"
        return 1
    fi
    require_file "$CONFIG_DAT" "missing_stage15_0_config_dat" || return 1
    CONFIG_CHECK_RUN_STATUS=1
    return 0
}

verify_config_dat() {
    status=0
    if [ "$STAGE15_0_ENABLE" = "0" ]; then
        require_key_value "$CONFIG_DAT" stage15_0_requested_flag 0 || status=1
        require_key_value "$CONFIG_DAT" stage15_0_disabled_by_default_status 1 || status=1
    fi
    require_key_value "$CONFIG_DAT" stage15_0_request_parse_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_structure_advance_enabled_flag 0 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_structure_advance_disabled_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_structure_advance_blocked_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_diagnostic_only_flag 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_diagnostic_only_default_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_diagnostic_only_enforced_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_require_stage14_closed_flag "$STAGE15_0_REQUIRE_STAGE14_CLOSED" || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_structure_state_allocation_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_structure_advance_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_bending_solve_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_tension_solve_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_fibre_position_update_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_fibre_velocity_update_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_wall_contact_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_multifibre_logic_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_fluid_field_access_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_fluid_field_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_rhs_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_pressure_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_projection_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_poisson_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_rk3_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_stage11_14_reimplementation_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_no_production_behavior_change_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_config_status 1 || status=1
    require_key_value "$CONFIG_DAT" stage15_0_check_status 1 || status=1
    if [ "$status" = "0" ]; then
        CONFIG_DEFAULT_SAFE_STATUS=1
        CONFIG_REQUEST_PARSE_STATUS=1
        CONFIG_STRUCTURE_ADVANCE_DISABLED_STATUS=1
        CONFIG_DIAGNOSTIC_ONLY_STATUS=1
        CONFIG_REQUIRE_STAGE14_CLOSED_STATUS=1
        CONFIG_NO_PRODUCTION_BEHAVIOR_CHANGE_STATUS=1
        CONFIG_NO_FLUID_RHS_CHANGE_STATUS=1
        CONFIG_STATUS=1
    fi
    return $status
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage15_0_requested_flag 1
stage15_0_build_status $BUILD_STATUS
stage15_0_stage14_11_prereq_status $STAGE14_11_PREREQ_STATUS
stage15_0_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_0_static_audit_status $STATIC_AUDIT_STATUS
stage15_0_stage15_source_static_status $STAGE15_SOURCE_STATIC_STATUS
stage15_0_no_structure_advance_static_status $NO_STRUCTURE_ADVANCE_STATIC_STATUS
stage15_0_no_bending_solve_static_status $NO_BENDING_SOLVE_STATIC_STATUS
stage15_0_no_tension_solve_static_status $NO_TENSION_SOLVE_STATIC_STATUS
stage15_0_no_fibre_position_update_static_status $NO_FIBRE_POSITION_UPDATE_STATIC_STATUS
stage15_0_no_fibre_velocity_update_static_status $NO_FIBRE_VELOCITY_UPDATE_STATIC_STATUS
stage15_0_no_wall_contact_static_status $NO_WALL_CONTACT_STATIC_STATUS
stage15_0_no_multifibre_static_status $NO_MULTIFIBRE_STATIC_STATUS
stage15_0_no_pressure_projection_poisson_static_status $NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS
stage15_0_no_rk3_channel_forcing_static_status $NO_RK3_CHANNEL_FORCING_STATIC_STATUS
stage15_0_no_production_ibm_static_status $NO_PRODUCTION_IBM_STATIC_STATUS
stage15_0_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_0_stage11_diagnostic_path_status $STAGE11_DIAGNOSTIC_PATH_STATUS
stage15_0_stage13_diagnostic_path_status $STAGE13_DIAGNOSTIC_PATH_STATUS
stage15_0_stage14_diagnostic_path_status $STAGE14_DIAGNOSTIC_PATH_STATUS
stage15_0_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_0_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_0_config_check_run_status $CONFIG_CHECK_RUN_STATUS
stage15_0_config_default_safe_status $CONFIG_DEFAULT_SAFE_STATUS
stage15_0_config_request_parse_status $CONFIG_REQUEST_PARSE_STATUS
stage15_0_config_structure_advance_disabled_status $CONFIG_STRUCTURE_ADVANCE_DISABLED_STATUS
stage15_0_config_diagnostic_only_status $CONFIG_DIAGNOSTIC_ONLY_STATUS
stage15_0_config_require_stage14_closed_status $CONFIG_REQUIRE_STAGE14_CLOSED_STATUS
stage15_0_config_no_production_behavior_change_status $CONFIG_NO_PRODUCTION_BEHAVIOR_CHANGE_STATUS
stage15_0_config_no_fluid_rhs_change_status $CONFIG_NO_FLUID_RHS_CHANGE_STATUS
stage15_0_config_status $CONFIG_STATUS
stage15_0_config_gate_status $final_status
stage15_0_enable_value $STAGE15_0_ENABLE
stage15_0_structure_advance_enable_value $STAGE15_0_STRUCTURE_ADVANCE_ENABLE
stage15_0_diagnostic_only_value $STAGE15_0_DIAGNOSTIC_ONLY
stage15_0_require_stage14_closed_value $STAGE15_0_REQUIRE_STAGE14_CLOSED
stage15_0_static_match_count $STATIC_MATCH_COUNT
EOF_DAT
}

run_static_audit || true
verify_stage14_closed || true

if [ "$BUILD_STATUS" = "1" ]; then
    build_target fibre_stage15_config_check || {
        BUILD_STATUS=0
        add_reason "build_failed_fibre_stage15_config_check"
    }
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_config_check; then
        verify_config_dat || true
    fi
else
    add_reason "stage15_0_config_check_skipped_due_to_build_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_11_PREREQ_STATUS" = "1" ] && [ "$STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$STAGE15_SOURCE_STATIC_STATUS" = "1" ] && \
   [ "$NO_STRUCTURE_ADVANCE_STATIC_STATUS" = "1" ] && [ "$NO_BENDING_SOLVE_STATIC_STATUS" = "1" ] && \
   [ "$NO_TENSION_SOLVE_STATIC_STATUS" = "1" ] && [ "$NO_FIBRE_POSITION_UPDATE_STATIC_STATUS" = "1" ] && \
   [ "$NO_FIBRE_VELOCITY_UPDATE_STATIC_STATUS" = "1" ] && [ "$NO_WALL_CONTACT_STATIC_STATUS" = "1" ] && \
   [ "$NO_MULTIFIBRE_STATIC_STATUS" = "1" ] && [ "$NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS" = "1" ] && \
   [ "$NO_RK3_CHANNEL_FORCING_STATIC_STATUS" = "1" ] && [ "$NO_PRODUCTION_IBM_STATIC_STATUS" = "1" ] && \
   [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && [ "$STAGE11_DIAGNOSTIC_PATH_STATUS" = "1" ] && \
   [ "$STAGE13_DIAGNOSTIC_PATH_STATUS" = "1" ] && [ "$STAGE14_DIAGNOSTIC_PATH_STATUS" = "1" ] && \
   [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && \
   [ "$CONFIG_CHECK_RUN_STATUS" = "1" ] && [ "$CONFIG_DEFAULT_SAFE_STATUS" = "1" ] && \
   [ "$CONFIG_REQUEST_PARSE_STATUS" = "1" ] && [ "$CONFIG_STRUCTURE_ADVANCE_DISABLED_STATUS" = "1" ] && \
   [ "$CONFIG_DIAGNOSTIC_ONLY_STATUS" = "1" ] && [ "$CONFIG_REQUIRE_STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$CONFIG_NO_PRODUCTION_BEHAVIOR_CHANGE_STATUS" = "1" ] && [ "$CONFIG_NO_FLUID_RHS_CHANGE_STATUS" = "1" ] && \
   [ "$CONFIG_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_output_dat 1
    echo 'STAGE 15.0 CONFIG VERDICT: PASS'
    echo 'STAGE 15.0 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 15.0 CONFIG VERDICT: FAIL'
echo 'STAGE 15.0 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_0_failure'
fi
exit 1
