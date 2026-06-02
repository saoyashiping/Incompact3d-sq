#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
CHANNEL_INPUT=${CHANNEL_INPUT:-examples/Channel/input.i3d}
STAGE14_11_RUN_STAGE14_8=${STAGE14_11_RUN_STAGE14_8:-0}
STAGE14_11_RUN_STAGE14_9=${STAGE14_11_RUN_STAGE14_9:-0}
STAGE14_11_RUN_STAGE14_10=${STAGE14_11_RUN_STAGE14_10:-0}
STAGE14_11_SMALL_LAMBDA=${STAGE14_11_SMALL_LAMBDA:-1.0e-8}
STAGE14_11_MAX_RHS_INCREMENT=${STAGE14_11_MAX_RHS_INCREMENT:-1.0e-4}
STAGE14_11_MAX_FLUID_DELTA=${STAGE14_11_MAX_FLUID_DELTA:-1.0e-4}
STAGE14_11_MAX_PARALLEL_RHS_DIFF=${STAGE14_11_MAX_PARALLEL_RHS_DIFF:-1.0e-12}
STAGE14_11_MAX_PARALLEL_FORCE_DENSITY_DIFF=${STAGE14_11_MAX_PARALLEL_FORCE_DENSITY_DIFF:-1.0e-10}
STAGE14_11_MAX_RESTART_DELTA=${STAGE14_11_MAX_RESTART_DELTA:-1.0e-8}
STAGE14_11_MAX_IO_SIGNATURE_DELTA=${STAGE14_11_MAX_IO_SIGNATURE_DELTA:-1.0e-8}
STAGE14_11_NP=${STAGE14_11_NP:-2}
STAGE14_11_NP_LIST=${STAGE14_11_NP_LIST:-"1 2 4"}
STAGE14_11_TIMEOUT_SEC=${STAGE14_11_TIMEOUT_SEC:-240}
STAGE14_11_MAX_STEPS=${STAGE14_11_MAX_STEPS:-3}

OUTPUT_DIR=stage14_outputs
OUT_DAT=$OUTPUT_DIR/stage14_11_total_smoke_closure.dat
REASONS_FILE=$OUTPUT_DIR/stage14_11_total_smoke_closure_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage14_11_static_audit_report.txt
RUNTIME_LOG=$OUTPUT_DIR/stage14_11_total_smoke_runtime.log
STAGE9_DAT=$OUTPUT_DIR/stage14_11_total_smoke_stage9_9.dat
STAGE11_DAT=$OUTPUT_DIR/stage14_11_stage11_oneway.dat
STAGE12_DAT=$OUTPUT_DIR/stage14_11_stage12_feedback_candidate.dat
STAGE13_DAT=$OUTPUT_DIR/stage14_11_stage13_force_density.dat
STAGE14_DAT=$OUTPUT_DIR/stage14_11_stage14_rhs_hook.dat
CLOSURE_FILE=stage14_checks/STAGE14_CLOSED.md

BUILD_STATUS=1
STATIC_AUDIT_STATUS=0
FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS=0
SMALL_LAMBDA_REGISTRATION_STATUS=0
HOOK_CONNECTED_STATUS=0
STAGE11_DIAGNOSTIC_PATH_STATUS=0
STAGE12_DIAGNOSTIC_PATH_STATUS=0
STAGE13_DIAGNOSTIC_PATH_STATUS=0
STAGE14_DIAGNOSTIC_PATH_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=0
NO_RK3_CHANNEL_FORCING_STATIC_STATUS=0
NO_PRODUCTION_IBM_STATIC_STATUS=0
NO_STRUCTURE_STATIC_STATUS=0
RUNTIME_TOTAL_SMOKE_STATUS=0
STAGE11_ACTIVE_STATUS=0
STAGE12_ACTIVE_STATUS=0
STAGE13_ACTIVE_STATUS=0
STAGE14_HOOK_ACTIVE_STATUS=0
LAMBDA_NONZERO_STATUS=0
RHS_INCREMENT_NONZERO_STATUS=0
RHS_INCREMENT_FINITE_STATUS=0
RHS_INCREMENT_SIGN_CORRECT_STATUS=0
RHS_INCREMENT_BOUNDED_STATUS=0
FLUID_RESPONSE_BOUNDED_STATUS=0
NO_NAN_INF_STATUS=0
REQUIRED_DIAGNOSTICS_STATUS=0
STAGE14_8_EVIDENCE_STATUS=0
STAGE14_9_EVIDENCE_STATUS=0
STAGE14_10_EVIDENCE_STATUS=0
NP124_EVIDENCE_STATUS=0
IO_RESTART_EVIDENCE_STATUS=0
CONTAMINATION_AUDIT_EVIDENCE_STATUS=0
NO_PRESSURE_MODIFICATION_STATUS=0
NO_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_MODIFICATION_STATUS=0
NO_CHANNEL_FORCING_MODIFICATION_STATUS=0
NO_PRODUCTION_IBM_FORCING_STATUS=0
NO_FEEDBACK_APPLICATION_STATUS=0
NO_TWOWAY_FORCE_STATUS=0
NO_STRUCTURE_ADVANCE_STATUS=0
NO_BENDING_SOLVE_STATUS=0
NO_TENSION_SOLVE_STATUS=0
NO_FIBRE_POSITION_UPDATE_STATUS=0
NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS=0
NO_WALL_CONTACT_STATUS=0
CLOSURE_FILE_GENERATED_STATUS=0

MAX_RHS_INCREMENT=0.0
MAX_FLUID_SIGNATURE_DELTA=0.0
MAX_PARALLEL_RHS_DIFF=0.0
MAX_PARALLEL_FORCE_DENSITY_DIFF=0.0
MAX_RESTART_DELTA=0.0
MAX_IO_SIGNATURE_DELTA=0.0
STATIC_MATCH_COUNT=0

mkdir -p stage14_outputs stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$RUNTIME_LOG" "$STAGE9_DAT" "$STAGE11_DAT" "$STAGE12_DAT" "$STAGE13_DAT" "$STAGE14_DAT"
rm -f "$CLOSURE_FILE"

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

xcompact3d_exe() {
    for exe in "$BUILD_DIR/bin/xcompact3d" "$BUILD_DIR/src/xcompact3d" "$BUILD_DIR/xcompact3d"; do
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

max_value() {
    first=$1
    second=$2
    awk -v a="$first" -v b="$second" 'BEGIN { print ((a+0.0) > (b+0.0)) ? a : b }'
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

require_key_exists_finite() {
    file=$1
    key=$2
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    is_finite_number "$value" || {
        add_reason "${key}_not_finite_in_${file}"
        return 1
    }
    return 0
}

require_real_gt() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value+0.0) > (limit+0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_gt_${limit}_got_${value}_in_${file}"
        return 1
    }
    return 0
}

require_real_le() {
    file=$1
    key=$2
    limit=$3
    value=$(get_value "$file" "$key" 2>/dev/null) || {
        add_reason "missing_${key}_in_${file}"
        return 1
    }
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value+0.0) <= (limit+0.0)) exit 0; exit 1 }' || {
        add_reason "${key}_expected_le_${limit}_got_${value}_in_${file}"
        return 1
    }
    return 0
}

scan_active_use_call_forbidden() {
    file=$1
    awk '
      BEGIN { bad=0 }
      /^[[:space:]]*!/ { next }
      {
        line=tolower($0)
        if (line ~ /^[[:space:]]*(use|call)[[:space:]]+/) {
          if (line ~ /(ibm|immersed|structure|advance|bending|bend|tension|wall|contact|position_update|velocity_update|fsi)/) {
            print FILENAME ":" FNR ":" $0
            bad=1
          }
        }
      }
      END { exit(bad ? 1 : 0) }
    ' "$file" >> "$STATIC_REPORT"
}

verify_small_lambda_value() {
    is_finite_number "$STAGE14_11_SMALL_LAMBDA" || {
        add_reason "stage14_11_small_lambda_not_finite"
        return 1
    }
    awk -v value="$STAGE14_11_SMALL_LAMBDA" 'BEGIN { v=value+0.0; if (v>0.0 && v<=1.0e-4) exit 0; exit 1 }' || {
        add_reason "stage14_11_small_lambda_not_small_positive"
        return 1
    }
    return 0
}

verify_np_values() {
    status=0
    awk -v value="$STAGE14_11_NP" 'BEGIN { v=value+0; if (v==1 || v==2 || v==4) exit 0; exit 1 }' || {
        add_reason "stage14_11_np_must_be_1_2_or_4"
        status=1
    }
    seen1=0; seen2=0; seen4=0
    for np in $STAGE14_11_NP_LIST; do
        case "$np" in
            1) seen1=1 ;;
            2) seen2=1 ;;
            4) seen4=1 ;;
            *) add_reason "stage14_11_np_list_contains_unsupported_${np}"; status=1 ;;
        esac
    done
    if [ "$seen1" != "1" ] || [ "$seen2" != "1" ] || [ "$seen4" != "1" ]; then
        add_reason "stage14_11_np_list_must_include_1_2_4"
        status=1
    fi
    return $status
}

run_static_audit() {
    status=0
    static_note "Stage 14.11 static regression audit"

    if grep -nE 'stage14_get_injection_gain[[:space:]]*\([[:space:]]*\)[[:space:]]*==[[:space:]]*0\.0(_[[:alnum:]_]+)?' src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "forbidden_stage14_lambda_zero_registration_gate_found"
        status=1
    else
        FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS=1
        SMALL_LAMBDA_REGISTRATION_STATUS=1
    fi

    if grep -q 'stage14_production_rhs_injection_apply' src/xcompact3d.f90 && \
       grep -q 'stage14_requested' src/xcompact3d.f90 && \
       grep -q 'stage14_rhs_injection_enabled' src/xcompact3d.f90 && \
       grep -q 'stage14_require_stage13' src/xcompact3d.f90 && \
       grep -q 'stage13_force_density_reg' src/xcompact3d.f90; then
        HOOK_CONNECTED_STATUS=1
    else
        add_reason "stage14_rhs_hook_connection_or_dependencies_missing"
        status=1
    fi

    if grep -q 'fibre_stage11_5_production_oneway_hook.dat' src/fibre_stage11_production_oneway_hook.f90 && \
       grep -q 'stage11_5_production_oneway_hook_status' src/fibre_stage11_production_oneway_hook.f90; then
        STAGE11_DIAGNOSTIC_PATH_STATUS=1
    else
        add_reason "stage11_5_production_oneway_diagnostic_path_missing"
        status=1
    fi

    if grep -q 'fibre_stage12_6_production_feedback_candidate.dat' src/fibre_stage12_production_feedback_candidate.f90 && \
       grep -q 'stage12_6_production_feedback_candidate_status' src/fibre_stage12_production_feedback_candidate.f90; then
        STAGE12_DIAGNOSTIC_PATH_STATUS=1
    else
        add_reason "stage12_6_production_feedback_candidate_diagnostic_path_missing"
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

    static_matches=0
    for file in src/poisson.f90 src/navier.f90 src/time_integrators.f90 src/derive.f90 src/Case-Channel.f90; do
        if [ -f "$file" ]; then
            count=$(grep -E 'fibre_stage14|stage14_' "$file" 2>/dev/null | wc -l | awk '{print $1}')
            static_matches=$((static_matches + count))
            if [ "$count" -gt 0 ]; then
                grep -nE 'fibre_stage14|stage14_' "$file" >> "$STATIC_REPORT" 2>/dev/null || true
            fi
        fi
    done
    STATIC_MATCH_COUNT=$static_matches
    if [ "$STATIC_MATCH_COUNT" -eq 0 ]; then
        NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=1
        NO_RK3_CHANNEL_FORCING_STATIC_STATUS=1
    else
        add_reason "stage14_static_matches_in_pressure_projection_poisson_rk3_channel_files_${STATIC_MATCH_COUNT}"
        status=1
    fi

    if scan_active_use_call_forbidden src/fibre_stage14_production_rhs_injection.f90; then
        NO_PRODUCTION_IBM_STATIC_STATUS=1
        NO_STRUCTURE_STATIC_STATUS=1
    else
        add_reason "stage14_forbidden_ibm_or_structure_use_call_found"
        status=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

prepare_input() {
    input_file="$OUTPUT_DIR/stage14_11_total_smoke_input_np${STAGE14_11_NP}.i3d"
    awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); print line }' \
        "$CHANNEL_INPUT" > "$input_file"
    echo "$input_file"
}

run_stage_chain() {
    lambda=$1
    exe=$2
    input=$3
    timeout "$STAGE14_11_TIMEOUT_SEC" env \
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
        X3D_STAGE14_RHS_INJECTION=1 \
        X3D_STAGE14_INJECTION_GAIN="$lambda" \
        X3D_STAGE14_MAX_STEPS=3 \
        X3D_STAGE14_REQUIRE_STAGE13=1 \
        X3D_STAGE14_DIAGNOSTIC_ONLY=1 \
        STAGE9_SKIP_PREREQS=1 \
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
        X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
        X3D_STAGE9_9_MAX_STEPS="$STAGE14_11_MAX_STEPS" \
        "$MPIEXEC" $MPIEXEC_FLAGS -np "$STAGE14_11_NP" "$exe" "$input"
}

run_total_smoke() {
    exe=$(xcompact3d_exe) || {
        add_reason "missing_xcompact3d_executable"
        return 1
    }
    input=$(prepare_input)
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat \
          stage11_outputs/fibre_stage11_5_production_oneway_hook.dat \
          stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat \
          stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat \
          stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
    if ! run_stage_chain "$STAGE14_11_SMALL_LAMBDA" "$exe" "$input" > "$RUNTIME_LOG" 2>&1; then
        add_reason "stage14_11_total_smoke_runtime_failed"
        if [ -f "$RUNTIME_LOG" ]; then
            tail -n 120 "$RUNTIME_LOG"
        fi
        return 1
    fi
    if ! grep -q 'STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS' "$RUNTIME_LOG"; then
        add_reason "stage14_11_runtime_missing_stage9_9_pass_line"
        return 1
    fi
    require_file stage9_outputs/fibre_stage9_9_parallel_consistency.dat "missing_stage9_9_parallel_consistency_dat" || return 1
    require_file stage11_outputs/fibre_stage11_5_production_oneway_hook.dat "missing_stage11_5_production_oneway_hook_diagnostics" || return 1
    require_file stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat "missing_stage12_6_production_feedback_candidate_diagnostics" || return 1
    require_file stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat "missing_stage13_6_production_force_density_candidate_diagnostics" || return 1
    require_file stage14_outputs/fibre_stage14_5_production_rhs_hook.dat "missing_stage14_5_production_rhs_hook_diagnostics" || return 1
    cp stage9_outputs/fibre_stage9_9_parallel_consistency.dat "$STAGE9_DAT"
    cp stage11_outputs/fibre_stage11_5_production_oneway_hook.dat "$STAGE11_DAT"
    cp stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat "$STAGE12_DAT"
    cp stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat "$STAGE13_DAT"
    cp stage14_outputs/fibre_stage14_5_production_rhs_hook.dat "$STAGE14_DAT"
    RUNTIME_TOTAL_SMOKE_STATUS=1
    REQUIRED_DIAGNOSTICS_STATUS=1
    return 0
}

verify_stage9_dat() {
    status=0
    require_key_value "$STAGE9_DAT" stage9_9_parallel_consistency_local_status 1 || status=1
    require_key_value "$STAGE9_DAT" stage9_9_decomposition_invariant_initial_state_status 1 || status=1
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        require_key_exists_finite "$STAGE9_DAT" "$metric" || status=1
    done
    return $status
}

verify_stage11_dat() {
    status=0
    require_key_value "$STAGE11_DAT" stage11_5_requested_flag 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_hook_initialized_status 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_hook_sample_called_status 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_no_ibm_spreading_status 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_no_feedback_force_status 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE11_DAT" stage11_5_production_oneway_hook_status 1 || status=1
    if [ "$status" = "0" ]; then
        STAGE11_ACTIVE_STATUS=1
    fi
    return $status
}

verify_stage12_dat() {
    status=0
    require_key_value "$STAGE12_DAT" stage12_6_requested_flag 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_hook_initialized_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_hook_sample_called_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_force_candidate_computed_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_force_candidate_finite_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_no_rhs_injection_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_no_ibm_spreading_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_no_feedback_application_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE12_DAT" stage12_6_production_feedback_candidate_status 1 || status=1
    if [ "$status" = "0" ]; then
        STAGE12_ACTIVE_STATUS=1
    fi
    return $status
}

verify_stage13_dat() {
    status=0
    require_key_value "$STAGE13_DAT" stage13_6_requested_flag 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_hook_initialized_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_hook_sample_called_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_force_density_candidate_computed_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_force_density_candidate_finite_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_force_density_norm_finite_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_integrated_force_finite_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_spreading_input_sign_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_wrong_sign_rejection_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_no_rhs_injection_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_no_feedback_application_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE13_DAT" stage13_6_production_force_density_candidate_status 1 || status=1
    require_key_exists_finite "$STAGE13_DAT" stage13_6_force_density_l2 || status=1
    require_key_exists_finite "$STAGE13_DAT" stage13_6_max_abs_force_density || status=1
    if [ "$status" = "0" ]; then
        STAGE13_ACTIVE_STATUS=1
    fi
    return $status
}

verify_stage14_dat() {
    status=0
    require_key_value "$STAGE14_DAT" stage14_5_requested_flag 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_lambda_zero_status 0 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_stage13_dependency_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_stage13_candidate_required_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_rhs_arrays_available_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_rhs_increment_computed_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_feedback_application_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE14_DAT" stage14_5_production_rhs_hook_status 1 || status=1
    require_key_exists_finite "$STAGE14_DAT" stage14_5_rhs_increment_l2 || status=1
    require_key_exists_finite "$STAGE14_DAT" stage14_5_rhs_increment_max_abs || status=1
    require_key_exists_finite "$STAGE14_DAT" stage14_5_rhs_signature_delta_l2 || status=1
    require_real_gt "$STAGE14_DAT" stage14_5_rhs_increment_l2 0.0 || status=1
    require_real_gt "$STAGE14_DAT" stage14_5_rhs_increment_max_abs 0.0 || status=1
    require_real_le "$STAGE14_DAT" stage14_5_rhs_increment_max_abs "$STAGE14_11_MAX_RHS_INCREMENT" || status=1
    require_real_le "$STAGE14_DAT" stage14_5_rhs_signature_delta_l2 "$STAGE14_11_MAX_FLUID_DELTA" || status=1
    max_abs=$(get_value "$STAGE14_DAT" stage14_5_rhs_increment_max_abs 2>/dev/null || echo 0.0)
    fluid_delta=$(get_value "$STAGE14_DAT" stage14_5_rhs_signature_delta_l2 2>/dev/null || echo 0.0)
    MAX_RHS_INCREMENT=$(max_value "$MAX_RHS_INCREMENT" "$max_abs")
    MAX_FLUID_SIGNATURE_DELTA=$(max_value "$MAX_FLUID_SIGNATURE_DELTA" "$fluid_delta")
    if [ "$status" = "0" ]; then
        STAGE14_HOOK_ACTIVE_STATUS=1
        LAMBDA_NONZERO_STATUS=1
        RHS_INCREMENT_NONZERO_STATUS=1
        RHS_INCREMENT_FINITE_STATUS=1
        RHS_INCREMENT_SIGN_CORRECT_STATUS=1
        RHS_INCREMENT_BOUNDED_STATUS=1
        FLUID_RESPONSE_BOUNDED_STATUS=1
        NO_PRESSURE_MODIFICATION_STATUS=1
        NO_PROJECTION_MODIFICATION_STATUS=1
        NO_POISSON_MODIFICATION_STATUS=1
        NO_RK3_MODIFICATION_STATUS=1
        NO_CHANNEL_FORCING_MODIFICATION_STATUS=1
        NO_PRODUCTION_IBM_FORCING_STATUS=1
        NO_FEEDBACK_APPLICATION_STATUS=1
        NO_TWOWAY_FORCE_STATUS=1
        NO_STRUCTURE_ADVANCE_STATUS=1
        NO_BENDING_SOLVE_STATUS=1
        NO_TENSION_SOLVE_STATUS=1
        NO_FIBRE_POSITION_UPDATE_STATUS=1
        NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS=1
        NO_WALL_CONTACT_STATUS=1
    fi
    return $status
}

verify_no_nan_inf() {
    status=0
    for file in "$STAGE9_DAT" "$STAGE11_DAT" "$STAGE12_DAT" "$STAGE13_DAT" "$STAGE14_DAT" "$RUNTIME_LOG"; do
        if [ -f "$file" ] && grep -Ei '(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)' "$file" >/dev/null 2>&1; then
            add_reason "nan_or_inf_detected_in_${file}"
            status=1
        fi
    done
    if [ "$status" = "0" ]; then
        NO_NAN_INF_STATUS=1
    fi
    return $status
}

verify_prior_stage14_8() {
    dat=stage14_outputs/stage14_8_parallel_small_lambda_response.dat
    if [ "$STAGE14_11_RUN_STAGE14_8" = "1" ]; then
        STAGE14_8_SMALL_LAMBDA="$STAGE14_11_SMALL_LAMBDA" \
        STAGE14_8_MAX_RHS_INCREMENT="$STAGE14_11_MAX_RHS_INCREMENT" \
        STAGE14_8_MAX_FLUID_DELTA="$STAGE14_11_MAX_FLUID_DELTA" \
            bash stage14_checks/run_stage14_8_parallel_small_lambda_response.sh > "$OUTPUT_DIR/stage14_11_stage14_8_prereq.log" 2>&1
        if [ $? -ne 0 ] || ! grep 'STAGE 14.8 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage14_11_stage14_8_prereq.log" >/dev/null 2>&1; then
            add_reason "stage14_8_prerequisite_failed"
            return 1
        fi
    fi
    require_file "$dat" "missing_stage14_8_parallel_small_lambda_response_dat" || return 1
    status=0
    require_key_value "$dat" stage14_8_parallel_small_lambda_response_status 1 || status=1
    require_key_value "$dat" stage14_8_normalized_rhs_parallel_consistency_status 1 || status=1
    require_key_value "$dat" stage14_8_force_density_parallel_consistency_status 1 || status=1
    require_key_value "$dat" stage14_8_final_signature_parallel_consistency_status 1 || status=1
    require_key_value "$dat" stage14_8_no_structure_advance_status 1 || status=1
    require_real_le "$dat" stage14_8_max_normalized_rhs_l2_delta "$STAGE14_11_MAX_PARALLEL_RHS_DIFF" || status=1
    require_real_le "$dat" stage14_8_max_force_density_delta "$STAGE14_11_MAX_PARALLEL_FORCE_DENSITY_DIFF" || status=1
    rhs_delta=$(get_value "$dat" stage14_8_max_normalized_rhs_l2_delta 2>/dev/null || echo 0.0)
    fd_delta=$(get_value "$dat" stage14_8_max_force_density_delta 2>/dev/null || echo 0.0)
    MAX_PARALLEL_RHS_DIFF=$(max_value "$MAX_PARALLEL_RHS_DIFF" "$rhs_delta")
    MAX_PARALLEL_FORCE_DENSITY_DIFF=$(max_value "$MAX_PARALLEL_FORCE_DENSITY_DIFF" "$fd_delta")
    if [ "$status" = "0" ]; then
        STAGE14_8_EVIDENCE_STATUS=1
        NP124_EVIDENCE_STATUS=1
    fi
    return $status
}

verify_prior_stage14_9() {
    dat=stage14_outputs/stage14_9_io_restart_stats_visu_rhs_injection.dat
    if [ "$STAGE14_11_RUN_STAGE14_9" = "1" ]; then
        STAGE14_9_SMALL_LAMBDA="$STAGE14_11_SMALL_LAMBDA" \
        STAGE14_9_MAX_RHS_INCREMENT="$STAGE14_11_MAX_RHS_INCREMENT" \
        STAGE14_9_MAX_FLUID_DELTA="$STAGE14_11_MAX_FLUID_DELTA" \
        STAGE14_9_MAX_RESTART_DELTA="$STAGE14_11_MAX_RESTART_DELTA" \
        STAGE14_9_MAX_IO_SIGNATURE_DELTA="$STAGE14_11_MAX_IO_SIGNATURE_DELTA" \
        STAGE14_9_NP="$STAGE14_11_NP" \
            bash stage14_checks/run_stage14_9_io_restart_stats_visu_rhs_injection.sh > "$OUTPUT_DIR/stage14_11_stage14_9_prereq.log" 2>&1
        if [ $? -ne 0 ] || ! grep 'STAGE 14.9 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage14_11_stage14_9_prereq.log" >/dev/null 2>&1; then
            add_reason "stage14_9_prerequisite_failed"
            return 1
        fi
    fi
    require_file "$dat" "missing_stage14_9_io_restart_stats_visu_rhs_injection_dat" || return 1
    status=0
    require_key_value "$dat" stage14_9_io_restart_stats_visu_rhs_injection_status 1 || status=1
    require_key_value "$dat" stage14_9_stats_compatibility_status 1 || status=1
    require_key_value "$dat" stage14_9_visu_compatibility_status 1 || status=1
    require_key_value "$dat" stage14_9_coarse_io_compatibility_status 1 || status=1
    require_key_value "$dat" stage14_9_restart_compatibility_status 1 || status=1
    require_key_value "$dat" stage14_9_no_structure_advance_status 1 || status=1
    require_real_le "$dat" stage14_9_max_restart_delta "$STAGE14_11_MAX_RESTART_DELTA" || status=1
    require_real_le "$dat" stage14_9_max_io_signature_delta "$STAGE14_11_MAX_IO_SIGNATURE_DELTA" || status=1
    restart_delta=$(get_value "$dat" stage14_9_max_restart_delta 2>/dev/null || echo 0.0)
    io_delta=$(get_value "$dat" stage14_9_max_io_signature_delta 2>/dev/null || echo 0.0)
    MAX_RESTART_DELTA=$(max_value "$MAX_RESTART_DELTA" "$restart_delta")
    MAX_IO_SIGNATURE_DELTA=$(max_value "$MAX_IO_SIGNATURE_DELTA" "$io_delta")
    if [ "$status" = "0" ]; then
        STAGE14_9_EVIDENCE_STATUS=1
        IO_RESTART_EVIDENCE_STATUS=1
    fi
    return $status
}

verify_prior_stage14_10() {
    dat=stage14_outputs/stage14_10_rhs_ibm_structure_contamination_audit.dat
    if [ "$STAGE14_11_RUN_STAGE14_10" = "1" ]; then
        STAGE14_10_SMALL_LAMBDA="$STAGE14_11_SMALL_LAMBDA" \
        STAGE14_10_MAX_RHS_INCREMENT="$STAGE14_11_MAX_RHS_INCREMENT" \
        STAGE14_10_MAX_FLUID_DELTA="$STAGE14_11_MAX_FLUID_DELTA" \
        STAGE14_10_NP="$STAGE14_11_NP" \
            bash stage14_checks/run_stage14_10_rhs_ibm_structure_contamination_audit.sh > "$OUTPUT_DIR/stage14_11_stage14_10_prereq.log" 2>&1
        if [ $? -ne 0 ] || ! grep 'STAGE 14.10 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage14_11_stage14_10_prereq.log" >/dev/null 2>&1; then
            add_reason "stage14_10_prerequisite_failed"
            return 1
        fi
    fi
    require_file "$dat" "missing_stage14_10_rhs_ibm_structure_contamination_audit_dat" || return 1
    status=0
    require_key_value "$dat" stage14_10_rhs_ibm_structure_contamination_audit_status 1 || status=1
    require_key_value "$dat" stage14_10_forbidden_lambda_gate_absent_status 1 || status=1
    require_key_value "$dat" stage14_10_rank0_diagnostic_status 1 || status=1
    require_key_value "$dat" stage14_10_stage13_np_sampling_repair_status 1 || status=1
    require_key_value "$dat" stage14_10_no_pressure_projection_contamination_status 1 || status=1
    require_key_value "$dat" stage14_10_no_poisson_modification_status 1 || status=1
    require_key_value "$dat" stage14_10_no_rk3_channel_forcing_contamination_status 1 || status=1
    require_key_value "$dat" stage14_10_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$dat" stage14_10_no_structure_advance_status 1 || status=1
    require_key_value "$dat" stage14_10_no_bending_solve_status 1 || status=1
    require_key_value "$dat" stage14_10_no_tension_solve_status 1 || status=1
    require_key_value "$dat" stage14_10_no_fibre_position_update_status 1 || status=1
    require_key_value "$dat" stage14_10_no_fibre_velocity_structural_update_status 1 || status=1
    require_key_value "$dat" stage14_10_no_wall_contact_status 1 || status=1
    if [ "$status" = "0" ]; then
        STAGE14_10_EVIDENCE_STATUS=1
        CONTAMINATION_AUDIT_EVIDENCE_STATUS=1
    fi
    return $status
}

write_closure_file() {
    cat > "$CLOSURE_FILE" <<EOF_CLOSED
# Stage 14 Closed

Stage 14 is closed only because the Stage 14.11 total smoke and closure gate passed.

## Objective

Stage 14 validated the controlled RHS-injection evidence path for the flexible-fibre/channel-DNS FSI project without enabling fibre-structure advance or new production physics.

## Controlled RHS formula

\`\`\`
RHS_new = RHS_old + lambda_14 * f_i_cand
\`\`\`

Here \`f_i_cand\` is the Stage 13 Eulerian force-density candidate with the already-audited fibre-on-fluid sign convention.

## Completed Stage 14 substeps

- Stage 14.0: configuration gate.
- Stage 14.1: RHS accumulator gate.
- Stage 14.2: RHS addition formula gate.
- Stage 14.3: RHS sign/scaling audit.
- Stage 14.4: RK timing audit.
- Stage 14.5: production RHS hook diagnostic gate.
- Stage 14.6: lambda=0 invariance gate.
- Stage 14.7: np=1 small-lambda controlled response gate.
- Stage 14.8: np=1/2/4 small-lambda parallel response gate.
- Stage 14.9: restart/statistics/visualization/coarse-I/O compatibility gate.
- Stage 14.10: RHS/IBM/structure contamination audit gate.
- Stage 14.11: total smoke and closure gate.

## Closure evidence

- Small lambda used by Stage 14.11: \`$STAGE14_11_SMALL_LAMBDA\`.
- Runtime np used by Stage 14.11 total smoke: \`$STAGE14_11_NP\`.
- Required np consistency evidence list: \`$STAGE14_11_NP_LIST\`.
- Maximum RHS increment observed by Stage 14.11 total smoke: \`$MAX_RHS_INCREMENT\`.
- Maximum Stage 14 total-smoke fluid signature delta: \`$MAX_FLUID_SIGNATURE_DELTA\`.
- Maximum inherited Stage 14.8 normalized-RHS parallel delta: \`$MAX_PARALLEL_RHS_DIFF\`.
- Maximum inherited Stage 14.8 force-density parallel delta: \`$MAX_PARALLEL_FORCE_DENSITY_DIFF\`.
- Maximum inherited Stage 14.9 restart delta: \`$MAX_RESTART_DELTA\`.
- Maximum inherited Stage 14.9 I/O signature delta: \`$MAX_IO_SIGNATURE_DELTA\`.

## Boundaries preserved

- No pressure/projection/Poisson modification is approved by Stage 14.
- No RK3/channel-forcing modification is approved by Stage 14.
- No production IBM forcing outside the Stage 13/14 diagnostic-to-RHS chain is approved by Stage 14.
- No fibre-structure advance is active yet.
- No bending solve, tension solve, fibre position update, structural velocity update, or wall/contact logic is active yet.

## Next recommended stage

Stage 15: production flexible-fibre structure advance connection.
EOF_CLOSED
    CLOSURE_FILE_GENERATED_STATUS=1
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage14_11_requested_flag 1
stage14_11_build_status $BUILD_STATUS
stage14_11_static_audit_status $STATIC_AUDIT_STATUS
stage14_11_forbidden_lambda_gate_absent_status $FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS
stage14_11_small_lambda_registration_status $SMALL_LAMBDA_REGISTRATION_STATUS
stage14_11_hook_connected_status $HOOK_CONNECTED_STATUS
stage14_11_stage11_diagnostic_path_status $STAGE11_DIAGNOSTIC_PATH_STATUS
stage14_11_stage12_diagnostic_path_status $STAGE12_DIAGNOSTIC_PATH_STATUS
stage14_11_stage13_diagnostic_path_status $STAGE13_DIAGNOSTIC_PATH_STATUS
stage14_11_stage14_diagnostic_path_status $STAGE14_DIAGNOSTIC_PATH_STATUS
stage14_11_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage14_11_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage14_11_no_pressure_projection_poisson_static_status $NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS
stage14_11_no_rk3_channel_forcing_static_status $NO_RK3_CHANNEL_FORCING_STATIC_STATUS
stage14_11_no_production_ibm_static_status $NO_PRODUCTION_IBM_STATIC_STATUS
stage14_11_no_structure_static_status $NO_STRUCTURE_STATIC_STATUS
stage14_11_runtime_total_smoke_status $RUNTIME_TOTAL_SMOKE_STATUS
stage14_11_required_diagnostics_status $REQUIRED_DIAGNOSTICS_STATUS
stage14_11_stage11_oneway_active_status $STAGE11_ACTIVE_STATUS
stage14_11_stage12_feedback_candidate_active_status $STAGE12_ACTIVE_STATUS
stage14_11_stage13_force_density_active_status $STAGE13_ACTIVE_STATUS
stage14_11_stage14_hook_active_status $STAGE14_HOOK_ACTIVE_STATUS
stage14_11_lambda_nonzero_status $LAMBDA_NONZERO_STATUS
stage14_11_rhs_increment_nonzero_status $RHS_INCREMENT_NONZERO_STATUS
stage14_11_rhs_increment_finite_status $RHS_INCREMENT_FINITE_STATUS
stage14_11_rhs_increment_sign_correct_status $RHS_INCREMENT_SIGN_CORRECT_STATUS
stage14_11_rhs_increment_bounded_status $RHS_INCREMENT_BOUNDED_STATUS
stage14_11_fluid_response_bounded_status $FLUID_RESPONSE_BOUNDED_STATUS
stage14_11_no_nan_inf_status $NO_NAN_INF_STATUS
stage14_11_stage14_8_evidence_status $STAGE14_8_EVIDENCE_STATUS
stage14_11_stage14_9_evidence_status $STAGE14_9_EVIDENCE_STATUS
stage14_11_stage14_10_evidence_status $STAGE14_10_EVIDENCE_STATUS
stage14_11_np124_evidence_status $NP124_EVIDENCE_STATUS
stage14_11_io_restart_evidence_status $IO_RESTART_EVIDENCE_STATUS
stage14_11_contamination_audit_evidence_status $CONTAMINATION_AUDIT_EVIDENCE_STATUS
stage14_11_no_pressure_modification_status $NO_PRESSURE_MODIFICATION_STATUS
stage14_11_no_projection_modification_status $NO_PROJECTION_MODIFICATION_STATUS
stage14_11_no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
stage14_11_no_rk3_modification_status $NO_RK3_MODIFICATION_STATUS
stage14_11_no_channel_forcing_modification_status $NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_11_no_production_ibm_forcing_status $NO_PRODUCTION_IBM_FORCING_STATUS
stage14_11_no_feedback_application_status $NO_FEEDBACK_APPLICATION_STATUS
stage14_11_no_twoway_force_status $NO_TWOWAY_FORCE_STATUS
stage14_11_no_structure_advance_status $NO_STRUCTURE_ADVANCE_STATUS
stage14_11_no_bending_solve_status $NO_BENDING_SOLVE_STATUS
stage14_11_no_tension_solve_status $NO_TENSION_SOLVE_STATUS
stage14_11_no_fibre_position_update_status $NO_FIBRE_POSITION_UPDATE_STATUS
stage14_11_no_fibre_velocity_structural_update_status $NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS
stage14_11_no_wall_contact_status $NO_WALL_CONTACT_STATUS
stage14_11_closure_file_generated_status $CLOSURE_FILE_GENERATED_STATUS
stage14_11_total_smoke_closure_status $final_status
stage14_11_lambda_14 $STAGE14_11_SMALL_LAMBDA
stage14_11_np $STAGE14_11_NP
stage14_11_np_list $STAGE14_11_NP_LIST
stage14_11_max_rhs_increment $MAX_RHS_INCREMENT
stage14_11_max_fluid_signature_delta $MAX_FLUID_SIGNATURE_DELTA
stage14_11_max_parallel_rhs_diff $MAX_PARALLEL_RHS_DIFF
stage14_11_max_parallel_force_density_diff $MAX_PARALLEL_FORCE_DENSITY_DIFF
stage14_11_max_restart_delta $MAX_RESTART_DELTA
stage14_11_max_io_signature_delta $MAX_IO_SIGNATURE_DELTA
stage14_11_max_rhs_increment_limit $STAGE14_11_MAX_RHS_INCREMENT
stage14_11_max_fluid_delta_limit $STAGE14_11_MAX_FLUID_DELTA
stage14_11_max_parallel_rhs_diff_limit $STAGE14_11_MAX_PARALLEL_RHS_DIFF
stage14_11_max_parallel_force_density_diff_limit $STAGE14_11_MAX_PARALLEL_FORCE_DENSITY_DIFF
stage14_11_max_restart_delta_limit $STAGE14_11_MAX_RESTART_DELTA
stage14_11_max_io_signature_delta_limit $STAGE14_11_MAX_IO_SIGNATURE_DELTA
stage14_11_static_match_count $STATIC_MATCH_COUNT
EOF_DAT
}

if ! verify_small_lambda_value; then
    BUILD_STATUS=0
fi
if ! verify_np_values; then
    BUILD_STATUS=0
fi

run_static_audit || true

if [ "$BUILD_STATUS" = "1" ]; then
    build_target xcompact3d || {
        BUILD_STATUS=0
        add_reason "build_failed_xcompact3d"
    }
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if run_total_smoke; then
        verify_stage9_dat || true
        verify_stage11_dat || true
        verify_stage12_dat || true
        verify_stage13_dat || true
        verify_stage14_dat || true
        verify_no_nan_inf || true
    fi
else
    add_reason "stage14_11_runtime_skipped_due_to_build_failure"
fi

verify_prior_stage14_8 || true
verify_prior_stage14_9 || true
verify_prior_stage14_10 || true

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && \
   [ "$FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && [ "$SMALL_LAMBDA_REGISTRATION_STATUS" = "1" ] && \
   [ "$HOOK_CONNECTED_STATUS" = "1" ] && [ "$STAGE11_DIAGNOSTIC_PATH_STATUS" = "1" ] && \
   [ "$STAGE12_DIAGNOSTIC_PATH_STATUS" = "1" ] && [ "$STAGE13_DIAGNOSTIC_PATH_STATUS" = "1" ] && \
   [ "$STAGE14_DIAGNOSTIC_PATH_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS" = "1" ] && \
   [ "$NO_RK3_CHANNEL_FORCING_STATIC_STATUS" = "1" ] && [ "$NO_PRODUCTION_IBM_STATIC_STATUS" = "1" ] && \
   [ "$NO_STRUCTURE_STATIC_STATUS" = "1" ] && [ "$RUNTIME_TOTAL_SMOKE_STATUS" = "1" ] && \
   [ "$REQUIRED_DIAGNOSTICS_STATUS" = "1" ] && [ "$STAGE11_ACTIVE_STATUS" = "1" ] && \
   [ "$STAGE12_ACTIVE_STATUS" = "1" ] && [ "$STAGE13_ACTIVE_STATUS" = "1" ] && \
   [ "$STAGE14_HOOK_ACTIVE_STATUS" = "1" ] && [ "$LAMBDA_NONZERO_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_NONZERO_STATUS" = "1" ] && [ "$RHS_INCREMENT_FINITE_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_SIGN_CORRECT_STATUS" = "1" ] && [ "$RHS_INCREMENT_BOUNDED_STATUS" = "1" ] && \
   [ "$FLUID_RESPONSE_BOUNDED_STATUS" = "1" ] && [ "$NO_NAN_INF_STATUS" = "1" ] && \
   [ "$STAGE14_8_EVIDENCE_STATUS" = "1" ] && [ "$STAGE14_9_EVIDENCE_STATUS" = "1" ] && \
   [ "$STAGE14_10_EVIDENCE_STATUS" = "1" ] && [ "$NP124_EVIDENCE_STATUS" = "1" ] && \
   [ "$IO_RESTART_EVIDENCE_STATUS" = "1" ] && [ "$CONTAMINATION_AUDIT_EVIDENCE_STATUS" = "1" ] && \
   [ "$NO_PRESSURE_MODIFICATION_STATUS" = "1" ] && [ "$NO_PROJECTION_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_POISSON_MODIFICATION_STATUS" = "1" ] && [ "$NO_RK3_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_CHANNEL_FORCING_MODIFICATION_STATUS" = "1" ] && [ "$NO_PRODUCTION_IBM_FORCING_STATUS" = "1" ] && \
   [ "$NO_FEEDBACK_APPLICATION_STATUS" = "1" ] && [ "$NO_TWOWAY_FORCE_STATUS" = "1" ] && \
   [ "$NO_STRUCTURE_ADVANCE_STATUS" = "1" ] && [ "$NO_BENDING_SOLVE_STATUS" = "1" ] && \
   [ "$NO_TENSION_SOLVE_STATUS" = "1" ] && [ "$NO_FIBRE_POSITION_UPDATE_STATUS" = "1" ] && \
   [ "$NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS" = "1" ] && [ "$NO_WALL_CONTACT_STATUS" = "1" ]; then
    write_closure_file
    write_output_dat 1
    echo 'STAGE 14.11 TOTAL SMOKE VERDICT: PASS'
    echo 'STAGE 14.11 FINAL VERDICT: PASS'
    echo 'STAGE14_CLOSED=YES'
    rm -f "$REASONS_FILE"
    exit 0
fi

rm -f "$CLOSURE_FILE"
write_output_dat 0
echo 'STAGE 14.11 TOTAL SMOKE VERDICT: FAIL'
echo 'STAGE 14.11 FINAL VERDICT: FAIL'
echo 'STAGE14_CLOSED=NO'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_11_failure'
fi
exit 1
