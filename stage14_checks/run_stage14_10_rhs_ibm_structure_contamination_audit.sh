#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
CHANNEL_INPUT=${CHANNEL_INPUT:-examples/Channel/input.i3d}
STAGE14_10_RUN_STAGE14_9=${STAGE14_10_RUN_STAGE14_9:-0}
STAGE14_10_SMALL_LAMBDA=${STAGE14_10_SMALL_LAMBDA:-1.0e-8}
STAGE14_10_MAX_RHS_INCREMENT=${STAGE14_10_MAX_RHS_INCREMENT:-1.0e-4}
STAGE14_10_MAX_FLUID_DELTA=${STAGE14_10_MAX_FLUID_DELTA:-1.0e-4}
STAGE14_10_MAX_STATIC_MATCHES=${STAGE14_10_MAX_STATIC_MATCHES:-0}
STAGE14_10_NP=${STAGE14_10_NP:-2}
STAGE14_10_TIMEOUT_SEC=${STAGE14_10_TIMEOUT_SEC:-240}
STAGE14_10_MAX_STEPS=${STAGE14_10_MAX_STEPS:-3}

OUTPUT_DIR=stage14_outputs
OUT_DAT=$OUTPUT_DIR/stage14_10_rhs_ibm_structure_contamination_audit.dat
REASONS_FILE=$OUTPUT_DIR/stage14_10_rhs_ibm_structure_contamination_audit_reasons.tmp
RUNTIME_LOG_SMALL=$OUTPUT_DIR/stage14_10_small_lambda_runtime.log
RUNTIME_LOG_ZERO=$OUTPUT_DIR/stage14_10_lambda0_runtime.log
STAGE9_DAT_SMALL=$OUTPUT_DIR/stage14_10_small_lambda_stage9_9.dat
STAGE9_DAT_ZERO=$OUTPUT_DIR/stage14_10_lambda0_stage9_9.dat
STAGE11_DAT_SMALL=$OUTPUT_DIR/stage14_10_small_lambda_stage11_oneway.dat
STAGE12_DAT_SMALL=$OUTPUT_DIR/stage14_10_small_lambda_stage12_feedback_candidate.dat
STAGE13_DAT_SMALL=$OUTPUT_DIR/stage14_10_small_lambda_stage13_force_density.dat
STAGE14_DAT_SMALL=$OUTPUT_DIR/stage14_10_small_lambda_rhs_hook.dat
STAGE14_DAT_ZERO=$OUTPUT_DIR/stage14_10_lambda0_rhs_hook.dat
STATIC_REPORT=$OUTPUT_DIR/stage14_10_static_audit_report.txt

BUILD_STATUS=1
STAGE14_9_PREREQ_STATUS=1
STATIC_AUDIT_STATUS=0
FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS=0
HOOK_CONNECTED_STATUS=0
NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=0
NO_RK3_CHANNEL_FORCING_STATIC_STATUS=0
NO_PRODUCTION_IBM_STATIC_STATUS=0
NO_STRUCTURE_STATIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_NP_SAMPLING_REPAIR_STATUS=0
RUNTIME_LAMBDA0_STATUS=0
RUNTIME_SMALL_LAMBDA_STATUS=0
STAGE11_ACTIVE_STATUS=0
STAGE12_ACTIVE_STATUS=0
STAGE13_ACTIVE_STATUS=0
STAGE14_HOOK_ACTIVE_STATUS=0
LAMBDA_ZERO_NO_CONTAMINATION_STATUS=0
LAMBDA_NONZERO_STATUS=0
RHS_INCREMENT_NONZERO_STATUS=0
RHS_INCREMENT_FINITE_STATUS=0
RHS_INCREMENT_SIGN_CORRECT_STATUS=0
RHS_INCREMENT_BOUNDED_STATUS=0
FLUID_RESPONSE_BOUNDED_STATUS=0
NO_NAN_INF_STATUS=0
NO_PRESSURE_PROJECTION_CONTAMINATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_CHANNEL_FORCING_CONTAMINATION_STATUS=0
NO_PRODUCTION_IBM_FORCING_STATUS=0
NO_FEEDBACK_APPLICATION_STATUS=0
NO_TWOWAY_FORCE_STATUS=0
NO_STRUCTURE_ADVANCE_STATUS=0
NO_BENDING_SOLVE_STATUS=0
NO_TENSION_SOLVE_STATUS=0
NO_FIBRE_POSITION_UPDATE_STATUS=0
NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS=0
NO_WALL_CONTACT_STATUS=0

MAX_RHS_INCREMENT=0.0
MAX_FLUID_SIGNATURE_DELTA=0.0
STATIC_MATCH_COUNT=0

mkdir -p stage14_outputs stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$RUNTIME_LOG_SMALL" "$RUNTIME_LOG_ZERO" "$STAGE9_DAT_SMALL" "$STAGE9_DAT_ZERO" \
      "$STAGE11_DAT_SMALL" "$STAGE12_DAT_SMALL" "$STAGE13_DAT_SMALL" "$STAGE14_DAT_SMALL" "$STAGE14_DAT_ZERO"

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

abs_value() {
    value=$1
    awk -v value="$value" 'BEGIN { v=value+0.0; if (v<0.0) v=-v; print v }'
}

max_value() {
    first=$1
    second=$2
    awk -v a="$first" -v b="$second" 'BEGIN { print ((a+0.0) > (b+0.0)) ? a : b }'
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

verify_small_lambda_value() {
    is_finite_number "$STAGE14_10_SMALL_LAMBDA" || {
        add_reason "stage14_10_small_lambda_not_finite"
        return 1
    }
    awk -v value="$STAGE14_10_SMALL_LAMBDA" 'BEGIN { v=value+0.0; if (v>0.0 && v<=1.0e-4) exit 0; exit 1 }' || {
        add_reason "stage14_10_small_lambda_not_small_positive"
        return 1
    }
    return 0
}

verify_np_value() {
    case "$STAGE14_10_NP" in
        1|2|4) return 0 ;;
        *) add_reason "stage14_10_np_must_be_1_2_or_4"; return 1 ;;
    esac
}

scan_active_use_call_forbidden() {
    label=$1
    shift
    status=0
    for file in "$@"; do
        [ -f "$file" ] || continue
        awk -v file="$file" -v label="$label" '
          function trim(s) { sub(/^[[:space:]]+/, "", s); sub(/[[:space:]]+$/, "", s); return s }
          BEGIN {
            n = split("poisson projection pressure_correct correct_pressure rk3 channel_forcing production_ibm ibm_forcing apply_ibm legacy_ibm structure_advance advance_structure fibre_structure bending tension wall_contact contact structural_velocity fibre_position_update", forbidden, " ")
            bad = 0
          }
          {
            raw = $0
            sub(/!.*/, "", raw)
            stmt = trim(tolower(raw))
            if (stmt !~ /^(use|call)[[:space:]]/) next
            for (i = 1; i <= n; i++) {
              if (index(stmt, forbidden[i]) > 0) {
                printf("%s active forbidden token %s in %s:%d: %s\n", label, forbidden[i], file, FNR, stmt)
                bad = 1
              }
            }
          }
          END { exit(bad ? 1 : 0) }
        ' "$file" >> "$STATIC_REPORT" || status=1
    done
    return $status
}

run_static_audit() {
    status=0
    static_note "Stage 14.10 static audit"

    if grep -nE 'stage14_get_injection_gain[[:space:]]*\([[:space:]]*\)[[:space:]]*==[[:space:]]*0\.0(_[[:alnum:]_]+)?' src/xcompact3d.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "forbidden_stage14_lambda_zero_registration_gate_found"
        status=1
    else
        FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS=1
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
    if [ "$STATIC_MATCH_COUNT" -le "$STAGE14_10_MAX_STATIC_MATCHES" ]; then
        NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS=1
        NO_RK3_CHANNEL_FORCING_STATIC_STATUS=1
    else
        add_reason "stage14_static_matches_in_pressure_projection_poisson_rk3_channel_files_${STATIC_MATCH_COUNT}"
        status=1
    fi

    if scan_active_use_call_forbidden stage14_forbidden_use_call src/fibre_stage14_production_rhs_injection.f90; then
        NO_PRODUCTION_IBM_STATIC_STATUS=1
        NO_STRUCTURE_STATIC_STATUS=1
    else
        add_reason "stage14_forbidden_ibm_or_structure_use_call_found"
        status=1
    fi

    if grep -q 'rank0_write_allowed' src/fibre_stage13_production_force_density_candidate.f90 && \
       { grep -q 'rank0_write_allowed' src/fibre_stage14_production_rhs_injection.f90 || \
         grep -q 'MPI_COMM_RANK' src/fibre_stage14_production_rhs_injection.f90 || \
         grep -q 'nrank' src/fibre_stage14_production_rhs_injection.f90; }; then
        RANK0_DIAGNOSTIC_STATUS=1
    else
        add_reason "rank0_safe_stage13_or_stage14_diagnostic_guard_missing"
        status=1
    fi

    if grep -nE 'i0[[:space:]]*=[[:space:]]*\(lbound\(ux,[[:space:]]*1\)[[:space:]]*\+[[:space:]]*ubound\(ux,[[:space:]]*1\)\)[[:space:]]*/[[:space:]]*2|j0[[:space:]]*=[[:space:]]*\(lbound\(ux,[[:space:]]*2\)[[:space:]]*\+[[:space:]]*ubound\(ux,[[:space:]]*2\)\)[[:space:]]*/[[:space:]]*2|k0[[:space:]]*=[[:space:]]*\(lbound\(ux,[[:space:]]*3\)[[:space:]]*\+[[:space:]]*ubound\(ux,[[:space:]]*3\)\)[[:space:]]*/[[:space:]]*2' src/fibre_stage13_production_force_density_candidate.f90 >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage13_force_density_local_subdomain_center_sampling_detected"
        status=1
    else
        STAGE13_NP_SAMPLING_REPAIR_STATUS=1
    fi

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

prepare_input() {
    label=$1
    output=$OUTPUT_DIR/stage14_10_${label}_input_np${STAGE14_10_NP}.i3d
    awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); print line }' "$CHANNEL_INPUT" > "$output"
    echo "$output"
}

run_stage_chain() {
    timeout "$STAGE14_10_TIMEOUT_SEC" env \
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
        X3D_STAGE14_INJECTION_GAIN="$1" \
        X3D_STAGE14_MAX_STEPS=3 \
        X3D_STAGE14_REQUIRE_STAGE13=1 \
        X3D_STAGE14_DIAGNOSTIC_ONLY=1 \
        STAGE9_SKIP_PREREQS=1 \
        X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
        X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
        X3D_STAGE9_9_MAX_STEPS="$STAGE14_10_MAX_STEPS" \
        "$2" $MPIEXEC_FLAGS -np "$STAGE14_10_NP" "$3" "$4"
}

copy_case_diagnostics() {
    label=$1
    stage9_copy=$2
    stage14_copy=$3
    if [ ! -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat ]; then
        add_reason "missing_stage9_9_parallel_consistency_dat_${label}"
        return 1
    fi
    cp stage9_outputs/fibre_stage9_9_parallel_consistency.dat "$stage9_copy"
    if [ ! -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics_${label}"
        return 1
    fi
    cp stage14_outputs/fibre_stage14_5_production_rhs_hook.dat "$stage14_copy"
    return 0
}

run_case() {
    label=$1
    lambda=$2
    log=$3
    stage9_copy=$4
    stage14_copy=$5
    exe=$(xcompact3d_exe) || {
        add_reason "missing_xcompact3d_executable_${label}"
        return 1
    }
    input=$(prepare_input "$label")
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency.dat \
          stage11_outputs/fibre_stage11_5_production_oneway_hook.dat \
          stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat \
          stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat \
          stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
    run_stage_chain "$lambda" "$MPIEXEC" "$exe" "$input" > "$log" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ]; then
        add_reason "stage14_10_${label}_runtime_failed"
        tail -n 120 "$log"
        return 1
    fi
    if ! grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log" >/dev/null 2>&1 && \
       ! grep 'STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS' "$log" >/dev/null 2>&1; then
        add_reason "missing_stage9_9_pass_verdict_${label}"
        return 1
    fi
    copy_case_diagnostics "$label" "$stage9_copy" "$stage14_copy" || return 1
    if [ "$label" = "small_lambda" ]; then
        if [ -f stage11_outputs/fibre_stage11_5_production_oneway_hook.dat ]; then
            cp stage11_outputs/fibre_stage11_5_production_oneway_hook.dat "$STAGE11_DAT_SMALL"
        else
            add_reason "missing_stage11_5_production_oneway_hook_diagnostics_small_lambda"
            return 1
        fi
        if [ -f stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat ]; then
            cp stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat "$STAGE12_DAT_SMALL"
        else
            add_reason "missing_stage12_6_production_feedback_candidate_diagnostics_small_lambda"
            return 1
        fi
        if [ -f stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat ]; then
            cp stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat "$STAGE13_DAT_SMALL"
        else
            add_reason "missing_stage13_6_force_density_candidate_diagnostics_small_lambda"
            return 1
        fi
    fi
    return 0
}

verify_stage11_stage12_stage13() {
    status=0
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_requested_flag 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_hook_sample_called_status 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_sample_performed_status 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_sampled_velocity_finite_status 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_no_ibm_spreading_status 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_no_feedback_force_status 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE11_DAT_SMALL" stage11_5_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_requested_flag 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_hook_sample_called_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_force_candidate_computed_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_force_candidate_finite_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_no_rhs_injection_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_no_ibm_spreading_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_no_feedback_application_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE12_DAT_SMALL" stage12_6_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_hook_initialized_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_hook_sample_called_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_force_density_candidate_computed_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_force_density_candidate_finite_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_force_density_norm_finite_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_spreading_input_sign_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_wrong_sign_rejection_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_no_feedback_application_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE13_DAT_SMALL" stage13_6_no_structure_advance_status 1 || status=1
    if [ "$status" = "0" ]; then
        STAGE11_ACTIVE_STATUS=1
        STAGE12_ACTIVE_STATUS=1
        STAGE13_ACTIVE_STATUS=1
    fi
    return $status
}

verify_stage14_zero() {
    status=0
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_requested_flag 1 || status=1
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_lambda_zero_status 1 || status=1
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_rhs_increment_zero_status 1 || status=1
    require_key_value "$STAGE14_DAT_ZERO" stage14_5_rhs_unchanged_status 1 || status=1
    require_real_le "$STAGE14_DAT_ZERO" stage14_5_rhs_increment_l2 1.0e-12 || status=1
    require_real_le "$STAGE14_DAT_ZERO" stage14_5_rhs_increment_max_abs 1.0e-12 || status=1
    if [ "$status" = "0" ]; then
        LAMBDA_ZERO_NO_CONTAMINATION_STATUS=1
    fi
    return $status
}

verify_stage14_small() {
    status=0
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_requested_flag 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_stage13_dependency_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_stage13_candidate_required_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_rhs_arrays_available_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_rhs_increment_computed_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_feedback_application_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_twoway_force_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_no_structure_advance_status 1 || status=1
    require_key_value "$STAGE14_DAT_SMALL" stage14_5_production_rhs_hook_status 1 || status=1
    lambda_zero=$(get_value "$STAGE14_DAT_SMALL" stage14_5_lambda_zero_status 2>/dev/null || echo missing)
    if [ "$lambda_zero" != "0" ]; then
        add_reason "stage14_10_small_lambda_not_recorded_as_nonzero"
        status=1
    fi
    gain=$(get_value "$STAGE14_DAT_SMALL" stage14_5_injection_gain 2>/dev/null || echo missing)
    is_finite_number "$gain" || { add_reason "stage14_10_injection_gain_not_finite"; status=1; }
    awk -v got="$gain" -v want="$STAGE14_10_SMALL_LAMBDA" 'BEGIN { d=(got+0.0)-(want+0.0); if (d<0.0) d=-d; tol=1.0e-14 + 1.0e-8*((want<0.0)?-want:want); if (d<=tol) exit 0; exit 1 }' || {
        add_reason "stage14_10_injection_gain_mismatch"
        status=1
    }
    l2=$(get_value "$STAGE14_DAT_SMALL" stage14_5_rhs_increment_l2 2>/dev/null || echo missing)
    max_abs=$(get_value "$STAGE14_DAT_SMALL" stage14_5_rhs_increment_max_abs 2>/dev/null || echo missing)
    is_finite_number "$l2" || { add_reason "stage14_10_rhs_increment_l2_not_finite"; status=1; }
    is_finite_number "$max_abs" || { add_reason "stage14_10_rhs_increment_max_abs_not_finite"; status=1; }
    require_real_gt "$STAGE14_DAT_SMALL" stage14_5_rhs_increment_l2 0.0 || status=1
    require_real_gt "$STAGE14_DAT_SMALL" stage14_5_rhs_increment_max_abs 0.0 || status=1
    require_real_le "$STAGE14_DAT_SMALL" stage14_5_rhs_increment_max_abs "$STAGE14_10_MAX_RHS_INCREMENT" || status=1
    MAX_RHS_INCREMENT=$(max_value "$MAX_RHS_INCREMENT" "$max_abs")
    if [ "$status" = "0" ]; then
        STAGE14_HOOK_ACTIVE_STATUS=1
        LAMBDA_NONZERO_STATUS=1
        RHS_INCREMENT_NONZERO_STATUS=1
        RHS_INCREMENT_FINITE_STATUS=1
        RHS_INCREMENT_SIGN_CORRECT_STATUS=1
        RHS_INCREMENT_BOUNDED_STATUS=1
        NO_PRESSURE_PROJECTION_CONTAMINATION_STATUS=1
        NO_POISSON_MODIFICATION_STATUS=1
        NO_RK3_CHANNEL_FORCING_CONTAMINATION_STATUS=1
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

compare_fluid_signatures() {
    status=0
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        zero=$(get_value "$STAGE9_DAT_ZERO" "$metric" 2>/dev/null || echo missing)
        small=$(get_value "$STAGE9_DAT_SMALL" "$metric" 2>/dev/null || echo missing)
        if ! is_finite_number "$zero" || ! is_finite_number "$small"; then
            add_reason "stage14_10_fluid_signature_not_finite_${metric}"
            status=1
            continue
        fi
        delta=$(awk -v a="$zero" -v b="$small" 'BEGIN { d=(b+0.0)-(a+0.0); if (d<0.0) d=-d; printf "%.17e", d }')
        MAX_FLUID_SIGNATURE_DELTA=$(max_value "$MAX_FLUID_SIGNATURE_DELTA" "$delta")
        awk -v d="$delta" -v limit="$STAGE14_10_MAX_FLUID_DELTA" 'BEGIN { if ((d+0.0) <= (limit+0.0)) exit 0; exit 1 }' || {
            add_reason "stage14_10_fluid_signature_delta_${metric}_exceeds_bound"
            status=1
        }
    done
    if [ "$status" = "0" ]; then
        FLUID_RESPONSE_BOUNDED_STATUS=1
        NO_NAN_INF_STATUS=1
    fi
    return $status
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage14_10_requested_flag 1
stage14_10_build_status $BUILD_STATUS
stage14_10_stage14_9_prereq_status $STAGE14_9_PREREQ_STATUS
stage14_10_static_audit_status $STATIC_AUDIT_STATUS
stage14_10_forbidden_lambda_gate_absent_status $FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS
stage14_10_hook_connected_status $HOOK_CONNECTED_STATUS
stage14_10_no_pressure_projection_poisson_static_status $NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS
stage14_10_no_rk3_channel_forcing_static_status $NO_RK3_CHANNEL_FORCING_STATIC_STATUS
stage14_10_no_production_ibm_static_status $NO_PRODUCTION_IBM_STATIC_STATUS
stage14_10_no_structure_static_status $NO_STRUCTURE_STATIC_STATUS
stage14_10_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage14_10_stage13_np_sampling_repair_status $STAGE13_NP_SAMPLING_REPAIR_STATUS
stage14_10_runtime_lambda0_status $RUNTIME_LAMBDA0_STATUS
stage14_10_runtime_small_lambda_status $RUNTIME_SMALL_LAMBDA_STATUS
stage14_10_stage11_oneway_active_status $STAGE11_ACTIVE_STATUS
stage14_10_stage12_feedback_candidate_active_status $STAGE12_ACTIVE_STATUS
stage14_10_stage13_force_density_active_status $STAGE13_ACTIVE_STATUS
stage14_10_stage14_hook_active_status $STAGE14_HOOK_ACTIVE_STATUS
stage14_10_lambda_zero_no_contamination_status $LAMBDA_ZERO_NO_CONTAMINATION_STATUS
stage14_10_lambda_nonzero_status $LAMBDA_NONZERO_STATUS
stage14_10_rhs_increment_nonzero_status $RHS_INCREMENT_NONZERO_STATUS
stage14_10_rhs_increment_finite_status $RHS_INCREMENT_FINITE_STATUS
stage14_10_rhs_increment_sign_correct_status $RHS_INCREMENT_SIGN_CORRECT_STATUS
stage14_10_rhs_increment_bounded_status $RHS_INCREMENT_BOUNDED_STATUS
stage14_10_fluid_response_bounded_status $FLUID_RESPONSE_BOUNDED_STATUS
stage14_10_no_nan_inf_status $NO_NAN_INF_STATUS
stage14_10_no_pressure_projection_contamination_status $NO_PRESSURE_PROJECTION_CONTAMINATION_STATUS
stage14_10_no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
stage14_10_no_rk3_channel_forcing_contamination_status $NO_RK3_CHANNEL_FORCING_CONTAMINATION_STATUS
stage14_10_no_production_ibm_forcing_status $NO_PRODUCTION_IBM_FORCING_STATUS
stage14_10_no_feedback_application_status $NO_FEEDBACK_APPLICATION_STATUS
stage14_10_no_twoway_force_status $NO_TWOWAY_FORCE_STATUS
stage14_10_no_structure_advance_status $NO_STRUCTURE_ADVANCE_STATUS
stage14_10_no_bending_solve_status $NO_BENDING_SOLVE_STATUS
stage14_10_no_tension_solve_status $NO_TENSION_SOLVE_STATUS
stage14_10_no_fibre_position_update_status $NO_FIBRE_POSITION_UPDATE_STATUS
stage14_10_no_fibre_velocity_structural_update_status $NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS
stage14_10_no_wall_contact_status $NO_WALL_CONTACT_STATUS
stage14_10_rhs_ibm_structure_contamination_audit_status $final_status
stage14_10_lambda_14 $STAGE14_10_SMALL_LAMBDA
stage14_10_np $STAGE14_10_NP
stage14_10_max_rhs_increment $MAX_RHS_INCREMENT
stage14_10_max_fluid_signature_delta $MAX_FLUID_SIGNATURE_DELTA
stage14_10_static_match_count $STATIC_MATCH_COUNT
stage14_10_max_rhs_increment_limit $STAGE14_10_MAX_RHS_INCREMENT
stage14_10_max_fluid_delta_limit $STAGE14_10_MAX_FLUID_DELTA
stage14_10_max_static_matches_limit $STAGE14_10_MAX_STATIC_MATCHES
EOF_DAT
}

if ! verify_small_lambda_value; then
    BUILD_STATUS=0
fi
if ! verify_np_value; then
    BUILD_STATUS=0
fi

run_static_audit || true

if [ "$BUILD_STATUS" = "1" ]; then
    build_target xcompact3d || {
        BUILD_STATUS=0
        add_reason "build_failed_xcompact3d"
    }
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_10_RUN_STAGE14_9" = "1" ]; then
    STAGE14_9_SMALL_LAMBDA="$STAGE14_10_SMALL_LAMBDA" \
    STAGE14_9_MAX_RHS_INCREMENT="$STAGE14_10_MAX_RHS_INCREMENT" \
    STAGE14_9_MAX_FLUID_DELTA="$STAGE14_10_MAX_FLUID_DELTA" \
    STAGE14_9_NP="$STAGE14_10_NP" \
        bash stage14_checks/run_stage14_9_io_restart_stats_visu_rhs_injection.sh > "$OUTPUT_DIR/stage14_10_stage14_9_prereq.log" 2>&1
    if [ $? -ne 0 ] || ! grep 'STAGE 14.9 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage14_10_stage14_9_prereq.log" >/dev/null 2>&1; then
        STAGE14_9_PREREQ_STATUS=0
        add_reason "optional_stage14_9_prerequisite_failed"
    fi
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_9_PREREQ_STATUS" = "1" ]; then
    if run_case lambda0 0.0 "$RUNTIME_LOG_ZERO" "$STAGE9_DAT_ZERO" "$STAGE14_DAT_ZERO"; then
        RUNTIME_LAMBDA0_STATUS=1
        verify_stage14_zero || true
    fi
    if run_case small_lambda "$STAGE14_10_SMALL_LAMBDA" "$RUNTIME_LOG_SMALL" "$STAGE9_DAT_SMALL" "$STAGE14_DAT_SMALL"; then
        RUNTIME_SMALL_LAMBDA_STATUS=1
        verify_stage11_stage12_stage13 || true
        verify_stage14_small || true
    fi
    if [ "$RUNTIME_LAMBDA0_STATUS" = "1" ] && [ "$RUNTIME_SMALL_LAMBDA_STATUS" = "1" ]; then
        compare_fluid_signatures || true
    fi
else
    add_reason "stage14_10_runtime_skipped_due_to_build_or_prereq_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_9_PREREQ_STATUS" = "1" ] && \
   [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$FORBIDDEN_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && \
   [ "$HOOK_CONNECTED_STATUS" = "1" ] && [ "$NO_PRESSURE_PROJECTION_POISSON_STATIC_STATUS" = "1" ] && \
   [ "$NO_RK3_CHANNEL_FORCING_STATIC_STATUS" = "1" ] && [ "$NO_PRODUCTION_IBM_STATIC_STATUS" = "1" ] && \
   [ "$NO_STRUCTURE_STATIC_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_NP_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$RUNTIME_LAMBDA0_STATUS" = "1" ] && \
   [ "$RUNTIME_SMALL_LAMBDA_STATUS" = "1" ] && [ "$STAGE11_ACTIVE_STATUS" = "1" ] && \
   [ "$STAGE12_ACTIVE_STATUS" = "1" ] && [ "$STAGE13_ACTIVE_STATUS" = "1" ] && \
   [ "$STAGE14_HOOK_ACTIVE_STATUS" = "1" ] && [ "$LAMBDA_ZERO_NO_CONTAMINATION_STATUS" = "1" ] && \
   [ "$LAMBDA_NONZERO_STATUS" = "1" ] && [ "$RHS_INCREMENT_NONZERO_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_FINITE_STATUS" = "1" ] && [ "$RHS_INCREMENT_SIGN_CORRECT_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_BOUNDED_STATUS" = "1" ] && [ "$FLUID_RESPONSE_BOUNDED_STATUS" = "1" ] && \
   [ "$NO_NAN_INF_STATUS" = "1" ] && [ "$NO_PRESSURE_PROJECTION_CONTAMINATION_STATUS" = "1" ] && \
   [ "$NO_POISSON_MODIFICATION_STATUS" = "1" ] && [ "$NO_RK3_CHANNEL_FORCING_CONTAMINATION_STATUS" = "1" ] && \
   [ "$NO_PRODUCTION_IBM_FORCING_STATUS" = "1" ] && [ "$NO_FEEDBACK_APPLICATION_STATUS" = "1" ] && \
   [ "$NO_TWOWAY_FORCE_STATUS" = "1" ] && [ "$NO_STRUCTURE_ADVANCE_STATUS" = "1" ] && \
   [ "$NO_BENDING_SOLVE_STATUS" = "1" ] && [ "$NO_TENSION_SOLVE_STATUS" = "1" ] && \
   [ "$NO_FIBRE_POSITION_UPDATE_STATUS" = "1" ] && [ "$NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS" = "1" ] && \
   [ "$NO_WALL_CONTACT_STATUS" = "1" ]; then
    write_output_dat 1
    echo 'STAGE 14.10 RHS/IBM/STRUCTURE CONTAMINATION AUDIT VERDICT: PASS'
    echo 'STAGE 14.10 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 14.10 RHS/IBM/STRUCTURE CONTAMINATION AUDIT VERDICT: FAIL'
echo 'STAGE 14.10 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_10_failure'
fi
exit 1
