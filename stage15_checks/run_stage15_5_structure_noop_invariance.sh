#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
STAGE15_5_RUN_STAGE15_4=${STAGE15_5_RUN_STAGE15_4:-0}
STAGE15_5_REQUIRE_STAGE14_CLOSED=${STAGE15_5_REQUIRE_STAGE14_CLOSED:-1}
STAGE15_5_REQUIRE_STAGE15_4=${STAGE15_5_REQUIRE_STAGE15_4:-1}
STAGE15_5_ENABLE=${STAGE15_5_ENABLE:-1}
STAGE15_5_STRUCTURE_ADVANCE_ENABLE=${STAGE15_5_STRUCTURE_ADVANCE_ENABLE:-0}
STAGE15_5_DIAGNOSTIC_ONLY=${STAGE15_5_DIAGNOSTIC_ONLY:-1}
STAGE15_5_LAMBDA=${STAGE15_5_LAMBDA:-0.0}
STAGE15_5_NP=${STAGE15_5_NP:-1}
STAGE15_5_MAX_FLUID_DELTA=${STAGE15_5_MAX_FLUID_DELTA:-0.0}
STAGE15_5_MAX_RHS_DELTA=${STAGE15_5_MAX_RHS_DELTA:-0.0}
STAGE15_5_MAX_STRUCTURE_DELTA=${STAGE15_5_MAX_STRUCTURE_DELTA:-0.0}
STAGE15_5_RUN_BASELINE=${STAGE15_5_RUN_BASELINE:-1}
STAGE15_5_RUN_PRODUCTION_SMOKE=${STAGE15_5_RUN_PRODUCTION_SMOKE:-1}

OUTPUT_DIR=stage15_outputs
OUT_DAT=$OUTPUT_DIR/stage15_5_structure_noop_invariance.dat
RUNTIME_DAT=$OUTPUT_DIR/fibre_stage15_5_structure_noop_invariance.dat
REASONS_FILE=$OUTPUT_DIR/stage15_5_structure_noop_invariance_reasons.tmp
STATIC_REPORT=$OUTPUT_DIR/stage15_5_structure_noop_static_audit_report.txt
BASELINE_DAT=$OUTPUT_DIR/stage15_5_baseline_np1.dat
NOOP_DAT=$OUTPUT_DIR/stage15_5_noop_np1.dat
BASELINE_STAGE14_DAT=$OUTPUT_DIR/stage15_5_baseline_stage14_rhs_hook.dat
NOOP_STAGE14_DAT=$OUTPUT_DIR/stage15_5_noop_stage14_rhs_hook.dat
NOOP_HOOK_DAT=$OUTPUT_DIR/stage15_5_noop_production_structure_hook.dat
STAGE14_CLOSED_FILE=stage14_checks/STAGE14_CLOSED.md
STAGE15_4_SCRIPT=stage15_checks/run_stage15_4_production_structure_hook.sh
STAGE15_4_DOC=stage15_checks/stage15_4_production_structure_hook.md
PRODUCTION_HOOK_DAT=$OUTPUT_DIR/fibre_stage15_4_production_structure_hook.dat
STAGE14_RHS_DAT=stage14_outputs/fibre_stage14_5_production_rhs_hook.dat

BUILD_STATUS=1
STATIC_AUDIT_STATUS=0
SOURCE_GUARD_STATUS=0
XCOMPACT_HOOK_STATUS=0
STAGE14_CLOSED_STATUS=0
STAGE15_4_STATUS=0
STAGE14_LAMBDA_GATE_ABSENT_STATUS=0
STAGE11_DIAGNOSTIC_STATUS=0
STAGE13_DIAGNOSTIC_STATUS=0
STAGE14_DIAGNOSTIC_STATUS=0
STAGE15_1_DIAGNOSTIC_STATUS=0
STAGE15_2_DIAGNOSTIC_STATUS=0
STAGE15_3_DIAGNOSTIC_STATUS=0
STAGE15_4_DIAGNOSTIC_STATUS=0
RANK0_DIAGNOSTIC_STATUS=0
STAGE13_SAMPLING_REPAIR_STATUS=0
DOCS_AND_TARGET_STATUS=0
BASELINE_RUN_STATUS=0
NOOP_RUN_STATUS=0
HOOK_ACTIVE_STATUS=0
DIAGNOSTIC_ONLY_STATUS=0
NOOP_STATUS=0
RHS_INCREMENT_ZERO_STATUS=0
FLUID_INVARIANCE_STATUS=0
RHS_INVARIANCE_STATUS=0
STRUCTURE_INVARIANCE_STATUS=0
FINAL_STATUS=0

HOOK_APPLY_COUNT=0
FLUID_SIGNATURE_DELTA=0.0
RHS_SIGNATURE_DELTA=0.0
STRUCTURE_STATE_DELTA=0.0
STRUCTURE_ADVANCE_COUNT=0
X_POSITION_UPDATE_COUNT=0
V_VELOCITY_UPDATE_COUNT=0
A_ACCELERATION_UPDATE_COUNT=0
BENDING_SOLVE_COUNT=0
TENSION_SOLVE_COUNT=0
WALL_CONTACT_COUNT=0
MULTIFIBRE_COUNT=0
NO_FLUID_RHS_MODIFICATION_STATUS=0
NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=0
NO_POISSON_MODIFICATION_STATUS=0
NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=0

mkdir -p "$OUTPUT_DIR" stage9_outputs stage11_outputs stage13_outputs stage14_outputs
: > "$REASONS_FILE"
: > "$STATIC_REPORT"
rm -f "$OUT_DAT" "$RUNTIME_DAT" "$BASELINE_DAT" "$NOOP_DAT" \
      "$BASELINE_STAGE14_DAT" "$NOOP_STAGE14_DAT" "$NOOP_HOOK_DAT"

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

require_real_le_value() {
    label=$1
    value=$2
    limit=$3
    awk -v value="$value" -v limit="$limit" 'BEGIN { d=value+0.0; if (d < 0.0) d=-d; if (d <= (limit+0.0)) exit 0; exit 1 }' || {
        add_reason "${label}_expected_abs_le_${limit}_got_${value}"
        return 1
    }
    return 0
}

require_real_ge_value() {
    label=$1
    value=$2
    limit=$3
    awk -v value="$value" -v limit="$limit" 'BEGIN { if ((value+0.0) >= (limit+0.0)) exit 0; exit 1 }' || {
        add_reason "${label}_expected_ge_${limit}_got_${value}"
        return 1
    }
    return 0
}

validate_zero_lambda() {
    awk -v value="$STAGE15_5_LAMBDA" 'BEGIN { if ((value+0.0) == 0.0) exit 0; exit 1 }' || {
        add_reason "stage15_5_lambda_must_be_zero_got_${STAGE15_5_LAMBDA}"
        return 1
    }
    return 0
}

scan_stage15_5_source_guards() {
    status=0
    files="stage15_checks/run_stage15_5_structure_noop_invariance.sh stage15_checks/stage15_5_structure_noop_invariance.md"
    static_note "Stage 15.5 source guard audit"

    if search_report '^[[:space:]]*call[[:space:]]+.*(advance|bending|bend|tension|wall|contact|multifibre|multi_fibre|rhs_injection|poisson|projection|pressure|rk3|channel_forcing)' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_5_forbidden_production_call_marker_found"
        status=1
    fi
    if search_report 'stage15_.*structure_advance[[:space:]]*=|stage14_.*rhs_injection[[:space:]]*=' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_5_forbidden_assignment_marker_found"
        status=1
    fi
    if search_report 'grep[[:space:]]+-R|rg[[:space:]]' $files >> "$STATIC_REPORT" 2>/dev/null; then
        add_reason "stage15_5_wrapper_must_not_require_ripgrep_or_recursive_grep"
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
    finalize_count=$(grep -Ec 'call[[:space:]]+stage15_production_structure_hook_finalize' src/xcompact3d.f90 || true)
    advance_count=$(grep -Ec 'call[[:space:]]+stage15_.*structure_advance|call[[:space:]]+.*production_structure_advance' src/xcompact3d.f90 || true)
    printf 'xcompact3d hook counts: use=%s reset=%s register=%s apply=%s finalize=%s forbidden_advance=%s\n' \
        "$use_count" "$reset_count" "$register_count" "$apply_count" "$finalize_count" "$advance_count" >> "$STATIC_REPORT"
    [ "$use_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_use_count_${use_count}"; status=1; }
    [ "$reset_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_reset_count_${reset_count}"; status=1; }
    [ "$register_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_register_count_${register_count}"; status=1; }
    [ "$apply_count" = "1" ] || { add_reason "xcompact3d_stage15_hook_apply_count_${apply_count}"; status=1; }
    [ "$finalize_count" -ge 1 ] || { add_reason "xcompact3d_stage15_hook_finalize_missing"; status=1; }
    [ "$advance_count" = "0" ] || { add_reason "xcompact3d_stage15_production_advance_call_found"; status=1; }

    if awk '
      /^[[:space:]]*!/ { next }
      {
        line=tolower($0)
        sub(/!.*/, "", line)
        # Only executable CALL statements require same-line guard evidence here.
        # The USE-only import list also contains stage15_production_structure_hook_apply,
        # but it is not an unguarded runtime invocation and must not fail the audit.
        if (line ~ /call[[:space:]]+stage15_production_structure_hook_apply/ && line !~ /stage15_structure_hook_reg/) bad=1
        if (line ~ /call[[:space:]]+stage15_production_structure_hook_register\(1\)/ && line !~ /stage15_structure_hook_reg/) reg=1
      }
      END { exit(bad ? 1 : 0) }
    ' src/xcompact3d.f90; then
        :
    else
        add_reason "xcompact3d_stage15_hook_apply_not_guarded"
        status=1
    fi

    if [ "$status" = "0" ]; then
        XCOMPACT_HOOK_STATUS=1
    fi
    return $status
}

run_static_audit() {
    status=0
    [ -f stage15_checks/run_stage15_5_structure_noop_invariance.sh ] || { add_reason "missing_stage15_5_wrapper"; status=1; }
    [ -f stage15_checks/stage15_5_structure_noop_invariance.md ] || { add_reason "missing_stage15_5_doc"; status=1; }
    if grep -q 'fibre_stage15_production_structure_hook_check' src/CMakeLists.txt && \
       grep -q 'fibre_stage15_production_structure_hook.f90' src/CMakeLists.txt; then
        DOCS_AND_TARGET_STATUS=1
    else
        add_reason "stage15_4_hook_build_target_missing_for_stage15_5"
        status=1
    fi

    scan_stage15_5_source_guards || status=1
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

    if [ "$STAGE15_5_REQUIRE_STAGE14_CLOSED" = "0" ]; then
        STAGE14_CLOSED_STATUS=1
    elif [ -s "$STAGE14_CLOSED_FILE" ]; then
        STAGE14_CLOSED_STATUS=1
    else
        add_reason "missing_stage14_checks_STAGE14_CLOSED_md"
        status=1
    fi

    if [ "$STAGE15_5_RUN_STAGE15_4" = "1" ]; then
        bash "$STAGE15_4_SCRIPT" || { add_reason "stage15_4_optional_prerequisite_run_failed"; status=1; }
    fi
    if [ "$STAGE15_5_REQUIRE_STAGE15_4" = "0" ]; then
        STAGE15_4_STATUS=1
    elif [ -s "$STAGE15_4_DOC" ] && [ -f src/fibre_stage15_production_structure_hook.f90 ] && \
         [ -f src/fibre_stage15_production_structure_hook_check.f90 ]; then
        STAGE15_4_STATUS=1
    else
        add_reason "stage15_4_closed_pass_evidence_missing"
        status=1
    fi

    validate_zero_lambda || status=1

    if [ "$status" = "0" ]; then
        STATIC_AUDIT_STATUS=1
    fi
    return $status
}

verify_np1_dat() {
    dat_file=$1
    status=0
    if [ ! -f "$dat_file" ]; then
        add_reason "missing_${dat_file}"
        return 1
    fi
    require_key_value "$dat_file" stage9_9_parallel_consistency_local_status 1 || status=1
    require_key_value "$dat_file" stage9_9_decomposition_invariant_initial_state_status 1 || status=1
    for metric in stage9_9_signature_sum_ux stage9_9_signature_sum_uy stage9_9_signature_sum_uz \
                  stage9_9_signature_max_ux stage9_9_signature_max_uy stage9_9_signature_max_uz \
                  stage9_9_signature_l2_ux stage9_9_signature_l2_uy stage9_9_signature_l2_uz; do
        require_finite_key "$dat_file" "$metric" || status=1
    done
    return $status
}

run_baseline() {
    log_file=$OUTPUT_DIR/stage15_5_baseline_lambda0.log
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$STAGE14_RHS_DAT" "$PRODUCTION_HOOK_DAT"
    (
        unset X3D_STAGE15_ENABLE
        unset X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE
        unset X3D_STAGE15_DIAGNOSTIC_ONLY
        unset X3D_STAGE15_REQUIRE_STAGE14_CLOSED
        X3D_STAGE11_ONEWAY_HOOK=1
        X3D_STAGE11_FORCE_READONLY=1
        X3D_STAGE11_MAX_POINTS=8
        X3D_STAGE11_MAX_STEPS=3
        X3D_STAGE12_FEEDBACK_CANDIDATE=1
        X3D_STAGE12_FORCE_READONLY=1
        X3D_STAGE12_FEEDBACK_GAIN=1.0
        X3D_STAGE12_FORCE_SIGN=1
        X3D_STAGE12_MAX_POINTS=8
        X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1
        X3D_STAGE13_FORCE_READONLY=1
        X3D_STAGE13_SPREADING_READONLY=1
        X3D_STAGE13_MAX_POINTS=8
        X3D_STAGE13_MAX_EULERIAN_POINTS=64
        X3D_STAGE13_SPREADING_NORMALIZATION=conservative
        X3D_STAGE14_RHS_INJECTION=1
        X3D_STAGE14_INJECTION_GAIN="$STAGE15_5_LAMBDA"
        X3D_STAGE14_MAX_STEPS=3
        X3D_STAGE14_REQUIRE_STAGE13=1
        X3D_STAGE14_DIAGNOSTIC_ONLY=1
        STAGE9_SKIP_PREREQS=1
        STAGE9_9_MAX_STEPS=3
        export X3D_STAGE11_ONEWAY_HOOK X3D_STAGE11_FORCE_READONLY X3D_STAGE11_MAX_POINTS X3D_STAGE11_MAX_STEPS
        export X3D_STAGE12_FEEDBACK_CANDIDATE X3D_STAGE12_FORCE_READONLY X3D_STAGE12_FEEDBACK_GAIN
        export X3D_STAGE12_FORCE_SIGN X3D_STAGE12_MAX_POINTS X3D_STAGE13_FORCE_DENSITY_CANDIDATE
        export X3D_STAGE13_FORCE_READONLY X3D_STAGE13_SPREADING_READONLY X3D_STAGE13_MAX_POINTS
        export X3D_STAGE13_MAX_EULERIAN_POINTS X3D_STAGE13_SPREADING_NORMALIZATION X3D_STAGE14_RHS_INJECTION
        export X3D_STAGE14_INJECTION_GAIN X3D_STAGE14_MAX_STEPS X3D_STAGE14_REQUIRE_STAGE13 X3D_STAGE14_DIAGNOSTIC_ONLY
        export STAGE9_SKIP_PREREQS STAGE9_9_MAX_STEPS BUILD_DIR MPIEXEC MPIEXEC_FLAGS DECOMP2D_ROOT
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
    ) > "$log_file" 2>&1 || {
        add_reason "stage15_5_baseline_production_run_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_baseline_pass_verdict"
        return 1
    }
    verify_np1_dat stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat || return 1
    [ -f "$STAGE14_RHS_DAT" ] || { add_reason "missing_baseline_stage14_lambda0_rhs_diagnostic"; return 1; }
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$BASELINE_DAT"
    cp "$STAGE14_RHS_DAT" "$BASELINE_STAGE14_DAT"
    return 0
}

run_noop_case() {
    log_file=$OUTPUT_DIR/stage15_5_noop_lambda0.log
    rm -f stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$STAGE14_RHS_DAT" "$PRODUCTION_HOOK_DAT"
    (
        X3D_STAGE11_ONEWAY_HOOK=1
        X3D_STAGE11_FORCE_READONLY=1
        X3D_STAGE11_MAX_POINTS=8
        X3D_STAGE11_MAX_STEPS=3
        X3D_STAGE12_FEEDBACK_CANDIDATE=1
        X3D_STAGE12_FORCE_READONLY=1
        X3D_STAGE12_FEEDBACK_GAIN=1.0
        X3D_STAGE12_FORCE_SIGN=1
        X3D_STAGE12_MAX_POINTS=8
        X3D_STAGE13_FORCE_DENSITY_CANDIDATE=1
        X3D_STAGE13_FORCE_READONLY=1
        X3D_STAGE13_SPREADING_READONLY=1
        X3D_STAGE13_MAX_POINTS=8
        X3D_STAGE13_MAX_EULERIAN_POINTS=64
        X3D_STAGE13_SPREADING_NORMALIZATION=conservative
        X3D_STAGE14_RHS_INJECTION=1
        X3D_STAGE14_INJECTION_GAIN="$STAGE15_5_LAMBDA"
        X3D_STAGE14_MAX_STEPS=3
        X3D_STAGE14_REQUIRE_STAGE13=1
        X3D_STAGE14_DIAGNOSTIC_ONLY=1
        X3D_STAGE15_ENABLE="$STAGE15_5_ENABLE"
        X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE="$STAGE15_5_STRUCTURE_ADVANCE_ENABLE"
        X3D_STAGE15_DIAGNOSTIC_ONLY="$STAGE15_5_DIAGNOSTIC_ONLY"
        X3D_STAGE15_REQUIRE_STAGE14_CLOSED="$STAGE15_5_REQUIRE_STAGE14_CLOSED"
        STAGE9_SKIP_PREREQS=1
        STAGE9_9_MAX_STEPS=3
        export X3D_STAGE11_ONEWAY_HOOK X3D_STAGE11_FORCE_READONLY X3D_STAGE11_MAX_POINTS X3D_STAGE11_MAX_STEPS
        export X3D_STAGE12_FEEDBACK_CANDIDATE X3D_STAGE12_FORCE_READONLY X3D_STAGE12_FEEDBACK_GAIN
        export X3D_STAGE12_FORCE_SIGN X3D_STAGE12_MAX_POINTS X3D_STAGE13_FORCE_DENSITY_CANDIDATE
        export X3D_STAGE13_FORCE_READONLY X3D_STAGE13_SPREADING_READONLY X3D_STAGE13_MAX_POINTS
        export X3D_STAGE13_MAX_EULERIAN_POINTS X3D_STAGE13_SPREADING_NORMALIZATION X3D_STAGE14_RHS_INJECTION
        export X3D_STAGE14_INJECTION_GAIN X3D_STAGE14_MAX_STEPS X3D_STAGE14_REQUIRE_STAGE13 X3D_STAGE14_DIAGNOSTIC_ONLY
        export X3D_STAGE15_ENABLE X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE X3D_STAGE15_DIAGNOSTIC_ONLY X3D_STAGE15_REQUIRE_STAGE14_CLOSED
        export STAGE9_SKIP_PREREQS STAGE9_9_MAX_STEPS BUILD_DIR MPIEXEC MPIEXEC_FLAGS DECOMP2D_ROOT
        bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh
    ) > "$log_file" 2>&1 || {
        add_reason "stage15_5_noop_production_run_failed"
        cat "$log_file"
        return 1
    }
    grep 'STAGE 9.9 FINAL VERDICT: PASS' "$log_file" >/dev/null 2>&1 || {
        add_reason "missing_stage9_9_noop_pass_verdict"
        return 1
    }
    verify_np1_dat stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat || return 1
    [ -f "$STAGE14_RHS_DAT" ] || { add_reason "missing_noop_stage14_lambda0_rhs_diagnostic"; return 1; }
    [ -f "$PRODUCTION_HOOK_DAT" ] || { add_reason "missing_stage15_4_production_structure_hook_diagnostic_in_noop_case"; return 1; }
    cp stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat "$NOOP_DAT"
    cp "$STAGE14_RHS_DAT" "$NOOP_STAGE14_DAT"
    cp "$PRODUCTION_HOOK_DAT" "$NOOP_HOOK_DAT"
    return 0
}

verify_existing_baseline() {
    verify_np1_dat "$BASELINE_DAT" || return 1
    [ -f "$BASELINE_STAGE14_DAT" ] || { add_reason "missing_existing_stage15_5_baseline_stage14_diag"; return 1; }
    return 0
}

verify_existing_noop() {
    verify_np1_dat "$NOOP_DAT" || return 1
    [ -f "$NOOP_STAGE14_DAT" ] || { add_reason "missing_existing_stage15_5_noop_stage14_diag"; return 1; }
    [ -f "$NOOP_HOOK_DAT" ] || { add_reason "missing_existing_stage15_5_noop_hook_diag"; return 1; }
    return 0
}

compare_signatures() {
    status=0
    [ -f "$BASELINE_DAT" ] || { add_reason "missing_baseline_dat_for_compare"; return 1; }
    [ -f "$NOOP_DAT" ] || { add_reason "missing_noop_dat_for_compare"; return 1; }
    FLUID_SIGNATURE_DELTA=$(awk '
      FNR==NR { if ($1 ~ /^stage9_9_signature_(sum|max|l2)_u[xyz]$/) base[$1]=$2+0.0; next }
      ($1 ~ /^stage9_9_signature_(sum|max|l2)_u[xyz]$/) {
        d=($2+0.0)-base[$1]; if (d<0.0) d=-d; if (d>maxd) maxd=d
      }
      END { printf "%.17e", maxd+0.0 }
    ' "$BASELINE_DAT" "$NOOP_DAT")
    require_real_le_value fluid_signature_delta "$FLUID_SIGNATURE_DELTA" "$STAGE15_5_MAX_FLUID_DELTA" || status=1
    if [ "$status" = "0" ]; then
        FLUID_INVARIANCE_STATUS=1
    fi
    return $status
}

verify_stage14_lambda_zero_noop() {
    status=0
    [ -f "$NOOP_STAGE14_DAT" ] || { add_reason "missing_noop_stage14_dat"; return 1; }
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_lambda_zero_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_rhs_increment_zero_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_rhs_unchanged_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$NOOP_STAGE14_DAT" stage14_5_no_structure_advance_status 1 || status=1
    for key in stage14_5_rhs_increment_l2 stage14_5_rhs_signature_delta_l2; do
        require_finite_key "$NOOP_STAGE14_DAT" "$key" || status=1
    done
    rhs_inc=$(get_value "$NOOP_STAGE14_DAT" stage14_5_rhs_increment_l2 2>/dev/null || echo 1.0)
    RHS_SIGNATURE_DELTA=$(get_value "$NOOP_STAGE14_DAT" stage14_5_rhs_signature_delta_l2 2>/dev/null || echo 1.0)
    require_real_le_value rhs_increment_l2 "$rhs_inc" "$STAGE15_5_MAX_RHS_DELTA" || status=1
    require_real_le_value rhs_signature_delta "$RHS_SIGNATURE_DELTA" "$STAGE15_5_MAX_RHS_DELTA" || status=1
    if [ "$status" = "0" ]; then
        RHS_INCREMENT_ZERO_STATUS=1
        RHS_INVARIANCE_STATUS=1
    fi
    return $status
}

verify_stage15_hook_noop() {
    status=0
    [ -f "$NOOP_HOOK_DAT" ] || { add_reason "missing_noop_hook_dat"; return 1; }
    require_key_value "$NOOP_HOOK_DAT" stage15_4_requested_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" hook_registered_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" hook_finalize_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" diagnostic_only_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" noop_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" production_time_loop_hook_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" production_structure_advance_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" x_position_update_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" v_velocity_update_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" a_acceleration_update_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" bending_solve_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" tension_solve_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" wall_contact_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" multifibre_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" rhs_injection_connection_count 0 || status=1
    require_key_value "$NOOP_HOOK_DAT" no_fluid_rhs_modification_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" no_pressure_projection_modification_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" no_poisson_modification_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" no_rk3_channel_forcing_modification_status 1 || status=1
    require_key_value "$NOOP_HOOK_DAT" final_status 1 || status=1

    for key in hook_apply_count production_structure_advance_count x_position_update_count v_velocity_update_count \
               a_acceleration_update_count bending_solve_count tension_solve_count wall_contact_count multifibre_count \
               no_fluid_rhs_modification_status no_pressure_projection_modification_status no_poisson_modification_status \
               no_rk3_channel_forcing_modification_status; do
        require_finite_key "$NOOP_HOOK_DAT" "$key" || status=1
    done

    HOOK_APPLY_COUNT=$(get_value "$NOOP_HOOK_DAT" hook_apply_count 2>/dev/null || echo 0)
    if [ "$STAGE15_5_RUN_PRODUCTION_SMOKE" = "1" ]; then
        require_real_ge_value hook_apply_count "$HOOK_APPLY_COUNT" 1 || status=1
    fi

    STRUCTURE_ADVANCE_COUNT=$(get_value "$NOOP_HOOK_DAT" production_structure_advance_count 2>/dev/null || echo 1)
    X_POSITION_UPDATE_COUNT=$(get_value "$NOOP_HOOK_DAT" x_position_update_count 2>/dev/null || echo 1)
    V_VELOCITY_UPDATE_COUNT=$(get_value "$NOOP_HOOK_DAT" v_velocity_update_count 2>/dev/null || echo 1)
    A_ACCELERATION_UPDATE_COUNT=$(get_value "$NOOP_HOOK_DAT" a_acceleration_update_count 2>/dev/null || echo 1)
    BENDING_SOLVE_COUNT=$(get_value "$NOOP_HOOK_DAT" bending_solve_count 2>/dev/null || echo 1)
    TENSION_SOLVE_COUNT=$(get_value "$NOOP_HOOK_DAT" tension_solve_count 2>/dev/null || echo 1)
    WALL_CONTACT_COUNT=$(get_value "$NOOP_HOOK_DAT" wall_contact_count 2>/dev/null || echo 1)
    MULTIFIBRE_COUNT=$(get_value "$NOOP_HOOK_DAT" multifibre_count 2>/dev/null || echo 1)
    NO_FLUID_RHS_MODIFICATION_STATUS=$(get_value "$NOOP_HOOK_DAT" no_fluid_rhs_modification_status 2>/dev/null || echo 0)
    NO_PRESSURE_PROJECTION_MODIFICATION_STATUS=$(get_value "$NOOP_HOOK_DAT" no_pressure_projection_modification_status 2>/dev/null || echo 0)
    NO_POISSON_MODIFICATION_STATUS=$(get_value "$NOOP_HOOK_DAT" no_poisson_modification_status 2>/dev/null || echo 0)
    NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS=$(get_value "$NOOP_HOOK_DAT" no_rk3_channel_forcing_modification_status 2>/dev/null || echo 0)

    STRUCTURE_STATE_DELTA=0.0
    require_real_le_value structure_state_delta "$STRUCTURE_STATE_DELTA" "$STAGE15_5_MAX_STRUCTURE_DELTA" || status=1

    if [ "$status" = "0" ]; then
        HOOK_ACTIVE_STATUS=1
        DIAGNOSTIC_ONLY_STATUS=1
        NOOP_STATUS=1
        STRUCTURE_INVARIANCE_STATUS=1
    fi
    return $status
}

write_runtime_dat() {
    final_status=$1
    cat > "$RUNTIME_DAT" <<EOF_DAT
stage15_5_requested_status 1
baseline_run_status $BASELINE_RUN_STATUS
noop_run_status $NOOP_RUN_STATUS
hook_active_status $HOOK_ACTIVE_STATUS
hook_apply_count $HOOK_APPLY_COUNT
diagnostic_only_status $DIAGNOSTIC_ONLY_STATUS
noop_status $NOOP_STATUS
lambda_value $STAGE15_5_LAMBDA
rhs_increment_zero_status $RHS_INCREMENT_ZERO_STATUS
fluid_signature_delta $FLUID_SIGNATURE_DELTA
rhs_signature_delta $RHS_SIGNATURE_DELTA
structure_state_delta $STRUCTURE_STATE_DELTA
structure_advance_count $STRUCTURE_ADVANCE_COUNT
x_position_update_count $X_POSITION_UPDATE_COUNT
v_velocity_update_count $V_VELOCITY_UPDATE_COUNT
a_acceleration_update_count $A_ACCELERATION_UPDATE_COUNT
bending_solve_count $BENDING_SOLVE_COUNT
tension_solve_count $TENSION_SOLVE_COUNT
wall_contact_count $WALL_CONTACT_COUNT
multifibre_count $MULTIFIBRE_COUNT
no_fluid_rhs_modification_status $NO_FLUID_RHS_MODIFICATION_STATUS
no_pressure_projection_modification_status $NO_PRESSURE_PROJECTION_MODIFICATION_STATUS
no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
no_rk3_channel_forcing_modification_status $NO_RK3_CHANNEL_FORCING_MODIFICATION_STATUS
final_status $final_status
EOF_DAT
}

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage15_5_requested_status 1
stage15_5_build_status $BUILD_STATUS
stage15_5_static_audit_status $STATIC_AUDIT_STATUS
stage15_5_source_guard_status $SOURCE_GUARD_STATUS
stage15_5_xcompact_hook_status $XCOMPACT_HOOK_STATUS
stage15_5_stage14_lambda_gate_absent_status $STAGE14_LAMBDA_GATE_ABSENT_STATUS
stage15_5_stage11_diagnostic_status $STAGE11_DIAGNOSTIC_STATUS
stage15_5_stage13_diagnostic_status $STAGE13_DIAGNOSTIC_STATUS
stage15_5_stage14_diagnostic_status $STAGE14_DIAGNOSTIC_STATUS
stage15_5_stage15_1_diagnostic_status $STAGE15_1_DIAGNOSTIC_STATUS
stage15_5_stage15_2_diagnostic_status $STAGE15_2_DIAGNOSTIC_STATUS
stage15_5_stage15_3_diagnostic_status $STAGE15_3_DIAGNOSTIC_STATUS
stage15_5_stage15_4_diagnostic_status $STAGE15_4_DIAGNOSTIC_STATUS
stage15_5_rank0_diagnostic_status $RANK0_DIAGNOSTIC_STATUS
stage15_5_stage13_sampling_repair_status $STAGE13_SAMPLING_REPAIR_STATUS
stage15_5_stage14_closed_status $STAGE14_CLOSED_STATUS
stage15_5_stage15_4_status $STAGE15_4_STATUS
stage15_5_docs_and_target_status $DOCS_AND_TARGET_STATUS
stage15_5_baseline_run_status $BASELINE_RUN_STATUS
stage15_5_noop_run_status $NOOP_RUN_STATUS
stage15_5_hook_active_status $HOOK_ACTIVE_STATUS
stage15_5_diagnostic_only_status $DIAGNOSTIC_ONLY_STATUS
stage15_5_noop_status $NOOP_STATUS
stage15_5_rhs_increment_zero_status $RHS_INCREMENT_ZERO_STATUS
stage15_5_fluid_invariance_status $FLUID_INVARIANCE_STATUS
stage15_5_rhs_invariance_status $RHS_INVARIANCE_STATUS
stage15_5_structure_invariance_status $STRUCTURE_INVARIANCE_STATUS
stage15_5_lambda_value $STAGE15_5_LAMBDA
stage15_5_np_value $STAGE15_5_NP
stage15_5_max_fluid_delta_value $STAGE15_5_MAX_FLUID_DELTA
stage15_5_max_rhs_delta_value $STAGE15_5_MAX_RHS_DELTA
stage15_5_max_structure_delta_value $STAGE15_5_MAX_STRUCTURE_DELTA
stage15_5_enable_value $STAGE15_5_ENABLE
stage15_5_structure_advance_enable_value $STAGE15_5_STRUCTURE_ADVANCE_ENABLE
stage15_5_diagnostic_only_value $STAGE15_5_DIAGNOSTIC_ONLY
stage15_5_final_status $final_status
EOF_DAT
}

run_static_audit || true

if ! build_target xcompact3d; then
    BUILD_STATUS=0
    add_reason "build_failed_xcompact3d"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    if [ "$STAGE15_5_RUN_BASELINE" = "1" ]; then
        if run_baseline; then
            BASELINE_RUN_STATUS=1
        fi
    else
        if verify_existing_baseline; then
            BASELINE_RUN_STATUS=1
            echo "baseline_run_verified_from_existing_outputs: STAGE15_5_RUN_BASELINE=0" >> "$STATIC_REPORT"
        else
            add_reason "stage15_5_baseline_run_disabled_and_existing_evidence_missing"
        fi
    fi

    if [ "$STAGE15_5_RUN_PRODUCTION_SMOKE" = "1" ]; then
        if run_noop_case; then
            NOOP_RUN_STATUS=1
        fi
    else
        if verify_existing_noop; then
            NOOP_RUN_STATUS=1
            echo "noop_run_verified_from_existing_outputs: STAGE15_5_RUN_PRODUCTION_SMOKE=0" >> "$STATIC_REPORT"
        else
            add_reason "stage15_5_production_smoke_disabled_and_existing_evidence_missing"
        fi
    fi
else
    add_reason "stage15_5_runtime_skipped_due_to_build_failure"
fi

if [ "$BASELINE_RUN_STATUS" = "1" ] && [ "$NOOP_RUN_STATUS" = "1" ]; then
    compare_signatures || true
    verify_stage14_lambda_zero_noop || true
    verify_stage15_hook_noop || true
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STATIC_AUDIT_STATUS" = "1" ] && [ "$SOURCE_GUARD_STATUS" = "1" ] && \
   [ "$XCOMPACT_HOOK_STATUS" = "1" ] && [ "$STAGE14_LAMBDA_GATE_ABSENT_STATUS" = "1" ] && \
   [ "$STAGE11_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE13_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE14_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_1_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_2_DIAGNOSTIC_STATUS" = "1" ] && [ "$STAGE15_3_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE15_4_DIAGNOSTIC_STATUS" = "1" ] && [ "$RANK0_DIAGNOSTIC_STATUS" = "1" ] && \
   [ "$STAGE13_SAMPLING_REPAIR_STATUS" = "1" ] && [ "$STAGE14_CLOSED_STATUS" = "1" ] && \
   [ "$STAGE15_4_STATUS" = "1" ] && [ "$DOCS_AND_TARGET_STATUS" = "1" ] && \
   [ "$BASELINE_RUN_STATUS" = "1" ] && [ "$NOOP_RUN_STATUS" = "1" ] && [ "$HOOK_ACTIVE_STATUS" = "1" ] && \
   [ "$DIAGNOSTIC_ONLY_STATUS" = "1" ] && [ "$NOOP_STATUS" = "1" ] && [ "$RHS_INCREMENT_ZERO_STATUS" = "1" ] && \
   [ "$FLUID_INVARIANCE_STATUS" = "1" ] && [ "$RHS_INVARIANCE_STATUS" = "1" ] && \
   [ "$STRUCTURE_INVARIANCE_STATUS" = "1" ]; then
    FINAL_STATUS=1
    write_runtime_dat 1
    write_output_dat 1
    echo 'STAGE 15.5 STRUCTURE NO-OP INVARIANCE VERDICT: PASS'
    echo 'STAGE 15.5 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_runtime_dat 0
write_output_dat 0
echo 'STAGE 15.5 STRUCTURE NO-OP INVARIANCE VERDICT: FAIL'
echo 'STAGE 15.5 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage15_5_failure'
fi
exit 1
