#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage12_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
STAGE12_DIAG="$OUTPUT_DIR/fibre_stage12_6_production_feedback_candidate.dat"
GATE_DAT="$OUTPUT_DIR/stage12_10_rhs_spreading_structure_contamination_audit.dat"
BUILD_LOG="$OUTPUT_DIR/stage12_10_build.log"
STATIC_LOG="$OUTPUT_DIR/stage12_10_static_audit.log"
STATIC_RESULTS="$OUTPUT_DIR/stage12_10_static_audit_results.tmp"
RUNTIME_LOG="$OUTPUT_DIR/stage12_10_stage9_9_runtime_smoke.log"
HOOK_FILE=src/fibre_stage12_production_feedback_candidate.f90
XCOMPACT_FILE=src/xcompact3d.f90
REASONS=""

ensure_build_dir() {
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        cmake -S . -B "$BUILD_DIR"
    fi
}

add_reason() {
    if [ -z "$REASONS" ]; then
        REASONS="$1"
    else
        REASONS="$REASONS; $1"
    fi
}

get_value() {
    awk -v key="$1" '$1 == key { print $2; found=1 } END { if (!found) exit 1 }' "$2"
}

check_diag_key() {
    key=$1
    expected=$2
    file=$3
    value=$(get_value "$key" "$file" 2>/dev/null)
    if [ "$value" = "$expected" ]; then
        return 0
    fi
    add_reason "$key"
    return 1
}

write_gate() {
    {
        echo "stage12_10_requested_flag $requested_flag"
        echo "stage12_10_build_status $build_status"
        echo "stage12_10_static_audit_status $static_audit_status"
        echo "stage12_10_velocity_intent_readonly_status $velocity_intent_readonly_status"
        echo "stage12_10_no_velocity_write_static_status $no_velocity_write_static_status"
        echo "stage12_10_no_pressure_write_static_status $no_pressure_write_static_status"
        echo "stage12_10_no_rhs_write_static_status $no_rhs_write_static_status"
        echo "stage12_10_no_eulerian_force_density_static_status $no_eulerian_force_density_static_status"
        echo "stage12_10_no_rhs_injection_static_status $no_rhs_injection_static_status"
        echo "stage12_10_no_ibm_spreading_static_status $no_ibm_spreading_static_status"
        echo "stage12_10_no_feedback_application_static_status $no_feedback_application_static_status"
        echo "stage12_10_no_twoway_force_static_status $no_twoway_force_static_status"
        echo "stage12_10_no_structure_advance_static_status $no_structure_advance_static_status"
        echo "stage12_10_runtime_smoke_status $runtime_smoke_status"
        echo "stage12_10_hook_active_status $hook_active_status"
        echo "stage12_10_force_candidate_computed_status $force_candidate_computed_status"
        echo "stage12_10_force_candidate_finite_status $force_candidate_finite_status"
        echo "stage12_10_force_norm_finite_status $force_norm_finite_status"
        echo "stage12_10_power_diagnostics_finite_status $power_diagnostics_finite_status"
        echo "stage12_10_no_field_modification_status $no_field_modification_status"
        echo "stage12_10_no_rhs_modification_status $no_rhs_modification_status"
        echo "stage12_10_no_eulerian_force_density_status $no_eulerian_force_density_status"
        echo "stage12_10_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage12_10_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage12_10_no_feedback_application_status $no_feedback_application_status"
        echo "stage12_10_no_twoway_force_status $no_twoway_force_status"
        echo "stage12_10_no_structure_advance_status $no_structure_advance_status"
        echo "stage12_10_rhs_spreading_structure_contamination_audit_status $rhs_spreading_structure_contamination_audit_status"
    } > "$GATE_DAT"
}

run_static_audit() {
    : > "$STATIC_LOG"
    : > "$STATIC_RESULTS"
    awk -v hook_file="$HOOK_FILE" -v xcompact_file="$XCOMPACT_FILE" '
        BEGIN {
            velocity_intent = 1
            velocity_write = 1
            pressure_write = 1
            rhs_write = 1
            eulerian_force_density = 1
            rhs_injection = 1
            ibm_spreading = 1
            feedback_application = 1
            twoway_force = 1
            structure_advance = 1
            hook_guard = 0
            guarded_calls = 1
            hook_call_count = 0
            offense_count = 0
        }
        function trim(s) {
            sub(/^[ \t\r\n]+/, "", s)
            sub(/[ \t\r\n]+$/, "", s)
            return s
        }
        function strip_fortran_comment(s,    i, c, out) {
            out = ""
            for (i = 1; i <= length(s); i++) {
                c = substr(s, i, 1)
                if (c == "!") break
                out = out c
            }
            return out
        }
        function has_field(s, f,    pattern) {
            pattern = "(^|[^a-z0-9_])" f "([^a-z0-9_]|$)"
            return (s ~ pattern)
        }
        function record(file, line, text, reason) {
            offense_count++
            print "OFFENSE", file, line, reason, text
        }
        function classify_forbidden(file, line, text, reason) {
            if (text ~ /(eulerian_force_density|force_density|fibre_ibm_force_buffer)/) {
                eulerian_force_density = 0
            }
            if (text ~ /(rhs_injection|stage12_[a-z0-9_]*rhs[a-z0-9_]*|fibre_stage5_rhs_wrapper|fibre_stage5_rhs_injection|fibre_stage7_rhs_candidate)/) {
                rhs_injection = 0
            }
            if (text ~ /(fibre_ibm_spreading|fibre_stage5_spreading_candidate|fibre_stage7_force_spreading|spread_force|force_spreading|spreading)/) {
                ibm_spreading = 0
            }
            if (text ~ /(fibre_ibm_feedback|fibre_stage8_feedback_candidate|feedback_application|apply_feedback)/) {
                feedback_application = 0
            }
            if (text ~ /(fibre_stage8_twoway_force_density|twoway|two_way)/) {
                twoway_force = 0
            }
            if (text ~ /(fibre_ibm_structure_coupling|fibre_structure_integrator|fibre_structure_advance|fibre_tension_solver|structure_advance|advance_structure)/) {
                structure_advance = 0
            }
            record(file, line, text, reason)
        }
        function inspect_use_call(file, line, text) {
            if (text ~ /^[ \t]*(use|call)[ \t]/) {
                if (text ~ /(fibre_ibm_spreading|fibre_ibm_feedback|fibre_ibm_force_buffer|fibre_ibm_structure_coupling|fibre_stage8_twoway_force_density|fibre_stage8_feedback_candidate|fibre_stage8_oneway_forcing|fibre_stage7_force_spreading|fibre_stage7_rhs_candidate|fibre_stage5_rhs_wrapper|fibre_stage5_spreading_candidate|fibre_stage5_rhs_injection|fibre_structure_integrator|fibre_structure_advance|fibre_tension_solver|rhs_injection|eulerian_force_density|force_density|spread_force|force_spreading|spreading|feedback_application|apply_feedback|twoway|two_way|structure_advance|advance_structure)/ || text ~ /stage12_[a-z0-9_]*rhs[a-z0-9_]*/) {
                    classify_forbidden(file, line, text, "forbidden_active_use_or_call")
                }
            }
        }
        function inspect_hook_intent(file, line, text) {
            if (file != hook_file) return
            if (text ~ /intent[ \t]*\([ \t]*(inout|out)[ \t]*\)/) {
                if (has_field(text, "ux") || has_field(text, "uy") || has_field(text, "uz")) {
                    velocity_intent = 0
                    record(file, line, text, "production_velocity_dummy_not_intent_in")
                }
                if (has_field(text, "pressure") || has_field(text, "pp3") || has_field(text, "phi") || has_field(text, "dpdyx1")) {
                    pressure_write = 0
                    record(file, line, text, "production_pressure_dummy_writable")
                }
                if (has_field(text, "rhs") || has_field(text, "gx") || has_field(text, "gy") || has_field(text, "gz") || has_field(text, "force") || has_field(text, "fpx") || has_field(text, "fpy") || has_field(text, "fpz") || has_field(text, "eulerian_force") || has_field(text, "force_density")) {
                    rhs_write = 0
                    record(file, line, text, "production_rhs_or_force_dummy_writable")
                }
            }
        }
        function inspect_hook_write(file, line, text,    lhs) {
            if (file != hook_file) return
            if (text !~ /=/ || text ~ /==|=>|<=|>=|\/=/) return
            lhs = text
            sub(/=.*/, "", lhs)
            lhs = trim(lhs)
            sub(/\(.*/, "", lhs)
            lhs = trim(lhs)
            if (lhs == "ux" || lhs == "uy" || lhs == "uz") {
                velocity_write = 0
                record(file, line, text, "production_velocity_assignment_lhs")
            }
            if (lhs == "pressure" || lhs == "pp3" || lhs == "phi" || lhs == "dpdyx1") {
                pressure_write = 0
                record(file, line, text, "production_pressure_assignment_lhs")
            }
            if (lhs == "rhs" || lhs == "gx" || lhs == "gy" || lhs == "gz" || lhs == "force" || lhs == "fpx" || lhs == "fpy" || lhs == "fpz" || lhs == "eulerian_force" || lhs == "force_density") {
                rhs_write = 0
                record(file, line, text, "production_rhs_or_force_assignment_lhs")
            }
        }
        function inspect_xcompact_guard(file, line, text) {
            if (file != xcompact_file) return
            if (text ~ /stage12_feedback_reg[ \t]*=/ && text ~ /stage12_requested[ \t]*\(/ && text ~ /stage12_readonly_mode[ \t]*\(/) {
                hook_guard = 1
            }
            if (text ~ /^[ \t]*(if[ \t]*\(.*\)[ \t]*)?call[ \t]+stage12_production_feedback_candidate_/) {
                hook_call_count++
                if (text !~ /^[ \t]*if[ \t]*\([^)]*stage12_feedback_reg[^)]*\)[ \t]*call[ \t]+stage12_production_feedback_candidate_/) {
                    guarded_calls = 0
                    record(file, line, text, "stage12_hook_call_not_inline_guarded")
                }
            }
        }
        {
            raw = $0
            active = tolower(trim(strip_fortran_comment(raw)))
            if (active == "") next
            inspect_use_call(FILENAME, FNR, active)
            inspect_hook_intent(FILENAME, FNR, active)
            inspect_hook_write(FILENAME, FNR, active)
            inspect_xcompact_guard(FILENAME, FNR, active)
        }
        END {
            if (hook_guard != 1) {
                guarded_calls = 0
                record(xcompact_file, 0, "stage12_feedback_reg guard missing", "stage12_guard_missing")
            }
            if (hook_call_count == 0) {
                guarded_calls = 0
                record(xcompact_file, 0, "stage12 production feedback candidate calls missing", "stage12_hook_calls_missing")
            }
            print "STATUS velocity_intent_readonly_status", velocity_intent
            print "STATUS no_velocity_write_static_status", velocity_write
            print "STATUS no_pressure_write_static_status", pressure_write
            print "STATUS no_rhs_write_static_status", rhs_write
            print "STATUS no_eulerian_force_density_static_status", eulerian_force_density
            print "STATUS no_rhs_injection_static_status", rhs_injection
            print "STATUS no_ibm_spreading_static_status", ibm_spreading
            print "STATUS no_feedback_application_static_status", feedback_application
            print "STATUS no_twoway_force_static_status", twoway_force
            print "STATUS no_structure_advance_static_status", structure_advance
            print "STATUS stage12_hook_guard_static_status", guarded_calls
            print "STATUS offense_count", offense_count
        }
    ' "$XCOMPACT_FILE" "$HOOK_FILE" > "$STATIC_RESULTS"

    while IFS= read -r line; do
        case $line in
            OFFENSE*)
                set -- $line
                file=$2
                number=$3
                reason=$4
                offending=$(printf '%s\n' "$line" | sed 's/^OFFENSE [^ ]* [^ ]* [^ ]* //')
                {
                    echo "file: $file"
                    echo "line number: $number"
                    echo "offending line: $offending"
                    echo "reason: $reason"
                } | tee -a "$STATIC_LOG"
                add_reason "static_${reason}"
                ;;
            "STATUS velocity_intent_readonly_status "*) velocity_intent_readonly_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_velocity_write_static_status "*) no_velocity_write_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_pressure_write_static_status "*) no_pressure_write_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_rhs_write_static_status "*) no_rhs_write_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_eulerian_force_density_static_status "*) no_eulerian_force_density_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_rhs_injection_static_status "*) no_rhs_injection_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_ibm_spreading_static_status "*) no_ibm_spreading_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_feedback_application_static_status "*) no_feedback_application_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_twoway_force_static_status "*) no_twoway_force_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_structure_advance_static_status "*) no_structure_advance_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS stage12_hook_guard_static_status "*) stage12_hook_guard_static_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS offense_count "*) static_offense_count=$(printf '%s\n' "$line" | awk '{print $3}') ;;
        esac
    done < "$STATIC_RESULTS"

    if [ "$velocity_intent_readonly_status" -eq 1 ] && \
       [ "$no_velocity_write_static_status" -eq 1 ] && \
       [ "$no_pressure_write_static_status" -eq 1 ] && \
       [ "$no_rhs_write_static_status" -eq 1 ] && \
       [ "$no_eulerian_force_density_static_status" -eq 1 ] && \
       [ "$no_rhs_injection_static_status" -eq 1 ] && \
       [ "$no_ibm_spreading_static_status" -eq 1 ] && \
       [ "$no_feedback_application_static_status" -eq 1 ] && \
       [ "$no_twoway_force_static_status" -eq 1 ] && \
       [ "$no_structure_advance_static_status" -eq 1 ] && \
       [ "$stage12_hook_guard_static_status" -eq 1 ] && \
       [ "$static_offense_count" -eq 0 ]; then
        static_audit_status=1
    else
        static_audit_status=0
    fi
}

verify_runtime_diag() {
    diag_ok=1
    requested_ok=0
    readonly_ok=0
    initialized_ok=0
    sample_called_ok=0

    if [ ! -f "$STAGE12_DIAG" ]; then
        add_reason "missing_stage12_6_feedback_candidate_diagnostics"
        return 1
    fi

    if check_diag_key stage12_6_requested_flag 1 "$STAGE12_DIAG"; then requested_ok=1; else diag_ok=0; fi
    if check_diag_key stage12_6_readonly_mode_status 1 "$STAGE12_DIAG"; then readonly_ok=1; else diag_ok=0; fi
    if check_diag_key stage12_6_hook_initialized_status 1 "$STAGE12_DIAG"; then initialized_ok=1; else diag_ok=0; fi
    if check_diag_key stage12_6_hook_sample_called_status 1 "$STAGE12_DIAG"; then sample_called_ok=1; else diag_ok=0; fi
    if check_diag_key stage12_6_sampled_velocity_available_status 1 "$STAGE12_DIAG"; then :; else diag_ok=0; fi
    if check_diag_key stage12_6_force_candidate_computed_status 1 "$STAGE12_DIAG"; then force_candidate_computed_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_force_candidate_finite_status 1 "$STAGE12_DIAG"; then force_candidate_finite_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_force_norm_finite_status 1 "$STAGE12_DIAG"; then force_norm_finite_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_power_diagnostics_finite_status 1 "$STAGE12_DIAG"; then power_diagnostics_finite_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_action_reaction_status 1 "$STAGE12_DIAG"; then :; else diag_ok=0; fi
    if check_diag_key stage12_6_pair_power_consistency_status 1 "$STAGE12_DIAG"; then :; else diag_ok=0; fi
    if check_diag_key stage12_6_field_modified_status 0 "$STAGE12_DIAG"; then no_field_modification_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_rhs_modified_status 0 "$STAGE12_DIAG"; then no_rhs_modification_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_no_eulerian_force_density_status 1 "$STAGE12_DIAG"; then no_eulerian_force_density_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_no_rhs_injection_status 1 "$STAGE12_DIAG"; then no_rhs_injection_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_no_ibm_spreading_status 1 "$STAGE12_DIAG"; then no_ibm_spreading_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_no_feedback_application_status 1 "$STAGE12_DIAG"; then no_feedback_application_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_no_twoway_force_status 1 "$STAGE12_DIAG"; then no_twoway_force_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_no_structure_advance_status 1 "$STAGE12_DIAG"; then no_structure_advance_status=1; else diag_ok=0; fi
    if check_diag_key stage12_6_production_feedback_candidate_status 1 "$STAGE12_DIAG"; then :; else diag_ok=0; fi

    if [ "$requested_ok" -eq 1 ] && [ "$readonly_ok" -eq 1 ] && \
       [ "$initialized_ok" -eq 1 ] && [ "$sample_called_ok" -eq 1 ]; then
        hook_active_status=1
    fi

    if [ "$diag_ok" -eq 1 ]; then
        return 0
    fi
    return 1
}

mkdir -p "$OUTPUT_DIR" "$STAGE11_OUTPUT_DIR" "$STAGE9_OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$RUNTIME_LOG"

requested_flag=1
build_status=0
static_audit_status=0
velocity_intent_readonly_status=0
no_velocity_write_static_status=0
no_pressure_write_static_status=0
no_rhs_write_static_status=0
no_eulerian_force_density_static_status=0
no_rhs_injection_static_status=0
no_ibm_spreading_static_status=0
no_feedback_application_static_status=0
no_twoway_force_static_status=0
no_structure_advance_static_status=0
stage12_hook_guard_static_status=0
static_offense_count=1
runtime_smoke_status=0
hook_active_status=0
force_candidate_computed_status=0
force_candidate_finite_status=0
force_norm_finite_status=0
power_diagnostics_finite_status=0
no_field_modification_status=0
no_rhs_modification_status=0
no_eulerian_force_density_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
rhs_spreading_structure_contamination_audit_status=0

if [ "${STAGE12_10_RUN_STAGE12_9:-0}" = "1" ]; then
    if ! bash stage12_checks/run_stage12_9_io_restart_stats_visu_force_candidate.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage12_9_prerequisite_failed"
    fi
fi

if ensure_build_dir >> "$BUILD_LOG" 2>&1; then
    build_ok=1
else
    build_ok=0
    add_reason "cmake_configure_failed"
fi

if [ "$build_ok" -eq 1 ]; then
    for target in \
        xcompact3d \
        fibre_stage11_config_check \
        fibre_stage11_lagrangian_state_check \
        fibre_stage11_grid_metadata_check \
        fibre_stage11_oneway_interpolation_check \
        fibre_stage11_controlled_interpolation_check \
        fibre_stage11_production_oneway_hook_check \
        fibre_stage12_config_check \
        fibre_stage12_force_buffer_check \
        fibre_stage12_prescribed_velocity_check \
        fibre_stage12_feedback_formula_check \
        fibre_stage12_sign_convention_audit_check \
        fibre_stage12_power_diagnostics_check \
        fibre_stage12_production_feedback_candidate_check
    do
        if ! cmake --build "$BUILD_DIR" --target "$target" -j >> "$BUILD_LOG" 2>&1; then
            build_ok=0
            add_reason "build_failed_$target"
        fi
    done
fi
build_status=$build_ok

run_static_audit

if [ "$build_ok" -eq 1 ]; then
    rm -f "$STAGE12_DIAG"
    if X3D_STAGE11_ONEWAY_HOOK=1 \
       X3D_STAGE11_FORCE_READONLY=1 \
       X3D_STAGE11_MAX_POINTS=8 \
       X3D_STAGE11_MAX_STEPS=3 \
       X3D_STAGE12_FEEDBACK_CANDIDATE=1 \
       X3D_STAGE12_FORCE_READONLY=1 \
       X3D_STAGE12_FEEDBACK_GAIN=1.0 \
       X3D_STAGE12_FORCE_SIGN=1 \
       X3D_STAGE12_MAX_POINTS=8 \
       STAGE9_SKIP_PREREQS=1 \
       X3D_STAGE9_9_PARALLEL_CONSISTENCY=1 \
       X3D_STAGE9_9_DETERMINISTIC_INIT=1 \
       X3D_STAGE9_9_MAX_STEPS=3 \
       bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh > "$RUNTIME_LOG" 2>&1; then
        if grep 'STAGE 9.9 FINAL VERDICT: PASS' "$RUNTIME_LOG" >/dev/null 2>&1; then
            runtime_smoke_status=1
        else
            add_reason "stage9_9_final_verdict_not_pass"
        fi
    else
        add_reason "stage9_9_runtime_smoke_failed"
    fi

    if ! verify_runtime_diag; then
        :
    fi
fi

if [ "$build_status" -eq 1 ] && \
   [ "$static_audit_status" -eq 1 ] && \
   [ "$runtime_smoke_status" -eq 1 ] && \
   [ "$hook_active_status" -eq 1 ] && \
   [ "$force_candidate_computed_status" -eq 1 ] && \
   [ "$force_candidate_finite_status" -eq 1 ] && \
   [ "$force_norm_finite_status" -eq 1 ] && \
   [ "$power_diagnostics_finite_status" -eq 1 ] && \
   [ "$no_field_modification_status" -eq 1 ] && \
   [ "$no_rhs_modification_status" -eq 1 ] && \
   [ "$no_eulerian_force_density_status" -eq 1 ] && \
   [ "$no_rhs_injection_status" -eq 1 ] && \
   [ "$no_ibm_spreading_status" -eq 1 ] && \
   [ "$no_feedback_application_status" -eq 1 ] && \
   [ "$no_twoway_force_status" -eq 1 ] && \
   [ "$no_structure_advance_status" -eq 1 ]; then
    rhs_spreading_structure_contamination_audit_status=1
fi

write_gate

if [ "$rhs_spreading_structure_contamination_audit_status" -eq 1 ]; then
    echo "STAGE 12.10 RHS SPREADING STRUCTURE CONTAMINATION AUDIT VERDICT: PASS"
    echo "STAGE 12.10 FINAL VERDICT: PASS"
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage12_10_failure"
fi

echo "STAGE 12.10 RHS SPREADING STRUCTURE CONTAMINATION AUDIT VERDICT: FAIL"
echo "STAGE 12.10 FINAL VERDICT: FAIL"
echo "Reasons: $REASONS"
exit 1
