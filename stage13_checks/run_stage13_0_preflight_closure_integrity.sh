#!/bin/sh

BUILD_DIR=${BUILD_DIR:-build_stage9}
OUTPUT_DIR=stage13_outputs
STAGE12_OUTPUT_DIR=stage12_outputs
STAGE11_OUTPUT_DIR=stage11_outputs
STAGE9_OUTPUT_DIR=stage9_outputs
GATE_DAT="$OUTPUT_DIR/stage13_0_preflight_closure_integrity.dat"
BUILD_LOG="$OUTPUT_DIR/stage13_0_build.log"
STATIC_LOG="$OUTPUT_DIR/stage13_0_xcompact_static.log"
STATIC_RESULTS="$OUTPUT_DIR/stage13_0_xcompact_static.tmp"
XCOMPACT_FILE=src/xcompact3d.f90
STAGE12_CLOSED=stage12_checks/STAGE12_CLOSED.md
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

write_gate() {
    {
        echo "stage13_0_requested_flag $requested_flag"
        echo "stage13_0_build_status $build_status"
        echo "stage13_0_stage12_closed_file_status $stage12_closed_file_status"
        echo "stage13_0_normal_exit_stage11_finalize_status $normal_exit_stage11_finalize_status"
        echo "stage13_0_stage9_7_early_stop_stage11_finalize_status $stage9_7_early_stop_stage11_finalize_status"
        echo "stage13_0_stage9_8_early_stop_stage11_finalize_status $stage9_8_early_stop_stage11_finalize_status"
        echo "stage13_0_stage9_9_early_stop_stage11_finalize_status $stage9_9_early_stop_stage11_finalize_status"
        echo "stage13_0_stage10_finalize_preserved_status $stage10_finalize_preserved_status"
        echo "stage13_0_stage12_finalize_preserved_status $stage12_finalize_preserved_status"
        echo "stage13_0_no_eulerian_force_density_status $no_eulerian_force_density_status"
        echo "stage13_0_no_rhs_injection_status $no_rhs_injection_status"
        echo "stage13_0_no_ibm_spreading_status $no_ibm_spreading_status"
        echo "stage13_0_no_feedback_application_status $no_feedback_application_status"
        echo "stage13_0_no_twoway_force_status $no_twoway_force_status"
        echo "stage13_0_no_structure_advance_status $no_structure_advance_status"
        echo "stage13_0_preflight_closure_integrity_status $preflight_closure_integrity_status"
    } > "$GATE_DAT"
}

run_static_check() {
    : > "$STATIC_LOG"
    : > "$STATIC_RESULTS"
    awk '
        BEGIN {
            normal_stage11 = 0
            s97_stage11 = 0
            s98_stage11 = 0
            s99_stage11 = 0
            stage10 = 0
            stage12 = 0
            no_eulerian = 1
            no_rhs = 1
            no_ibm = 1
            no_feedback = 1
            no_twoway = 1
            no_structure = 1
            offense_count = 0
            context = ""
        }
        function trim(s) { sub(/^[ \t\r\n]+/, "", s); sub(/[ \t\r\n]+$/, "", s); return s }
        function active_line(s,    i,c,out) {
            out = ""
            for (i = 1; i <= length(s); i++) {
                c = substr(s, i, 1)
                if (c == "!") break
                out = out c
            }
            return tolower(trim(out))
        }
        function offense(line, reason) {
            offense_count++
            print "OFFENSE", FNR, reason, line
        }
        function inspect_no_physics(line) {
            if (line ~ /^[ \t]*(use|call)[ \t]/ || line ~ /=/) {
                if (line ~ /(eulerian_force_density|force_density)/) { no_eulerian = 0; offense(line, "eulerian_force_density") }
                if (line ~ /(rhs_injection|stage13_[a-z0-9_]*rhs|rhs[ \t]*<-|rhs[ \t]*=)/) { no_rhs = 0; offense(line, "rhs_injection") }
                if (line ~ /(ibm_spreading|spread_force|force_spreading|spreading)/) { no_ibm = 0; offense(line, "ibm_spreading") }
                if (line ~ /(feedback_application|apply_feedback)/) { no_feedback = 0; offense(line, "feedback_application") }
                if (line ~ /(twoway|two_way)/) { no_twoway = 0; offense(line, "twoway_force") }
                if (line ~ /(structure_advance|advance_structure|fibre_structure_integrator|fibre_tension_solver)/) { no_structure = 0; offense(line, "structure_advance") }
            }
        }
        {
            line = active_line($0)
            if (line == "") next
            inspect_no_physics(line)
            if (line ~ /stage9_8_stop_now/) context = "stage9_8"
            if (line ~ /stage9_9_stop_now/) context = "stage9_9"
            if (line ~ /stage9_7_stop_now/) context = "stage9_7"
            if (line ~ /if[ \t]*\(stage10_reg\)[ \t]*call[ \t]+stage10_hook_finalize/) stage10 = 1
            if (line ~ /if[ \t]*\(stage12_feedback_reg\)[ \t]*call[ \t]+stage12_production_feedback_candidate_finalize/) stage12 = 1
            if (line ~ /if[ \t]*\(stage11_oneway_reg\)[ \t]*call[ \t]+stage11_production_oneway_finalize/) {
                if (context == "stage9_8") s98_stage11 = 1
                else if (context == "stage9_9") s99_stage11 = 1
                else if (context == "stage9_7") s97_stage11 = 1
                else normal_stage11 = 1
            }
            if (line ~ /^call[ \t]+finalise_xcompact3d/) context = ""
        }
        END {
            print "STATUS normal", normal_stage11
            print "STATUS stage9_7", s97_stage11
            print "STATUS stage9_8", s98_stage11
            print "STATUS stage9_9", s99_stage11
            print "STATUS stage10", stage10
            print "STATUS stage12", stage12
            print "STATUS no_eulerian", no_eulerian
            print "STATUS no_rhs", no_rhs
            print "STATUS no_ibm", no_ibm
            print "STATUS no_feedback", no_feedback
            print "STATUS no_twoway", no_twoway
            print "STATUS no_structure", no_structure
            print "STATUS offense_count", offense_count
        }
    ' "$XCOMPACT_FILE" > "$STATIC_RESULTS"

    while IFS= read -r line; do
        case $line in
            OFFENSE*)
                echo "$line" >> "$STATIC_LOG"
                set -- $line
                add_reason "static_$3"
                ;;
            "STATUS normal "*) normal_exit_stage11_finalize_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS stage9_7 "*) stage9_7_early_stop_stage11_finalize_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS stage9_8 "*) stage9_8_early_stop_stage11_finalize_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS stage9_9 "*) stage9_9_early_stop_stage11_finalize_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS stage10 "*) stage10_finalize_preserved_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS stage12 "*) stage12_finalize_preserved_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_eulerian "*) no_eulerian_force_density_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_rhs "*) no_rhs_injection_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_ibm "*) no_ibm_spreading_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_feedback "*) no_feedback_application_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_twoway "*) no_twoway_force_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
            "STATUS no_structure "*) no_structure_advance_status=$(printf '%s\n' "$line" | awk '{print $3}') ;;
        esac
    done < "$STATIC_RESULTS"

    if [ "$normal_exit_stage11_finalize_status" -ne 1 ]; then add_reason "normal_exit_missing_stage11_finalize"; fi
    if [ "$stage9_7_early_stop_stage11_finalize_status" -ne 1 ]; then add_reason "stage9_7_early_stop_missing_stage11_finalize"; fi
    if [ "$stage9_8_early_stop_stage11_finalize_status" -ne 1 ]; then add_reason "stage9_8_early_stop_missing_stage11_finalize"; fi
    if [ "$stage9_9_early_stop_stage11_finalize_status" -ne 1 ]; then add_reason "stage9_9_early_stop_missing_stage11_finalize"; fi
    if [ "$stage10_finalize_preserved_status" -ne 1 ]; then add_reason "stage10_finalize_missing"; fi
    if [ "$stage12_finalize_preserved_status" -ne 1 ]; then add_reason "stage12_finalize_missing"; fi
}

mkdir -p "$OUTPUT_DIR" "$STAGE12_OUTPUT_DIR" "$STAGE11_OUTPUT_DIR" "$STAGE9_OUTPUT_DIR"
: > "$BUILD_LOG"
: > "$STATIC_LOG"

requested_flag=1
build_status=0
stage12_closed_file_status=0
normal_exit_stage11_finalize_status=0
stage9_7_early_stop_stage11_finalize_status=0
stage9_8_early_stop_stage11_finalize_status=0
stage9_9_early_stop_stage11_finalize_status=0
stage10_finalize_preserved_status=0
stage12_finalize_preserved_status=0
no_eulerian_force_density_status=0
no_rhs_injection_status=0
no_ibm_spreading_status=0
no_feedback_application_status=0
no_twoway_force_status=0
no_structure_advance_status=0
preflight_closure_integrity_status=0

if [ -f "$STAGE12_CLOSED" ]; then
    stage12_closed_file_status=1
else
    add_reason "missing_stage12_closed_md"
fi

if [ "${STAGE13_0_RUN_STAGE12_CLOSURE:-0}" = "1" ]; then
    if ! bash stage12_checks/run_stage12_11_total_smoke.sh >> "$BUILD_LOG" 2>&1; then
        add_reason "optional_stage12_11_closure_failed"
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

run_static_check

if [ "$build_status" -eq 1 ] && \
   [ "$stage12_closed_file_status" -eq 1 ] && \
   [ "$normal_exit_stage11_finalize_status" -eq 1 ] && \
   [ "$stage9_7_early_stop_stage11_finalize_status" -eq 1 ] && \
   [ "$stage9_8_early_stop_stage11_finalize_status" -eq 1 ] && \
   [ "$stage9_9_early_stop_stage11_finalize_status" -eq 1 ] && \
   [ "$stage10_finalize_preserved_status" -eq 1 ] && \
   [ "$stage12_finalize_preserved_status" -eq 1 ] && \
   [ "$no_eulerian_force_density_status" -eq 1 ] && \
   [ "$no_rhs_injection_status" -eq 1 ] && \
   [ "$no_ibm_spreading_status" -eq 1 ] && \
   [ "$no_feedback_application_status" -eq 1 ] && \
   [ "$no_twoway_force_status" -eq 1 ] && \
   [ "$no_structure_advance_status" -eq 1 ]; then
    preflight_closure_integrity_status=1
fi

write_gate

if [ "$preflight_closure_integrity_status" -eq 1 ]; then
    echo "STAGE 13.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: PASS"
    echo "STAGE 13.0 FINAL VERDICT: PASS"
    exit 0
fi

if [ -z "$REASONS" ]; then
    REASONS="unknown_stage13_0_failure"
fi

echo "STAGE 13.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: FAIL"
echo "STAGE 13.0 FINAL VERDICT: FAIL"
echo "Reasons: $REASONS"
exit 1
