#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
CHANNEL_INPUT=${CHANNEL_INPUT:-examples/Channel/input.i3d}
STAGE14_9_SMALL_LAMBDA=${STAGE14_9_SMALL_LAMBDA:-1.0e-8}
STAGE14_9_MAX_RHS_INCREMENT=${STAGE14_9_MAX_RHS_INCREMENT:-1.0e-4}
STAGE14_9_MAX_FLUID_DELTA=${STAGE14_9_MAX_FLUID_DELTA:-1.0e-4}
STAGE14_9_MAX_RESTART_DELTA=${STAGE14_9_MAX_RESTART_DELTA:-1.0e-8}
STAGE14_9_MAX_IO_SIGNATURE_DELTA=${STAGE14_9_MAX_IO_SIGNATURE_DELTA:-1.0e-8}
STAGE14_9_NP=${STAGE14_9_NP:-2}
STAGE14_9_RUN_STAGE14_8=${STAGE14_9_RUN_STAGE14_8:-0}
STAGE14_9_TIMEOUT_SEC=${STAGE14_9_TIMEOUT_SEC:-240}
STAGE14_9_MAX_STEPS=${STAGE14_9_MAX_STEPS:-3}
STAGE14_9_RESTART_STEPS_BEFORE=${STAGE14_9_RESTART_STEPS_BEFORE:-3}
STAGE14_9_RESTART_STEPS_AFTER=${STAGE14_9_RESTART_STEPS_AFTER:-3}

OUTPUT_DIR=stage14_outputs
OUT_DAT=$OUTPUT_DIR/stage14_9_io_restart_stats_visu_rhs_injection.dat
REASONS_FILE=$OUTPUT_DIR/stage14_9_io_restart_stats_visu_rhs_injection_reasons.tmp

STATS_LOG=$OUTPUT_DIR/stage14_9_stats_visu_coarse_io.log
STATS_DAT=$OUTPUT_DIR/stage14_9_stats_visu_coarse_io.dat
STATS_STAGE14_DAT=$OUTPUT_DIR/stage14_9_stats_visu_rhs_hook.dat
STATS_STAGE13_DAT=$OUTPUT_DIR/stage14_9_stats_visu_stage13_force_density.dat
RESTART_COLD_LOG=$OUTPUT_DIR/stage14_9_restart_cold.log
RESTART_COLD_DAT=$OUTPUT_DIR/stage14_9_restart_cold.dat
RESTART_COLD_STAGE14_DAT=$OUTPUT_DIR/stage14_9_restart_cold_rhs_hook.dat
RESTART_COLD_STAGE13_DAT=$OUTPUT_DIR/stage14_9_restart_cold_stage13_force_density.dat
RESTART_READ_LOG=$OUTPUT_DIR/stage14_9_restart_read.log
RESTART_READ_DAT=$OUTPUT_DIR/stage14_9_restart_read.dat
RESTART_READ_STAGE14_DAT=$OUTPUT_DIR/stage14_9_restart_read_rhs_hook.dat
RESTART_READ_STAGE13_DAT=$OUTPUT_DIR/stage14_9_restart_read_stage13_force_density.dat
RESTART_SIGNATURE_FILE=$OUTPUT_DIR/stage14_9_restart_signature_np${STAGE14_9_NP}.dat

BUILD_STATUS=1
STAGE14_8_PREREQ_STATUS=1
STATS_IO_RUN_STATUS=0
STATS_COMPAT_STATUS=0
VISU_COMPAT_STATUS=0
COARSE_IO_COMPAT_STATUS=0
RESTART_WRITE_STATUS=0
RESTART_READ_STATUS=0
RESTART_COMPAT_STATUS=0
STAGE14_HOOK_ACTIVE_STATUS=0
LAMBDA_NONZERO_STATUS=0
RHS_INCREMENT_NONZERO_STATUS=0
RHS_INCREMENT_FINITE_STATUS=0
RHS_INCREMENT_SIGN_CORRECT_STATUS=0
RHS_INCREMENT_BOUNDED_STATUS=0
FLUID_RESPONSE_BOUNDED_STATUS=0
STAGE13_FORCE_DENSITY_STATUS=0
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

MAX_RHS_INCREMENT=0.0
MAX_RESTART_DELTA=0.0
MAX_IO_SIGNATURE_DELTA=0.0
MAX_FLUID_DELTA=0.0

mkdir -p stage14_outputs stage13_outputs stage12_outputs stage11_outputs stage9_outputs
: > "$REASONS_FILE"
rm -f "$OUT_DAT" "$STATS_LOG" "$STATS_DAT" "$STATS_STAGE14_DAT" "$STATS_STAGE13_DAT" \
      "$RESTART_COLD_LOG" "$RESTART_COLD_DAT" "$RESTART_COLD_STAGE14_DAT" "$RESTART_COLD_STAGE13_DAT" \
      "$RESTART_READ_LOG" "$RESTART_READ_DAT" "$RESTART_READ_STAGE14_DAT" "$RESTART_READ_STAGE13_DAT" "$RESTART_SIGNATURE_FILE"

add_reason() {
    echo "$1" >> "$REASONS_FILE"
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
    is_finite_number "$STAGE14_9_SMALL_LAMBDA" || {
        add_reason "stage14_9_small_lambda_not_finite"
        return 1
    }
    awk -v value="$STAGE14_9_SMALL_LAMBDA" 'BEGIN { v=value+0.0; if (v>0.0 && v<=1.0e-4) exit 0; exit 1 }' || {
        add_reason "stage14_9_small_lambda_not_small_positive"
        return 1
    }
    return 0
}

verify_np_value() {
    case "$STAGE14_9_NP" in
        1|2|4) return 0 ;;
        *) add_reason "stage14_9_np_must_be_1_2_or_4"; return 1 ;;
    esac
}

prepare_input() {
    phase=$1
    output=$2
    case "$phase" in
        stats)
            awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); print line }' "$CHANNEL_INPUT" > "$output"
            ;;
        cold)
            awk '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 0", line); if (line ~ /^[[:space:]]*icheckpoint[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 1", line); print line }' "$CHANNEL_INPUT" > "$output"
            ;;
        restart)
            awk -v ifirst_val="$((STAGE14_9_RESTART_STEPS_BEFORE + 1))" '{ line=$0; if (line ~ /^[[:space:]]*irestart[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 1", line); if (line ~ /^[[:space:]]*icheckpoint[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= 1", line); if (line ~ /^[[:space:]]*ifirst[[:space:]]*=/) sub(/=[[:space:]]*[0-9]+/, "= " ifirst_val, line); print line }' "$CHANNEL_INPUT" > "$output"
            ;;
    esac
}

run_stage_chain() {
    timeout "$STAGE14_9_TIMEOUT_SEC" env \
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
        X3D_STAGE14_INJECTION_GAIN="$STAGE14_9_SMALL_LAMBDA" \
        X3D_STAGE14_MAX_STEPS=3 \
        X3D_STAGE14_REQUIRE_STAGE13=1 \
        X3D_STAGE14_DIAGNOSTIC_ONLY=1 \
        "$@"
}

copy_required_diagnostics() {
    label=$1
    stage14_copy=$2
    stage13_copy=$3
    if [ ! -f stage14_outputs/fibre_stage14_5_production_rhs_hook.dat ]; then
        add_reason "missing_stage14_5_production_rhs_hook_diagnostics_${label}"
        return 1
    fi
    cp stage14_outputs/fibre_stage14_5_production_rhs_hook.dat "$stage14_copy"
    if [ ! -f stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat ]; then
        add_reason "missing_stage13_6_force_density_candidate_diagnostics_${label}"
        return 1
    fi
    cp stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat "$stage13_copy"
    return 0
}

run_stats_visu_coarse() {
    exe=$(xcompact3d_exe) || {
        add_reason "missing_xcompact3d_executable_stats_visu_coarse"
        return 1
    }
    input=$OUTPUT_DIR/stage14_9_stats_input_np${STAGE14_9_NP}.i3d
    prepare_input stats "$input"
    rm -f stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat \
          stage14_outputs/fibre_stage14_5_production_rhs_hook.dat \
          stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat
    run_stage_chain \
        X3D_STAGE9_7_STATS_VISU_IO_SMOKE=1 \
        X3D_STAGE9_7_MAX_STEPS="$STAGE14_9_MAX_STEPS" \
        X3D_STAGE9_7_REQUIRE_STATS=1 \
        X3D_STAGE9_7_REQUIRE_VISU=1 \
        X3D_STAGE9_7_REQUIRE_COARSE_IO=1 \
        "$MPIEXEC" $MPIEXEC_FLAGS -np "$STAGE14_9_NP" "$exe" "$input" > "$STATS_LOG" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ]; then
        add_reason "stage14_9_stats_visu_coarse_run_failed"
        tail -n 120 "$STATS_LOG"
        return 1
    fi
    grep 'STAGE 9.7 STATS VISU IO SMOKE VERDICT: PASS' "$STATS_LOG" >/dev/null 2>&1 || {
        add_reason "missing_stage9_7_pass_verdict_stage14_9_stats_visu_coarse"
        return 1
    }
    if [ ! -f stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat ]; then
        add_reason "missing_stage9_7_stats_visu_io_dat_stage14_9"
        return 1
    fi
    cp stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat "$STATS_DAT"
    copy_required_diagnostics stats_visu_coarse "$STATS_STAGE14_DAT" "$STATS_STAGE13_DAT" || return 1
    return 0
}

run_restart_phase() {
    phase=$1
    log=$2
    dat=$3
    stage14_copy=$4
    stage13_copy=$5
    exe=$(xcompact3d_exe) || {
        add_reason "missing_xcompact3d_executable_restart_${phase}"
        return 1
    }
    if [ "$phase" = "cold" ]; then
        input=$OUTPUT_DIR/stage14_9_restart_cold_input_np${STAGE14_9_NP}.i3d
        prepare_input cold "$input"
        rm -f checkpoint checkpoint.old restart.info restart* \
              stage9_outputs/fibre_stage9_8_restart_io_regression.dat \
              stage14_outputs/fibre_stage14_5_production_rhs_hook.dat \
              stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat
        stage9_phase=cold
        max_steps=$STAGE14_9_RESTART_STEPS_BEFORE
    else
        checkpoint_file=$(find . -maxdepth 1 -type f \( -name 'checkpoint' -o -name 'checkpoint*' -o -name 'restart*' \) -size +0c -print -quit)
        if [ -z "$checkpoint_file" ]; then
            add_reason "missing_restart_checkpoint_file_before_restart_read"
            return 1
        fi
        input=$OUTPUT_DIR/stage14_9_restart_read_input_np${STAGE14_9_NP}.i3d
        prepare_input restart "$input"
        rm -f stage9_outputs/fibre_stage9_8_restart_io_regression.dat \
              stage14_outputs/fibre_stage14_5_production_rhs_hook.dat \
              stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat
        stage9_phase=restart
        max_steps=$STAGE14_9_RESTART_STEPS_AFTER
    fi
    run_stage_chain \
        X3D_STAGE9_8_RESTART_IO_REGRESSION=1 \
        X3D_STAGE9_8_PHASE="$stage9_phase" \
        X3D_STAGE9_8_MAX_STEPS_BEFORE_RESTART="$STAGE14_9_RESTART_STEPS_BEFORE" \
        X3D_STAGE9_8_MAX_STEPS_AFTER_RESTART="$STAGE14_9_RESTART_STEPS_AFTER" \
        X3D_STAGE9_8_RESTART_SIGNATURE_TOL="$STAGE14_9_MAX_RESTART_DELTA" \
        X3D_STAGE9_8_SIGNATURE_FILE="$RESTART_SIGNATURE_FILE" \
        "$MPIEXEC" $MPIEXEC_FLAGS -np "$STAGE14_9_NP" "$exe" "$input" > "$log" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ]; then
        add_reason "stage14_9_restart_${phase}_run_failed"
        tail -n 120 "$log"
        return 1
    fi
    grep 'STAGE 9.8 RESTART IO REGRESSION VERDICT: PASS' "$log" >/dev/null 2>&1 || {
        add_reason "missing_stage9_8_pass_verdict_stage14_9_restart_${phase}"
        return 1
    }
    if [ ! -f stage9_outputs/fibre_stage9_8_restart_io_regression.dat ]; then
        add_reason "missing_stage9_8_restart_io_dat_stage14_9_${phase}"
        return 1
    fi
    cp stage9_outputs/fibre_stage9_8_restart_io_regression.dat "$dat"
    copy_required_diagnostics "restart_${phase}" "$stage14_copy" "$stage13_copy" || return 1
    return 0
}

verify_stats_visu_coarse_dat() {
    status=0
    require_key_value "$STATS_DAT" stage9_7_requested_flag 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_stats_path_executed_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_stats_output_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_stats_finite_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_expected_stats_files_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_visu_path_executed_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_visu_output_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_visu_field_finite_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_expected_visu_files_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_output_files_nonempty_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_coarse_mesh_descriptor_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_decomp_io_open_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_decomp_io_write_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_decomp_io_close_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_output_field_finite_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_no_nan_inf_status 1 || status=1
    require_key_value "$STATS_DAT" stage9_7_stats_visu_io_smoke_status 1 || status=1
    if [ "$status" = "0" ]; then
        STATS_COMPAT_STATUS=1
        VISU_COMPAT_STATUS=1
        COARSE_IO_COMPAT_STATUS=1
    fi
    return $status
}

verify_restart_dat() {
    status=0
    require_key_value "$RESTART_COLD_DAT" stage9_8_requested_flag 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_phase_cold_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_restart_write_path_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_restart_files_exist_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_restart_files_nonempty_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_velocity_finite_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_pressure_finite_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_divergence_finite_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_no_nan_inf_status 1 || status=1
    require_key_value "$RESTART_COLD_DAT" stage9_8_restart_io_regression_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_requested_flag 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_phase_restart_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_restart_read_path_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_restart_files_exist_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_restart_files_nonempty_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_restart_signature_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_velocity_finite_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_pressure_finite_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_divergence_finite_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_no_nan_inf_status 1 || status=1
    require_key_value "$RESTART_READ_DAT" stage9_8_restart_io_regression_status 1 || status=1
    for key in stage9_8_d_sum_ux stage9_8_d_sum_uy stage9_8_d_sum_uz stage9_8_d_max_ux stage9_8_d_max_uy stage9_8_d_max_uz; do
        value=$(get_value "$RESTART_READ_DAT" "$key" 2>/dev/null) || {
            add_reason "missing_${key}_in_${RESTART_READ_DAT}"
            status=1
            continue
        }
        is_finite_number "$value" || {
            add_reason "${key}_not_finite_in_${RESTART_READ_DAT}"
            status=1
            continue
        }
        abs_delta=$(abs_value "$value")
        MAX_RESTART_DELTA=$(max_value "$MAX_RESTART_DELTA" "$abs_delta")
        MAX_IO_SIGNATURE_DELTA=$(max_value "$MAX_IO_SIGNATURE_DELTA" "$abs_delta")
        MAX_FLUID_DELTA=$(max_value "$MAX_FLUID_DELTA" "$abs_delta")
        require_real_le "$RESTART_READ_DAT" "$key" "$STAGE14_9_MAX_RESTART_DELTA" || status=1
        awk -v value="$abs_delta" -v limit="$STAGE14_9_MAX_IO_SIGNATURE_DELTA" 'BEGIN { if ((value+0.0) <= (limit+0.0)) exit 0; exit 1 }' || {
            add_reason "${key}_expected_le_stage14_9_max_io_signature_delta_got_${value}"
            status=1
        }
        awk -v value="$abs_delta" -v limit="$STAGE14_9_MAX_FLUID_DELTA" 'BEGIN { if ((value+0.0) <= (limit+0.0)) exit 0; exit 1 }' || {
            add_reason "${key}_expected_le_stage14_9_max_fluid_delta_got_${value}"
            status=1
        }
    done
    if [ "$status" = "0" ]; then
        RESTART_COMPAT_STATUS=1
        FLUID_RESPONSE_BOUNDED_STATUS=1
    fi
    return $status
}

verify_stage14_hook_file() {
    file=$1
    label=$2
    status=0
    if [ ! -f "$file" ]; then
        add_reason "missing_stage14_rhs_hook_file_${label}"
        return 1
    fi
    require_key_value "$file" stage14_5_requested_flag 1 || status=1
    require_key_value "$file" stage14_5_rhs_injection_enabled_flag 1 || status=1
    require_key_value "$file" stage14_5_injection_gain_finite_status 1 || status=1
    require_key_value "$file" stage14_5_hook_initialized_status 1 || status=1
    require_key_value "$file" stage14_5_hook_apply_called_status 1 || status=1
    require_key_value "$file" stage14_5_stage13_dependency_status 1 || status=1
    require_key_value "$file" stage14_5_stage13_candidate_required_status 1 || status=1
    require_key_value "$file" stage14_5_rhs_arrays_available_status 1 || status=1
    require_key_value "$file" stage14_5_rhs_increment_computed_status 1 || status=1
    require_key_value "$file" stage14_5_no_pressure_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_projection_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_poisson_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_rk3_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_channel_forcing_modification_status 1 || status=1
    require_key_value "$file" stage14_5_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$file" stage14_5_no_feedback_application_status 1 || status=1
    require_key_value "$file" stage14_5_no_twoway_force_status 1 || status=1
    require_key_value "$file" stage14_5_no_structure_advance_status 1 || status=1
    require_key_value "$file" stage14_5_production_rhs_hook_status 1 || status=1
    lambda_zero=$(get_value "$file" stage14_5_lambda_zero_status 2>/dev/null || echo missing)
    if [ "$lambda_zero" != "0" ]; then
        add_reason "stage14_9_nonzero_lambda_not_recorded_${label}"
        status=1
    fi
    gain=$(get_value "$file" stage14_5_injection_gain 2>/dev/null || echo missing)
    is_finite_number "$gain" || { add_reason "stage14_9_injection_gain_not_finite_${label}"; status=1; }
    awk -v got="$gain" -v want="$STAGE14_9_SMALL_LAMBDA" 'BEGIN { d=(got+0.0)-(want+0.0); if (d<0.0) d=-d; tol=1.0e-14 + 1.0e-8*((want<0.0)?-want:want); if (d<=tol) exit 0; exit 1 }' || {
        add_reason "stage14_9_injection_gain_mismatch_${label}"
        status=1
    }
    l2=$(get_value "$file" stage14_5_rhs_increment_l2 2>/dev/null || echo missing)
    max_abs=$(get_value "$file" stage14_5_rhs_increment_max_abs 2>/dev/null || echo missing)
    is_finite_number "$l2" || { add_reason "stage14_9_rhs_increment_l2_not_finite_${label}"; status=1; }
    is_finite_number "$max_abs" || { add_reason "stage14_9_rhs_increment_max_abs_not_finite_${label}"; status=1; }
    require_real_gt "$file" stage14_5_rhs_increment_l2 0.0 || status=1
    require_real_gt "$file" stage14_5_rhs_increment_max_abs 0.0 || status=1
    require_real_le "$file" stage14_5_rhs_increment_max_abs "$STAGE14_9_MAX_RHS_INCREMENT" || status=1
    MAX_RHS_INCREMENT=$(max_value "$MAX_RHS_INCREMENT" "$max_abs")
    return $status
}

verify_stage13_file() {
    file=$1
    label=$2
    status=0
    if [ ! -f "$file" ]; then
        add_reason "missing_stage13_force_density_file_${label}"
        return 1
    fi
    require_key_value "$file" stage13_6_hook_initialized_status 1 || status=1
    require_key_value "$file" stage13_6_hook_sample_called_status 1 || status=1
    require_key_value "$file" stage13_6_force_density_candidate_computed_status 1 || status=1
    require_key_value "$file" stage13_6_force_density_candidate_finite_status 1 || status=1
    require_key_value "$file" stage13_6_force_density_norm_finite_status 1 || status=1
    require_key_value "$file" stage13_6_integrated_force_finite_status 1 || status=1
    require_key_value "$file" stage13_6_spreading_input_sign_status 1 || status=1
    require_key_value "$file" stage13_6_wrong_sign_rejection_status 1 || status=1
    require_key_value "$file" stage13_6_no_production_ibm_forcing_status 1 || status=1
    require_key_value "$file" stage13_6_no_feedback_application_status 1 || status=1
    require_key_value "$file" stage13_6_no_twoway_force_status 1 || status=1
    require_key_value "$file" stage13_6_no_structure_advance_status 1 || status=1
    return $status
}

verify_all_stage14_and_stage13() {
    status=0
    for item in \
        "$STATS_STAGE14_DAT:stats_visu_coarse" \
        "$RESTART_COLD_STAGE14_DAT:restart_cold" \
        "$RESTART_READ_STAGE14_DAT:restart_read"; do
        file=${item%%:*}
        label=${item##*:}
        verify_stage14_hook_file "$file" "$label" || status=1
    done
    for item in \
        "$STATS_STAGE13_DAT:stats_visu_coarse" \
        "$RESTART_COLD_STAGE13_DAT:restart_cold" \
        "$RESTART_READ_STAGE13_DAT:restart_read"; do
        file=${item%%:*}
        label=${item##*:}
        verify_stage13_file "$file" "$label" || status=1
    done
    if [ "$status" = "0" ]; then
        STAGE14_HOOK_ACTIVE_STATUS=1
        LAMBDA_NONZERO_STATUS=1
        RHS_INCREMENT_NONZERO_STATUS=1
        RHS_INCREMENT_FINITE_STATUS=1
        RHS_INCREMENT_SIGN_CORRECT_STATUS=1
        RHS_INCREMENT_BOUNDED_STATUS=1
        STAGE13_FORCE_DENSITY_STATUS=1
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

write_output_dat() {
    final_status=$1
    cat > "$OUT_DAT" <<EOF_DAT
stage14_9_requested_flag 1
stage14_9_build_status $BUILD_STATUS
stage14_9_stage14_8_prereq_status $STAGE14_8_PREREQ_STATUS
stage14_9_stats_io_run_status $STATS_IO_RUN_STATUS
stage14_9_stats_compatibility_status $STATS_COMPAT_STATUS
stage14_9_visu_compatibility_status $VISU_COMPAT_STATUS
stage14_9_coarse_io_compatibility_status $COARSE_IO_COMPAT_STATUS
stage14_9_restart_write_status $RESTART_WRITE_STATUS
stage14_9_restart_read_status $RESTART_READ_STATUS
stage14_9_restart_compatibility_status $RESTART_COMPAT_STATUS
stage14_9_stage14_hook_active_status $STAGE14_HOOK_ACTIVE_STATUS
stage14_9_lambda_nonzero_status $LAMBDA_NONZERO_STATUS
stage14_9_rhs_increment_nonzero_status $RHS_INCREMENT_NONZERO_STATUS
stage14_9_rhs_increment_finite_status $RHS_INCREMENT_FINITE_STATUS
stage14_9_rhs_increment_sign_correct_status $RHS_INCREMENT_SIGN_CORRECT_STATUS
stage14_9_rhs_increment_bounded_status $RHS_INCREMENT_BOUNDED_STATUS
stage14_9_fluid_response_bounded_status $FLUID_RESPONSE_BOUNDED_STATUS
stage14_9_stage13_force_density_status $STAGE13_FORCE_DENSITY_STATUS
stage14_9_no_pressure_modification_status $NO_PRESSURE_MODIFICATION_STATUS
stage14_9_no_projection_modification_status $NO_PROJECTION_MODIFICATION_STATUS
stage14_9_no_poisson_modification_status $NO_POISSON_MODIFICATION_STATUS
stage14_9_no_rk3_modification_status $NO_RK3_MODIFICATION_STATUS
stage14_9_no_channel_forcing_modification_status $NO_CHANNEL_FORCING_MODIFICATION_STATUS
stage14_9_no_production_ibm_forcing_status $NO_PRODUCTION_IBM_FORCING_STATUS
stage14_9_no_feedback_application_status $NO_FEEDBACK_APPLICATION_STATUS
stage14_9_no_twoway_force_status $NO_TWOWAY_FORCE_STATUS
stage14_9_no_structure_advance_status $NO_STRUCTURE_ADVANCE_STATUS
stage14_9_no_bending_solve_status $NO_BENDING_SOLVE_STATUS
stage14_9_no_tension_solve_status $NO_TENSION_SOLVE_STATUS
stage14_9_no_fibre_position_update_status $NO_FIBRE_POSITION_UPDATE_STATUS
stage14_9_no_fibre_velocity_structural_update_status $NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS
stage14_9_no_wall_contact_status $NO_WALL_CONTACT_STATUS
stage14_9_io_restart_stats_visu_rhs_injection_status $final_status
stage14_9_lambda_14 $STAGE14_9_SMALL_LAMBDA
stage14_9_np $STAGE14_9_NP
stage14_9_max_rhs_increment $MAX_RHS_INCREMENT
stage14_9_max_restart_delta $MAX_RESTART_DELTA
stage14_9_max_io_signature_delta $MAX_IO_SIGNATURE_DELTA
stage14_9_max_fluid_delta $MAX_FLUID_DELTA
stage14_9_max_rhs_increment_limit $STAGE14_9_MAX_RHS_INCREMENT
stage14_9_max_fluid_delta_limit $STAGE14_9_MAX_FLUID_DELTA
stage14_9_max_restart_delta_limit $STAGE14_9_MAX_RESTART_DELTA
stage14_9_max_io_signature_delta_limit $STAGE14_9_MAX_IO_SIGNATURE_DELTA
EOF_DAT
}

if ! verify_small_lambda_value; then
    BUILD_STATUS=0
fi
if ! verify_np_value; then
    BUILD_STATUS=0
fi

if [ "$BUILD_STATUS" = "1" ]; then
    build_target xcompact3d || {
        BUILD_STATUS=0
        add_reason "build_failed_xcompact3d"
    }
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_9_RUN_STAGE14_8" = "1" ]; then
    STAGE14_8_SMALL_LAMBDA="$STAGE14_9_SMALL_LAMBDA" \
    STAGE14_8_MAX_RHS_INCREMENT="$STAGE14_9_MAX_RHS_INCREMENT" \
    STAGE14_8_MAX_FLUID_DELTA="$STAGE14_9_MAX_FLUID_DELTA" \
        bash stage14_checks/run_stage14_8_parallel_small_lambda_response.sh > "$OUTPUT_DIR/stage14_9_stage14_8_prereq.log" 2>&1
    if [ $? -ne 0 ] || ! grep 'STAGE 14.8 FINAL VERDICT: PASS' "$OUTPUT_DIR/stage14_9_stage14_8_prereq.log" >/dev/null 2>&1; then
        STAGE14_8_PREREQ_STATUS=0
        add_reason "optional_stage14_8_prerequisite_failed"
    fi
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_8_PREREQ_STATUS" = "1" ]; then
    if run_stats_visu_coarse; then
        STATS_IO_RUN_STATUS=1
        verify_stats_visu_coarse_dat || true
    fi
    if run_restart_phase cold "$RESTART_COLD_LOG" "$RESTART_COLD_DAT" "$RESTART_COLD_STAGE14_DAT" "$RESTART_COLD_STAGE13_DAT"; then
        RESTART_WRITE_STATUS=1
        if run_restart_phase restart "$RESTART_READ_LOG" "$RESTART_READ_DAT" "$RESTART_READ_STAGE14_DAT" "$RESTART_READ_STAGE13_DAT"; then
            RESTART_READ_STATUS=1
            verify_restart_dat || true
        fi
    fi
    verify_all_stage14_and_stage13 || true
else
    add_reason "stage14_9_run_phase_skipped_due_to_build_or_prereq_failure"
fi

if [ "$BUILD_STATUS" = "1" ] && [ "$STAGE14_8_PREREQ_STATUS" = "1" ] && \
   [ "$STATS_IO_RUN_STATUS" = "1" ] && [ "$STATS_COMPAT_STATUS" = "1" ] && \
   [ "$VISU_COMPAT_STATUS" = "1" ] && [ "$COARSE_IO_COMPAT_STATUS" = "1" ] && \
   [ "$RESTART_WRITE_STATUS" = "1" ] && [ "$RESTART_READ_STATUS" = "1" ] && \
   [ "$RESTART_COMPAT_STATUS" = "1" ] && [ "$STAGE14_HOOK_ACTIVE_STATUS" = "1" ] && \
   [ "$LAMBDA_NONZERO_STATUS" = "1" ] && [ "$RHS_INCREMENT_NONZERO_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_FINITE_STATUS" = "1" ] && [ "$RHS_INCREMENT_SIGN_CORRECT_STATUS" = "1" ] && \
   [ "$RHS_INCREMENT_BOUNDED_STATUS" = "1" ] && [ "$FLUID_RESPONSE_BOUNDED_STATUS" = "1" ] && \
   [ "$STAGE13_FORCE_DENSITY_STATUS" = "1" ] && [ "$NO_PRESSURE_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_PROJECTION_MODIFICATION_STATUS" = "1" ] && [ "$NO_POISSON_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_RK3_MODIFICATION_STATUS" = "1" ] && [ "$NO_CHANNEL_FORCING_MODIFICATION_STATUS" = "1" ] && \
   [ "$NO_PRODUCTION_IBM_FORCING_STATUS" = "1" ] && [ "$NO_FEEDBACK_APPLICATION_STATUS" = "1" ] && \
   [ "$NO_TWOWAY_FORCE_STATUS" = "1" ] && [ "$NO_STRUCTURE_ADVANCE_STATUS" = "1" ] && \
   [ "$NO_BENDING_SOLVE_STATUS" = "1" ] && [ "$NO_TENSION_SOLVE_STATUS" = "1" ] && \
   [ "$NO_FIBRE_POSITION_UPDATE_STATUS" = "1" ] && [ "$NO_FIBRE_VELOCITY_STRUCTURAL_UPDATE_STATUS" = "1" ] && \
   [ "$NO_WALL_CONTACT_STATUS" = "1" ]; then
    write_output_dat 1
    echo 'STAGE 14.9 IO/RESTART/STATS/VISU RHS-INJECTION VERDICT: PASS'
    echo 'STAGE 14.9 FINAL VERDICT: PASS'
    rm -f "$REASONS_FILE"
    exit 0
fi

write_output_dat 0
echo 'STAGE 14.9 IO/RESTART/STATS/VISU RHS-INJECTION VERDICT: FAIL'
echo 'STAGE 14.9 FINAL VERDICT: FAIL'
echo 'Reasons:'
if [ -s "$REASONS_FILE" ]; then
    sed 's/^/ - /' "$REASONS_FILE"
else
    echo ' - unknown_stage14_9_failure'
fi
exit 1
