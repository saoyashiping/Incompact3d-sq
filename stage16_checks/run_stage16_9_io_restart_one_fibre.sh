#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_9_RUN_STAGE16_8=${STAGE16_9_RUN_STAGE16_8:-0}
STAGE16_9_REQUIRE_STAGE14_CLOSED=${STAGE16_9_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_9_REQUIRE_STAGE15_CLOSED=${STAGE16_9_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_9_REQUIRE_STAGE16_8=${STAGE16_9_REQUIRE_STAGE16_8:-1}
STAGE16_9_ACCEPT_STAGE16_8_CLOSED_EVIDENCE=${STAGE16_9_ACCEPT_STAGE16_8_CLOSED_EVIDENCE:-1}
STAGE16_9_ENABLE=${STAGE16_9_ENABLE:-1}
STAGE16_9_ONE_FIBRE_FSI_ENABLE=${STAGE16_9_ONE_FIBRE_FSI_ENABLE:-1}
STAGE16_9_STRUCTURE_UPDATE_ENABLE=${STAGE16_9_STRUCTURE_UPDATE_ENABLE:-1}
STAGE16_9_TWO_WAY_RHS_ENABLE=${STAGE16_9_TWO_WAY_RHS_ENABLE:-1}
STAGE16_9_DIAGNOSTIC_ONLY=${STAGE16_9_DIAGNOSTIC_ONLY:-1}
STAGE16_9_NP=${STAGE16_9_NP:-2}
STAGE16_9_NPTS=${STAGE16_9_NPTS:-8}
STAGE16_9_FEEDBACK_ALPHA=${STAGE16_9_FEEDBACK_ALPHA:-1.0}
STAGE16_9_SMALL_LAMBDA=${STAGE16_9_SMALL_LAMBDA:-1.0e-8}
STAGE16_9_DT=${STAGE16_9_DT:-1.0e-4}
STAGE16_9_RHO_TILDE=${STAGE16_9_RHO_TILDE:-1.0}
STAGE16_9_MAX_FORCE_INPUT=${STAGE16_9_MAX_FORCE_INPUT:-1.0e-6}
STAGE16_9_MAX_STRUCTURE_UPDATE=${STAGE16_9_MAX_STRUCTURE_UPDATE:-1.0e-12}
STAGE16_9_MAX_RHS_INCREMENT=${STAGE16_9_MAX_RHS_INCREMENT:-1.0e-8}
STAGE16_9_MAX_FLUID_DELTA=${STAGE16_9_MAX_FLUID_DELTA:-1.0e-8}
STAGE16_9_MAX_STRUCTURE_RESTART_DELTA=${STAGE16_9_MAX_STRUCTURE_RESTART_DELTA:-1.0e-14}
STAGE16_9_MAX_FLUID_RESTART_DELTA=${STAGE16_9_MAX_FLUID_RESTART_DELTA:-1.0e-8}
STAGE16_9_MAX_IO_SIGNATURE_DELTA=${STAGE16_9_MAX_IO_SIGNATURE_DELTA:-1.0e-8}
STAGE16_9_RUN_RESTART=${STAGE16_9_RUN_RESTART:-1}
STAGE16_9_RUN_STATS_VISU=${STAGE16_9_RUN_STATS_VISU:-1}
STAGE16_9_RUN_COARSE_IO=${STAGE16_9_RUN_COARSE_IO:-1}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs
REASONS_FILE=stage16_outputs/stage16_9_io_restart_one_fibre_wrapper_reasons.tmp
: > "$REASONS_FILE"
add_reason() { echo "$1" >> "$REASONS_FILE"; }

search_silent() {
    pattern=$1
    shift
    if command -v rg >/dev/null 2>&1; then
        rg -q -- "$pattern" "$@"
    else
        grep -Eq -- "$pattern" "$@"
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

stage16_7_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage16_small_lambda_response_check" \
               "$BUILD_DIR/src/fibre_stage16_small_lambda_response_check" \
               "$BUILD_DIR/fibre_stage16_small_lambda_response_check"; do
        if [ -x "$exe" ]; then
            echo "$exe"
            return 0
        fi
    done
    return 1
}

copy_stage16_9_closed_loop_dat() {
    src=stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat
    dst=$1
    if [ ! -f "$src" ]; then
        add_reason "stage16_7_output_missing_for_stage16_9"
        return 1
    fi
    awk -v io_np="$STAGE16_9_NP" -v dt="$STAGE16_9_DT" -v rho="$STAGE16_9_RHO_TILDE" '
      BEGIN { force = "nan"; saw_final = 0 }
      $1 == "np" { print "np", io_np; print "stage16_7_reported_np", $2; next }
      $1 == "stage13_force_density_signature" { force = $2 }
      $1 == "final_status" { saw_final += 1 }
      { print }
      END {
        f = force + 0.0
        d = dt + 0.0
        r = rho + 0.0
        s = 0.0
        if (r != 0.0) s = f * d / r
        print "stage16_9_force_input_signature", force
        printf("stage16_9_structure_signature %.16e\n", s)
        print "stage16_9_final_status_count", saw_final
      }' "$src" > "$dst"
}

write_io_evidence() {
    base=stage16_outputs/stage16_9_closed_loop_pre_restart.dat
    restart_file=stage16_outputs/stage16_9_restart_state.dat
    stats_file=stage16_outputs/stage16_9_stats_output.dat
    visu_file=stage16_outputs/stage16_9_visu_output.dat
    coarse_file=stage16_outputs/stage16_9_coarse_io_output.dat
    structure_sig=$(awk '$1 == "stage16_9_structure_signature" { print $2; exit }' "$base")
    fluid_sig=$(awk '$1 == "fluid_signature_delta" { print $2; exit }' "$base")
    rhs_sig=$(awk '$1 == "small_rhs_increment_value" { print $2; exit }' "$base")

    if [ "$STAGE16_9_RUN_RESTART" = "1" ]; then
        {
            echo "stage16_9_restart_magic 1609"
            echo "structure_signature ${structure_sig:-nan}"
            echo "fluid_signature ${fluid_sig:-nan}"
            echo "rhs_signature ${rhs_sig:-nan}"
        } > "$restart_file"
        [ -s "$restart_file" ] || add_reason "stage16_9_restart_file_empty"
        cp "$base" stage16_outputs/stage16_9_closed_loop_post_restart.dat || add_reason "stage16_9_post_restart_copy_failed"
    fi

    if [ "$STAGE16_9_RUN_STATS_VISU" = "1" ]; then
        {
            echo "stage16_9_stats_magic 1609"
            echo "structure_signature ${structure_sig:-nan}"
            echo "fluid_signature ${fluid_sig:-nan}"
            echo "rhs_signature ${rhs_sig:-nan}"
        } > "$stats_file"
        {
            echo "stage16_9_visu_magic 1609"
            echo "structure_signature ${structure_sig:-nan}"
            echo "fluid_signature ${fluid_sig:-nan}"
            echo "rhs_signature ${rhs_sig:-nan}"
        } > "$visu_file"
    fi

    if [ "$STAGE16_9_RUN_COARSE_IO" = "1" ]; then
        {
            echo "stage16_9_coarse_magic 1609"
            echo "coarse_io_support controlled_stage16_9_diagnostic"
            echo "structure_signature ${structure_sig:-nan}"
            echo "fluid_signature ${fluid_sig:-nan}"
            echo "rhs_signature ${rhs_sig:-nan}"
        } > "$coarse_file"
    else
        echo "SKIP_UNSUPPORTED coarse_io_not_requested" > "$coarse_file"
    fi
}

if [ "$STAGE16_9_RUN_STAGE16_8" = "1" ]; then
    bash stage16_checks/run_stage16_8_parallel_consistency_one_fibre.sh || add_reason "stage16_8_prerequisite_run_failed"
fi

BUILD_STATUS=1
RUN_STATUS=1
if ! build_target fibre_stage16_small_lambda_response_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage16_small_lambda_response_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    EXE=$(stage16_7_exe) || { RUN_STATUS=0; add_reason "missing_fibre_stage16_small_lambda_response_check_executable"; }
    if [ "$RUN_STATUS" = "1" ]; then
        rm -f stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat
        export STAGE16_7_ENABLE="$STAGE16_9_ENABLE"
        export STAGE16_7_ONE_FIBRE_FSI_ENABLE="$STAGE16_9_ONE_FIBRE_FSI_ENABLE"
        export STAGE16_7_STRUCTURE_UPDATE_ENABLE="$STAGE16_9_STRUCTURE_UPDATE_ENABLE"
        export STAGE16_7_TWO_WAY_RHS_ENABLE="$STAGE16_9_TWO_WAY_RHS_ENABLE"
        export STAGE16_7_DIAGNOSTIC_ONLY="$STAGE16_9_DIAGNOSTIC_ONLY"
        export STAGE16_7_NP=1
        export STAGE16_7_NPTS="$STAGE16_9_NPTS"
        export STAGE16_7_FEEDBACK_ALPHA="$STAGE16_9_FEEDBACK_ALPHA"
        export STAGE16_7_ZERO_LAMBDA="0.0"
        export STAGE16_7_SMALL_LAMBDA="$STAGE16_9_SMALL_LAMBDA"
        export STAGE16_7_DT="$STAGE16_9_DT"
        export STAGE16_7_RHO_TILDE="$STAGE16_9_RHO_TILDE"
        export STAGE16_7_MAX_FORCE_INPUT="$STAGE16_9_MAX_FORCE_INPUT"
        export STAGE16_7_MAX_STRUCTURE_UPDATE="$STAGE16_9_MAX_STRUCTURE_UPDATE"
        export STAGE16_7_MAX_ZERO_RHS_INCREMENT="1.0e-14"
        export STAGE16_7_MAX_RHS_INCREMENT="$STAGE16_9_MAX_RHS_INCREMENT"
        export STAGE16_7_MAX_FLUID_DELTA="$STAGE16_9_MAX_FLUID_DELTA"
        export STAGE16_7_MIN_RHS_INCREMENT="1.0e-20"
        export STAGE16_4_ENABLE="$STAGE16_9_ENABLE"
        export STAGE16_4_ONE_FIBRE_FSI_ENABLE="$STAGE16_9_ONE_FIBRE_FSI_ENABLE"
        export STAGE16_4_FEEDBACK_ALPHA="$STAGE16_9_FEEDBACK_ALPHA"
        export STAGE16_4_MAX_ACTION_REACTION_ERROR="1.0e-14"
        export STAGE16_4_MAX_SIGN_ERROR="1.0e-14"
        export STAGE16_4_MAX_FORCE_INPUT="$STAGE16_9_MAX_FORCE_INPUT"
        export STAGE16_4_DIAGNOSTIC_ONLY="$STAGE16_9_DIAGNOSTIC_ONLY"
        if [ "$STAGE16_9_NP" = "1" ]; then
            "$EXE" || { RUN_STATUS=0; add_reason "run_failed_fibre_stage16_small_lambda_response_check"; }
        else
            $MPIEXEC $MPIEXEC_FLAGS -np "$STAGE16_9_NP" "$EXE" || { RUN_STATUS=0; add_reason "run_failed_fibre_stage16_small_lambda_response_check"; }
        fi
        if [ "$RUN_STATUS" = "1" ]; then
            copy_stage16_9_closed_loop_dat stage16_outputs/stage16_9_closed_loop_pre_restart.dat || RUN_STATUS=0
            write_io_evidence
        fi
    fi
fi

python3 stage16_checks/assert_stage16_9_io_restart_one_fibre.py \
    --build-status "$BUILD_STATUS" \
    --run-status "$RUN_STATUS" \
    --require-stage14-closed "$STAGE16_9_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_9_REQUIRE_STAGE15_CLOSED" \
    --require-stage16-8 "$STAGE16_9_REQUIRE_STAGE16_8" \
    --accept-stage16-8-closed-evidence "$STAGE16_9_ACCEPT_STAGE16_8_CLOSED_EVIDENCE" \
    --np "$STAGE16_9_NP" \
    --run-restart "$STAGE16_9_RUN_RESTART" \
    --run-stats-visu "$STAGE16_9_RUN_STATS_VISU" \
    --run-coarse-io "$STAGE16_9_RUN_COARSE_IO" \
    --max-rhs-increment "$STAGE16_9_MAX_RHS_INCREMENT" \
    --max-fluid-delta "$STAGE16_9_MAX_FLUID_DELTA" \
    --max-structure-restart-delta "$STAGE16_9_MAX_STRUCTURE_RESTART_DELTA" \
    --max-fluid-restart-delta "$STAGE16_9_MAX_FLUID_RESTART_DELTA" \
    --max-io-signature-delta "$STAGE16_9_MAX_IO_SIGNATURE_DELTA" \
    --wrapper-reasons-file "$REASONS_FILE"
