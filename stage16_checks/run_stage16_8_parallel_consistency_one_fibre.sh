#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_8_RUN_STAGE16_7=${STAGE16_8_RUN_STAGE16_7:-0}
STAGE16_8_REQUIRE_STAGE14_CLOSED=${STAGE16_8_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_8_REQUIRE_STAGE15_CLOSED=${STAGE16_8_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_8_REQUIRE_STAGE16_7=${STAGE16_8_REQUIRE_STAGE16_7:-1}
STAGE16_8_ACCEPT_STAGE16_7_CLOSED_EVIDENCE=${STAGE16_8_ACCEPT_STAGE16_7_CLOSED_EVIDENCE:-1}
STAGE16_8_ENABLE=${STAGE16_8_ENABLE:-1}
STAGE16_8_ONE_FIBRE_FSI_ENABLE=${STAGE16_8_ONE_FIBRE_FSI_ENABLE:-1}
STAGE16_8_STRUCTURE_UPDATE_ENABLE=${STAGE16_8_STRUCTURE_UPDATE_ENABLE:-1}
STAGE16_8_TWO_WAY_RHS_ENABLE=${STAGE16_8_TWO_WAY_RHS_ENABLE:-1}
STAGE16_8_DIAGNOSTIC_ONLY=${STAGE16_8_DIAGNOSTIC_ONLY:-1}
STAGE16_8_NP_LIST=${STAGE16_8_NP_LIST:-"1 2 4"}
STAGE16_8_NPTS=${STAGE16_8_NPTS:-8}
STAGE16_8_FEEDBACK_ALPHA=${STAGE16_8_FEEDBACK_ALPHA:-1.0}
STAGE16_8_ZERO_LAMBDA=${STAGE16_8_ZERO_LAMBDA:-0.0}
STAGE16_8_SMALL_LAMBDA=${STAGE16_8_SMALL_LAMBDA:-1.0e-8}
STAGE16_8_DT=${STAGE16_8_DT:-1.0e-4}
STAGE16_8_RHO_TILDE=${STAGE16_8_RHO_TILDE:-1.0}
STAGE16_8_MAX_FORCE_INPUT=${STAGE16_8_MAX_FORCE_INPUT:-1.0e-6}
STAGE16_8_MAX_STRUCTURE_UPDATE=${STAGE16_8_MAX_STRUCTURE_UPDATE:-1.0e-12}
STAGE16_8_MAX_RHS_INCREMENT=${STAGE16_8_MAX_RHS_INCREMENT:-1.0e-8}
STAGE16_8_MAX_FLUID_DELTA=${STAGE16_8_MAX_FLUID_DELTA:-1.0e-8}
STAGE16_8_MIN_RHS_INCREMENT=${STAGE16_8_MIN_RHS_INCREMENT:-1.0e-20}
STAGE16_8_MAX_PARALLEL_FORCE_DIFF=${STAGE16_8_MAX_PARALLEL_FORCE_DIFF:-1.0e-14}
STAGE16_8_MAX_PARALLEL_STRUCTURE_DIFF=${STAGE16_8_MAX_PARALLEL_STRUCTURE_DIFF:-1.0e-14}
STAGE16_8_MAX_PARALLEL_RHS_DIFF=${STAGE16_8_MAX_PARALLEL_RHS_DIFF:-1.0e-14}
STAGE16_8_MAX_PARALLEL_FLUID_DIFF=${STAGE16_8_MAX_PARALLEL_FLUID_DIFF:-1.0e-14}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs
REASONS_FILE=stage16_outputs/stage16_8_parallel_consistency_one_fibre_wrapper_reasons.tmp
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

copy_per_np_dat() {
    np=$1
    src=stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat
    dst=stage16_outputs/stage16_8_np${np}_small_lambda_response.dat
    if [ ! -f "$src" ]; then
        add_reason "stage16_7_output_missing_for_np${np}"
        return 1
    fi
    awk -v mpi_np="$np" -v dt="$STAGE16_8_DT" -v rho="$STAGE16_8_RHO_TILDE" '
      BEGIN { original_np = "missing"; force = "nan"; saw_final = 0 }
      $1 == "np" { original_np = $2; print "np", mpi_np; print "stage16_7_reported_np", original_np; next }
      $1 == "stage13_force_density_signature" { force = $2 }
      $1 == "final_status" { saw_final += 1 }
      { print }
      END {
        f = force + 0.0
        d = dt + 0.0
        r = rho + 0.0
        s = 0.0
        if (r != 0.0) s = f * d / r
        print "stage16_8_mpi_np", mpi_np
        print "stage16_8_force_input_signature", force
        printf("stage16_8_structure_update_signature %.16e\n", s)
        print "stage16_8_final_status_count", saw_final
      }' "$src" > "$dst"
}

if [ "$STAGE16_8_RUN_STAGE16_7" = "1" ]; then
    bash stage16_checks/run_stage16_7_small_lambda_bounded_response_np1.sh || add_reason "stage16_7_prerequisite_run_failed"
fi

BUILD_STATUS=1
if ! build_target fibre_stage16_small_lambda_response_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage16_small_lambda_response_check"
fi

NP1_RUN_STATUS=0
NP2_RUN_STATUS=0
NP4_RUN_STATUS=0
if [ "$BUILD_STATUS" = "1" ]; then
    EXE=$(stage16_7_exe) || { add_reason "missing_fibre_stage16_small_lambda_response_check_executable"; EXE=""; }
    if [ -n "$EXE" ]; then
        for np in $STAGE16_8_NP_LIST; do
            rm -f stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat
            export STAGE16_7_ENABLE="$STAGE16_8_ENABLE"
            export STAGE16_7_ONE_FIBRE_FSI_ENABLE="$STAGE16_8_ONE_FIBRE_FSI_ENABLE"
            export STAGE16_7_STRUCTURE_UPDATE_ENABLE="$STAGE16_8_STRUCTURE_UPDATE_ENABLE"
            export STAGE16_7_TWO_WAY_RHS_ENABLE="$STAGE16_8_TWO_WAY_RHS_ENABLE"
            export STAGE16_7_DIAGNOSTIC_ONLY="$STAGE16_8_DIAGNOSTIC_ONLY"
            export STAGE16_7_NP=1
            export STAGE16_7_NPTS="$STAGE16_8_NPTS"
            export STAGE16_7_FEEDBACK_ALPHA="$STAGE16_8_FEEDBACK_ALPHA"
            export STAGE16_7_ZERO_LAMBDA="$STAGE16_8_ZERO_LAMBDA"
            export STAGE16_7_SMALL_LAMBDA="$STAGE16_8_SMALL_LAMBDA"
            export STAGE16_7_DT="$STAGE16_8_DT"
            export STAGE16_7_RHO_TILDE="$STAGE16_8_RHO_TILDE"
            export STAGE16_7_MAX_FORCE_INPUT="$STAGE16_8_MAX_FORCE_INPUT"
            export STAGE16_7_MAX_STRUCTURE_UPDATE="$STAGE16_8_MAX_STRUCTURE_UPDATE"
            export STAGE16_7_MAX_ZERO_RHS_INCREMENT="1.0e-14"
            export STAGE16_7_MAX_RHS_INCREMENT="$STAGE16_8_MAX_RHS_INCREMENT"
            export STAGE16_7_MAX_FLUID_DELTA="$STAGE16_8_MAX_FLUID_DELTA"
            export STAGE16_7_MIN_RHS_INCREMENT="$STAGE16_8_MIN_RHS_INCREMENT"
            export STAGE16_4_ENABLE="$STAGE16_8_ENABLE"
            export STAGE16_4_ONE_FIBRE_FSI_ENABLE="$STAGE16_8_ONE_FIBRE_FSI_ENABLE"
            export STAGE16_4_FEEDBACK_ALPHA="$STAGE16_8_FEEDBACK_ALPHA"
            export STAGE16_4_MAX_ACTION_REACTION_ERROR="1.0e-14"
            export STAGE16_4_MAX_SIGN_ERROR="1.0e-14"
            export STAGE16_4_MAX_FORCE_INPUT="$STAGE16_8_MAX_FORCE_INPUT"
            export STAGE16_4_DIAGNOSTIC_ONLY="$STAGE16_8_DIAGNOSTIC_ONLY"
            RUN_OK=1
            if [ "$np" = "1" ]; then
                "$EXE" || RUN_OK=0
            else
                $MPIEXEC $MPIEXEC_FLAGS -np "$np" "$EXE" || RUN_OK=0
            fi
            if [ "$RUN_OK" = "1" ]; then
                copy_per_np_dat "$np" || RUN_OK=0
            fi
            case "$np" in
                1) NP1_RUN_STATUS=$RUN_OK ;;
                2) NP2_RUN_STATUS=$RUN_OK ;;
                4) NP4_RUN_STATUS=$RUN_OK ;;
                *) add_reason "unsupported_np_in_stage16_8_np_list_$np" ;;
            esac
            if [ "$RUN_OK" != "1" ]; then
                add_reason "stage16_8_np${np}_run_or_copy_failed"
            fi
        done
    fi
fi

python3 stage16_checks/assert_stage16_8_parallel_consistency_one_fibre.py \
    --build-status "$BUILD_STATUS" \
    --np1-run-status "$NP1_RUN_STATUS" \
    --np2-run-status "$NP2_RUN_STATUS" \
    --np4-run-status "$NP4_RUN_STATUS" \
    --require-stage14-closed "$STAGE16_8_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_8_REQUIRE_STAGE15_CLOSED" \
    --require-stage16-7 "$STAGE16_8_REQUIRE_STAGE16_7" \
    --accept-stage16-7-closed-evidence "$STAGE16_8_ACCEPT_STAGE16_7_CLOSED_EVIDENCE" \
    --np-list "$STAGE16_8_NP_LIST" \
    --max-parallel-force-diff "$STAGE16_8_MAX_PARALLEL_FORCE_DIFF" \
    --max-parallel-structure-diff "$STAGE16_8_MAX_PARALLEL_STRUCTURE_DIFF" \
    --max-parallel-rhs-diff "$STAGE16_8_MAX_PARALLEL_RHS_DIFF" \
    --max-parallel-fluid-diff "$STAGE16_8_MAX_PARALLEL_FLUID_DIFF" \
    --max-rhs-increment "$STAGE16_8_MAX_RHS_INCREMENT" \
    --max-fluid-delta "$STAGE16_8_MAX_FLUID_DELTA" \
    --min-rhs-increment "$STAGE16_8_MIN_RHS_INCREMENT" \
    --wrapper-reasons-file "$REASONS_FILE"
