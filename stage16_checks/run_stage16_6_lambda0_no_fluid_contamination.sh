#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_6_RUN_STAGE16_5=${STAGE16_6_RUN_STAGE16_5:-0}
STAGE16_6_REQUIRE_STAGE14_CLOSED=${STAGE16_6_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_6_REQUIRE_STAGE15_CLOSED=${STAGE16_6_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_6_REQUIRE_STAGE16_5=${STAGE16_6_REQUIRE_STAGE16_5:-1}
STAGE16_6_ACCEPT_STAGE16_5_CLOSED_EVIDENCE=${STAGE16_6_ACCEPT_STAGE16_5_CLOSED_EVIDENCE:-1}
STAGE16_6_ENABLE=${STAGE16_6_ENABLE:-1}
STAGE16_6_ONE_FIBRE_FSI_ENABLE=${STAGE16_6_ONE_FIBRE_FSI_ENABLE:-1}
STAGE16_6_STRUCTURE_UPDATE_ENABLE=${STAGE16_6_STRUCTURE_UPDATE_ENABLE:-1}
STAGE16_6_TWO_WAY_RHS_ENABLE=${STAGE16_6_TWO_WAY_RHS_ENABLE:-1}
STAGE16_6_DIAGNOSTIC_ONLY=${STAGE16_6_DIAGNOSTIC_ONLY:-1}
STAGE16_6_NP=${STAGE16_6_NP:-1}
STAGE16_6_NPTS=${STAGE16_6_NPTS:-8}
STAGE16_6_FEEDBACK_ALPHA=${STAGE16_6_FEEDBACK_ALPHA:-1.0}
STAGE16_6_LAMBDA=${STAGE16_6_LAMBDA:-0.0}
STAGE16_6_DT=${STAGE16_6_DT:-1.0e-4}
STAGE16_6_RHO_TILDE=${STAGE16_6_RHO_TILDE:-1.0}
STAGE16_6_MAX_FORCE_INPUT=${STAGE16_6_MAX_FORCE_INPUT:-1.0e-6}
STAGE16_6_MAX_STRUCTURE_UPDATE=${STAGE16_6_MAX_STRUCTURE_UPDATE:-1.0e-12}
STAGE16_6_MAX_RHS_INCREMENT=${STAGE16_6_MAX_RHS_INCREMENT:-1.0e-14}
STAGE16_6_MAX_FLUID_DELTA=${STAGE16_6_MAX_FLUID_DELTA:-1.0e-14}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs
REASONS_FILE=stage16_outputs/stage16_6_lambda0_no_fluid_contamination_wrapper_reasons.tmp
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

stage16_6_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage16_lambda0_no_contamination_check" \
               "$BUILD_DIR/src/fibre_stage16_lambda0_no_contamination_check" \
               "$BUILD_DIR/fibre_stage16_lambda0_no_contamination_check"; do
        if [ -x "$exe" ]; then
            echo "$exe"
            return 0
        fi
    done
    return 1
}

if [ "$STAGE16_6_RUN_STAGE16_5" = "1" ]; then
    bash stage16_checks/run_stage16_5_closed_loop_dryrun_np1.sh || add_reason "stage16_5_prerequisite_run_failed"
fi

BUILD_STATUS=1
RUN_STATUS=1
if ! build_target fibre_stage16_lambda0_no_contamination_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage16_lambda0_no_contamination_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    EXE=$(stage16_6_exe) || { RUN_STATUS=0; add_reason "missing_fibre_stage16_lambda0_no_contamination_check_executable"; }
    if [ "$RUN_STATUS" = "1" ]; then
        export STAGE16_6_ENABLE STAGE16_6_ONE_FIBRE_FSI_ENABLE STAGE16_6_STRUCTURE_UPDATE_ENABLE
        export STAGE16_6_TWO_WAY_RHS_ENABLE STAGE16_6_DIAGNOSTIC_ONLY STAGE16_6_NP STAGE16_6_NPTS
        export STAGE16_6_FEEDBACK_ALPHA STAGE16_6_LAMBDA STAGE16_6_DT STAGE16_6_RHO_TILDE
        export STAGE16_6_MAX_FORCE_INPUT STAGE16_6_MAX_STRUCTURE_UPDATE STAGE16_6_MAX_RHS_INCREMENT
        export STAGE16_6_MAX_FLUID_DELTA
        export STAGE16_4_ENABLE="$STAGE16_6_ENABLE"
        export STAGE16_4_ONE_FIBRE_FSI_ENABLE="$STAGE16_6_ONE_FIBRE_FSI_ENABLE"
        export STAGE16_4_FEEDBACK_ALPHA="$STAGE16_6_FEEDBACK_ALPHA"
        export STAGE16_4_MAX_ACTION_REACTION_ERROR="1.0e-14"
        export STAGE16_4_MAX_SIGN_ERROR="1.0e-14"
        export STAGE16_4_MAX_FORCE_INPUT="$STAGE16_6_MAX_FORCE_INPUT"
        export STAGE16_4_DIAGNOSTIC_ONLY="$STAGE16_6_DIAGNOSTIC_ONLY"
        "$EXE" || { RUN_STATUS=0; add_reason "run_failed_fibre_stage16_lambda0_no_contamination_check"; }
    fi
fi

python3 stage16_checks/assert_stage16_6_lambda0_no_fluid_contamination.py \
    --build-status "$BUILD_STATUS" \
    --run-status "$RUN_STATUS" \
    --require-stage14-closed "$STAGE16_6_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_6_REQUIRE_STAGE15_CLOSED" \
    --require-stage16-5 "$STAGE16_6_REQUIRE_STAGE16_5" \
    --accept-stage16-5-closed-evidence "$STAGE16_6_ACCEPT_STAGE16_5_CLOSED_EVIDENCE" \
    --wrapper-reasons-file "$REASONS_FILE"
