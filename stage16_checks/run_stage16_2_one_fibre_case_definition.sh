#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_2_RUN_STAGE16_1=${STAGE16_2_RUN_STAGE16_1:-0}
STAGE16_2_REQUIRE_STAGE14_CLOSED=${STAGE16_2_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_2_REQUIRE_STAGE15_CLOSED=${STAGE16_2_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_2_REQUIRE_STAGE16_1=${STAGE16_2_REQUIRE_STAGE16_1:-1}
STAGE16_2_ENABLE=${STAGE16_2_ENABLE:-1}
STAGE16_2_ONE_FIBRE_FSI_ENABLE=${STAGE16_2_ONE_FIBRE_FSI_ENABLE:-1}
STAGE16_2_NPTS=${STAGE16_2_NPTS:-8}
STAGE16_2_MIN_WALL_CLEARANCE=${STAGE16_2_MIN_WALL_CLEARANCE:-1.0e-2}
STAGE16_2_MIN_POINT_SPACING=${STAGE16_2_MIN_POINT_SPACING:-1.0e-6}
STAGE16_2_MAX_INITIAL_VELOCITY=${STAGE16_2_MAX_INITIAL_VELOCITY:-1.0e-8}
STAGE16_2_MAX_INITIAL_ACCELERATION=${STAGE16_2_MAX_INITIAL_ACCELERATION:-1.0e-8}
STAGE16_2_MAX_STRUCTURE_UPDATE=${STAGE16_2_MAX_STRUCTURE_UPDATE:-1.0e-12}
STAGE16_2_MAX_RHS_INCREMENT=${STAGE16_2_MAX_RHS_INCREMENT:-1.0e-8}
STAGE16_2_DIAGNOSTIC_ONLY=${STAGE16_2_DIAGNOSTIC_ONLY:-1}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs
REASONS_FILE=stage16_outputs/stage16_2_one_fibre_case_definition_wrapper_reasons.tmp
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

stage16_2_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage16_one_fibre_case_check" \
               "$BUILD_DIR/src/fibre_stage16_one_fibre_case_check" \
               "$BUILD_DIR/fibre_stage16_one_fibre_case_check"; do
        if [ -x "$exe" ]; then
            echo "$exe"
            return 0
        fi
    done
    return 1
}

if [ "$STAGE16_2_RUN_STAGE16_1" = "1" ]; then
    bash stage16_checks/run_stage16_1_config.sh || add_reason "stage16_1_prerequisite_run_failed"
fi

BUILD_STATUS=1
RUN_STATUS=1
if ! build_target fibre_stage16_one_fibre_case_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage16_one_fibre_case_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    EXE=$(stage16_2_exe) || { RUN_STATUS=0; add_reason "missing_fibre_stage16_one_fibre_case_check_executable"; }
    if [ "$RUN_STATUS" = "1" ]; then
        export STAGE16_2_ENABLE STAGE16_2_ONE_FIBRE_FSI_ENABLE STAGE16_2_NPTS
        export STAGE16_2_MIN_WALL_CLEARANCE STAGE16_2_MIN_POINT_SPACING
        export STAGE16_2_MAX_INITIAL_VELOCITY STAGE16_2_MAX_INITIAL_ACCELERATION
        export STAGE16_2_MAX_STRUCTURE_UPDATE STAGE16_2_MAX_RHS_INCREMENT STAGE16_2_DIAGNOSTIC_ONLY
        "$EXE" || { RUN_STATUS=0; add_reason "run_failed_fibre_stage16_one_fibre_case_check"; }
    fi
fi

python3 stage16_checks/assert_stage16_2_one_fibre_case_definition.py \
    --build-status "$BUILD_STATUS" \
    --run-status "$RUN_STATUS" \
    --require-stage14-closed "$STAGE16_2_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_2_REQUIRE_STAGE15_CLOSED" \
    --require-stage16-1 "$STAGE16_2_REQUIRE_STAGE16_1" \
    --wrapper-reasons-file "$REASONS_FILE"
