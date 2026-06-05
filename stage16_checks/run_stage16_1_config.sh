#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_1_RUN_STAGE16_0=${STAGE16_1_RUN_STAGE16_0:-0}
STAGE16_1_REQUIRE_STAGE14_CLOSED=${STAGE16_1_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_1_REQUIRE_STAGE15_CLOSED=${STAGE16_1_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_1_REQUIRE_STAGE16_0=${STAGE16_1_REQUIRE_STAGE16_0:-1}
STAGE16_1_ENABLE=${STAGE16_1_ENABLE:-0}
STAGE16_1_ONE_FIBRE_FSI_ENABLE=${STAGE16_1_ONE_FIBRE_FSI_ENABLE:-0}
STAGE16_1_STRUCTURE_ADVANCE_ENABLE=${STAGE16_1_STRUCTURE_ADVANCE_ENABLE:-0}
STAGE16_1_TWO_WAY_RHS_ENABLE=${STAGE16_1_TWO_WAY_RHS_ENABLE:-0}
STAGE16_1_DIAGNOSTIC_ONLY=${STAGE16_1_DIAGNOSTIC_ONLY:-1}
STAGE16_1_FEEDBACK_ALPHA=${STAGE16_1_FEEDBACK_ALPHA:-1.0}
STAGE16_1_LAMBDA=${STAGE16_1_LAMBDA:-0.0}
STAGE16_1_MAX_STRUCTURE_UPDATE=${STAGE16_1_MAX_STRUCTURE_UPDATE:-1.0e-12}
STAGE16_1_MAX_RHS_INCREMENT=${STAGE16_1_MAX_RHS_INCREMENT:-1.0e-8}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs
REASONS_FILE=stage16_outputs/stage16_1_config_wrapper_reasons.tmp
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

stage16_1_exe() {
    for exe in "$BUILD_DIR/bin/fibre_stage16_config_check" \
               "$BUILD_DIR/src/fibre_stage16_config_check" \
               "$BUILD_DIR/fibre_stage16_config_check"; do
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

if [ "$STAGE16_1_RUN_STAGE16_0" = "1" ]; then
    bash stage16_checks/run_stage16_0_preflight_closure_integrity.sh || add_reason "stage16_0_prerequisite_run_failed"
fi

BUILD_STATUS=1
RUN_STATUS=1
if ! build_target fibre_stage16_config_check; then
    BUILD_STATUS=0
    add_reason "build_failed_fibre_stage16_config_check"
fi

if [ "$BUILD_STATUS" = "1" ]; then
    EXE=$(stage16_1_exe) || { RUN_STATUS=0; add_reason "missing_fibre_stage16_config_check_executable"; }
    if [ "$RUN_STATUS" = "1" ]; then
        export STAGE16_1_ENABLE STAGE16_1_ONE_FIBRE_FSI_ENABLE STAGE16_1_STRUCTURE_ADVANCE_ENABLE
        export STAGE16_1_TWO_WAY_RHS_ENABLE STAGE16_1_DIAGNOSTIC_ONLY STAGE16_1_FEEDBACK_ALPHA
        export STAGE16_1_LAMBDA STAGE16_1_MAX_STRUCTURE_UPDATE STAGE16_1_MAX_RHS_INCREMENT
        export STAGE16_1_REQUIRE_STAGE15_CLOSED
        "$EXE" || { RUN_STATUS=0; add_reason "run_failed_fibre_stage16_config_check"; }
    fi
fi

python3 stage16_checks/assert_stage16_1_config.py \
    --build-status "$BUILD_STATUS" \
    --run-status "$RUN_STATUS" \
    --require-stage14-closed "$STAGE16_1_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_1_REQUIRE_STAGE15_CLOSED" \
    --require-stage16-0 "$STAGE16_1_REQUIRE_STAGE16_0" \
    --wrapper-reasons-file "$REASONS_FILE"
