#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_0_REQUIRE_STAGE14_CLOSED=${STAGE16_0_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_0_REQUIRE_STAGE15_CLOSED=${STAGE16_0_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE=${STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE:-1}
STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING=${STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING:-0}
STAGE16_0_ENABLE=${STAGE16_0_ENABLE:-1}
STAGE16_0_DIAGNOSTIC_ONLY=${STAGE16_0_DIAGNOSTIC_ONLY:-1}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs

add_reason_file=stage16_outputs/stage16_0_preflight_closure_integrity_wrapper_reasons.tmp
: > "$add_reason_file"
add_wrapper_reason() { echo "$1" >> "$add_reason_file"; }

# Ripgrep is optional. Keep a grep fallback so this Stage 16.0 wrapper works on
# hosts where rg is not installed.
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

if [ "$STAGE16_0_ENABLE" != "1" ]; then
    add_wrapper_reason "stage16_0_enable_not_set_to_1"
fi
if [ "$STAGE16_0_DIAGNOSTIC_ONLY" != "1" ]; then
    add_wrapper_reason "stage16_0_diagnostic_only_not_set_to_1"
fi

if [ ! -s stage15_checks/STAGE15_CLOSED.md ] && [ "$STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING" = "1" ]; then
    if [ -x stage15_checks/run_stage15_11_total_smoke_closure.sh ] || [ -f stage15_checks/run_stage15_11_total_smoke_closure.sh ]; then
        bash stage15_checks/run_stage15_11_total_smoke_closure.sh || add_wrapper_reason "stage15_11_auto_run_failed"
    else
        add_wrapper_reason "stage15_11_auto_run_requested_but_script_missing"
    fi
fi

if ! ensure_build_dir; then
    add_wrapper_reason "stage16_0_build_dir_configure_failed"
fi

python3 stage16_checks/assert_stage16_0_preflight_closure_integrity.py \
    --require-stage14-closed "$STAGE16_0_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_0_REQUIRE_STAGE15_CLOSED" \
    --accept-closed-stage15-evidence "$STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE" \
    --enable "$STAGE16_0_ENABLE" \
    --diagnostic-only "$STAGE16_0_DIAGNOSTIC_ONLY" \
    --wrapper-reasons-file "$add_reason_file"
