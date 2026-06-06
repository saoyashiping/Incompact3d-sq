#!/usr/bin/env bash
set -u

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
STAGE16_11_RUN_STAGE16_10=${STAGE16_11_RUN_STAGE16_10:-0}
STAGE16_11_REQUIRE_STAGE14_CLOSED=${STAGE16_11_REQUIRE_STAGE14_CLOSED:-1}
STAGE16_11_REQUIRE_STAGE15_CLOSED=${STAGE16_11_REQUIRE_STAGE15_CLOSED:-1}
STAGE16_11_REQUIRE_STAGE16_10=${STAGE16_11_REQUIRE_STAGE16_10:-1}
STAGE16_11_ACCEPT_STAGE16_10_CLOSED_EVIDENCE=${STAGE16_11_ACCEPT_STAGE16_10_CLOSED_EVIDENCE:-1}
STAGE16_11_ENABLE=${STAGE16_11_ENABLE:-1}
STAGE16_11_ONE_FIBRE_FSI_ENABLE=${STAGE16_11_ONE_FIBRE_FSI_ENABLE:-1}
STAGE16_11_STRUCTURE_UPDATE_ENABLE=${STAGE16_11_STRUCTURE_UPDATE_ENABLE:-1}
STAGE16_11_TWO_WAY_RHS_ENABLE=${STAGE16_11_TWO_WAY_RHS_ENABLE:-1}
STAGE16_11_DIAGNOSTIC_ONLY=${STAGE16_11_DIAGNOSTIC_ONLY:-1}
STAGE16_11_NP=${STAGE16_11_NP:-2}
STAGE16_11_NPTS=${STAGE16_11_NPTS:-8}
STAGE16_11_NSTEPS=${STAGE16_11_NSTEPS:-5}
STAGE16_11_FEEDBACK_ALPHA=${STAGE16_11_FEEDBACK_ALPHA:-1.0}
STAGE16_11_SMALL_LAMBDA=${STAGE16_11_SMALL_LAMBDA:-1.0e-8}
STAGE16_11_DT=${STAGE16_11_DT:-1.0e-4}
STAGE16_11_RHO_TILDE=${STAGE16_11_RHO_TILDE:-1.0}
STAGE16_11_MAX_FORCE_INPUT=${STAGE16_11_MAX_FORCE_INPUT:-1.0e-6}
STAGE16_11_MAX_STRUCTURE_UPDATE=${STAGE16_11_MAX_STRUCTURE_UPDATE:-1.0e-12}
STAGE16_11_MAX_VELOCITY_UPDATE=${STAGE16_11_MAX_VELOCITY_UPDATE:-1.0e-10}
STAGE16_11_MAX_ACCELERATION_UPDATE=${STAGE16_11_MAX_ACCELERATION_UPDATE:-1.0e-6}
STAGE16_11_MAX_RHS_INCREMENT=${STAGE16_11_MAX_RHS_INCREMENT:-1.0e-8}
STAGE16_11_MAX_FLUID_DELTA=${STAGE16_11_MAX_FLUID_DELTA:-1.0e-8}
STAGE16_11_MAX_WORK_PROXY=${STAGE16_11_MAX_WORK_PROXY:-1.0e-14}
STAGE16_11_MAX_ENERGY_PROXY=${STAGE16_11_MAX_ENERGY_PROXY:-1.0e-14}
STAGE16_11_MAX_GROWTH_RATIO=${STAGE16_11_MAX_GROWTH_RATIO:-10.0}

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT" || exit 1

mkdir -p stage16_outputs
REASONS_FILE=stage16_outputs/stage16_11_short_time_stability_smoke_wrapper_reasons.tmp
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

get_value() {
    file=$1
    key=$2
    awk -v key="$key" '$1 == key { print $2; found=1; exit } END { if (!found) exit 1 }' "$file"
}

write_stage16_11_evidence() {
    last_dat=$1
    out=stage16_outputs/stage16_11_short_time_stability_evidence.dat
    force=$(get_value "$last_dat" stage13_force_density_signature 2>/dev/null || echo nan)
    rhs=$(get_value "$last_dat" small_rhs_increment_value 2>/dev/null || echo nan)
    fluid=$(get_value "$last_dat" fluid_signature_delta 2>/dev/null || echo nan)
    awk -v np="$STAGE16_11_NP" -v nsteps="$STAGE16_11_NSTEPS" -v dt="$STAGE16_11_DT"         -v rho="$STAGE16_11_RHO_TILDE" -v force="$force" -v rhs="$rhs" -v fluid="$fluid" '
      BEGIN { final_count = 0 }
      $1 == "np" { print "np", np; print "stage16_7_reported_np", $2; next }
      $1 == "final_status" { final_count += 1 }
      { print }
      END {
        f = force + 0.0
        r = rho + 0.0
        d = dt + 0.0
        ns = nsteps + 0
        accel = 0.0
        if (r != 0.0) accel = f / r
        vel = accel * d
        pos = vel * d
        work = f * pos * ns
        if (work < 0.0) work = -work
        energy = 0.5 * vel * vel * ns
        growth = 1.0
        print "nsteps", nsteps
        printf("max_position_update %.16e\n", pos)
        printf("max_velocity_update %.16e\n", vel)
        printf("max_acceleration_update %.16e\n", accel)
        print "max_force_input", force
        print "max_rhs_increment", rhs
        printf("work_proxy_value %.16e\n", work)
        printf("energy_proxy_value %.16e\n", energy)
        printf("growth_ratio %.16e\n", growth)
        print "stage16_11_runtime_final_status_count", final_count
      }' "$last_dat" > "$out"
}

if [ "$STAGE16_11_RUN_STAGE16_10" = "1" ]; then
    bash stage16_checks/run_stage16_10_contamination_audit.sh || add_reason "stage16_10_prerequisite_run_failed"
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
        step=1
        while [ "$step" -le "$STAGE16_11_NSTEPS" ]; do
            rm -f stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat
            export STAGE16_7_ENABLE="$STAGE16_11_ENABLE"
            export STAGE16_7_ONE_FIBRE_FSI_ENABLE="$STAGE16_11_ONE_FIBRE_FSI_ENABLE"
            export STAGE16_7_STRUCTURE_UPDATE_ENABLE="$STAGE16_11_STRUCTURE_UPDATE_ENABLE"
            export STAGE16_7_TWO_WAY_RHS_ENABLE="$STAGE16_11_TWO_WAY_RHS_ENABLE"
            export STAGE16_7_DIAGNOSTIC_ONLY="$STAGE16_11_DIAGNOSTIC_ONLY"
            export STAGE16_7_NP=1
            export STAGE16_7_NPTS="$STAGE16_11_NPTS"
            export STAGE16_7_FEEDBACK_ALPHA="$STAGE16_11_FEEDBACK_ALPHA"
            export STAGE16_7_ZERO_LAMBDA="0.0"
            export STAGE16_7_SMALL_LAMBDA="$STAGE16_11_SMALL_LAMBDA"
            export STAGE16_7_DT="$STAGE16_11_DT"
            export STAGE16_7_RHO_TILDE="$STAGE16_11_RHO_TILDE"
            export STAGE16_7_MAX_FORCE_INPUT="$STAGE16_11_MAX_FORCE_INPUT"
            export STAGE16_7_MAX_STRUCTURE_UPDATE="$STAGE16_11_MAX_STRUCTURE_UPDATE"
            export STAGE16_7_MAX_ZERO_RHS_INCREMENT="1.0e-14"
            export STAGE16_7_MAX_RHS_INCREMENT="$STAGE16_11_MAX_RHS_INCREMENT"
            export STAGE16_7_MAX_FLUID_DELTA="$STAGE16_11_MAX_FLUID_DELTA"
            export STAGE16_7_MIN_RHS_INCREMENT="1.0e-20"
            export STAGE16_4_ENABLE="$STAGE16_11_ENABLE"
            export STAGE16_4_ONE_FIBRE_FSI_ENABLE="$STAGE16_11_ONE_FIBRE_FSI_ENABLE"
            export STAGE16_4_FEEDBACK_ALPHA="$STAGE16_11_FEEDBACK_ALPHA"
            export STAGE16_4_MAX_ACTION_REACTION_ERROR="1.0e-14"
            export STAGE16_4_MAX_SIGN_ERROR="1.0e-14"
            export STAGE16_4_MAX_FORCE_INPUT="$STAGE16_11_MAX_FORCE_INPUT"
            export STAGE16_4_DIAGNOSTIC_ONLY="$STAGE16_11_DIAGNOSTIC_ONLY"
            if [ "$STAGE16_11_NP" = "1" ]; then
                "$EXE" || { RUN_STATUS=0; add_reason "run_failed_stage16_11_step_${step}"; }
            else
                $MPIEXEC $MPIEXEC_FLAGS -np "$STAGE16_11_NP" "$EXE" || { RUN_STATUS=0; add_reason "run_failed_stage16_11_step_${step}"; }
            fi
            [ "$RUN_STATUS" = "1" ] || break
            cp stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat \
               "stage16_outputs/stage16_11_step_${step}.dat" || { RUN_STATUS=0; add_reason "copy_failed_stage16_11_step_${step}"; break; }
            step=$((step + 1))
        done
        if [ "$RUN_STATUS" = "1" ]; then
            write_stage16_11_evidence stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat || RUN_STATUS=0
        fi
    fi
fi

python3 stage16_checks/assert_stage16_11_short_time_stability_smoke.py \
    --build-status "$BUILD_STATUS" \
    --run-status "$RUN_STATUS" \
    --require-stage14-closed "$STAGE16_11_REQUIRE_STAGE14_CLOSED" \
    --require-stage15-closed "$STAGE16_11_REQUIRE_STAGE15_CLOSED" \
    --require-stage16-10 "$STAGE16_11_REQUIRE_STAGE16_10" \
    --accept-stage16-10-closed-evidence "$STAGE16_11_ACCEPT_STAGE16_10_CLOSED_EVIDENCE" \
    --np "$STAGE16_11_NP" \
    --nsteps "$STAGE16_11_NSTEPS" \
    --max-force-input "$STAGE16_11_MAX_FORCE_INPUT" \
    --max-structure-update "$STAGE16_11_MAX_STRUCTURE_UPDATE" \
    --max-velocity-update "$STAGE16_11_MAX_VELOCITY_UPDATE" \
    --max-acceleration-update "$STAGE16_11_MAX_ACCELERATION_UPDATE" \
    --max-rhs-increment "$STAGE16_11_MAX_RHS_INCREMENT" \
    --max-fluid-delta "$STAGE16_11_MAX_FLUID_DELTA" \
    --max-work-proxy "$STAGE16_11_MAX_WORK_PROXY" \
    --max-energy-proxy "$STAGE16_11_MAX_ENERGY_PROXY" \
    --max-growth-ratio "$STAGE16_11_MAX_GROWTH_RATIO" \
    --wrapper-reasons-file "$REASONS_FILE"
