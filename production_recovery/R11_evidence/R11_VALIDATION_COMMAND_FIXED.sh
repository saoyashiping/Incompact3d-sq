#!/usr/bin/env bash
set -u

ROOT_DIR="$(pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/R11_evidence"
BUILD_DIR="$ROOT_DIR/build_r11_mpi_consistency"
DECOMP_ROOT="/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4"
HOOK_LOG="$ROOT_DIR/production_recovery/R11_RUN_LOG_hook_check.txt"
BUILD_LOG="$ROOT_DIR/production_recovery/R11_BUILD_LOG.txt"
X3D_EXE="$BUILD_DIR/bin/xcompact3d"
HOOK_EXE="$BUILD_DIR/bin/fibre_prod_main_hook_check"

mkdir -p "$EVIDENCE_DIR"
rm -rf "$BUILD_DIR"

write_blocked_run_logs() {
  local reason="$1"
  for np in 1 2 4; do
    cat > "$ROOT_DIR/production_recovery/R11_RUN_LOG_lambda0_np${np}.txt" <<RUNEOF
# R11 lambda=0 np=${np} run log
BLOCKED: ${reason}
FIBRE_PROD_ENABLE=1
FIBRE_PROD_LAMBDA=0
FIBRE_PROD_DIAGNOSTICS=1
RUNEOF
    cat > "$ROOT_DIR/production_recovery/R11_RUN_LOG_smalllambda_np${np}.txt" <<RUNEOF
# R11 small-lambda np=${np} run log
BLOCKED: ${reason}
FIBRE_PROD_ENABLE=1
FIBRE_PROD_LAMBDA=1.0e-8
FIBRE_PROD_DIAGNOSTICS=1
RUNEOF
    cat > "$EVIDENCE_DIR/R11_LAMBDA0_NP${np}_AUDIT.txt" <<AUDITEOF
Result: BLOCKED
Reason: ${reason}
Required PASS evidence was not produced for lambda=0 np=${np}.
AUDITEOF
    cat > "$EVIDENCE_DIR/R11_SMALLLAMBDA_NP${np}_AUDIT.txt" <<AUDITEOF
Result: BLOCKED
Reason: ${reason}
Required PASS evidence was not produced for small-lambda np=${np}.
AUDITEOF
    cat > "$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_lambda0_np${np}.txt" <<DIAGEOF
Result: BLOCKED
enabled=UNKNOWN
lambda_fsi=0
modified_cells=UNKNOWN
last_status=UNKNOWN
failed_calls=UNKNOWN
before_finite=UNKNOWN
after_finite=UNKNOWN
no_contamination=UNKNOWN
DIAGEOF
    cat > "$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_smalllambda_np${np}.txt" <<DIAGEOF
Result: BLOCKED
enabled=UNKNOWN
lambda_fsi=1.0e-8
modified_cells=UNKNOWN
last_status=UNKNOWN
failed_calls=UNKNOWN
before_finite=UNKNOWN
after_finite=UNKNOWN
small_lambda_response=UNKNOWN
DIAGEOF
  done
}

write_mpi_audit_blocked() {
  local reason="$1"
  cat > "$EVIDENCE_DIR/R11_MPI_CONSISTENCY_AUDIT.md" <<AUDITEOF
# R11 MPI Consistency Audit

Result: BLOCKED

## Reason

${reason}

## np=1/2/4 lambda=0 result

BLOCKED for np=1, np=2, and np=4.

## np=1/2/4 small-lambda result

BLOCKED for np=1, np=2, and np=4.

## Hook diagnostics

No real hook diagnostics were produced because the build/run workflow did not complete.

## NaN/Inf status

No runtime NaN/Inf evidence exists because no xcompact3d MPI run completed.

## R12 status

R12 was not entered.

## Evidence boundary

R11 evidence is limited to np=1/2/4 consistency for the controlled R10 hook.  BLOCKED does not prove PASS or FAIL and does not certify production DNS-FSI final closure.
AUDITEOF
}

{
  echo "# R11 build log"
  echo "ROOT_DIR=$ROOT_DIR"
  echo "BUILD_DIR=$BUILD_DIR"
  echo "DECOMP_ROOT=$DECOMP_ROOT"
  echo
  echo "## Configure command"
  echo "cmake -S . -B \"$BUILD_DIR\" -DCMAKE_Fortran_COMPILER=mpif90 -DDECOMP2D_ROOT=\"$DECOMP_ROOT\" -DCMAKE_PREFIX_PATH=\"$DECOMP_ROOT\""
} > "$BUILD_LOG"

cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DDECOMP2D_ROOT="$DECOMP_ROOT" \
  -DCMAKE_PREFIX_PATH="$DECOMP_ROOT" \
  2>&1 | tee -a "$BUILD_LOG"
configure_status=${PIPESTATUS[0]}

{
  echo
  echo "## Build hook-check command"
  echo "cmake --build \"$BUILD_DIR\" --target fibre_prod_main_hook_check -j 2"
} >> "$BUILD_LOG"
cmake --build "$BUILD_DIR" --target fibre_prod_main_hook_check -j 2 \
  2>&1 | tee -a "$BUILD_LOG"
hook_build_status=${PIPESTATUS[0]}

{
  echo "# R11 hook-check run log"
  echo "## Run command"
  echo "\"$HOOK_EXE\""
} > "$HOOK_LOG"
if [ -x "$HOOK_EXE" ]; then
  "$HOOK_EXE" 2>&1 | tee -a "$HOOK_LOG"
  hook_run_status=${PIPESTATUS[0]}
else
  echo "BLOCKED: $HOOK_EXE was not produced." | tee -a "$HOOK_LOG"
  hook_run_status=127
fi

{
  echo
  echo "## Build xcompact3d command"
  echo "cmake --build \"$BUILD_DIR\" --target xcompact3d -j 2"
} >> "$BUILD_LOG"
cmake --build "$BUILD_DIR" --target xcompact3d -j 2 \
  2>&1 | tee -a "$BUILD_LOG"
xcompact_build_status=${PIPESTATUS[0]}

{
  echo
  echo "## R11 command statuses"
  echo "configure_status=$configure_status"
  echo "hook_build_status=$hook_build_status"
  echo "hook_run_status=$hook_run_status"
  echo "xcompact_build_status=$xcompact_build_status"
} >> "$BUILD_LOG"

if [ "$configure_status" -ne 0 ] || [ "$hook_build_status" -ne 0 ] || [ "$hook_run_status" -ne 0 ] || [ "$xcompact_build_status" -ne 0 ] || [ ! -x "$X3D_EXE" ]; then
  reason="R11 build/run prerequisites were unavailable or incomplete; configure_status=$configure_status, hook_build_status=$hook_build_status, hook_run_status=$hook_run_status, xcompact_build_status=$xcompact_build_status."
  write_blocked_run_logs "$reason"
  write_mpi_audit_blocked "$reason"
  exit 0
fi

input_source=""
if [ -f "$ROOT_DIR/input.i3d" ]; then
  input_source="$ROOT_DIR/input.i3d"
elif [ -f "$ROOT_DIR/tests/Channel/input_test_x.i3d" ]; then
  input_source="$ROOT_DIR/tests/Channel/input_test_x.i3d"
elif [ -f "$ROOT_DIR/examples/Channel/input.i3d" ]; then
  input_source="$ROOT_DIR/examples/Channel/input.i3d"
else
  reason="No input.i3d, tests/Channel/input_test_x.i3d, or examples/Channel/input.i3d was available."
  write_blocked_run_logs "$reason"
  write_mpi_audit_blocked "$reason"
  exit 0
fi

for mode in lambda0 smalllambda; do
  for np in 1 2 4; do
    run_dir="$EVIDENCE_DIR/r11_run_${mode}_np${np}"
    mkdir -p "$run_dir"
    cp "$input_source" "$run_dir/input.i3d"
    if [ "$mode" = "lambda0" ]; then
      log="$ROOT_DIR/production_recovery/R11_RUN_LOG_lambda0_np${np}.txt"
      lambda="0"
    else
      log="$ROOT_DIR/production_recovery/R11_RUN_LOG_smalllambda_np${np}.txt"
      lambda="1.0e-8"
    fi
    (
      cd "$run_dir" || exit 99
      rm -f R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt \
            R10_MAIN_HOOK_DIAGNOSTICS_lambda0.txt R10_MAIN_HOOK_DIAGNOSTICS_smalllambda.txt \
            R10_MAIN_HOOK_DIAGNOSTICS.txt
      FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA="$lambda" FIBRE_PROD_DIAGNOSTICS=1 \
        FIBRE_PROD_DIAGNOSTICS_DIR="$EVIDENCE_DIR" \
        mpirun --oversubscribe -np "$np" "$X3D_EXE" 2>&1
    ) | tee "$log"
  done
done
