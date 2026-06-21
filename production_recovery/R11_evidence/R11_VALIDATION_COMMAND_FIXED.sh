#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR" || exit 2

EVIDENCE_DIR="$ROOT_DIR/production_recovery/R11_evidence"
BUILD_DIR="$ROOT_DIR/build_r11_mpi_consistency"
DECOMP_ROOT="/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4"
HOOK_LOG="$ROOT_DIR/production_recovery/R11_RUN_LOG_hook_check.txt"
BUILD_LOG="$ROOT_DIR/production_recovery/R11_BUILD_LOG.txt"
PASS_FAIL="$ROOT_DIR/production_recovery/R11_PASS_FAIL.md"
X3D_EXE="$BUILD_DIR/bin/xcompact3d"
HOOK_EXE="$BUILD_DIR/bin/fibre_prod_main_hook_check"

mkdir -p "$EVIDENCE_DIR"
rm -rf "$BUILD_DIR"

# Remove stale R10/R11 audit outputs before the new run.  R11 must never pass or fail
# from old BLOCKED/PASS files left by a previous attempt.
rm -f "$EVIDENCE_DIR"/R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt \
      "$EVIDENCE_DIR"/R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt \
      "$EVIDENCE_DIR"/R10_MAIN_HOOK_DIAGNOSTICS.txt \
      "$EVIDENCE_DIR"/R10_MAIN_HOOK_DIAGNOSTICS_lambda0.txt \
      "$EVIDENCE_DIR"/R10_MAIN_HOOK_DIAGNOSTICS_smalllambda.txt \
      "$EVIDENCE_DIR"/R11_LAMBDA0_NP*_AUDIT.txt \
      "$EVIDENCE_DIR"/R11_SMALLLAMBDA_NP*_AUDIT.txt \
      "$EVIDENCE_DIR"/R11_MAIN_HOOK_DIAGNOSTICS_lambda0_np*.txt \
      "$EVIDENCE_DIR"/R11_MAIN_HOOK_DIAGNOSTICS_smalllambda_np*.txt \
      "$EVIDENCE_DIR"/R11_MPI_CONSISTENCY_AUDIT.md

write_blocked_run_logs() {
  local reason="$1"
  for np in 1 2 4; do
    cat > "$ROOT_DIR/production_recovery/R11_RUN_LOG_lambda0_np${np}.txt" <<RUNEOF
# R11 lambda=0 np=${np} run log
BLOCKED: ${reason}
RUNEOF
    cat > "$ROOT_DIR/production_recovery/R11_RUN_LOG_smalllambda_np${np}.txt" <<RUNEOF
# R11 small-lambda np=${np} run log
BLOCKED: ${reason}
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

## R12 status

R12 was not entered.
AUDITEOF
  cat > "$PASS_FAIL" <<PFE
# Production Recovery R11 PASS/FAIL

Result: BLOCKED

## Reason

${reason}
PFE
}

write_fail_file() {
  local reason="$1"
  cat > "$PASS_FAIL" <<PFE
# Production Recovery R11 PASS/FAIL

Result: FAIL

## Reason

${reason}

## Boundary

R11 FAIL means at least one required np=1/2/4 controlled-hook consistency criterion failed. R12 was not entered.
PFE
}

write_pass_file() {
  cat > "$PASS_FAIL" <<'PFE'
# Production Recovery R11 PASS/FAIL

Result: PASS

## Meaning

R11 PASS means the R10 controlled main-loop hook passed np=1/2/4 MPI smoke consistency checks for both lambda=0 and small lambda.

## Boundary

R11 PASS does not mean R12 paper-level validation matrix PASS.
R11 PASS does not mean production DNS-FSI final closure.
PFE
}

copy_or_fail() {
  local src="$1"
  local dst="$2"
  local label="$3"
  if [ -f "$src" ]; then
    cp "$src" "$dst"
  else
    cat > "$dst" <<EOF_MISS
Result: FAIL
Reason: missing generated file ${label}: ${src}
EOF_MISS
  fi
}

run_case() {
  local mode="$1"
  local np="$2"
  local lambda="$3"
  local input_source="$4"
  local run_dir="$EVIDENCE_DIR/r11_run_${mode}_np${np}"
  local log="$ROOT_DIR/production_recovery/R11_RUN_LOG_${mode}_np${np}.txt"

  rm -rf "$run_dir"
  mkdir -p "$run_dir"
  cp "$input_source" "$run_dir/input.i3d"

  (
    cd "$run_dir" || exit 99
    rm -f R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt \
          R10_MAIN_HOOK_DIAGNOSTICS.txt R10_MAIN_HOOK_DIAGNOSTICS_lambda0.txt \
          R10_MAIN_HOOK_DIAGNOSTICS_smalllambda.txt
    FIBRE_PROD_ENABLE=1 \
    FIBRE_PROD_LAMBDA="$lambda" \
    FIBRE_PROD_DIAGNOSTICS=1 \
    FIBRE_PROD_DIAGNOSTICS_DIR="$run_dir" \
    mpirun --oversubscribe -np "$np" "$X3D_EXE"
  ) 2>&1 | tee "$log"
  local run_status=${PIPESTATUS[0]}
  echo "Return status: $run_status" >> "$log"

  if [ "$mode" = "lambda0" ]; then
    copy_or_fail "$run_dir/R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt" \
      "$EVIDENCE_DIR/R11_LAMBDA0_NP${np}_AUDIT.txt" "lambda0 np=${np} audit"
    copy_or_fail "$run_dir/R10_MAIN_HOOK_DIAGNOSTICS_lambda0.txt" \
      "$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_lambda0_np${np}.txt" "lambda0 np=${np} diagnostics"
  else
    copy_or_fail "$run_dir/R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt" \
      "$EVIDENCE_DIR/R11_SMALLLAMBDA_NP${np}_AUDIT.txt" "smalllambda np=${np} audit"
    copy_or_fail "$run_dir/R10_MAIN_HOOK_DIAGNOSTICS_smalllambda.txt" \
      "$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_smalllambda_np${np}.txt" "smalllambda np=${np} diagnostics"
  fi

  return "$run_status"
}

pass_grep() { grep -q "$1" "$2" 2>/dev/null; }
value_positive() {
  local key="$1"; local file="$2"
  awk -F= -v key="$key" '$1==key {gsub(/^[ \t]+|[ \t]+$/, "", $2); if (($2+0) > 0) ok=1} END {exit ok?0:1}' "$file" 2>/dev/null
}

diag_lambda0_ok() {
  local file="$1"
  pass_grep '^enabled=T' "$file" && \
  pass_grep '^modified_cells=0$' "$file" && \
  pass_grep '^last_status=0$' "$file" && \
  pass_grep '^failed_calls=0$' "$file" && \
  pass_grep '^before_finite=T' "$file" && \
  pass_grep '^after_finite=T' "$file" && \
  pass_grep '^no_contamination=T' "$file"
}

diag_small_ok() {
  local file="$1"
  pass_grep '^enabled=T' "$file" && \
  value_positive 'modified_cells' "$file" && \
  pass_grep '^last_status=0$' "$file" && \
  pass_grep '^failed_calls=0$' "$file" && \
  pass_grep '^before_finite=T' "$file" && \
  pass_grep '^after_finite=T' "$file" && \
  pass_grep '^small_lambda_response=T' "$file"
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

# Execute all six MPI smoke runs. Do not stop early; collect complete evidence.
declare -A run_statuses
for np in 1 2 4; do
  run_case lambda0 "$np" 0 "$input_source"; run_statuses["lambda0_np${np}"]=$?
done
for np in 1 2 4; do
  run_case smalllambda "$np" 1.0e-8 "$input_source"; run_statuses["smalllambda_np${np}"]=$?
done

all_pass=1
hook_check_pass=0
pass_grep 'R10_FIBRE_PROD_MAIN_HOOK_CHECK PASS' "$HOOK_LOG" && hook_check_pass=1
[ "$hook_check_pass" -eq 1 ] || all_pass=0

for np in 1 2 4; do
  log0="$ROOT_DIR/production_recovery/R11_RUN_LOG_lambda0_np${np}.txt"
  logs="$ROOT_DIR/production_recovery/R11_RUN_LOG_smalllambda_np${np}.txt"
  a0="$EVIDENCE_DIR/R11_LAMBDA0_NP${np}_AUDIT.txt"
  as="$EVIDENCE_DIR/R11_SMALLLAMBDA_NP${np}_AUDIT.txt"
  d0="$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_lambda0_np${np}.txt"
  ds="$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_smalllambda_np${np}.txt"

  [ "${run_statuses[lambda0_np${np}]}" -eq 0 ] || all_pass=0
  [ "${run_statuses[smalllambda_np${np}]}" -eq 0 ] || all_pass=0
  pass_grep 'Good job! Xcompact3d finished successfully!' "$log0" || all_pass=0
  pass_grep 'Good job! Xcompact3d finished successfully!' "$logs" || all_pass=0
  pass_grep '^Result: PASS' "$a0" || all_pass=0
  pass_grep '^Result: PASS' "$as" || all_pass=0
  diag_lambda0_ok "$d0" || all_pass=0
  diag_small_ok "$ds" || all_pass=0
done

{
  echo "# R11 MPI Consistency Audit"
  echo
  if [ "$all_pass" -eq 1 ]; then
    echo "Result: PASS"
  else
    echo "Result: FAIL"
  fi
  echo
  echo "## Hook check"
  echo
  echo "hook_check_pass=$hook_check_pass"
  echo
  echo "## np=1/2/4 lambda=0 result"
  echo
  for np in 1 2 4; do
    echo "### np=$np"
    echo "run_status=${run_statuses[lambda0_np${np}]}"
    grep -E '^(Result:|hook_calls=|enabled=|lambda_fsi=|modified_cells=|last_status=|failed_calls=|before_finite=|after_finite=|no_contamination=)' \
      "$EVIDENCE_DIR/R11_LAMBDA0_NP${np}_AUDIT.txt" "$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_lambda0_np${np}.txt" 2>/dev/null || true
    echo
  done
  echo "## np=1/2/4 small-lambda result"
  echo
  for np in 1 2 4; do
    echo "### np=$np"
    echo "run_status=${run_statuses[smalllambda_np${np}]}"
    grep -E '^(Result:|hook_calls=|enabled=|lambda_fsi=|modified_cells=|last_status=|failed_calls=|before_finite=|after_finite=|small_lambda_response=)' \
      "$EVIDENCE_DIR/R11_SMALLLAMBDA_NP${np}_AUDIT.txt" "$EVIDENCE_DIR/R11_MAIN_HOOK_DIAGNOSTICS_smalllambda_np${np}.txt" 2>/dev/null || true
    echo
  done
  echo "## NaN/Inf status"
  echo
  echo "before_finite/after_finite are checked in all diagnostics."
  echo
  echo "## R12 status"
  echo
  echo "R12 was not entered."
  echo
  echo "## Evidence boundary"
  echo
  echo "R11 evidence is limited to np=1/2/4 consistency for the controlled R10 hook."
} > "$EVIDENCE_DIR/R11_MPI_CONSISTENCY_AUDIT.md"

if [ "$all_pass" -eq 1 ]; then
  write_pass_file
else
  write_fail_file "One or more R11 MPI consistency criteria failed. See production_recovery/R11_evidence/R11_MPI_CONSISTENCY_AUDIT.md."
fi

cat "$PASS_FAIL"
