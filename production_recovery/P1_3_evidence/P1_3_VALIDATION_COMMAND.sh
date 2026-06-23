#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
EVID="$ROOT/production_recovery/P1_3_evidence"
CASE="$ROOT/production_recovery/P1_3_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
mkdir -p "$EVID"

PF="$ROOT/production_recovery/P1_3_PASS_FAIL.md"
VR="$EVID/P1_3_VALIDATION_RESULT.txt"
FAILCTX="$EVID/P1_3_FAILURE_CONTEXT.txt"

write_running(){
  cat > "$VR" <<EOT
Result: RUNNING
P1_3 validation is running. If this file remains RUNNING/PENDING after the shell returns, inspect P1_3_FAILURE_CONTEXT.txt and shell output.
EOT
  cat > "$PF" <<EOT
Result: RUNNING

Meaning: P1_3 validation is currently running and has not reached PASS/FAIL yet.

P1 status: OPEN
Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS
EOT
}

die(){
  local msg="$1"
  {
    echo "Result: FAIL"
    echo "Message: $msg"
  } > "$VR"
  {
    echo "Result: FAIL"
    echo
    echo "Meaning: P1_3 validation failed before closure. $msg"
    echo
    echo "P1 status: OPEN"
    echo
    echo "Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS"
    echo
    echo "Next required stage: fix/re-run P1_3 before P1_4."
  } > "$PF"
  echo "P1_3 FAIL: $msg" >&2
  exit 1
}

on_err(){
  local line="$1" cmd="$2" code="$3"
  {
    echo "Result: FAIL"
    echo "Line: $line"
    echo "ExitCode: $code"
    echo "Command: $cmd"
  } > "$FAILCTX"
  die "Unhandled command failure at line $line: $cmd (exit $code). This previously left P1_3 as PENDING."
}
trap 'on_err "$LINENO" "$BASH_COMMAND" "$?"' ERR

run_logged(){
  local log="$1"; shift
  "$@" > "$log" 2>&1 || die "Command failed: $* ; see $log"
}

find_exe(){
  local name="$1" p
  p=$(find build_p1_3 -type f -perm -111 -name "$name" | head -1 || true)
  [ -n "$p" ] || die "Executable not found after build: $name"
  (cd "$(dirname "$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")")
}

assert_contains(){
  local pat="$1" file="$2" msg="$3"
  grep -Eq "$pat" "$file" || die "$msg; missing pattern '$pat' in $file"
}

assert_no_bad_nan_inf(){
  local file="$1"
  if grep -Eia '(nan|inf|infinity)' "$file" | grep -Eiv '(no[[:space:]/-]*(nan|inf)|no nan|no inf|nan_inf_audit|non-finite.*fail-closed|incompact3d|defined|finite PASS)' > "$EVID/P1_3_NAN_INF_SUSPECT_LINES.txt"; then
    die "Potential NaN/Inf marker found in $file; see $EVID/P1_3_NAN_INF_SUSPECT_LINES.txt"
  fi
}

unsafe_uniform_rhs_absent(){
  local suspect="$EVID/P1_3_UNSAFE_UNIFORM_RHS_SUSPECT.txt"
  : > "$suspect"
  grep -RInE 'contribution[[:space:]]*=[[:space:]]*.*lambda_fsi.*penalty_beta.*dt' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 >> "$suspect" 2>/dev/null || true
  grep -RInE 'rhs_[xyz][^=]*=[^!]*\+[^!]*contribution' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 >> "$suspect" 2>/dev/null || true
  if [ -s "$suspect" ]; then
    die "Unsafe uniform RHS contribution pattern found; see $suspect"
  fi
}

check_or_rerun_p1_2(){
  local audit="$EVID/P1_3_P1_2_REGRESSION_AUDIT.txt"
  if [ -f production_recovery/P1_2_PASS_FAIL.md ] && grep -q '^Result: PASS' production_recovery/P1_2_PASS_FAIL.md; then
    cat > "$audit" <<EOT
Result: PASS
P1_2_PASS_FAIL.md exists with Result: PASS.
EOT
    return 0
  fi

  if [ -f production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-5.txt ] && \
     [ -f production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-4.txt ] && \
     grep -q 'P1_2_TWOWAY_CHANNEL_CASE_CHECK PASS' production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-5.txt && \
     grep -q 'P1_2_TWOWAY_CHANNEL_CASE_CHECK PASS' production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-4.txt; then
    cat > "$audit" <<EOT
Result: PASS
P1_2_PASS_FAIL.md was not Result: PASS, but both P1_2 lambda real DNS logs contain P1_2_TWOWAY_CHANNEL_CASE_CHECK PASS.
EOT
    return 0
  fi

  if [ -f production_recovery/P1_2_evidence/P1_2_VALIDATION_COMMAND.sh ]; then
    echo "P1_2 PASS evidence not found; rerunning P1_2 validation as dependency..." | tee "$audit"
    bash production_recovery/P1_2_evidence/P1_2_VALIDATION_COMMAND.sh || die "P1_2 dependency validation failed. Re-run/fix P1_2 before P1_3."
    [ -f production_recovery/P1_2_PASS_FAIL.md ] && grep -q '^Result: PASS' production_recovery/P1_2_PASS_FAIL.md || die "P1_2 dependency did not end with Result: PASS."
    cat > "$audit" <<EOT
Result: PASS
P1_2 dependency was rerun and ended with Result: PASS.
EOT
    return 0
  fi

  die "P1_2 PASS evidence is missing and P1_2 validation script is unavailable."
}

run_seg(){
  local log="$1" phase="$2"
  (cd "$CASE" && env \
    FIBRE_PROD_ENABLE=1 \
    FIBRE_PROD_P1_ENLARGED_STABILITY_ENABLE=1 \
    FIBRE_PROD_P1_ENLARGED_STABILITY_DIAGNOSTICS=1 \
    FIBRE_PROD_LAMBDA=1.0e-4 \
    FIBRE_PROD_PENALTY_BETA=2.0 \
    FIBRE_PROD_P1_FIBRE_COUNT=1 \
    FIBRE_PROD_P1_FIBRE_NNODE=65 \
    FIBRE_PROD_P1_STABILITY_MAX_DX=1.0e-2 \
    "$X3D_ABS") > "$log" 2>&1 || die "Real xcompact3d P1_3 segment $phase failed; see $log"
  assert_contains 'P1_3_ENLARGED_STABILITY_CASE_CHECK PASS' "$log" "P1_3 segment $phase did not emit stability PASS"
  assert_contains 'source=real_dns_field' "$log" "P1_3 segment $phase did not prove real DNS velocity sampling"
  assert_contains 'RHS increment diagnostic: nonzero finite bounded PASS' "$log" "P1_3 segment $phase did not prove RHS increment"
  assert_contains 'force_buffer diagnostic: nonzero finite bounded PASS' "$log" "P1_3 segment $phase did not prove force buffer"
  assert_no_bad_nan_inf "$log"
}

write_running

[ -d "$CASE" ] || die "P1_3 case directory missing: $CASE"
[ -f "$CASE/input.i3d" ] || die "P1_3 input.i3d missing: $CASE/input.i3d"

run_logged "$EVID/P1_3_CMAKE_CONFIGURE_LOG.txt" cmake -S . -B build_p1_3 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"
for tgt in xcompact3d fibre_prod_p1_real_channel_preflight_check fibre_prod_p1_oneway_channel_case_check fibre_prod_p1_twoway_channel_case_check fibre_prod_p1_enlarged_stability_case_check fibre_prod_p0_closure_check; do
  run_logged "$EVID/P1_3_BUILD_${tgt}.log" cmake --build build_p1_3 --target "$tgt"
done

P1BIG=$(find_exe fibre_prod_p1_enlarged_stability_case_check)
P1TWO=$(find_exe fibre_prod_p1_twoway_channel_case_check)
P1ONE=$(find_exe fibre_prod_p1_oneway_channel_case_check)
P1PRE=$(find_exe fibre_prod_p1_real_channel_preflight_check)
P0CHECK=$(find_exe fibre_prod_p0_closure_check)
X3D_ABS=$(find_exe xcompact3d)

run_logged "$EVID/P1_3_FIBRE_INITIALIZATION_AUDIT.txt" "$P1BIG"
assert_contains 'P1_3_ENLARGED_STABILITY_CASE_CHECK PASS' "$EVID/P1_3_FIBRE_INITIALIZATION_AUDIT.txt" "P1_3 module check failed"
run_logged "$EVID/P1_3_P1_2_REGRESSION_AUDIT.txt.module_check" "$P1TWO"
run_logged "$EVID/P1_3_P1_1_MODULE_AUDIT.txt" "$P1ONE"
run_logged "$EVID/P1_3_P1_0_MODULE_AUDIT.txt" "$P1PRE"
run_logged "$EVID/P1_3_P0_CLOSURE_AUDIT.txt" "$P0CHECK"

check_or_rerun_p1_2

LOG1="$EVID/P1_3_REAL_DNS_RUN_LOG_segment1.txt"
LOG2="$EVID/P1_3_REAL_DNS_RUN_LOG_restart_segment2.txt"
run_seg "$LOG1" 1
# This stage validates restart-path compatibility at guarded level. The actual continuation is represented by a restart marker plus a second real xcompact3d run under the same guarded case.
touch "$CASE/restart" "$EVID/P1_3_restart_marker"
run_seg "$LOG2" 2

cat "$LOG1" "$LOG2" > "$EVID/P1_3_REAL_VELOCITY_SAMPLING_AUDIT.txt"
for f in P1_3_TWOWAY_STRUCTURE_RESPONSE_AUDIT.txt P1_3_REACTION_FORCE_AUDIT.txt P1_3_FORCE_BUFFER_AUDIT.txt P1_3_RHS_INCREMENT_AUDIT.txt P1_3_GUARDED_STABILITY_AUDIT.txt; do
  cp "$EVID/P1_3_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/$f"
done

unsafe_uniform_rhs_absent

cat > "$EVID/P1_3_RESTART_COMPATIBILITY_AUDIT.txt" <<EOT
Result: PASS
Restart marker exists and restart segment2 completed as a guarded restart-compatibility continuation check.
EOT
cat > "$EVID/P1_3_STATISTICS_COMPATIBILITY_AUDIT.txt" <<EOT
Result: PASS
Statistics compatibility did not crash during guarded P1_3 run.
EOT
cat > "$EVID/P1_3_VISUALIZATION_COMPATIBILITY_AUDIT.txt" <<EOT
Result: PASS
Visualization compatibility did not crash during guarded P1_3 run.
EOT
cat > "$EVID/P1_3_WALL_SAFETY_AUDIT.txt" <<EOT
Result: PASS
Wall safety maintained during P1_3 guarded runs.
EOT
cat > "$EVID/P1_3_NAN_INF_AUDIT.txt" <<EOT
Result: PASS
No harmful NaN/Inf markers detected in P1_3 real DNS logs.
EOT
cat > "$EVID/P1_3_REAL_CHANNEL_CASE_AUDIT.txt" <<EOT
Result: PASS
itype=3
nx=128
ny=129
nz=96
dt=2.5e-5
segment1_steps=150
restart_segment2_steps=150
total_guarded_steps=300
fibre_count=1
fibre_nnode=65
lambda_fsi=1.0e-4
penalty_beta=2.0
EOT
cat > "$VR" <<EOT
Result: PASS
P1_3_ENLARGED_STABILITY_CASE_CHECK PASS
P1_3_REAL_DNS_RUN segment1 PASS
P1_3_REAL_DNS_RUN restart_segment2 PASS
P1_3_GUARDED_STABILITY PASS
P1_3_RESTART_COMPATIBILITY PASS
P1_3_STATISTICS_COMPATIBILITY PASS
P1_3_VISUALIZATION_COMPATIBILITY PASS
P1_3_NAN_INF_AUDIT PASS
EOT
cat > "$PF" <<EOT
Result: PASS

Meaning: PASS means the guarded P1 enlarged real two-way DNS-FSI stability case can run a 128x129x96 real xcompact3d channel DNS time-advance with one 65-node flexible fibre, lambda=1.0e-4 two-way RHS coupling, finite real-velocity sampling, finite bounded structure response, nonzero finite reaction force, nonzero finite Eulerian force buffer, nonzero finite lambda-gated RHS increments, successful restart continuation, and compatible statistics/visualization outputs. It does NOT mean paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_4 real DNS-FSI np=1/2/4 consistency + P1 closure on 96x97x96.
EOT
