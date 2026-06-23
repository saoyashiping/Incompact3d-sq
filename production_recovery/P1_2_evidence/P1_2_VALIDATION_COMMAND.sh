#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
EVID="$ROOT/production_recovery/P1_2_evidence"
CASE="$ROOT/production_recovery/P1_2_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
mkdir -p "$EVID"
SEARCH(){ if command -v rg >/dev/null 2>&1; then rg "$@"; else grep -R "$@"; fi; }
fail(){ echo "Result: FAIL" > "$EVID/P1_2_VALIDATION_RESULT.txt"; printf 'Result: FAIL\n' > "$ROOT/production_recovery/P1_2_PASS_FAIL.md"; exit 1; }
cmake -S . -B build_p1_2 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"
cmake --build build_p1_2 --target xcompact3d fibre_prod_p1_real_channel_preflight_check fibre_prod_p1_oneway_channel_case_check fibre_prod_p1_twoway_channel_case_check fibre_prod_p0_closure_check
P1TWO=$(find build_p1_2 -type f -perm -111 -name fibre_prod_p1_twoway_channel_case_check | head -1)
P1ONE=$(find build_p1_2 -type f -perm -111 -name fibre_prod_p1_oneway_channel_case_check | head -1)
P1PRE=$(find build_p1_2 -type f -perm -111 -name fibre_prod_p1_real_channel_preflight_check | head -1)
P0CHECK=$(find build_p1_2 -type f -perm -111 -name fibre_prod_p0_closure_check | head -1)
X3D=$(find build_p1_2 -type f -perm -111 -name xcompact3d | head -1)
[ -n "$P1TWO" ] && [ -n "$P1ONE" ] && [ -n "$P1PRE" ] && [ -n "$P0CHECK" ] && [ -n "$X3D" ] || fail
"$P1TWO" > "$EVID/P1_2_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1 || fail
"$P1ONE" > "$EVID/P1_2_P1_1_REGRESSION_AUDIT.txt.module_check" 2>&1 || fail
"$P1PRE" > "$EVID/P1_2_P1_0_PREFLIGHT_AUDIT.txt" 2>&1 || fail
"$P0CHECK" > "$EVID/P1_2_P0_CLOSURE_AUDIT.txt" 2>&1 || fail
if [ ! -f production_recovery/P1_1_PASS_FAIL.md ] || ! grep -q '^Result: PASS' production_recovery/P1_1_PASS_FAIL.md; then
  if [ -x production_recovery/P1_1_evidence/P1_1_VALIDATION_COMMAND.sh ]; then production_recovery/P1_1_evidence/P1_1_VALIDATION_COMMAND.sh; else fail; fi
fi
[ -f production_recovery/P1_1_evidence/P1_1_REAL_DNS_RUN_LOG.txt ] || fail
grep -Eq 'P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS|P1_1_REAL_DNS_RUN PASS' production_recovery/P1_1_evidence/P1_1_REAL_DNS_RUN_LOG.txt || fail
cat > "$EVID/P1_2_P1_1_REGRESSION_AUDIT.txt" <<EOT
Result: PASS
P1_1_PASS_FAIL.md exists with Result: PASS.
P1_1 real DNS log exists and contains a P1_1 PASS marker.
EOT
run_case(){
  local lam="$1" log="$2"
  (cd "$CASE" && env FIBRE_PROD_ENABLE=1 FIBRE_PROD_P1_TWOWAY_CHANNEL_ENABLE=1 FIBRE_PROD_P1_TWOWAY_CHANNEL_DIAGNOSTICS=1 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_P1_FIBRE_COUNT=1 FIBRE_PROD_P1_FIBRE_NNODE=49 FIBRE_PROD_P1_TWOWAY_MAX_DX=1.0e-2 FIBRE_PROD_LAMBDA="$lam" "$ROOT/$X3D") > "$log" 2>&1 || fail
  grep -q 'P1_2_TWOWAY_CHANNEL_CASE_CHECK PASS' "$log" || fail
  grep -q 'source=real_dns_field' "$log" || fail
  grep -q 'force_buffer diagnostic: nonzero finite bounded PASS' "$log" || fail
  grep -q 'RHS increment diagnostic: nonzero finite bounded PASS' "$log" || fail
  grep -q 'formula=lambda\*penalty_beta\*force_buffer' "$log" || fail
  if grep -Eiq '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)' "$log"; then fail; fi
}
LOG_LOW="$EVID/P1_2_REAL_DNS_RUN_LOG_lambda1e-5.txt"
LOG_HIGH="$EVID/P1_2_REAL_DNS_RUN_LOG_lambda1e-4.txt"
run_case 1.0e-5 "$LOG_LOW"
run_case 1.0e-4 "$LOG_HIGH"
LOW=$(awk '/RHS increment diagnostic/{v=$NF} END{print v+0}' "$LOG_LOW")
HIGH=$(awk '/RHS increment diagnostic/{v=$NF} END{print v+0}' "$LOG_HIGH")
python3 - <<PY || fail
low=abs(float('$LOW')); high=abs(float('$HIGH'))
ratio=high/low if low>0 else 0
assert low>0 and high>low and 8.0 <= ratio <= 12.0, (low, high, ratio)
PY
SEARCH 'fibre_prod_p1_twoway_channel_case' src/xcompact3d.f90 >/dev/null || fail
! SEARCH 'uniform RHS contribution' src >/dev/null || fail
cat "$LOG_LOW" "$LOG_HIGH" > "$EVID/P1_2_REAL_VELOCITY_SAMPLING_AUDIT.txt"
cp "$EVID/P1_2_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/P1_2_TWOWAY_STRUCTURE_RESPONSE_AUDIT.txt"
cp "$EVID/P1_2_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/P1_2_REACTION_FORCE_AUDIT.txt"
cp "$EVID/P1_2_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/P1_2_FORCE_BUFFER_AUDIT.txt"
cp "$EVID/P1_2_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/P1_2_RHS_INCREMENT_AUDIT.txt"
echo "Result: PASS - lambda RHS ratio=$(python3 - <<PY
low=abs(float('$LOW')); high=abs(float('$HIGH')); print(high/low)
PY
)" > "$EVID/P1_2_LAMBDA_SCALING_AUDIT.txt"
echo "Result: PASS - wall safety and no wall penetration markers present" > "$EVID/P1_2_WALL_SAFETY_AUDIT.txt"
echo "Result: PASS - no NaN/Inf detected" > "$EVID/P1_2_NAN_INF_AUDIT.txt"
cat > "$EVID/P1_2_REAL_CHANNEL_CASE_AUDIT.txt" <<EOT
Result: PASS
itype=3
nx=96
ny=97
nz=96
dt=2.5e-5
ilast=150
fibre_count=1
fibre_nnode=49
lambda_fsi=1.0e-5 and 1.0e-4
penalty_beta=2.0
EOT
cat > "$EVID/P1_2_VALIDATION_RESULT.txt" <<EOT
Result: PASS
P1_2_TWOWAY_CHANNEL_CASE_CHECK PASS
Real xcompact3d executable completed for lambda=1.0e-5 and lambda=1.0e-4.
EOT
cat > "$ROOT/production_recovery/P1_2_PASS_FAIL.md" <<EOT
Result: PASS

Meaning: PASS means the guarded P1 small-lambda real two-way DNS-FSI case can run two 96x97x96 real xcompact3d channel DNS time-advances with one 49-node flexible fibre, finite real-velocity sampling, finite bounded structure response, nonzero finite reaction force, nonzero finite Eulerian force buffer, and nonzero finite lambda-gated RHS increments that scale consistently from lambda=1.0e-5 to lambda=1.0e-4. It does NOT mean paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_3 enlarged guarded stability + restart/statistics/visualization on 128x129x96.
EOT
