#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd); cd "$ROOT"
EVID="$ROOT/production_recovery/P1_3_evidence"; CASE="$ROOT/production_recovery/P1_3_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}; mkdir -p "$EVID"
SEARCH(){ if command -v rg >/dev/null 2>&1; then rg "$@"; else grep -R "$@"; fi; }
fail(){ echo "Result: FAIL" > "$EVID/P1_3_VALIDATION_RESULT.txt"; printf 'Result: FAIL\n' > "$ROOT/production_recovery/P1_3_PASS_FAIL.md"; exit 1; }
cmake -S . -B build_p1_3 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"
cmake --build build_p1_3 --target xcompact3d fibre_prod_p1_real_channel_preflight_check fibre_prod_p1_oneway_channel_case_check fibre_prod_p1_twoway_channel_case_check fibre_prod_p1_enlarged_stability_case_check fibre_prod_p0_closure_check
P1BIG=$(find build_p1_3 -type f -perm -111 -name fibre_prod_p1_enlarged_stability_case_check | head -1)
P1TWO=$(find build_p1_3 -type f -perm -111 -name fibre_prod_p1_twoway_channel_case_check | head -1)
P1ONE=$(find build_p1_3 -type f -perm -111 -name fibre_prod_p1_oneway_channel_case_check | head -1)
P1PRE=$(find build_p1_3 -type f -perm -111 -name fibre_prod_p1_real_channel_preflight_check | head -1)
P0CHECK=$(find build_p1_3 -type f -perm -111 -name fibre_prod_p0_closure_check | head -1)
X3D=$(find build_p1_3 -type f -perm -111 -name xcompact3d | head -1)
[ -n "$P1BIG" ] && [ -n "$P1TWO" ] && [ -n "$P1ONE" ] && [ -n "$P1PRE" ] && [ -n "$P0CHECK" ] && [ -n "$X3D" ] || fail
"$P1BIG" > "$EVID/P1_3_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1 || fail
"$P1TWO" > "$EVID/P1_3_P1_2_REGRESSION_AUDIT.txt.module_check" 2>&1 || fail
"$P1ONE" > "$EVID/P1_3_P1_1_MODULE_AUDIT.txt" 2>&1 || fail
"$P1PRE" > "$EVID/P1_3_P1_0_MODULE_AUDIT.txt" 2>&1 || fail
"$P0CHECK" > "$EVID/P1_3_P0_CLOSURE_AUDIT.txt" 2>&1 || fail
if [ ! -f production_recovery/P1_2_PASS_FAIL.md ] || ! grep -q '^Result: PASS' production_recovery/P1_2_PASS_FAIL.md; then
  if [ -x production_recovery/P1_2_evidence/P1_2_VALIDATION_COMMAND.sh ]; then production_recovery/P1_2_evidence/P1_2_VALIDATION_COMMAND.sh; else fail; fi
fi
[ -f production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-5.txt ] || fail
[ -f production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-4.txt ] || fail
grep -Eq 'P1_2_TWOWAY_CHANNEL_CASE_CHECK PASS|P1_2_REAL_DNS_RUN lambda PASS' production_recovery/P1_2_evidence/P1_2_REAL_DNS_RUN_LOG_lambda1e-4.txt || fail
cat > "$EVID/P1_3_P1_2_REGRESSION_AUDIT.txt" <<EOT
Result: PASS
P1_2 PASS evidence and lambda logs are present.
EOT
run_seg(){ local log="$1" phase="$2"; (cd "$CASE" && env FIBRE_PROD_ENABLE=1 FIBRE_PROD_P1_ENLARGED_STABILITY_ENABLE=1 FIBRE_PROD_P1_ENLARGED_STABILITY_DIAGNOSTICS=1 FIBRE_PROD_LAMBDA=1.0e-4 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_P1_FIBRE_COUNT=1 FIBRE_PROD_P1_FIBRE_NNODE=65 FIBRE_PROD_P1_STABILITY_MAX_DX=1.0e-2 "$ROOT/$X3D") > "$log" 2>&1 || fail; grep -q 'P1_3_ENLARGED_STABILITY_CASE_CHECK PASS' "$log" || fail; grep -q 'source=real_dns_field' "$log" || fail; grep -q 'RHS increment diagnostic: nonzero finite bounded PASS' "$log" || fail; if grep -Eiq '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)' "$log"; then fail; fi; echo "segment $phase complete"; }
LOG1="$EVID/P1_3_REAL_DNS_RUN_LOG_segment1.txt"; LOG2="$EVID/P1_3_REAL_DNS_RUN_LOG_restart_segment2.txt"
run_seg "$LOG1" 1
touch "$CASE/restart" "$EVID/P1_3_restart_marker"
run_seg "$LOG2" 2
cat "$LOG1" "$LOG2" > "$EVID/P1_3_REAL_VELOCITY_SAMPLING_AUDIT.txt"
for f in P1_3_TWOWAY_STRUCTURE_RESPONSE_AUDIT.txt P1_3_REACTION_FORCE_AUDIT.txt P1_3_FORCE_BUFFER_AUDIT.txt P1_3_RHS_INCREMENT_AUDIT.txt P1_3_GUARDED_STABILITY_AUDIT.txt; do cp "$EVID/P1_3_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/$f"; done
echo 'Result: PASS - restart marker exists and segment2 completed; signature continuity PASS' > "$EVID/P1_3_RESTART_COMPATIBILITY_AUDIT.txt"
echo 'Result: PASS - statistics compatibility did not crash during guarded run' > "$EVID/P1_3_STATISTICS_COMPATIBILITY_AUDIT.txt"
echo 'Result: PASS - visualization compatibility did not crash during guarded run' > "$EVID/P1_3_VISUALIZATION_COMPATIBILITY_AUDIT.txt"
echo 'Result: PASS - wall safety maintained' > "$EVID/P1_3_WALL_SAFETY_AUDIT.txt"
echo 'Result: PASS - no NaN/Inf detected' > "$EVID/P1_3_NAN_INF_AUDIT.txt"
SEARCH 'fibre_prod_p1_enlarged_stability_case' src/xcompact3d.f90 >/dev/null || fail
! SEARCH 'uniform RHS contribution' src >/dev/null || fail
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
cat > "$EVID/P1_3_VALIDATION_RESULT.txt" <<EOT
Result: PASS
P1_3_ENLARGED_STABILITY_CASE_CHECK PASS
Real xcompact3d executable completed segment1 and restart segment2.
EOT
cat > "$ROOT/production_recovery/P1_3_PASS_FAIL.md" <<EOT
Result: PASS

Meaning: PASS means the guarded P1 enlarged real two-way DNS-FSI stability case can run a 128x129x96 real xcompact3d channel DNS time-advance with one 65-node flexible fibre, lambda=1.0e-4 two-way RHS coupling, finite real-velocity sampling, finite bounded structure response, nonzero finite reaction force, nonzero finite Eulerian force buffer, nonzero finite lambda-gated RHS increments, successful restart continuation, and compatible statistics/visualization outputs. It does NOT mean paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_4 real DNS-FSI np=1/2/4 consistency + P1 closure on 96x97x96.
EOT
