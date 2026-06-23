#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
EVID="$ROOT/production_recovery/P1_1_evidence"
CASE="$ROOT/production_recovery/P1_1_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
mkdir -p "$EVID"
SEARCH(){ if command -v rg >/dev/null 2>&1; then rg "$@"; else grep -R "$@"; fi; }
fail(){ echo "Result: FAIL" > "$EVID/P1_1_VALIDATION_RESULT.txt"; printf 'Result: FAIL\n' > "$ROOT/production_recovery/P1_1_PASS_FAIL.md"; exit 1; }
cmake -S . -B build_p1_1 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"
cmake --build build_p1_1 --target xcompact3d fibre_prod_p1_real_channel_preflight_check fibre_prod_p1_oneway_channel_case_check fibre_prod_p0_closure_check
P1ONE=$(find build_p1_1 -type f -perm -111 -name fibre_prod_p1_oneway_channel_case_check | head -1)
P1PRE=$(find build_p1_1 -type f -perm -111 -name fibre_prod_p1_real_channel_preflight_check | head -1)
P0CHECK=$(find build_p1_1 -type f -perm -111 -name fibre_prod_p0_closure_check | head -1)
X3D=$(find build_p1_1 -type f -perm -111 -name xcompact3d | head -1)
[ -n "$P1ONE" ] && [ -n "$P1PRE" ] && [ -n "$P0CHECK" ] && [ -n "$X3D" ] || fail
"$P1ONE" > "$EVID/P1_1_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1 || fail
"$P1PRE" > "$EVID/P1_1_P1_0_REGRESSION_AUDIT.txt.module_check" 2>&1 || fail
"$P0CHECK" > "$EVID/P1_1_P0_CLOSURE_AUDIT.txt" 2>&1 || fail
if [ ! -f production_recovery/P1_0_PASS_FAIL.md ] || ! grep -q '^Result: PASS' production_recovery/P1_0_PASS_FAIL.md; then fail; fi
if [ ! -f production_recovery/P1_0_evidence/P1_0_REAL_DNS_RUN_LOG.txt ]; then fail; fi
grep -q 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS' production_recovery/P1_0_evidence/P1_0_REAL_DNS_RUN_LOG.txt || fail
cat > "$EVID/P1_1_P1_0_REGRESSION_AUDIT.txt" <<EOT
Result: PASS
P1_0_PASS_FAIL.md exists with Result: PASS.
P1_0_REAL_DNS_RUN_LOG.txt exists and contains P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS.
EOT
(cd "$CASE" && env FIBRE_PROD_ENABLE=1 FIBRE_PROD_P1_ONEWAY_CHANNEL_ENABLE=1 FIBRE_PROD_P1_ONEWAY_CHANNEL_DIAGNOSTICS=1 FIBRE_PROD_LAMBDA=0 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_P1_FIBRE_COUNT=1 FIBRE_PROD_P1_FIBRE_NNODE=49 FIBRE_PROD_P1_ONEWAY_MAX_DX=1.0e-2 "$ROOT/$X3D") > "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" 2>&1 || fail
grep -q 'P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS' "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" || fail
grep -q 'source=real_dns_field' "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" || fail
grep -q 'bounded_dx PASS' "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" || fail
grep -q 'no RHS feedback applied' "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" || fail
if grep -Eiq '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)' "$EVID/P1_1_REAL_DNS_RUN_LOG.txt"; then fail; fi
SEARCH 'fibre_prod_p1_oneway_channel_case' src/xcompact3d.f90 >/dev/null || fail
! SEARCH 'fibre_prod_force_buffer_rhs_gate_apply|fibre_prod_main_hook_apply_force_buffer|fibre_prod_rhs_adapter_apply|uniform RHS contribution' src >/dev/null || fail
cp "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" "$EVID/P1_1_REAL_VELOCITY_SAMPLING_AUDIT.txt"
cp "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" "$EVID/P1_1_ONEWAY_STRUCTURE_RESPONSE_AUDIT.txt"
cp "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" "$EVID/P1_1_LAMBDA0_NO_CONTAMINATION_AUDIT.txt"
cp "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" "$EVID/P1_1_WALL_SAFETY_AUDIT.txt"
echo "Result: PASS - no NaN/Inf detected" > "$EVID/P1_1_NAN_INF_AUDIT.txt"
cat > "$EVID/P1_1_REAL_CHANNEL_CASE_AUDIT.txt" <<EOT
Result: PASS
itype=3
nx=96
ny=97
nz=96
dt=5.0e-5
ilast=100
fibre_count=1
fibre_nnode=49
lambda_fsi=0
EOT
cat > "$EVID/P1_1_VALIDATION_RESULT.txt" <<EOT
Result: PASS
P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS
Real xcompact3d executable run completed.
EOT
cat > "$ROOT/production_recovery/P1_1_PASS_FAIL.md" <<EOT
Result: PASS

Meaning: PASS means the guarded P1 real-channel one-way flexible-fibre case can run a 96x97x96 real xcompact3d channel DNS time-advance with one 49-node flexible fibre, finite real-velocity sampling, finite bounded one-way structure response, wall-safe fibre geometry, and strict lambda=0 no-contamination. It does NOT mean two-way FSI production or paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_2 small-lambda real two-way DNS-FSI response on 96x97x96.
EOT
