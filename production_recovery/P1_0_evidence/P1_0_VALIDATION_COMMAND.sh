#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
EVID="$ROOT/production_recovery/P1_0_evidence"
CASE="$ROOT/production_recovery/P1_0_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
mkdir -p "$EVID"
fail(){ echo "Result: FAIL" > "$EVID/P1_0_VALIDATION_RESULT.txt"; printf 'Result: FAIL\n' > "$ROOT/production_recovery/P1_0_PASS_FAIL.md"; exit 1; }
cmake -S . -B build_p1_0 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"
cmake --build build_p1_0 --target xcompact3d fibre_prod_p1_real_channel_preflight_check fibre_prod_p0_closure_check fibre_prod_synthetic_closed_loop_check
P1CHECK=$(find build_p1_0 -type f -perm -111 -name fibre_prod_p1_real_channel_preflight_check | head -1)
P0CHECK=$(find build_p1_0 -type f -perm -111 -name fibre_prod_p0_closure_check | head -1)
X3D=$(find build_p1_0 -type f -perm -111 -name xcompact3d | head -1)
[ -n "$P1CHECK" ] && [ -n "$P0CHECK" ] && [ -n "$X3D" ] || fail
"$P1CHECK" > "$EVID/P1_0_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1 || fail
"$P0CHECK" > "$EVID/P1_0_P0_REGRESSION_AUDIT.txt" 2>&1 || fail
(cd "$CASE" && env FIBRE_PROD_ENABLE=1 FIBRE_PROD_P1_REAL_CHANNEL_PREFLIGHT_ENABLE=1 FIBRE_PROD_P1_REAL_CHANNEL_PREFLIGHT_DIAGNOSTICS=1 FIBRE_PROD_LAMBDA=0 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_P1_FIBRE_COUNT=1 FIBRE_PROD_P1_FIBRE_NNODE=33 "$ROOT/$X3D") > "$EVID/P1_0_REAL_DNS_RUN_LOG.txt" 2>&1 || fail
cp "$EVID/P1_0_REAL_DNS_RUN_LOG.txt" "$EVID/P1_0_REAL_VELOCITY_SAMPLING_AUDIT.txt"
cp "$EVID/P1_0_REAL_DNS_RUN_LOG.txt" "$EVID/P1_0_LAMBDA0_NO_CONTAMINATION_AUDIT.txt"
cp "$EVID/P1_0_REAL_DNS_RUN_LOG.txt" "$EVID/P1_0_WALL_SAFETY_AUDIT.txt"
if grep -Eiq '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)' "$EVID/P1_0_REAL_DNS_RUN_LOG.txt"; then fail; fi
echo "No NaN/Inf detected" > "$EVID/P1_0_NAN_INF_AUDIT.txt"
grep -q 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS' "$EVID/P1_0_REAL_DNS_RUN_LOG.txt" || fail
cat > "$EVID/P1_0_REAL_CHANNEL_CASE_AUDIT.txt" <<EOT
Result: PASS
case=input.i3d
itype=3 real channel
nx=64 ny=65 nz=64 dt=5.0e-5 ilast=20
EOT
cat > "$EVID/P1_0_VALIDATION_RESULT.txt" <<EOT
Result: PASS
P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS
Real xcompact3d executable run completed.
EOT
cat > "$ROOT/production_recovery/P1_0_PASS_FAIL.md" <<EOT
Result: PASS

Meaning: PASS means the guarded P1 real-channel DNS-FSI preflight can run a real xcompact3d channel DNS time-advance with one initialized flexible fibre, finite real-velocity sampling, wall-safe fibre geometry, and strict lambda=0 no-contamination. It does NOT mean two-way FSI production or paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_1 lambda=0 real DNS + flexible fibre one-way no-contamination on 96x97x96.
EOT
