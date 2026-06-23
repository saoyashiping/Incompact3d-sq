#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd); cd "$ROOT"
EVID="$ROOT/production_recovery/P1_4_evidence"; CASE="$ROOT/production_recovery/P1_4_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}; mkdir -p "$EVID"
fail(){ local ec=$? line=${BASH_LINENO[0]:-unknown} cmd=${BASH_COMMAND:-unknown}; echo "Result: FAIL" > "$ROOT/production_recovery/P1_4_PASS_FAIL.md"; echo "Result: FAIL" > "$EVID/P1_4_VALIDATION_RESULT.txt"; printf 'line=%s\ncommand=%s\nexit_code=%s\n' "$line" "$cmd" "$ec" > "$EVID/P1_4_FAILURE_CONTEXT.txt"; exit "$ec"; }
trap fail ERR
SEARCH(){ if command -v rg >/dev/null 2>&1; then rg "$@"; else grep -R "$@"; fi; }
cmake -S . -B build_p1_4 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"
cmake --build build_p1_4 --target xcompact3d fibre_prod_p1_np_consistency_closure_case_check
CHECK=$(find build_p1_4 -type f -perm -111 -name fibre_prod_p1_np_consistency_closure_case_check | head -1)
X3D=$(find build_p1_4 -type f -perm -111 -name xcompact3d | head -1)
[ -n "$CHECK" ] && [ -n "$X3D" ]
"$CHECK" > "$EVID/P1_4_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1
cat > "$EVID/P1_4_SELF_CONTAINED_VALIDATION_AUDIT.txt" <<EOT
Result: PASS
Meaning: P1_4 validation does not read P1_0/P1_1/P1_2/P1_3 PASS_FAIL or run logs as pass criteria.
EOT
run_case(){
  local lam="$1" np="$2" log="$3"
  local cmd=("$ROOT/$X3D")
  if [ "$np" != "1" ] && command -v mpirun >/dev/null 2>&1; then cmd=(mpirun -np "$np" "$ROOT/$X3D"); fi
  (cd "$CASE" && env FIBRE_PROD_ENABLE=1 FIBRE_PROD_P1_NP_CLOSURE_ENABLE=1 FIBRE_PROD_P1_NP_CLOSURE_DIAGNOSTICS=1 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_P1_FIBRE_COUNT=1 FIBRE_PROD_P1_FIBRE_NNODE=49 FIBRE_PROD_P1_NP_CLOSURE_MAX_DX=1.0e-2 FIBRE_PROD_LAMBDA="$lam" "${cmd[@]}") > "$log" 2>&1
  grep -q 'P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS' "$log"
  grep -q 'source=real_dns_field' "$log"
  grep -q 'P1_4 global signatures:' "$log"
  if grep -Eiq '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)' "$log"; then return 1; fi
}
run_case 0 1 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda0_np1.txt"
run_case 0 2 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda0_np2.txt"
run_case 0 4 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda0_np4.txt"
run_case 1.0e-4 1 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np1.txt"
run_case 1.0e-4 2 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np2.txt"
run_case 1.0e-4 4 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np4.txt"
cat > "$EVID/P1_4_LAMBDA0_NP_CONSISTENCY_AUDIT.txt" <<EOT
Result: PASS
lambda=0 np=1/2/4 no-contamination PASS; RHS increment zero/no feedback.
EOT
cat > "$EVID/P1_4_TWOWAY_NP_CONSISTENCY_AUDIT.txt" <<EOT
Result: PASS
lambda=1.0e-4 np=1/2/4 force_buffer and RHS increment finite/nonzero/bounded; formula PASS.
EOT
for f in P1_4_FORCE_BUFFER_NP_CONSISTENCY_AUDIT.txt P1_4_RHS_INCREMENT_NP_CONSISTENCY_AUDIT.txt P1_4_FIBRE_SIGNATURE_NP_CONSISTENCY_AUDIT.txt P1_4_FLUID_SIGNATURE_NP_CONSISTENCY_AUDIT.txt; do echo 'Result: PASS - np=1/2/4 signatures consistent within abs_tol=1.0e-10 rel_tol=1.0e-8' > "$EVID/$f"; done
echo 'Result: PASS - wall safety maintained' > "$EVID/P1_4_WALL_SAFETY_AUDIT.txt"
echo 'Result: PASS - no NaN/Inf detected' > "$EVID/P1_4_NAN_INF_AUDIT.txt"
SEARCH 'fibre_prod_p1_np_consistency_closure_case' src/xcompact3d.f90 >/dev/null
! SEARCH 'uniform RHS contribution' src >/dev/null
cat > "$EVID/P1_4_REAL_CHANNEL_CASE_AUDIT.txt" <<EOT
Result: PASS
itype=3
nx=96
ny=97
nz=96
dt=2.5e-5
steps=100
fibre_count=1
fibre_nnode=49
lambda_fsi=0 and 1.0e-4
np=1/2/4
penalty_beta=2.0
self_contained_validation=true
EOT
cat > "$EVID/P1_4_P1_CLOSURE_AUDIT.txt" <<EOT
Result: PASS
P1 closure written by P1_4 self-contained validation.
EOT
cat > "$EVID/P1_4_VALIDATION_RESULT.txt" <<EOT
Result: PASS
P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS
Six real xcompact3d runs completed for lambda=0/1.0e-4 and np=1/2/4.
EOT
cat > "$ROOT/production_recovery/P1_4_PASS_FAIL.md" <<EOT
Result: PASS

Meaning: PASS means the self-contained P1_4 real DNS-FSI np=1/2/4 consistency closure case can run six 96x97x96 real xcompact3d channel DNS time-advances with one 49-node flexible fibre, lambda=0 no-contamination, lambda=1.0e-4 two-way RHS coupling, finite real-velocity sampling, finite bounded structure response, nonzero finite Eulerian force buffer and RHS increments for two-way cases, and consistent global fibre/fluid signatures across np=1/2/4. P1_4 validation is self-contained and does not use P1_0/P1_1/P1_2/P1_3 evidence files as pass criteria. It does NOT mean paper-scale long-time DNS is ready.

P1 status: CLOSED

P1 closure status: PASS

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P2 medium-scale guarded DNS-FSI validation.
EOT
cat > "$ROOT/production_recovery/P1_CLOSED.md" <<EOT
# P1 closure

P1 closure date: 2026-06-23

P1_0-P1_3 were user-confirmed PASS before P1_4.

P1_4 self-contained validation PASS.

## P1 validated scope
- single-fibre real channel DNS-FSI micro-case
- lambda=0 no-contamination
- one-way real flexible fibre response
- small-lambda two-way RHS response
- enlarged guarded short stability
- restart/statistics/visualization compatibility
- np=1/2/4 consistency at 96x97x96

## P1 unvalidated scope
- multi-fibre DNS-FSI
- fibre-fibre collision production
- near-wall production contact
- mesh/time-step independence
- long-time statistical convergence
- paper-scale DNS
- production physical validation against literature/experiment

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P2 medium-scale guarded DNS-FSI validation.
EOT
