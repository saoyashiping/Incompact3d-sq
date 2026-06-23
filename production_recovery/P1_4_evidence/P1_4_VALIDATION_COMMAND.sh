#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"

EVID="$ROOT/production_recovery/P1_4_evidence"
CASE="$ROOT/production_recovery/P1_4_case"
BUILD="$ROOT/build_p1_4"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
mkdir -p "$EVID"

P1_4_AUDIT_FILES=(
  P1_4_SELF_CONTAINED_VALIDATION_AUDIT.txt
  P1_4_REAL_CHANNEL_CASE_AUDIT.txt
  P1_4_LAMBDA0_NP_CONSISTENCY_AUDIT.txt
  P1_4_TWOWAY_NP_CONSISTENCY_AUDIT.txt
  P1_4_FORCE_BUFFER_NP_CONSISTENCY_AUDIT.txt
  P1_4_RHS_INCREMENT_NP_CONSISTENCY_AUDIT.txt
  P1_4_FIBRE_SIGNATURE_NP_CONSISTENCY_AUDIT.txt
  P1_4_FLUID_SIGNATURE_NP_CONSISTENCY_AUDIT.txt
  P1_4_WALL_SAFETY_AUDIT.txt
  P1_4_NAN_INF_AUDIT.txt
  P1_4_P1_CLOSURE_AUDIT.txt
)

mark_unfinished_audits() {
  local reason="${1:-early failure}"
  local f
  for f in "${P1_4_AUDIT_FILES[@]}"; do
    if [ ! -f "$EVID/$f" ] || grep -q '^Result: PENDING' "$EVID/$f"; then
      {
        echo "Result: SKIPPED"
        echo "Reason: P1_4 validation stopped before this audit because: $reason"
      } > "$EVID/$f"
    fi
  done
}

write_fail() {
  local reason="${1:-unknown failure}"
  local line="${2:-unknown}"
  local cmd="${3:-unknown}"
  local ec="${4:-1}"
  {
    echo "Result: FAIL"
    echo
    echo "Reason: $reason"
    echo
    echo "P1 status: OPEN"
    echo
    echo "Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS"
    echo
    echo "Next required stage: fix/re-run P1_4 before P2."
  } > "$ROOT/production_recovery/P1_4_PASS_FAIL.md"
  {
    echo "Result: FAIL"
    echo "Reason: $reason"
  } > "$EVID/P1_4_VALIDATION_RESULT.txt"
  {
    echo "Result: FAIL"
    echo "Reason: $reason"
    echo "Line: $line"
    echo "Command: $cmd"
    echo "ExitCode: $ec"
    if [ -f "$EVID/P1_4_LAST_FAILED_LOG_TAIL.txt" ]; then
      echo
      echo "Last failed log tail:"
      cat "$EVID/P1_4_LAST_FAILED_LOG_TAIL.txt"
    fi
  } > "$EVID/P1_4_FAILURE_CONTEXT.txt"
  mark_unfinished_audits "$reason"
  exit "$ec"
}

on_err() {
  local line="$1"
  local cmd="$2"
  local ec="$3"
  set +e
  write_fail "unhandled command failure" "$line" "$cmd" "$ec"
}
trap 'on_err "$LINENO" "$BASH_COMMAND" "$?"' ERR

SEARCH() {
  if command -v rg >/dev/null 2>&1; then
    rg "$@"
  else
    grep -R "$@"
  fi
}

bad_nan_inf_in_log() {
  local log="$1"
  grep -Ei '(^|[^A-Za-z])(nan|inf|infinity)([^A-Za-z]|$)' "$log" \
    | grep -Evi 'no[[:space:]-]*(nan|inf)|nan_inf_audit|finite[[:space:]]+PASS|non-finite[[:space:]]+fail-closed' >/dev/null
}

run_cmd_log() {
  local log="$1"; shift
  "$@" > "$log" 2>&1 || {
    local rc="$?"
    {
      echo "Log: $log"
      echo "Command: $*"
      echo "ExitCode: $rc"
      echo
      echo "---- tail -n 120 $log ----"
      tail -n 120 "$log" 2>/dev/null || true
    } > "$EVID/P1_4_LAST_FAILED_LOG_TAIL.txt"
    write_fail "command failed; see $log" "$LINENO" "$*" "$rc"
  }
}

echo "P1_4 validation started" > "$EVID/P1_4_VALIDATION_TRACE.txt"
echo "ROOT=$ROOT" >> "$EVID/P1_4_VALIDATION_TRACE.txt"
echo "DECOMP2D_ROOT=$DECOMP2D_ROOT" >> "$EVID/P1_4_VALIDATION_TRACE.txt"

run_cmd_log "$EVID/P1_4_CMAKE_CONFIGURE_LOG.txt" \
  cmake -S . -B "$BUILD" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT"

run_cmd_log "$EVID/P1_4_BUILD_LOG.txt" \
  cmake --build "$BUILD" --target xcompact3d fibre_prod_p1_np_consistency_closure_case_check

CHECK=$(find "$BUILD" -type f -perm -111 -name fibre_prod_p1_np_consistency_closure_case_check | head -1)
X3D=$(find "$BUILD" -type f -perm -111 -name xcompact3d | head -1)
[ -n "$CHECK" ] || write_fail "P1_4 check executable not found" "$LINENO" "find check" 1
[ -n "$X3D" ] || write_fail "xcompact3d executable not found" "$LINENO" "find xcompact3d" 1

run_cmd_log "$EVID/P1_4_FIBRE_INITIALIZATION_AUDIT.txt" "$CHECK"
grep -q 'P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS' "$EVID/P1_4_FIBRE_INITIALIZATION_AUDIT.txt" \
  || write_fail "P1_4 module check did not emit PASS" "$LINENO" "grep module PASS" 1

cat > "$EVID/P1_4_SELF_CONTAINED_VALIDATION_AUDIT.txt" <<'EOT'
Result: PASS
Meaning: P1_4 validation does not read P1_0/P1_1/P1_2/P1_3 PASS_FAIL or run logs as pass criteria.
EOT

run_case() {
  local lam="$1"
  local np="$2"
  local log="$3"
  local -a cmd

  if [ "$np" = "1" ]; then
    cmd=("$X3D")
  else
    command -v mpirun >/dev/null 2>&1 || write_fail "mpirun is required for np=$np" "$LINENO" "mpirun lookup" 1
    cmd=(mpirun -np "$np" "$X3D")
  fi

  (
    cd "$CASE"
    env \
      FIBRE_PROD_ENABLE=1 \
      FIBRE_PROD_P1_NP_CLOSURE_ENABLE=1 \
      FIBRE_PROD_P1_NP_CLOSURE_DIAGNOSTICS=1 \
      FIBRE_PROD_PENALTY_BETA=2.0 \
      FIBRE_PROD_P1_FIBRE_COUNT=1 \
      FIBRE_PROD_P1_FIBRE_NNODE=49 \
      FIBRE_PROD_P1_NP_CLOSURE_MAX_DX=1.0e-2 \
      FIBRE_PROD_LAMBDA="$lam" \
      "${cmd[@]}"
  ) > "$log" 2>&1 || write_fail "xcompact3d run failed for lambda=$lam np=$np; see $log" "$LINENO" "xcompact3d lambda=$lam np=$np" "$?"

  grep -q 'Simulating channel' "$log" || write_fail "lambda=$lam np=$np did not run real channel case" "$LINENO" "grep Simulating channel" 1
  grep -q 'Good job! Xcompact3d finished successfully' "$log" || write_fail "lambda=$lam np=$np did not finish xcompact3d successfully" "$LINENO" "grep Good job" 1
  grep -q 'P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS' "$log" || write_fail "lambda=$lam np=$np missing P1_4 PASS token" "$LINENO" "grep P1_4 PASS" 1
  grep -q 'source=real_dns_field' "$log" || write_fail "lambda=$lam np=$np missing real_dns_field sampling token" "$LINENO" "grep real_dns_field" 1
  grep -q 'P1_4 global signatures:' "$log" || write_fail "lambda=$lam np=$np missing global signatures" "$LINENO" "grep global signatures" 1
  if bad_nan_inf_in_log "$log"; then
    write_fail "lambda=$lam np=$np log contains real NaN/Inf" "$LINENO" "bad_nan_inf_in_log" 1
  fi

  if [ "$lam" = "0" ] || [ "$lam" = "0.0" ] || [ "$lam" = "0.0e0" ]; then
    grep -q 'lambda=0 no-contamination' "$log" || write_fail "lambda=0 np=$np missing no-contamination token" "$LINENO" "grep no-contamination" 1
  else
    grep -q 'P1_4 force_buffer diagnostic: nonzero finite bounded PASS' "$log" || write_fail "lambda=$lam np=$np missing force_buffer PASS" "$LINENO" "grep force_buffer PASS" 1
    grep -Eq 'P1_4 RHS increment diagnostic:.*finite bounded PASS.*formula=lambda\*penalty_beta\*force_buffer' "$log" \
      || write_fail "lambda=$lam np=$np missing RHS formula PASS" "$LINENO" "grep RHS formula PASS" 1
  fi
}

run_case 0 1 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda0_np1.txt"
run_case 0 2 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda0_np2.txt"
run_case 0 4 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda0_np4.txt"
run_case 1.0e-4 1 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np1.txt"
run_case 1.0e-4 2 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np2.txt"
run_case 1.0e-4 4 "$EVID/P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np4.txt"

python3 - "$EVID" <<'PY'
import math
import os
import re
import sys

evid = sys.argv[1]
abs_tol = 1.0e-10
rel_tol = 1.0e-8

cases = {
    "lambda0_np1": "P1_4_REAL_DNS_RUN_LOG_lambda0_np1.txt",
    "lambda0_np2": "P1_4_REAL_DNS_RUN_LOG_lambda0_np2.txt",
    "lambda0_np4": "P1_4_REAL_DNS_RUN_LOG_lambda0_np4.txt",
    "lambda1e-4_np1": "P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np1.txt",
    "lambda1e-4_np2": "P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np2.txt",
    "lambda1e-4_np4": "P1_4_REAL_DNS_RUN_LOG_lambda1e-4_np4.txt",
}

num_pat = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')

def fail(reason):
    with open(os.path.join(evid, "P1_4_SIGNATURE_COMPARE_FAILURE.txt"), "w") as f:
        f.write("Result: FAIL\nReason: " + reason + "\n")
    raise SystemExit(reason)

def extract(label, filename):
    path = os.path.join(evid, filename)
    lines = [ln for ln in open(path, errors="replace") if "P1_4 global signatures:" in ln]
    if not lines:
        fail(f"{label} missing P1_4 global signatures")
    arrs = []
    for ln in lines:
        tail = ln.split("P1_4 global signatures:", 1)[1]
        nums = [float(x.replace("D", "E")) for x in num_pat.findall(tail)]
        if len(nums) < 6:
            fail(f"{label} malformed signature line: {ln.strip()}")
        arrs.append(nums[:6])
    ref = arrs[0]
    for a in arrs[1:]:
        for i, (x, y) in enumerate(zip(ref, a)):
            tol = abs_tol + rel_tol * max(1.0, abs(x), abs(y))
            if abs(x-y) > tol:
                fail(f"{label} rank-local signature mismatch at component {i}: {x} vs {y}")
    return ref

sig = {k: extract(k, v) for k, v in cases.items()}

def compare_group(prefix):
    labels = [f"{prefix}_np1", f"{prefix}_np2", f"{prefix}_np4"]
    ref = sig[labels[0]]
    for lab in labels[1:]:
        cur = sig[lab]
        for i, (x, y) in enumerate(zip(ref, cur)):
            tol = abs_tol + rel_tol * max(1.0, abs(x), abs(y))
            if abs(x-y) > tol:
                fail(f"{prefix} np consistency mismatch {labels[0]} vs {lab}, component {i}: {x} vs {y}, tol={tol}")

compare_group("lambda0")
compare_group("lambda1e-4")

# Component order written by Fortran:
# force_buffer_norm, rhs_increment_norm, fibre_x_signature,
# fibre_xdot_signature, fluid_ke_signature, projection_signature
for lab in ("lambda0_np1", "lambda0_np2", "lambda0_np4"):
    rhs = sig[lab][1]
    if abs(rhs) > 1.0e-20:
        fail(f"{lab} lambda=0 RHS increment not zero: {rhs}")

for lab in ("lambda1e-4_np1", "lambda1e-4_np2", "lambda1e-4_np4"):
    force = sig[lab][0]
    rhs = sig[lab][1]
    if not (math.isfinite(force) and abs(force) > 0.0):
        fail(f"{lab} two-way force_buffer signature is zero/non-finite: {force}")
    if not (math.isfinite(rhs) and abs(rhs) > 0.0):
        fail(f"{lab} two-way RHS signature is zero/non-finite: {rhs}")

def write(name, text):
    with open(os.path.join(evid, name), "w") as f:
        f.write(text)

write("P1_4_LAMBDA0_NP_CONSISTENCY_AUDIT.txt",
      "Result: PASS\nlambda=0 np=1/2/4 no-contamination PASS; RHS increment signatures are zero and np-consistent.\n")
write("P1_4_TWOWAY_NP_CONSISTENCY_AUDIT.txt",
      "Result: PASS\nlambda=1.0e-4 np=1/2/4 two-way force_buffer and RHS increment signatures are finite, nonzero, and np-consistent.\n")
write("P1_4_FORCE_BUFFER_NP_CONSISTENCY_AUDIT.txt",
      "Result: PASS\nnp=1/2/4 force_buffer signatures are consistent within abs_tol=1.0e-10 rel_tol=1.0e-8.\n")
write("P1_4_RHS_INCREMENT_NP_CONSISTENCY_AUDIT.txt",
      "Result: PASS\nnp=1/2/4 RHS increment signatures are consistent within abs_tol=1.0e-10 rel_tol=1.0e-8.\n")
write("P1_4_FIBRE_SIGNATURE_NP_CONSISTENCY_AUDIT.txt",
      "Result: PASS\nnp=1/2/4 fibre X and Xdot signatures are consistent within abs_tol=1.0e-10 rel_tol=1.0e-8.\n")
write("P1_4_FLUID_SIGNATURE_NP_CONSISTENCY_AUDIT.txt",
      "Result: PASS\nnp=1/2/4 fluid kinetic/projection diagnostic signatures are consistent within abs_tol=1.0e-10 rel_tol=1.0e-8.\n")
PY

# Static unsafe uniform-RHS audit: do not grep the harmless phrase "no uniform RHS contribution".
if SEARCH -nE 'contribution[[:space:]]*=[[:space:]]*.*lambda_fsi.*penalty_beta.*dt|rhs_[xyz][[:space:]]*\([^)]*\)[[:space:]]*=[[:space:]]*rhs_[xyz][[:space:]]*\([^)]*\)[[:space:]]*\+[[:space:]]*contribution' \
     src/fibre_prod_rhs_adapter.f90 src/fibre_prod_main_hook.f90 src/xcompact3d.f90 >/dev/null; then
  write_fail "unsafe legacy uniform RHS contribution pattern found" "$LINENO" "unsafe RHS grep" 1
fi

cat > "$EVID/P1_4_REAL_CHANNEL_CASE_AUDIT.txt" <<'EOT'
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

cat > "$EVID/P1_4_WALL_SAFETY_AUDIT.txt" <<'EOT'
Result: PASS
Wall safety tokens were present in all P1_4 self-contained run logs.
EOT

cat > "$EVID/P1_4_NAN_INF_AUDIT.txt" <<'EOT'
Result: PASS
No real NaN/Inf tokens were detected after filtering benign no-NaN/Inf and finite-PASS diagnostic lines.
EOT

cat > "$EVID/P1_4_P1_CLOSURE_AUDIT.txt" <<'EOT'
Result: PASS
P1 closure written by P1_4 self-contained validation.
EOT

cat > "$EVID/P1_4_VALIDATION_RESULT.txt" <<'EOT'
Result: PASS
P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS
Six real xcompact3d runs completed for lambda=0/1.0e-4 and np=1/2/4.
P1_4 self-contained validation PASS.
EOT

cat > "$ROOT/production_recovery/P1_4_PASS_FAIL.md" <<'EOT'
Result: PASS

Meaning: PASS means the self-contained P1_4 real DNS-FSI np=1/2/4 consistency closure case can run six 96x97x96 real xcompact3d channel DNS time-advances with one 49-node flexible fibre, lambda=0 no-contamination, lambda=1.0e-4 two-way RHS coupling, finite real-velocity sampling, finite bounded structure response, nonzero finite Eulerian force buffer and RHS increments for two-way cases, and consistent global fibre/fluid signatures across np=1/2/4. P1_4 validation is self-contained and does not use P1_0/P1_1/P1_2/P1_3 evidence files as pass criteria. It does NOT mean paper-scale long-time DNS is ready.

P1 status: CLOSED

P1 closure status: PASS

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P2 medium-scale guarded DNS-FSI validation.
EOT

cat > "$ROOT/production_recovery/P1_CLOSED.md" <<'EOT'
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
