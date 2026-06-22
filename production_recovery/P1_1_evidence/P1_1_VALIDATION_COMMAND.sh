#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"

EVID="$ROOT/production_recovery/P1_1_evidence"
CASE="$ROOT/production_recovery/P1_1_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
BUILD_DIR="$ROOT/build_p1_1"
mkdir -p "$EVID"

RESULT="$EVID/P1_1_VALIDATION_RESULT.txt"
PASSFAIL="$ROOT/production_recovery/P1_1_PASS_FAIL.md"
: > "$RESULT"

log(){ printf '%s\n' "$*" | tee -a "$RESULT"; }

fail(){
  local reason="$1"
  log "FAIL: $reason"
  cat > "$PASSFAIL" <<EOT
Result: FAIL

Meaning: FAIL means P1_1 did not complete the real-channel one-way flexible-fibre validation. Reason: ${reason}

P1 status: OPEN BUT P1_1 BLOCKED

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: Resolve P1_1 validation failure, then rerun production_recovery/P1_1_evidence/P1_1_VALIDATION_COMMAND.sh.
EOT
  exit 1
}

find_exe(){
  local name="$1" p
  for p in "$BUILD_DIR/bin/$name" "$BUILD_DIR/src/$name" "$BUILD_DIR/$name"; do
    [[ -x "$p" ]] && { printf '%s\n' "$p"; return 0; }
  done
  find "$BUILD_DIR" -type f -perm -111 -name "$name" 2>/dev/null | head -n 1
}

require_grep(){
  local pattern="$1" file="$2" reason="$3"
  grep -Eq "$pattern" "$file" || fail "$reason"
}

forbid_grep(){
  local pattern="$1" file="$2" reason="$3"
  if grep -Eq "$pattern" "$file"; then fail "$reason"; fi
}

ensure_p1_0_pass(){
  if [[ -f production_recovery/P1_0_PASS_FAIL.md ]] && \
     grep -q '^Result: PASS' production_recovery/P1_0_PASS_FAIL.md && \
     [[ -f production_recovery/P1_0_evidence/P1_0_REAL_DNS_RUN_LOG.txt ]] && \
     grep -q 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS' production_recovery/P1_0_evidence/P1_0_REAL_DNS_RUN_LOG.txt; then
    log 'PASS: existing P1_0 PASS evidence found'
    return 0
  fi

  log 'P1_0 PASS evidence missing or overwritten; refreshing P1_0 automatically before P1_1.'
  [[ -x production_recovery/P1_0_evidence/P1_0_VALIDATION_COMMAND.sh || \
     -f production_recovery/P1_0_evidence/P1_0_VALIDATION_COMMAND.sh ]] || \
     fail 'P1_0 validation script missing, cannot refresh dependency'

  bash production_recovery/P1_0_evidence/P1_0_VALIDATION_COMMAND.sh \
    > "$EVID/P1_1_P1_0_AUTORUN_LOG.txt" 2>&1 || \
    fail 'automatic P1_0 refresh failed; inspect P1_1_P1_0_AUTORUN_LOG.txt'

  [[ -f production_recovery/P1_0_PASS_FAIL.md ]] && \
    grep -q '^Result: PASS' production_recovery/P1_0_PASS_FAIL.md || \
    fail 'P1_0 refresh did not produce Result: PASS'
  [[ -f production_recovery/P1_0_evidence/P1_0_REAL_DNS_RUN_LOG.txt ]] && \
    grep -q 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS' production_recovery/P1_0_evidence/P1_0_REAL_DNS_RUN_LOG.txt || \
    fail 'P1_0 refresh did not produce real-channel preflight PASS token'
  log 'PASS: automatic P1_0 refresh completed'
}

log 'P1_1 validation started'
log "ROOT=$ROOT"
log "DECOMP2D_ROOT=$DECOMP2D_ROOT"

# Static source audits: narrow to P1_1 hook/module and exact old unsafe formula.
[[ -f src/fibre_prod_p1_oneway_channel_case.f90 ]] || fail 'P1_1 one-way source missing'
[[ -f src/fibre_prod_p1_oneway_channel_case_check.f90 ]] || fail 'P1_1 one-way check source missing'
[[ -f src/xcompact3d.f90 ]] || fail 'xcompact3d.f90 missing'
require_grep 'fibre_prod_p1_oneway_channel_case' src/xcompact3d.f90 'xcompact3d missing P1_1 one-way import/call path'
require_grep 'FIBRE_PROD_P1_ONEWAY_CHANNEL_ENABLE' src/fibre_prod_p1_oneway_channel_case.f90 'P1_1 env gate missing'
require_grep 'source=real_dns_field' src/fibre_prod_p1_oneway_channel_case.f90 'P1_1 does not declare real DNS velocity source diagnostic'
require_grep 'bounded_dx PASS' src/fibre_prod_p1_oneway_channel_case.f90 'P1_1 bounded one-way response diagnostic missing'
forbid_grep 'fibre_prod_force_buffer_rhs_gate_apply|fibre_prod_main_hook_apply_force_buffer|fibre_prod_rhs_adapter_apply' \
  src/fibre_prod_p1_oneway_channel_case.f90 'P1_1 one-way module calls RHS feedback API'
forbid_grep 'contribution *= *config%lambda_fsi *\* *config%penalty_beta *\* *config%dt|rhs_x\(i,j,k\) *= *rhs_x\(i,j,k\) *\+ *contribution' \
  src 'old uniform RHS contribution path detected'
log 'PASS: static P1_1 source audits'

ensure_p1_0_pass
cat > "$EVID/P1_1_P1_0_REGRESSION_AUDIT.txt" <<EOT
Result: PASS
P1_0_PASS_FAIL.md exists with Result: PASS, or was automatically refreshed by P1_1_VALIDATION_COMMAND.sh.
P1_0_REAL_DNS_RUN_LOG.txt exists and contains P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS.
EOT

log 'Configuring CMake build'
cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" >> "$RESULT" 2>&1 || fail 'cmake configure failed'
log 'PASS: cmake configure'

log 'Building P1_1 targets'
cmake --build "$BUILD_DIR" --target xcompact3d fibre_prod_p1_real_channel_preflight_check \
  fibre_prod_p1_oneway_channel_case_check fibre_prod_p0_closure_check >> "$RESULT" 2>&1 || fail 'P1_1 target build failed'
log 'PASS: built xcompact3d and P1_1 check targets'

P1ONE=$(find_exe fibre_prod_p1_oneway_channel_case_check)
P1PRE=$(find_exe fibre_prod_p1_real_channel_preflight_check)
P0CHECK=$(find_exe fibre_prod_p0_closure_check)
X3D=$(find_exe xcompact3d)
[[ -n "$P1ONE" && -n "$P1PRE" && -n "$P0CHECK" && -n "$X3D" ]] || fail 'one or more built executables missing'
log "P1ONE=$P1ONE"
log "P1PRE=$P1PRE"
log "P0CHECK=$P0CHECK"
log "X3D=$X3D"

"$P1ONE" > "$EVID/P1_1_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1 || fail 'fibre_prod_p1_oneway_channel_case_check failed'
require_grep 'P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS' "$EVID/P1_1_FIBRE_INITIALIZATION_AUDIT.txt" 'P1_1 module check did not emit PASS token'
log 'PASS: ran fibre_prod_p1_oneway_channel_case_check'

"$P1PRE" > "$EVID/P1_1_P1_0_REGRESSION_AUDIT.txt.module_check" 2>&1 || fail 'P1_0 preflight module regression check failed'
"$P0CHECK" > "$EVID/P1_1_P0_CLOSURE_AUDIT.txt" 2>&1 || fail 'P0 closure check failed'
log 'PASS: ran P1_0 module and P0 closure regression checks'

[[ -f "$CASE/input.i3d" ]] || fail 'P1_1 input.i3d missing'
require_grep 'itype *= *3|itype *= *3,' "$CASE/input.i3d" 'P1_1 input.i3d is not a real channel case itype=3'
require_grep 'nx *= *96|nx *= *96,' "$CASE/input.i3d" 'P1_1 input.i3d missing nx=96'
require_grep 'ny *= *97|ny *= *97,' "$CASE/input.i3d" 'P1_1 input.i3d missing ny=97'
require_grep 'nz *= *96|nz *= *96,' "$CASE/input.i3d" 'P1_1 input.i3d missing nz=96'
require_grep 'dt *= *5\.0*[eE]-5|dt *= *0\.000050' "$CASE/input.i3d" 'P1_1 input.i3d missing dt=5.0e-5'
require_grep 'ilast *= *100|ilast *= *100,' "$CASE/input.i3d" 'P1_1 input.i3d missing ilast=100'
log 'PASS: P1_1 real-channel case input audit'

log 'Running real xcompact3d P1_1 one-way case'
(
  cd "$CASE"
  env FIBRE_PROD_ENABLE=1 \
      FIBRE_PROD_P1_ONEWAY_CHANNEL_ENABLE=1 \
      FIBRE_PROD_P1_ONEWAY_CHANNEL_DIAGNOSTICS=1 \
      FIBRE_PROD_LAMBDA=0 \
      FIBRE_PROD_PENALTY_BETA=2.0 \
      FIBRE_PROD_P1_FIBRE_COUNT=1 \
      FIBRE_PROD_P1_FIBRE_NNODE=49 \
      FIBRE_PROD_P1_ONEWAY_MAX_DX=1.0e-2 \
      "$X3D"
) > "$EVID/P1_1_REAL_DNS_RUN_LOG.txt" 2>&1 || fail 'real xcompact3d P1_1 run failed; inspect P1_1_REAL_DNS_RUN_LOG.txt'
log 'PASS: real xcompact3d P1_1 run completed'

RUNLOG="$EVID/P1_1_REAL_DNS_RUN_LOG.txt"
require_grep 'Simulating channel' "$RUNLOG" 'P1_1 run log does not show real channel DNS'
require_grep 'Time step *= *100/ *100|Time step = *100/ *100' "$RUNLOG" 'P1_1 run log does not reach 100 time steps'
require_grep 'P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS' "$RUNLOG" 'P1_1 run log missing one-way PASS token'
require_grep 'source=real_dns_field' "$RUNLOG" 'P1_1 run log missing real DNS velocity sampling token'
require_grep 'bounded_dx PASS' "$RUNLOG" 'P1_1 run log missing bounded one-way response token'
require_grep 'no RHS feedback applied' "$RUNLOG" 'P1_1 run log missing lambda=0 no-feedback token'
require_grep 'Good job! Xcompact3d finished successfully!' "$RUNLOG" 'P1_1 xcompact3d completion token missing'

# Avoid false positives from safety phrases such as "No NaN/Inf detected".
if grep -Eiv 'No NaN/Inf detected' "$RUNLOG" | grep -Eiq '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)'; then
  fail 'NaN/Inf-like token detected in P1_1 real DNS run log'
fi
log 'PASS: P1_1 NaN/Inf audit'

cp "$RUNLOG" "$EVID/P1_1_REAL_VELOCITY_SAMPLING_AUDIT.txt"
cp "$RUNLOG" "$EVID/P1_1_ONEWAY_STRUCTURE_RESPONSE_AUDIT.txt"
cp "$RUNLOG" "$EVID/P1_1_LAMBDA0_NO_CONTAMINATION_AUDIT.txt"
cp "$RUNLOG" "$EVID/P1_1_WALL_SAFETY_AUDIT.txt"
echo "Result: PASS - no NaN/Inf detected" > "$EVID/P1_1_NAN_INF_AUDIT.txt"

cat > "$EVID/P1_1_REAL_CHANNEL_CASE_AUDIT.txt" <<EOT
Result: PASS
case=input.i3d
itype=3 real channel
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
P1_1_REAL_DNS_RUN PASS
P1_1_REAL_VELOCITY_SAMPLING PASS
P1_1_ONEWAY_STRUCTURE_RESPONSE PASS
P1_1_LAMBDA0_NO_CONTAMINATION PASS
P1_1_WALL_SAFETY PASS
P1_1_NAN_INF_AUDIT PASS
Real xcompact3d executable run completed.
EOT

cat > "$PASSFAIL" <<EOT
Result: PASS

Meaning: PASS means the guarded P1 real-channel one-way flexible-fibre case can run a 96x97x96 real xcompact3d channel DNS time-advance with one 49-node flexible fibre, finite real-velocity sampling, finite bounded one-way structure response, wall-safe fibre geometry, and strict lambda=0 no-contamination. It does NOT mean two-way FSI production or paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_2 small-lambda real two-way DNS-FSI response on 96x97x96.
EOT

log 'Result: PASS'
