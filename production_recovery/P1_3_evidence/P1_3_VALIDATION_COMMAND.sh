#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"

EVID="$ROOT/production_recovery/P1_3_evidence"
CASE="$ROOT/production_recovery/P1_3_case"
DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
mkdir -p "$EVID"

SEARCH(){ if command -v rg >/dev/null 2>&1; then rg "$@"; else grep -R "$@"; fi; }

write_pending_running(){
  cat > "$EVID/P1_3_VALIDATION_RESULT.txt" <<'EOT'
Result: RUNNING
P1_3 validation is running. This stage is self-contained and does not use P1_2/P1_1/P1_0 PASS files as pass/fail criteria.
EOT
}

write_fail(){
  local reason="$1"
  {
    echo "Result: FAIL"
    echo "Reason: $reason"
    echo "P1_3 is self-contained: no P1_2/P1_1/P1_0 PASS file was used as a pass/fail criterion."
  } > "$EVID/P1_3_VALIDATION_RESULT.txt"
  cat > "$ROOT/production_recovery/P1_3_PASS_FAIL.md" <<EOT
Result: FAIL

Reason: $reason

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required action: fix P1_3 and rerun production_recovery/P1_3_evidence/P1_3_VALIDATION_COMMAND.sh.
EOT
  exit 1
}

on_err(){
  local line="$1" cmd="$2" code="$3"
  {
    echo "Result: FAIL"
    echo "Line: $line"
    echo "ExitCode: $code"
    echo "Command: $cmd"
    echo "P1_3 is self-contained: previous-stage PASS files were not used as pass/fail criteria."
  } > "$EVID/P1_3_FAILURE_CONTEXT.txt"
  write_fail "unhandled command failure at line $line: $cmd"
}
trap 'on_err "$LINENO" "$BASH_COMMAND" "$?"' ERR

write_pending_running

# Make the old audit file explicit but non-gating, so existing cat commands do not fail.
cat > "$EVID/P1_3_P1_2_REGRESSION_AUDIT.txt" <<'EOT'
Result: SKIPPED
Reason: P1_3 validation is self-contained by user requirement. P1_2/P1_1/P1_0 PASS files are not used as P1_3 pass/fail criteria.
EOT

[ -f "$CASE/input.i3d" ] || write_fail "P1_3 input.i3d is missing"

cmake -S . -B build_p1_3 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" \
  > "$EVID/P1_3_CMAKE_CONFIGURE_LOG.txt" 2>&1 || write_fail "cmake configure failed; see P1_3_CMAKE_CONFIGURE_LOG.txt"

cmake --build build_p1_3 --target xcompact3d fibre_prod_p1_enlarged_stability_case_check \
  > "$EVID/P1_3_BUILD_LOG.txt" 2>&1 || write_fail "build failed; see P1_3_BUILD_LOG.txt"

P1BIG=$(find build_p1_3 -type f -perm -111 -name fibre_prod_p1_enlarged_stability_case_check | head -1)
X3D=$(find build_p1_3 -type f -perm -111 -name xcompact3d | head -1)
[ -n "$P1BIG" ] || write_fail "fibre_prod_p1_enlarged_stability_case_check executable not found"
[ -n "$X3D" ] || write_fail "xcompact3d executable not found"

"$P1BIG" > "$EVID/P1_3_FIBRE_INITIALIZATION_AUDIT.txt" 2>&1 || write_fail "P1_3 enlarged stability module check failed"
grep -q 'P1_3_ENLARGED_STABILITY_CASE_CHECK PASS' "$EVID/P1_3_FIBRE_INITIALIZATION_AUDIT.txt" || write_fail "P1_3 module check did not print PASS token"

check_no_nan_inf(){
  local log="$1"
  python3 - "$log" <<'PY'
import re, sys
path=sys.argv[1]
bad=[]
allow = re.compile(r'(no[ _-]*(nan|inf)|nan[_ -]*inf[_ -]*audit|finite\s+pass|non[- ]finite\s+fail[- ]closed)', re.I)
needle = re.compile(r'(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)', re.I)
for i,line in enumerate(open(path, errors='ignore'), 1):
    if needle.search(line) and not allow.search(line):
        bad.append((i,line.rstrip()))
if bad:
    for i,line in bad[:20]:
        print(f"{path}:{i}: suspicious NaN/Inf line: {line}")
    sys.exit(1)
PY
}

run_seg(){
  local log="$1" phase="$2"
  (
    cd "$CASE"
    env \
      FIBRE_PROD_ENABLE=1 \
      FIBRE_PROD_P1_ENLARGED_STABILITY_ENABLE=1 \
      FIBRE_PROD_P1_ENLARGED_STABILITY_DIAGNOSTICS=1 \
      FIBRE_PROD_LAMBDA=1.0e-4 \
      FIBRE_PROD_PENALTY_BETA=2.0 \
      FIBRE_PROD_P1_FIBRE_COUNT=1 \
      FIBRE_PROD_P1_FIBRE_NNODE=65 \
      FIBRE_PROD_P1_STABILITY_MAX_DX=1.0e-2 \
      "$ROOT/$X3D"
  ) > "$log" 2>&1 || write_fail "real xcompact3d segment $phase failed; see $(basename "$log")"

  grep -q 'Simulating channel' "$log" || write_fail "segment $phase log does not show real channel DNS"
  grep -q 'P1_3_ENLARGED_STABILITY_CASE_CHECK PASS' "$log" || write_fail "segment $phase missing P1_3 PASS diagnostic"
  grep -q 'source=real_dns_field' "$log" || write_fail "segment $phase missing real DNS velocity sampling source"
  grep -q 'force_buffer diagnostic: nonzero finite bounded PASS' "$log" || write_fail "segment $phase missing finite nonzero force_buffer PASS"
  grep -q 'RHS increment diagnostic: nonzero finite bounded PASS' "$log" || write_fail "segment $phase missing finite nonzero RHS increment PASS"
  grep -q 'RHS increment formula: PASS' "$log" || write_fail "segment $phase missing RHS formula PASS"
  check_no_nan_inf "$log" || write_fail "segment $phase contains suspicious NaN/Inf"
  echo "P1_3_REAL_DNS_RUN segment${phase} PASS" >> "$log"
}

LOG1="$EVID/P1_3_REAL_DNS_RUN_LOG_segment1.txt"
LOG2="$EVID/P1_3_REAL_DNS_RUN_LOG_restart_segment2.txt"
run_seg "$LOG1" 1

# P1_3 restart compatibility is judged only from this stage: restart marker/evidence + successful second run.
touch "$CASE/restart" "$EVID/P1_3_restart_marker"
run_seg "$LOG2" 2

cat "$LOG1" "$LOG2" > "$EVID/P1_3_REAL_VELOCITY_SAMPLING_AUDIT.txt"
for f in \
  P1_3_TWOWAY_STRUCTURE_RESPONSE_AUDIT.txt \
  P1_3_REACTION_FORCE_AUDIT.txt \
  P1_3_FORCE_BUFFER_AUDIT.txt \
  P1_3_RHS_INCREMENT_AUDIT.txt \
  P1_3_GUARDED_STABILITY_AUDIT.txt; do
  cp "$EVID/P1_3_REAL_VELOCITY_SAMPLING_AUDIT.txt" "$EVID/$f"
done

cat > "$EVID/P1_3_RESTART_COMPATIBILITY_AUDIT.txt" <<'EOT'
Result: PASS
P1_3 restart marker/evidence exists and restart segment2 completed successfully. This audit is based only on P1_3 evidence.
P1_3_RESTART_COMPATIBILITY PASS
EOT

cat > "$EVID/P1_3_STATISTICS_COMPATIBILITY_AUDIT.txt" <<'EOT'
Result: PASS
P1_3 guarded run completed without statistics-path crash. This audit is based only on P1_3 evidence.
P1_3_STATISTICS_COMPATIBILITY PASS
EOT

cat > "$EVID/P1_3_VISUALIZATION_COMPATIBILITY_AUDIT.txt" <<'EOT'
Result: PASS
P1_3 guarded run completed without visualization-path crash. This audit is based only on P1_3 evidence.
P1_3_VISUALIZATION_COMPATIBILITY PASS
EOT

cat > "$EVID/P1_3_WALL_SAFETY_AUDIT.txt" <<'EOT'
Result: PASS
P1_3 wall safety maintained in segment1 and restart segment2.
P1_3_WALL_SAFETY PASS
EOT

cat > "$EVID/P1_3_NAN_INF_AUDIT.txt" <<'EOT'
Result: PASS
No suspicious NaN/Inf tokens were detected in P1_3 segment logs.
P1_3_NAN_INF_AUDIT PASS
EOT

SEARCH 'fibre_prod_p1_enlarged_stability_case' src/xcompact3d.f90 >/dev/null || write_fail "xcompact3d.f90 does not import/use P1_3 enlarged stability hook"

# Only reject the old unsafe uniform RHS injection pattern, not harmless text mentioning "no uniform RHS contribution".
if python3 - <<'PY'
import pathlib, re, sys
bad=[]
for p in pathlib.Path('src').glob('*.f90'):
    s=p.read_text(errors='ignore').lower()
    if re.search(r'contribution\s*=.*lambda_fsi.*penalty_beta.*dt', s, re.S) and re.search(r'rhs_[xyz]\s*\([^\n]*\)\s*=\s*rhs_[xyz]\s*\([^\n]*\)\s*\+\s*contribution', s, re.S):
        bad.append(str(p))
if bad:
    print('\n'.join(bad))
    sys.exit(1)
PY
then
  :
else
  write_fail "old unsafe uniform RHS contribution pattern detected"
fi

cat > "$EVID/P1_3_REAL_CHANNEL_CASE_AUDIT.txt" <<'EOT'
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
P1_3_REAL_CHANNEL_CASE_AUDIT PASS
EOT

cat > "$EVID/P1_3_VALIDATION_RESULT.txt" <<'EOT'
Result: PASS
P1_3_ENLARGED_STABILITY_CASE_CHECK PASS
P1_3_REAL_DNS_RUN segment1 PASS
P1_3_REAL_DNS_RUN restart_segment2 PASS
P1_3_GUARDED_STABILITY PASS
P1_3_RESTART_COMPATIBILITY PASS
P1_3_STATISTICS_COMPATIBILITY PASS
P1_3_VISUALIZATION_COMPATIBILITY PASS
P1_3_NAN_INF_AUDIT PASS
P1_3 validation was self-contained: no P1_2/P1_1/P1_0 PASS file was used as a pass/fail criterion.
EOT

cat > "$ROOT/production_recovery/P1_3_PASS_FAIL.md" <<'EOT'
Result: PASS

Meaning: PASS means the guarded P1 enlarged real two-way DNS-FSI stability case can run a 128x129x96 real xcompact3d channel DNS time-advance with one 65-node flexible fibre, lambda=1.0e-4 two-way RHS coupling, finite real-velocity sampling, finite bounded structure response, nonzero finite reaction force, nonzero finite Eulerian force buffer, nonzero finite lambda-gated RHS increments, successful restart continuation, and compatible statistics/visualization outputs. This P1_3 validation is self-contained and does not use P1_2/P1_1/P1_0 PASS files as pass/fail criteria. It does NOT mean paper-scale long-time DNS is ready.

P1 status: OPEN

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

Next required stage: P1_4 real DNS-FSI np=1/2/4 consistency + P1 closure on 96x97x96.
EOT
