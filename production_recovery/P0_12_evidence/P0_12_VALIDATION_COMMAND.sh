#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/P0_12_evidence"
RESULT_FILE="$EVIDENCE_DIR/P0_12_VALIDATION_RESULT.txt"
PASS_FAIL="$ROOT_DIR/production_recovery/P0_12_PASS_FAIL.md"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
BUILD_DIR="$ROOT_DIR/build_p0_12"
: > "$RESULT_FILE"

log() { printf '%s\n' "$*" | tee -a "$RESULT_FILE"; }
failures=0
search() {
  if command -v rg >/dev/null 2>&1; then
    rg "$@"
  else
    pattern="$1"; shift
    files=()
    for path in "$@"; do
      if [[ -d "$path" ]]; then
        while IFS= read -r file; do files+=("$file"); done < <(find "$path" -type f)
      else
        files+=("$path")
      fi
    done
    [[ ${#files[@]} -gt 0 ]] || return 1
    grep -n -E "$pattern" "${files[@]}"
  fi
}
pass() { log "PASS: $*"; }
fail() { log "FAIL: $*"; failures=$((failures+1)); }

log "P0_12 validation started"
log "ROOT_DIR=$ROOT_DIR"
log "DECOMP2D_ROOT=$DECOMP2D_ROOT"
cd "$ROOT_DIR" || exit 1

[[ -f src/fibre_prod_synthetic_closed_loop.f90 ]] && pass "synthetic closed-loop source exists" || fail "synthetic closed-loop source missing"
[[ -f src/fibre_prod_synthetic_closed_loop_check.f90 ]] && pass "synthetic closed-loop check source exists" || fail "synthetic closed-loop check source missing"
search "fibre_prod_synthetic_closed_loop" src/xcompact3d.f90 >/dev/null && pass "xcompact3d imports synthetic closed-loop" || fail "xcompact3d import missing"
search "FIBRE_PROD_SYNTHETIC_CLOSED_LOOP_ENABLE" src/fibre_prod_synthetic_closed_loop.f90 src/xcompact3d.f90 >/dev/null && pass "P0_12 path is env gated" || fail "P0_12 env gate missing"
search "fibre_prod_synthetic_closed_loop_check" src/CMakeLists.txt >/dev/null && pass "P0_12 CMake target exists" || fail "P0_12 CMake target missing"
if search "contribution[[:space:]]*=[[:space:]]*.*lambda.*penalty_beta.*dt" src >/dev/null 2>&1; then
  fail "old uniform RHS contribution path found"
else
  pass "old uniform RHS contribution path absent"
fi
if search "pressure|projection|rk3|channel-forcing" src/fibre_prod_synthetic_closed_loop.f90 >/dev/null 2>&1; then
  fail "synthetic module references forbidden pressure/projection/RK3/channel-forcing path"
else
  pass "synthetic module has no pressure/projection/RK3/channel-forcing path"
fi

{
  echo "P0_12 synthetic closed-loop API audit"
  search "fibre_prod_synthetic_closed_loop_(init|run|get_signature|env_enabled)" src/fibre_prod_synthetic_closed_loop.f90 || true
} > "$EVIDENCE_DIR/P0_12_SYNTHETIC_CLOSED_LOOP_AUDIT.txt"
{
  echo "P0_12 lambda0 full path audit"
  search "lambda_fsi == 0.0_dp|max_abs_rhs_increment == 0.0_dp|lambda0_rhs_increment_norm" src/fibre_prod_synthetic_closed_loop.f90 src/fibre_prod_synthetic_closed_loop_check.f90 || true
} > "$EVIDENCE_DIR/P0_12_LAMBDA0_FULL_PATH_AUDIT.txt"
{
  echo "P0_12 small-lambda full path audit"
  search "1.0e-3_dp|small-lambda|force_buffer|rhs_increment" src/fibre_prod_synthetic_closed_loop.f90 src/fibre_prod_synthetic_closed_loop_check.f90 || true
} > "$EVIDENCE_DIR/P0_12_SMALL_LAMBDA_FULL_PATH_AUDIT.txt"
{
  echo "P0_12 xcompact3d synthetic call path audit"
  search "fibre_prod_synthetic_closed_loop" src/xcompact3d.f90 || true
} > "$EVIDENCE_DIR/P0_12_XCOMPACT3D_SYNTHETIC_CALL_PATH_AUDIT.txt"
{
  echo "P0_12 no uniform RHS audit"
  search "uniform|contribution|force_buffer|lambda_fsi.*penalty_beta" src/fibre_prod_synthetic_closed_loop.f90 src/fibre_prod_force_buffer_rhs_gate.f90 || true
} > "$EVIDENCE_DIR/P0_12_NO_UNIFORM_RHS_AUDIT.txt"

log "Configuring CMake build"
if cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" >> "$RESULT_FILE" 2>&1; then
  pass "cmake configure"
else
  fail "cmake configure failed"
fi

targets=(xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check fibre_prod_reaction_spreading_buffer_check fibre_prod_force_buffer_rhs_gate_check fibre_prod_synthetic_closed_loop_check)
if [[ $failures -eq 0 ]]; then
  for target in "${targets[@]}"; do
    if cmake --build "$BUILD_DIR" --target "$target" >> "$RESULT_FILE" 2>&1; then pass "built $target"; else fail "build failed $target"; fi
  done
fi

find_exe() {
  local target="$1"
  for candidate in "$BUILD_DIR/bin/$target" "$BUILD_DIR/src/$target" "$BUILD_DIR/$target"; do
    [[ -x "$candidate" ]] && { printf '%s\n' "$candidate"; return 0; }
  done
  find "$BUILD_DIR" -type f -name "$target" -perm -111 2>/dev/null | head -n 1
}
run_check() {
  local target="$1" exe
  exe="$(find_exe "$target")"
  [[ -n "$exe" ]] || { fail "executable not found $target"; return; }
  if [[ "$target" == "fibre_prod_force_buffer_rhs_path_check" ]]; then
    if FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-3 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0 "$exe" >> "$RESULT_FILE" 2>&1; then pass "ran $target"; else fail "run failed $target"; fi
  else
    if "$exe" >> "$RESULT_FILE" 2>&1; then pass "ran $target"; else fail "run failed $target"; fi
  fi
}

checks=(fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check fibre_prod_reaction_spreading_buffer_check fibre_prod_force_buffer_rhs_gate_check)
if [[ $failures -eq 0 ]]; then
  for check in "${checks[@]}"; do run_check "$check"; done
  exe="$(find_exe fibre_prod_synthetic_closed_loop_check)"
  if [[ -z "$exe" ]]; then
    fail "synthetic closed-loop executable not found"
  else
    for np in 1 2 4; do
      out="$EVIDENCE_DIR/P0_12_NP${np}_SIGNATURE.txt"
      if command -v mpirun >/dev/null 2>&1; then
        if mpirun -np "$np" "$exe" > "$out" 2>&1; then pass "ran synthetic check np=$np"; else fail "synthetic check failed np=$np"; fi
      elif [[ "$np" == "1" ]]; then
        if "$exe" > "$out" 2>&1; then pass "ran synthetic check np=1 without mpirun"; else fail "synthetic check failed np=1"; fi
      else
        echo "mpirun unavailable; np=$np required run failed" > "$out"
        fail "mpirun unavailable for np=$np"
      fi
      cat "$out" >> "$RESULT_FILE" 2>/dev/null || true
    done
  fi
fi


for np in 1 2 4; do
  sig_file="$EVIDENCE_DIR/P0_12_NP${np}_SIGNATURE.txt"
  if [[ ! -f "$sig_file" ]]; then
    {
      echo "P0_12_SIGNATURE np=${np} unavailable"
      echo "P0_12_SYNTHETIC_CLOSED_LOOP_CHECK FAIL: executable was not run"
    } > "$sig_file"
  fi
done
{
  echo "P0_12 np=1/2/4 consistency audit"
  if [[ -f "$EVIDENCE_DIR/P0_12_NP1_SIGNATURE.txt" && -f "$EVIDENCE_DIR/P0_12_NP2_SIGNATURE.txt" && -f "$EVIDENCE_DIR/P0_12_NP4_SIGNATURE.txt" ]]; then
    python3 - "$EVIDENCE_DIR" <<'PY'
import sys, re, math, pathlib
root=pathlib.Path(sys.argv[1])
keys=['sampled_norm','force_buffer_norm','rhs_increment_norm','lambda0_rhs_increment_norm','smalllambda_rhs_increment_norm','rhs_increment_max','conservation_scale_error']
def parse(path):
    vals={}
    for line in path.read_text(errors='ignore').splitlines():
        m=re.match(r'P0_12_SIGNATURE\s+([A-Za-z0-9_]+)=\s*([-+0-9.Ee]+)', line.strip())
        if m: vals[m.group(1)]=float(m.group(2))
    return vals
ref=parse(root/'P0_12_NP1_SIGNATURE.txt')
failed=False
for np in (2,4):
    cur=parse(root/f'P0_12_NP{np}_SIGNATURE.txt')
    for key in keys:
        a=ref.get(key); b=cur.get(key)
        if a is None or b is None or abs(a-b) > 1e-10 + 1e-10*max(abs(a), abs(b)):
            print(f'FAIL np={np} key={key} ref={a} cur={b}')
            failed=True
        else:
            print(f'PASS np={np} key={key} value={b:.16e}')
sys.exit(1 if failed else 0)
PY
  else
    echo "FAIL missing np signature files"
  fi
} > "$EVIDENCE_DIR/P0_12_NP124_CONSISTENCY_AUDIT.txt"
if [[ -f "$EVIDENCE_DIR/P0_12_NP124_CONSISTENCY_AUDIT.txt" ]] && grep -q '^FAIL' "$EVIDENCE_DIR/P0_12_NP124_CONSISTENCY_AUDIT.txt"; then
  fail "np=1/2/4 consistency audit failed"
elif [[ -f "$EVIDENCE_DIR/P0_12_NP124_CONSISTENCY_AUDIT.txt" ]]; then
  pass "np=1/2/4 consistency audit"
fi

required=(P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK P0_3_RUNTIME_BRIDGE_CHECK P0_4_VELOCITY_BRIDGE_CHECK P0_5_STATE_VELOCITY_ATTACHMENT_CHECK P0_6_HYDRO_INPUT_CANDIDATE_CHECK P0_7_STRUCTURE_INPUT_HANDOFF_CHECK P0_8_STRUCTURE_DRY_STEP_CHECK P0_9_STRUCTURE_COMMIT_GATE_CHECK P0_10_REACTION_SPREADING_BUFFER_CHECK P0_11_FORCE_BUFFER_RHS_GATE_CHECK P0_12_SYNTHETIC_CLOSED_LOOP_CHECK)
if [[ $failures -eq 0 ]]; then
  for token in "${required[@]}"; do
    if grep -q "$token" "$RESULT_FILE"; then pass "found token $token"; else fail "missing token $token"; fi
  done
fi

if [[ $failures -eq 0 ]]; then
  cat > "$PASS_FAIL" <<'PASSDOC'
Result: PASS

Meaning: PASS means the guarded production synthetic path from velocity sampling through fibre-state attachment, hydrodynamic input, structure dry-step/commit, reaction-force spreading, Eulerian force buffer, and lambda-gated RHS coupling can complete a deterministic end-to-end single-step closed-loop with strict lambda=0 no-contamination, bounded small-lambda RHS response, and np=1/2/4 signature consistency, while P0.2/P0.3/P0.4/P0.5/P0.6/P0.7/P0.8/P0.9/P0.10/P0.11 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0_13 no-fibre regression + restart/statistics/visualization + P0 closure/P1 entry gate.
PASSDOC
  log "Result: PASS"
  exit 0
else
  cat > "$PASS_FAIL" <<'FAILDOC'
Result: FAIL

Meaning: FAIL means one or more P0_12 static/build/run/np-consistency checks did not complete successfully. The synthetic closed-loop path must not be treated as production DNS-FSI readiness until this stage passes in an environment with the required toolchain and MPI runner.

Production-run status: STILL BLOCKED

Next required stage: Resolve P0_12 validation failures, then rerun P0_12 before P0_13.
FAILDOC
  log "Result: FAIL"
  exit 1
fi
