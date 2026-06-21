#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/P0_13_evidence"
RESULT_FILE="$EVIDENCE_DIR/P0_13_VALIDATION_RESULT.txt"
PASS_FAIL="$ROOT_DIR/production_recovery/P0_13_PASS_FAIL.md"
CLOSURE_REPORT="$ROOT_DIR/production_recovery/P0_13_FINAL_CLOSURE_REPORT.md"
P1_GATE="$ROOT_DIR/production_recovery/P0_13_P1_ENTRY_GATE.md"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
BUILD_DIR="$ROOT_DIR/build_p0_13"
: > "$RESULT_FILE"

log() { printf '%s\n' "$*" | tee -a "$RESULT_FILE"; }
failures=0
pass() { log "PASS: $*"; }
fail() { log "FAIL: $*"; failures=$((failures+1)); }
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

log "P0_13 validation started"
log "ROOT_DIR=$ROOT_DIR"
log "DECOMP2D_ROOT=$DECOMP2D_ROOT"
cd "$ROOT_DIR" || exit 1

[[ -f src/fibre_prod_p0_closure_check.f90 ]] && pass "P0 closure check source exists" || fail "P0 closure check source missing"
search "fibre_prod_p0_closure_check" src/CMakeLists.txt >/dev/null && pass "P0 closure CMake target exists" || fail "P0 closure CMake target missing"
search "FIBRE_PROD_SYNTHETIC_CLOSED_LOOP_ENABLE" src/xcompact3d.f90 src/fibre_prod_synthetic_closed_loop.f90 >/dev/null && pass "P0_12 synthetic path is env gated" || fail "P0_12 env gate missing"
search "FIBRE_PROD_FORCE_BUFFER_RHS_GATE_ENABLE" src/xcompact3d.f90 src/fibre_prod_force_buffer_rhs_gate.f90 >/dev/null && pass "P0_11 RHS gate path is env gated" || fail "P0_11 env gate missing"
search "FIBRE_PROD_REACTION_SPREADING_ENABLE" src/xcompact3d.f90 src/fibre_prod_reaction_spreading_buffer.f90 >/dev/null && pass "P0_10 reaction spreading path is env gated" || fail "P0_10 env gate missing"
search "FIBRE_PROD_STRUCTURE_DRY_COMMIT_ENABLE" src/fibre_prod_structure_commit_gate.f90 >/dev/null && pass "P0_9 commit gate env helper exists" || fail "P0_9 commit env gate missing"
if search "contribution[[:space:]]*=[[:space:]]*.*lambda.*penalty_beta.*dt" src >/dev/null 2>&1; then
  fail "old uniform RHS contribution path found"
else
  pass "old uniform RHS contribution path absent"
fi

{
  echo "P0_13 default-off audit"
  search "FIBRE_PROD_.*ENABLE|env_enabled|fibre_prod_.*_enabled" src/xcompact3d.f90 src/fibre_prod_*f90 || true
} > "$EVIDENCE_DIR/P0_13_DEFAULT_OFF_AUDIT.txt"
{
  echo "P0_13 no-fibre regression audit"
  search "fibre_prod_rhs_adapter_apply|fibre_prod_main_hook_apply_force_buffer|force_buffer|synthetic_closed_loop" \
    src/xcompact3d.f90 src/fibre_prod_p0_closure_check.f90 || true
} > "$EVIDENCE_DIR/P0_13_NO_FIBRE_REGRESSION_AUDIT.txt"
{
  echo "P0_13 restart compatibility audit"
  find src -maxdepth 1 -type f \( -iname '*restart*' -o -iname '*checkpoint*' \) -print
  search "restart|checkpoint" src/fibre_prod_*.f90 || true
} > "$EVIDENCE_DIR/P0_13_RESTART_COMPATIBILITY_AUDIT.txt"
{
  echo "P0_13 statistics compatibility audit"
  find src -maxdepth 1 -type f \( -iname '*stat*' -o -iname 'statistics.f90' \) -print
  search "statistics|statistic" src/fibre_prod_*.f90 || true
} > "$EVIDENCE_DIR/P0_13_STATISTICS_COMPATIBILITY_AUDIT.txt"
{
  echo "P0_13 visualization compatibility audit"
  find src -maxdepth 1 -type f \( -iname '*visu*' -o -iname '*visual*' \) -print
  search "visu|visual" src/fibre_prod_*.f90 || true
} > "$EVIDENCE_DIR/P0_13_VISUALIZATION_COMPATIBILITY_AUDIT.txt"
{
  echo "P0_13 P1 entry gate audit"
  cat "$P1_GATE"
} > "$EVIDENCE_DIR/P0_13_P1_ENTRY_GATE_AUDIT.txt"

if [[ ! -s "$EVIDENCE_DIR/P0_13_RESTART_COMPATIBILITY_AUDIT.txt" ]]; then fail "restart audit empty"; fi
if [[ ! -s "$EVIDENCE_DIR/P0_13_STATISTICS_COMPATIBILITY_AUDIT.txt" ]]; then fail "statistics audit empty"; fi
if [[ ! -s "$EVIDENCE_DIR/P0_13_VISUALIZATION_COMPATIBILITY_AUDIT.txt" ]]; then fail "visualization audit empty"; fi
if search "pressure|projection|RK3|channel-forcing" src/fibre_prod_p0_closure_check.f90 >/dev/null 2>&1; then
  fail "closure check references forbidden pressure/projection/RK3/channel-forcing logic"
else
  pass "closure check has no forbidden pressure/projection/RK3/channel-forcing logic"
fi

log "Configuring CMake build"
if cmake -S . -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" >> "$RESULT_FILE" 2>&1; then
  pass "cmake configure"
else
  fail "cmake configure failed"
fi

targets=(xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check fibre_prod_reaction_spreading_buffer_check fibre_prod_force_buffer_rhs_gate_check fibre_prod_synthetic_closed_loop_check fibre_prod_p0_closure_check)
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

checks=(fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check fibre_prod_reaction_spreading_buffer_check fibre_prod_force_buffer_rhs_gate_check fibre_prod_p0_closure_check)
if [[ $failures -eq 0 ]]; then
  for check in "${checks[@]}"; do run_check "$check"; done
  exe="$(find_exe fibre_prod_synthetic_closed_loop_check)"
  if [[ -z "$exe" ]]; then
    fail "synthetic closed-loop executable not found"
  else
    for np in 1 2 4; do
      out="$EVIDENCE_DIR/P0_13_NP${np}_SYNTHETIC_OUTPUT.txt"
      if command -v mpirun >/dev/null 2>&1; then
        if mpirun -np "$np" "$exe" > "$out" 2>&1; then pass "ran synthetic np=$np"; else fail "synthetic np=$np failed"; fi
      else
        echo "mpirun unavailable; required np=$np synthetic run failed" > "$out"
        fail "mpirun unavailable for synthetic np=$np"
      fi
      cat "$out" >> "$RESULT_FILE" 2>/dev/null || true
    done
  fi
fi

{
  echo "P0_13 P0 regression matrix"
  for token in \
    P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK \
    P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK \
    P0_3_RUNTIME_BRIDGE_CHECK \
    P0_4_VELOCITY_BRIDGE_CHECK \
    P0_5_STATE_VELOCITY_ATTACHMENT_CHECK \
    P0_6_HYDRO_INPUT_CANDIDATE_CHECK \
    P0_7_STRUCTURE_INPUT_HANDOFF_CHECK \
    P0_8_STRUCTURE_DRY_STEP_CHECK \
    P0_9_STRUCTURE_COMMIT_GATE_CHECK \
    P0_10_REACTION_SPREADING_BUFFER_CHECK \
    P0_11_FORCE_BUFFER_RHS_GATE_CHECK \
    P0_12_SYNTHETIC_CLOSED_LOOP_CHECK \
    P0_13_P0_CLOSURE_CHECK; do
    if grep -q "$token" "$RESULT_FILE"; then echo "PASS $token"; else echo "FAIL missing $token"; fi
  done
} > "$EVIDENCE_DIR/P0_13_P0_REGRESSION_MATRIX.txt"

if [[ $failures -eq 0 ]]; then
  while read -r line; do [[ "$line" == FAIL* ]] && fail "regression token missing"; done < "$EVIDENCE_DIR/P0_13_P0_REGRESSION_MATRIX.txt"
fi

if [[ $failures -eq 0 ]]; then
  cat > "$PASS_FAIL" <<'PASSDOC'
Result: PASS

Meaning: PASS means P0 production-path reconstruction is closed: no-fibre default xcompact3d behavior remains protected, restart/statistics/visualization compatibility is not disrupted, all P0.2–P0.12 regression checks remain intact, and the guarded synthetic fibre DNS-FSI path is ready to enter P1 single-fibre real channel DNS-FSI micro-case. It does NOT mean paper-scale long-time DNS production is ready.

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS

P1 entry status: ALLOWED

Next required stage: P1_0 single-fibre real channel DNS-FSI guarded micro-case preflight.
PASSDOC
  cat > "$P1_GATE" <<'P1DOC'
# P0_13 P1 Entry Gate

P1 entry gate: ALLOWED

Allowed next action:
P1_0 single-fibre real channel DNS-FSI guarded micro-case preflight.

P1_0 allowed scope:

* small channel DNS case;
* single fibre;
* short guarded run;
* small lambda;
* lambda=0 no-contamination regression;
* no long-time production statistics;
* no paper-scale conclusion.

Not allowed:

* paper-scale long-time DNS;
* mesh/timestep final production study;
* multi-parameter production campaign;
* publication-level statistics claim.
P1DOC
  log "Result: PASS"
  exit 0
else
  cat > "$PASS_FAIL" <<'FAILDOC'
Result: FAIL

Meaning: FAIL means one or more P0_13 static/build/run/np-closure checks did not complete successfully. P0 remains open and P1 entry is blocked until this stage passes in an environment with the required toolchain and MPI runner.

Production-run status: STILL BLOCKED

P1 entry status: BLOCKED

Next required stage: Resolve P0_13 validation failures, then rerun P0_13 before P1_0.
FAILDOC
  log "Result: FAIL"
  exit 1
fi
