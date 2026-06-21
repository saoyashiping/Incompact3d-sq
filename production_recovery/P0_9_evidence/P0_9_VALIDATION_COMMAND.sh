#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/P0_9_evidence"
RESULT_FILE="$EVIDENCE_DIR/P0_9_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="$ROOT_DIR/production_recovery/P0_9_PASS_FAIL.md"
BUILD_DIR="$ROOT_DIR/build_p0_9"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"

mkdir -p "$EVIDENCE_DIR"
: > "$RESULT_FILE"
fail=0

log() { printf '%s\n' "$*" | tee -a "$RESULT_FILE"; }
pass() { log "PASS: $*"; }
fail_check() { log "FAIL: $*"; fail=1; }

contains() {
  local file="$1" pattern="$2"
  [[ -f "$file" ]] && rg -q "$pattern" "$file"
}

not_contains() {
  local file="$1" pattern="$2"
  [[ -f "$file" ]] && ! rg -q "$pattern" "$file"
}

log "P0.9 validation started"
log "ROOT_DIR=$ROOT_DIR"
log "DECOMP2D_ROOT=$DECOMP2D_ROOT"

[[ -f "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" ]] && pass "commit gate source exists" || fail_check "commit gate source missing"
contains "$ROOT_DIR/src/fibre_prod_state.f90" "structure_u|fibre_prod_state_commit_structure_trial" && pass "state has structure velocity commit storage/helper" || fail_check "state commit storage/helper missing"
contains "$ROOT_DIR/src/xcompact3d.f90" "fibre_prod_structure_commit_gate" && pass "xcompact3d imports commit gate" || fail_check "xcompact3d commit gate import missing"
contains "$ROOT_DIR/src/xcompact3d.f90" "FIBRE_PROD_STRUCTURE_DRY_COMMIT_ENABLE|fibre_prod_structure_commit_gate_env_enabled" && pass "commit gate call path is env gated" || fail_check "commit gate env gate missing"
not_contains "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "fibre_prod_rhs_adapter_apply" && pass "commit gate does not call RHS adapter" || fail_check "commit gate calls RHS adapter"
not_contains "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "fibre_prod_main_hook_apply_force_buffer" && pass "commit gate does not call force-buffer main hook" || fail_check "commit gate calls force-buffer main hook"
not_contains "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "force_buffer%f[xyz]" && pass "commit gate does not write Eulerian force buffers" || fail_check "commit gate writes Eulerian force buffers"
not_contains "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "spread|spreading" && pass "commit gate has no IBM spreading path" || fail_check "commit gate contains spreading path"
not_contains "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "reaction" && pass "commit gate contains no reaction-force path" || fail_check "commit gate contains reaction-force path"
not_contains "$ROOT_DIR/src/fibre_prod_rhs_adapter.f90" "contribution *=|contribution = lambda_fsi|uniform contribution" && pass "old uniform RHS contribution path absent" || fail_check "old uniform RHS contribution path found"
contains "$ROOT_DIR/src/CMakeLists.txt" "fibre_prod_structure_commit_gate_check" && pass "P0.9 CMake target exists" || fail_check "P0.9 CMake target missing"

{
  echo "P0.9 structure commit gate source:"
  rg -n "commit_gate|reject_code|gate_enabled|commit_to_state" "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" || true
} > "$EVIDENCE_DIR/P0_9_STRUCTURE_COMMIT_GATE_AUDIT.txt"
{
  echo "P0.9 accept/reject evidence:"
  rg -n "reject|accepted|committed|nonfinite|bound|disabled" "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "$ROOT_DIR/src/fibre_prod_structure_commit_gate_check.f90" || true
} > "$EVIDENCE_DIR/P0_9_COMMIT_ACCEPT_REJECT_AUDIT.txt"
{
  echo "P0.9 xcompact3d call path evidence:"
  rg -n "structure_commit_gate|FIBRE_PROD_STRUCTURE_DRY_COMMIT_ENABLE" "$ROOT_DIR/src/xcompact3d.f90" "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" || true
} > "$EVIDENCE_DIR/P0_9_XCOMPACT3D_COMMIT_GATE_CALL_PATH_AUDIT.txt"
{
  echo "P0.9 no RHS feedback evidence:"
  rg -n "rhs_adapter|apply_force_buffer|force_buffer|RHS|rhs" "$ROOT_DIR/src/fibre_prod_structure_commit_gate.f90" "$ROOT_DIR/src/fibre_prod_structure_commit_gate_check.f90" || true
} > "$EVIDENCE_DIR/P0_9_NO_RHS_FEEDBACK_AUDIT.txt"

if [[ "$fail" -eq 0 ]]; then
  log "Configuring CMake build"
  if cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" >> "$RESULT_FILE" 2>&1; then
    for target in xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check \
                  fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check \
                  fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check \
                  fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check \
                  fibre_prod_structure_commit_gate_check; do
      log "Building target $target"
      cmake --build "$BUILD_DIR" --target "$target" >> "$RESULT_FILE" 2>&1 || { fail_check "build target $target"; break; }
    done
  else
    fail_check "cmake configure failed"
  fi
fi

if [[ "$fail" -eq 0 ]]; then
  for exe in fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check \
             fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check \
             fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check \
             fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check; do
    log "Running $exe"
    "$BUILD_DIR/src/$exe" >> "$RESULT_FILE" 2>&1 || { fail_check "run $exe"; break; }
  done
fi

if [[ "$fail" -eq 0 ]]; then
  required_tokens=(
    "P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS"
    "P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK PASS"
    "P0_3_RUNTIME_BRIDGE_CHECK PASS"
    "P0_4_VELOCITY_BRIDGE_CHECK PASS"
    "P0_5_STATE_VELOCITY_ATTACHMENT_CHECK PASS"
    "P0_6_HYDRO_INPUT_CANDIDATE_CHECK PASS"
    "P0_7_STRUCTURE_INPUT_HANDOFF_CHECK PASS"
    "P0_8_STRUCTURE_DRY_STEP_CHECK PASS"
    "P0_9_STRUCTURE_COMMIT_GATE_CHECK PASS"
  )
  for token in "${required_tokens[@]}"; do
    rg -q "$token" "$RESULT_FILE" && pass "found token $token" || fail_check "missing token $token"
  done
fi

if [[ "$fail" -eq 0 ]]; then
  cat > "$PASS_FAIL_FILE" <<'PASS_EOF'
Result: PASS

Meaning: PASS means a controlled structure dry-step commit gate can accept finite/bounded trial structure states and optionally commit them to the production fibre structure state without generating reaction force, modifying velocity fields, force buffers, or DNS RHS, while P0.2/P0.3/P0.4/P0.5/P0.6/P0.7/P0.8 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0.10 structure-to-fluid reaction-force candidate diagnostics, no Eulerian spreading and no RHS feedback.
PASS_EOF
  log "Result: PASS"
else
  cat > "$PASS_FAIL_FILE" <<'FAIL_EOF'
Result: FAIL

Meaning: FAIL means one or more P0.9 static/build/run checks did not complete successfully. The controlled commit gate must not be used for production DNS-FSI until this stage passes in an environment with the required toolchain.

Production-run status: STILL BLOCKED

Next required stage: Resolve P0.9 validation failures, then rerun P0.9 before P0.10.
FAIL_EOF
  log "Result: FAIL"
fi

cat "$PASS_FAIL_FILE" >> "$RESULT_FILE"
[[ "$fail" -eq 0 ]] || exit 1
