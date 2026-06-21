#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/P0_11_evidence"
RESULT_FILE="$EVIDENCE_DIR/P0_11_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="$ROOT_DIR/production_recovery/P0_11_PASS_FAIL.md"
BUILD_DIR="$ROOT_DIR/build_p0_11"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
fail=0
mkdir -p "$EVIDENCE_DIR"
: > "$RESULT_FILE"

log() { printf '%s\n' "$*" | tee -a "$RESULT_FILE"; }
pass() { log "PASS: $*"; }
fail_check() { log "FAIL: $*"; fail=1; }
search_file() { if command -v rg >/dev/null 2>&1; then rg -q "$1" "$2"; else grep -Eq "$1" "$2"; fi; }
not_search_file() { if command -v rg >/dev/null 2>&1; then ! rg -q "$1" "$2"; else ! grep -Eq "$1" "$2"; fi; }
append_matches() { if command -v rg >/dev/null 2>&1; then rg -n "$1" "${@:2}" || true; else grep -En "$1" "${@:2}" || true; fi; }
find_exe() {
  local target="$1"
  for path in "$BUILD_DIR/bin/$target" "$BUILD_DIR/src/$target" "$BUILD_DIR/$target"; do
    [[ -x "$path" ]] && { printf '%s\n' "$path"; return 0; }
  done
  find "$BUILD_DIR" -type f -name "$target" -perm -111 -print -quit 2>/dev/null
}

log "P0_11 validation started"
log "ROOT_DIR=$ROOT_DIR"
log "DECOMP2D_ROOT=$DECOMP2D_ROOT"

[[ -f "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" ]] && pass "RHS gate source exists" || fail_check "RHS gate source missing"
[[ -f "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate_check.f90" ]] && pass "RHS gate check source exists" || fail_check "RHS gate check source missing"
search_file "fibre_prod_force_buffer_rhs_gate" "$ROOT_DIR/src/xcompact3d.f90" && pass "xcompact3d imports RHS gate" || fail_check "xcompact3d RHS gate import missing"
search_file "FIBRE_PROD_FORCE_BUFFER_RHS_GATE_ENABLE|fibre_prod_force_buffer_rhs_gate_env_enabled" "$ROOT_DIR/src/xcompact3d.f90" && pass "P0_11 xcompact3d path is env gated" || fail_check "P0_11 env gate missing"
search_file "fibre_prod_main_hook_apply_force_buffer" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" && pass "RHS gate uses force-buffer main-hook path" || fail_check "RHS gate does not use force-buffer main-hook path"
not_search_file "contribution \*=|contribution = lambda_fsi|uniform contribution" "$ROOT_DIR/src/fibre_prod_rhs_adapter.f90" && pass "old uniform RHS contribution path absent" || fail_check "old uniform RHS contribution path found"
search_file "fibre_prod_force_buffer_rhs_gate_check" "$ROOT_DIR/src/CMakeLists.txt" && pass "P0_11 CMake target exists" || fail_check "P0_11 CMake target missing"
not_search_file "pressur|projection|channel-forcing" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" && pass "RHS gate module has no pressure/projection/channel-forcing path" || fail_check "RHS gate module touches pressure/projection/channel-forcing text"

append_matches "rhs_gate|force_buffer|main_hook|expected_scale|measured_scale_error" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" > "$EVIDENCE_DIR/P0_11_FORCE_BUFFER_RHS_GATE_AUDIT.txt"
append_matches "lambda_zero|lambda0|lambda_fsi = 0|lambda=0" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate_check.f90" > "$EVIDENCE_DIR/P0_11_LAMBDA0_NO_CONTAMINATION_AUDIT.txt"
append_matches "small-lambda|1.0e-3|bounded|max_abs_rhs_increment" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate_check.f90" > "$EVIDENCE_DIR/P0_11_SMALL_LAMBDA_RESPONSE_AUDIT.txt"
append_matches "linear|lambda_2|penalty|ratio|inc_beta" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate_check.f90" > "$EVIDENCE_DIR/P0_11_LINEAR_SCALING_AUDIT.txt"
append_matches "force_buffer_rhs_gate|FORCE_BUFFER_RHS_GATE" "$ROOT_DIR/src/xcompact3d.f90" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate.f90" > "$EVIDENCE_DIR/P0_11_XCOMPACT3D_RHS_GATE_CALL_PATH_AUDIT.txt"
append_matches "uniform|force_buffer|increment" "$ROOT_DIR/src/fibre_prod_force_buffer_rhs_gate_check.f90" "$ROOT_DIR/src/fibre_prod_rhs_adapter.f90" > "$EVIDENCE_DIR/P0_11_NO_UNIFORM_RHS_AUDIT.txt"

if [[ "$fail" -eq 0 ]]; then
  log "Configuring CMake build"
  if cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" >> "$RESULT_FILE" 2>&1; then
    targets=(xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check
             fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check
             fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check
             fibre_prod_reaction_spreading_buffer_check fibre_prod_force_buffer_rhs_gate_check)
    for target in "${targets[@]}"; do
      log "Building target $target"
      cmake --build "$BUILD_DIR" --target "$target" >> "$RESULT_FILE" 2>&1 || { fail_check "build target $target"; break; }
    done
  else
    fail_check "cmake configure failed"
  fi
fi

if [[ "$fail" -eq 0 ]]; then
  run_targets=(fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check
               fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check
               fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check
               fibre_prod_reaction_spreading_buffer_check fibre_prod_force_buffer_rhs_gate_check)
  for target in "${run_targets[@]}"; do
    exe="$(find_exe "$target")"
    [[ -n "$exe" ]] || { fail_check "executable not found $target"; break; }
    log "Running $target at $exe"
    if [[ "$target" == "fibre_prod_force_buffer_rhs_path_check" ]]; then
      FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-3 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0 "$exe" >> "$RESULT_FILE" 2>&1 || { fail_check "run $target"; break; }
    else
      "$exe" >> "$RESULT_FILE" 2>&1 || { fail_check "run $target"; break; }
    fi
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
    "P0_10_REACTION_SPREADING_BUFFER_CHECK PASS"
    "P0_11_FORCE_BUFFER_RHS_GATE_CHECK PASS"
  )
  for token in "${required_tokens[@]}"; do
    search_file "$token" "$RESULT_FILE" && pass "found token $token" || fail_check "missing token $token"
  done
fi

if [[ "$fail" -eq 0 ]]; then
  cat > "$PASS_FAIL_FILE" <<'PASS_EOF'
Result: PASS

Meaning: PASS means a nonzero Eulerian production force buffer generated by the reaction-force spreading path can be passed through a lambda-gated RHS coupling path, producing strict no-contamination at lambda=0 and bounded, linearly scaled RHS increments at small nonzero lambda, while P0.2/P0.3/P0.4/P0.5/P0.6/P0.7/P0.8/P0.9/P0.10 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0_12 end-to-end synthetic closed-loop + np=1/2/4 consistency.
PASS_EOF
  log "Result: PASS"
else
  cat > "$PASS_FAIL_FILE" <<'FAIL_EOF'
Result: FAIL

Meaning: FAIL means one or more P0_11 static/build/run checks did not complete successfully. The force-buffer RHS gate must not be used for production DNS-FSI until this stage passes in an environment with the required toolchain.

Production-run status: STILL BLOCKED

Next required stage: Resolve P0_11 validation failures, then rerun P0_11 before P0_12.
FAIL_EOF
  log "Result: FAIL"
fi

cat "$PASS_FAIL_FILE" >> "$RESULT_FILE"
[[ "$fail" -eq 0 ]] || exit 1
