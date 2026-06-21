#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/P0_10_evidence"
RESULT_FILE="$EVIDENCE_DIR/P0_10_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="$ROOT_DIR/production_recovery/P0_10_PASS_FAIL.md"
BUILD_DIR="$ROOT_DIR/build_p0_10"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
fail=0

mkdir -p "$EVIDENCE_DIR"
: > "$RESULT_FILE"

log() { printf '%s\n' "$*" | tee -a "$RESULT_FILE"; }
pass() { log "PASS: $*"; }
fail_check() { log "FAIL: $*"; fail=1; }
search_file() {
  local pattern="$1" file="$2"
  if command -v rg >/dev/null 2>&1; then
    rg -q "$pattern" "$file"
  else
    grep -Eq "$pattern" "$file"
  fi
}
not_search_file() {
  local pattern="$1" file="$2"
  if command -v rg >/dev/null 2>&1; then
    ! rg -q "$pattern" "$file"
  else
    ! grep -Eq "$pattern" "$file"
  fi
}
append_matches() {
  local pattern="$1"; shift
  if command -v rg >/dev/null 2>&1; then
    rg -n "$pattern" "$@" || true
  else
    grep -En "$pattern" "$@" || true
  fi
}
find_exe() {
  local target="$1"
  for path in "$BUILD_DIR/bin/$target" "$BUILD_DIR/src/$target" "$BUILD_DIR/$target"; do
    [[ -x "$path" ]] && { printf '%s\n' "$path"; return 0; }
  done
  find "$BUILD_DIR" -type f -name "$target" -perm -111 -print -quit 2>/dev/null
}

log "P0_10 validation started"
log "ROOT_DIR=$ROOT_DIR"
log "DECOMP2D_ROOT=$DECOMP2D_ROOT"

[[ -f "$ROOT_DIR/src/fibre_prod_reaction_force_candidate.f90" ]] && pass "reaction force candidate source exists" || fail_check "reaction force candidate source missing"
[[ -f "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" ]] && pass "reaction spreading source exists" || fail_check "reaction spreading source missing"
[[ -f "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer_check.f90" ]] && pass "reaction spreading check source exists" || fail_check "reaction spreading check source missing"
search_file "fibre_prod_reaction_force_candidate" "$ROOT_DIR/src/xcompact3d.f90" && pass "xcompact3d imports reaction candidate" || fail_check "xcompact3d reaction candidate import missing"
search_file "fibre_prod_reaction_spreading_buffer" "$ROOT_DIR/src/xcompact3d.f90" && pass "xcompact3d imports reaction spreading" || fail_check "xcompact3d reaction spreading import missing"
search_file "FIBRE_PROD_REACTION_SPREADING_ENABLE|fibre_prod_reaction_spreading_buffer_env_enabled" "$ROOT_DIR/src/xcompact3d.f90" && pass "P0_10 xcompact3d path is env gated" || fail_check "P0_10 env gate missing"
not_search_file "fibre_prod_rhs_adapter_apply" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" && pass "reaction spreading does not call RHS adapter" || fail_check "reaction spreading calls RHS adapter"
not_search_file "fibre_prod_main_hook_apply_force_buffer" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" && pass "reaction spreading does not call force-buffer main hook" || fail_check "reaction spreading calls force-buffer main hook"
not_search_file "rhs_x|rhs_y|rhs_z" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" && pass "reaction spreading module has no RHS arrays" || fail_check "reaction spreading module mentions RHS arrays"
search_file "force_buffer%fx|force_buffer%fy|force_buffer%fz|fibre_prod_spread_multiple_point_forces" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" && pass "reaction spreading writes/uses force buffer path" || fail_check "reaction spreading force buffer path missing"
not_search_file "contribution \*=|contribution = lambda_fsi|uniform contribution" "$ROOT_DIR/src/fibre_prod_rhs_adapter.f90" && pass "old uniform RHS contribution path absent" || fail_check "old uniform RHS contribution path found"
search_file "fibre_prod_reaction_spreading_buffer_check" "$ROOT_DIR/src/CMakeLists.txt" && pass "P0_10 CMake target exists" || fail_check "P0_10 CMake target missing"

append_matches "reaction_force|structure_input|net_reaction" "$ROOT_DIR/src/fibre_prod_reaction_force_candidate.f90" > "$EVIDENCE_DIR/P0_10_REACTION_FORCE_CANDIDATE_AUDIT.txt"
append_matches "force_buffer|has_spread|max_abs|sum_abs" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer_check.f90" > "$EVIDENCE_DIR/P0_10_EULERIAN_FORCE_BUFFER_AUDIT.txt"
append_matches "conservation|net_lagrangian|net_eulerian" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer_check.f90" > "$EVIDENCE_DIR/P0_10_SPREADING_CONSERVATION_AUDIT.txt"
append_matches "reaction_spreading|REACTION_SPREADING|REACTION_FORCE_CANDIDATE" "$ROOT_DIR/src/xcompact3d.f90" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" > "$EVIDENCE_DIR/P0_10_XCOMPACT3D_REACTION_SPREADING_CALL_PATH_AUDIT.txt"
append_matches "rhs|RHS|adapter|apply_force_buffer|force_buffer" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer.f90" "$ROOT_DIR/src/fibre_prod_reaction_spreading_buffer_check.f90" > "$EVIDENCE_DIR/P0_10_NO_RHS_FEEDBACK_AUDIT.txt"

if [[ "$fail" -eq 0 ]]; then
  log "Configuring CMake build"
  if cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT" -DDECOMP2D_ROOT="$DECOMP2D_ROOT" >> "$RESULT_FILE" 2>&1; then
    targets=(xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check
             fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check
             fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check fibre_prod_structure_commit_gate_check
             fibre_prod_reaction_spreading_buffer_check)
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
               fibre_prod_reaction_spreading_buffer_check)
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
  )
  for token in "${required_tokens[@]}"; do
    search_file "$token" "$RESULT_FILE" && pass "found token $token" || fail_check "missing token $token"
  done
fi

if [[ "$fail" -eq 0 ]]; then
  cat > "$PASS_FAIL_FILE" <<'PASS_EOF'
Result: PASS

Meaning: PASS means a structure-to-fluid reaction-force candidate can be constructed from the committed structure-side input and spread into an Eulerian production force buffer without modifying velocity fields or DNS RHS, while P0.2/P0.3/P0.4/P0.5/P0.6/P0.7/P0.8/P0.9 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0_11 force-buffer → RHS lambda-gated synthetic coupling.
PASS_EOF
  log "Result: PASS"
else
  cat > "$PASS_FAIL_FILE" <<'FAIL_EOF'
Result: FAIL

Meaning: FAIL means one or more P0_10 static/build/run checks did not complete successfully. The reaction spreading buffer path must not be used for production DNS-FSI until this stage passes in an environment with the required toolchain.

Production-run status: STILL BLOCKED

Next required stage: Resolve P0_10 validation failures, then rerun P0_10 before P0_11.
FAIL_EOF
  log "Result: FAIL"
fi

cat "$PASS_FAIL_FILE" >> "$RESULT_FILE"
[[ "$fail" -eq 0 ]] || exit 1
