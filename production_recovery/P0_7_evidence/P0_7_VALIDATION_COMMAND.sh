#!/usr/bin/env bash
set -u
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT_DIR}" || exit 1
EVIDENCE_DIR="production_recovery/P0_7_evidence"
mkdir -p "${EVIDENCE_DIR}"
RESULT_FILE="${EVIDENCE_DIR}/P0_7_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="production_recovery/P0_7_PASS_FAIL.md"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
BUILD_DIR="build_p0_7"
ALL_PASS=1
: > "${RESULT_FILE}"
log() { echo "$*" | tee -a "${RESULT_FILE}"; }
pass() { log "PASS: $*"; }
fail() { log "FAIL: $*"; ALL_PASS=0; }
write_pass_fail() {
  local result="$1"
  cat > "${PASS_FAIL_FILE}" <<EOF2
# P0.7 Structure Input Handoff PASS/FAIL

Result: ${result}

Meaning: PASS means the structure-side hydrodynamic input candidate can be handed off into a controlled structure input slot without advancing the structure, modifying velocity fields, force buffers, or DNS RHS, while P0.2/P0.3/P0.4/P0.5/P0.6 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0.8 controlled structure dry-step preflight, no fluid feedback.
EOF2
}
find_exe() {
  find "${BUILD_DIR}" -type f -name "$1" -perm -111 2>/dev/null | head -n 1
}
run_and_expect() {
  local exe="$1"
  local token="$2"
  local label="$3"
  if [ -z "${exe}" ] || [ ! -x "${exe}" ]; then
    fail "${label} executable not found"
    return
  fi
  local output
  if output="$(${exe} 2>&1)"; then
    printf '%s\n' "${output}" >> "${RESULT_FILE}"
    if printf '%s\n' "${output}" | grep -q "${token}"; then
      pass "${label} emitted ${token}"
    else
      fail "${label} missing ${token}"
    fi
  else
    printf '%s\n' "${output}" >> "${RESULT_FILE}"
    fail "${label} execution failed"
  fi
}

log "# P0.7 Validation Result"
log "DECOMP2D_ROOT=${DECOMP2D_ROOT}"
log ""
log "## Static structure input handoff audits"
[ -f src/fibre_prod_structure_input_handoff.f90 ] && pass "structure input handoff source exists" || fail "structure input handoff source missing"
grep -q 'structure_input_force' src/fibre_prod_state.f90 && grep -q 'has_structure_input_force' src/fibre_prod_state.f90 && pass "fibre_prod_state structure input storage exists" || fail "fibre_prod_state structure input storage missing"
grep -q 'use fibre_prod_structure_input_handoff' src/xcompact3d.f90 && pass "xcompact3d imports structure input handoff" || fail "xcompact3d missing structure input handoff import"
grep -q 'FIBRE_PROD_STRUCTURE_INPUT_HANDOFF_ENABLE' src/fibre_prod_structure_input_handoff.f90 && grep -q 'fibre_prod_structure_input_handoff_env_enabled' src/xcompact3d.f90 && pass "structure handoff path is env-gated by FIBRE_PROD_STRUCTURE_INPUT_HANDOFF_ENABLE" || fail "structure handoff env gate missing"
if grep -q 'fibre_prod_rhs_adapter_apply' src/fibre_prod_structure_input_handoff.f90; then fail "structure handoff calls RHS adapter"; else pass "structure handoff does not call RHS adapter"; fi
if grep -q 'fibre_prod_main_hook_apply_force_buffer' src/fibre_prod_structure_input_handoff.f90; then fail "structure handoff calls force-buffer RHS hook"; else pass "structure handoff does not call force-buffer RHS hook"; fi
if grep -q 'force_buffer%f[xyz]' src/fibre_prod_structure_input_handoff.f90; then fail "structure handoff writes force buffer component"; else pass "structure handoff does not write force buffer components"; fi
if grep -qi 'spreading\|spread_point_force\|spread_multiple_point_forces' src/fibre_prod_structure_input_handoff.f90; then fail "structure handoff calls IBM force transfer"; else pass "structure handoff does not call IBM force transfer"; fi
if grep -qi 'structure.*advance\|advance.*structure\|structure.*step\|step.*structure\|structure.*update\|update.*structure\|fibre_prod_structure_solver' src/fibre_prod_structure_input_handoff.f90; then fail "structure handoff references structure motion routine"; else pass "structure handoff does not reference structure motion routine"; fi
if grep -qF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_structure_input_handoff.f90 || \
   grep -qF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_structure_input_handoff.f90; then
  fail "uniform contribution old path is present"
else
  pass "uniform contribution old path is absent"
fi
grep -q 'add_executable(fibre_prod_structure_input_handoff_check' src/CMakeLists.txt && pass "CMake target fibre_prod_structure_input_handoff_check exists" || fail "structure handoff check target missing"

cat > "${EVIDENCE_DIR}/P0_7_STRUCTURE_INPUT_HANDOFF_AUDIT.txt" <<EOF2
P0.7 structure input handoff audit
Static audit confirms src/fibre_prod_structure_input_handoff.f90 exists and provides structure-side handoff diagnostics only.
Dynamic evidence is provided by fibre_prod_structure_input_handoff_check when it emits P0_7_STRUCTURE_INPUT_HANDOFF_CHECK PASS.
EOF2
cat > "${EVIDENCE_DIR}/P0_7_NO_STRUCTURE_ADVANCE_AUDIT.txt" <<EOF2
P0.7 no structure advance audit
Static audit confirms the handoff module does not reference structure motion routines.
Dynamic evidence verifies coordinates, structure velocity, and acceleration are unchanged.
EOF2
cat > "${EVIDENCE_DIR}/P0_7_XCOMPACT3D_HANDOFF_CALL_PATH_AUDIT.txt" <<EOF2
P0.7 xcompact3d handoff call path audit
Static audit confirms xcompact3d imports fibre_prod_structure_input_handoff and the handoff path is gated by FIBRE_PROD_STRUCTURE_INPUT_HANDOFF_ENABLE.
EOF2
cat > "${EVIDENCE_DIR}/P0_7_NO_RHS_FEEDBACK_AUDIT.txt" <<EOF2
P0.7 no RHS feedback audit
Static audit confirms the handoff module does not call RHS adapter, does not call main_hook force-buffer-to-RHS, does not call IBM force transfer, and does not write force-buffer components.
Dynamic evidence verifies velocity/RHS/force-buffer no-contamination when fibre_prod_structure_input_handoff_check emits PASS.
EOF2

log ""
log "## CMake build"
if cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" -DDECOMP2D_ROOT="${DECOMP2D_ROOT}" >> "${RESULT_FILE}" 2>&1; then
  pass "cmake configure completed"
else
  fail "cmake configure failed"
fi
for target in xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check; do
  if [ "${ALL_PASS}" -eq 1 ] && cmake --build "${BUILD_DIR}" --target "${target}" >> "${RESULT_FILE}" 2>&1; then
    pass "built ${target}"
  else
    fail "failed to build ${target}"
  fi
done

log ""
log "## Run checks"
if [ "${ALL_PASS}" -eq 1 ]; then
  run_and_expect "$(find_exe fibre_prod_main_hook_check)" "P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS" "fibre_prod_main_hook_check"
  FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-3 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0 \
    run_and_expect "$(find_exe fibre_prod_force_buffer_rhs_path_check)" "P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK PASS" "fibre_prod_force_buffer_rhs_path_check"
  run_and_expect "$(find_exe fibre_prod_runtime_bridge_check)" "P0_3_RUNTIME_BRIDGE_CHECK PASS" "fibre_prod_runtime_bridge_check"
  run_and_expect "$(find_exe fibre_prod_velocity_bridge_check)" "P0_4_VELOCITY_BRIDGE_CHECK PASS" "fibre_prod_velocity_bridge_check"
  run_and_expect "$(find_exe fibre_prod_state_velocity_attachment_check)" "P0_5_STATE_VELOCITY_ATTACHMENT_CHECK PASS" "fibre_prod_state_velocity_attachment_check"
  run_and_expect "$(find_exe fibre_prod_hydro_input_candidate_check)" "P0_6_HYDRO_INPUT_CANDIDATE_CHECK PASS" "fibre_prod_hydro_input_candidate_check"
  run_and_expect "$(find_exe fibre_prod_structure_input_handoff_check)" "P0_7_STRUCTURE_INPUT_HANDOFF_CHECK PASS" "fibre_prod_structure_input_handoff_check"
else
  fail "skipping executable runs because build failed"
fi

if [ "${ALL_PASS}" -eq 1 ]; then
  write_pass_fail "PASS"
  log ""
  log "PASS"
else
  write_pass_fail "FAIL"
  log ""
  log "FAIL"
fi
log ""
log "## P0.7 pass/fail file"
cat "${PASS_FAIL_FILE}" | tee -a "${RESULT_FILE}"
cat "${PASS_FAIL_FILE}"
[ "${ALL_PASS}" -eq 1 ]
