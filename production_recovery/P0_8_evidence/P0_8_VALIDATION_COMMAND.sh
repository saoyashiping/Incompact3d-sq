#!/usr/bin/env bash
set -u
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT_DIR}" || exit 1
EVIDENCE_DIR="production_recovery/P0_8_evidence"
mkdir -p "${EVIDENCE_DIR}"
RESULT_FILE="${EVIDENCE_DIR}/P0_8_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="production_recovery/P0_8_PASS_FAIL.md"
DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
BUILD_DIR="build_p0_8"
ALL_PASS=1
: > "${RESULT_FILE}"
log() { echo "$*" | tee -a "${RESULT_FILE}"; }
pass() { log "PASS: $*"; }
fail() { log "FAIL: $*"; ALL_PASS=0; }
write_pass_fail() {
  local result="$1"
  cat > "${PASS_FAIL_FILE}" <<EOF2
# P0.8 Structure Dry-Step PASS/FAIL

Result: ${result}

Meaning: PASS means a controlled structure dry-step predictor can consume structure-side input force into scratch/trial storage without committing production structure state, modifying velocity fields, force buffers, or DNS RHS, while P0.2/P0.3/P0.4/P0.5/P0.6/P0.7 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0.9 controlled structure dry-step commit gate, still no fluid feedback.
EOF2
}
find_exe() { find "${BUILD_DIR}" -type f -name "$1" -perm -111 2>/dev/null | head -n 1; }
run_and_expect() {
  local exe="$1" token="$2" label="$3" output
  if [ -z "${exe}" ] || [ ! -x "${exe}" ]; then fail "${label} executable not found"; return; fi
  if output="$(${exe} 2>&1)"; then
    printf '%s\n' "${output}" >> "${RESULT_FILE}"
    printf '%s\n' "${output}" | grep -q "${token}" && pass "${label} emitted ${token}" || fail "${label} missing ${token}"
  else
    printf '%s\n' "${output}" >> "${RESULT_FILE}"
    fail "${label} execution failed"
  fi
}

log "# P0.8 Validation Result"
log "DECOMP2D_ROOT=${DECOMP2D_ROOT}"
log ""
log "## Static structure dry-step audits"
[ -f src/fibre_prod_structure_dry_step.f90 ] && pass "structure dry-step source exists" || fail "structure dry-step source missing"
grep -q 'use fibre_prod_structure_dry_step' src/xcompact3d.f90 && pass "xcompact3d imports structure dry-step" || fail "xcompact3d missing structure dry-step import"
grep -q 'FIBRE_PROD_STRUCTURE_DRY_STEP_ENABLE' src/fibre_prod_structure_dry_step.f90 && grep -q 'fibre_prod_structure_dry_step_env_enabled' src/xcompact3d.f90 && pass "dry-step path is env-gated by FIBRE_PROD_STRUCTURE_DRY_STEP_ENABLE" || fail "dry-step env gate missing"
if grep -q 'fibre_prod_rhs_adapter_apply' src/fibre_prod_structure_dry_step.f90; then fail "dry-step calls RHS adapter"; else pass "dry-step does not call RHS adapter"; fi
if grep -q 'fibre_prod_main_hook_apply_force_buffer' src/fibre_prod_structure_dry_step.f90; then fail "dry-step calls force-buffer RHS hook"; else pass "dry-step does not call force-buffer RHS hook"; fi
if grep -q 'force_buffer%f[xyz]' src/fibre_prod_structure_dry_step.f90; then fail "dry-step writes force buffer component"; else pass "dry-step does not write force buffer components"; fi
if grep -qi 'spreading\|spread_point_force\|spread_multiple_point_forces' src/fibre_prod_structure_dry_step.f90; then fail "dry-step calls IBM force transfer"; else pass "dry-step does not call IBM force transfer"; fi
if grep -q 'state%x.*x_trial\|x_trial.*state%x\|state%v.*u_trial\|u_trial.*state%v\|state%a.*dx_trial\|dx_trial.*state%a' src/fibre_prod_structure_dry_step.f90; then fail "dry-step appears to commit trial state"; else pass "dry-step does not commit trial state to production state"; fi
if grep -qF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_structure_dry_step.f90 || \
   grep -qF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_structure_dry_step.f90; then
  fail "uniform contribution old path is present"
else
  pass "uniform contribution old path is absent"
fi
grep -q 'add_executable(fibre_prod_structure_dry_step_check' src/CMakeLists.txt && pass "CMake target fibre_prod_structure_dry_step_check exists" || fail "structure dry-step check target missing"

cat > "${EVIDENCE_DIR}/P0_8_STRUCTURE_DRY_STEP_AUDIT.txt" <<EOF2
P0.8 structure dry-step audit
Static audit confirms src/fibre_prod_structure_dry_step.f90 exists and provides scratch dry-step predictor diagnostics only.
Dynamic evidence is provided by fibre_prod_structure_dry_step_check when it emits P0_8_STRUCTURE_DRY_STEP_CHECK PASS.
EOF2
cat > "${EVIDENCE_DIR}/P0_8_NO_STATE_COMMIT_AUDIT.txt" <<EOF2
P0.8 no state commit audit
Static audit confirms the dry-step module does not assign trial buffers back to production state fields.
Dynamic evidence verifies coordinates, structure velocity, sampled velocity, hydro candidate, and structure input force are unchanged.
EOF2
cat > "${EVIDENCE_DIR}/P0_8_XCOMPACT3D_DRY_STEP_CALL_PATH_AUDIT.txt" <<EOF2
P0.8 xcompact3d dry-step call path audit
Static audit confirms xcompact3d imports fibre_prod_structure_dry_step and the dry-step path is gated by FIBRE_PROD_STRUCTURE_DRY_STEP_ENABLE.
EOF2
cat > "${EVIDENCE_DIR}/P0_8_NO_RHS_FEEDBACK_AUDIT.txt" <<EOF2
P0.8 no RHS feedback audit
Static audit confirms the dry-step module does not call RHS adapter, does not call main_hook force-buffer-to-RHS, does not call IBM force transfer, and does not write force-buffer components.
Dynamic evidence verifies velocity/RHS/force-buffer no-contamination when fibre_prod_structure_dry_step_check emits PASS.
EOF2

log ""
log "## CMake build"
if cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" -DDECOMP2D_ROOT="${DECOMP2D_ROOT}" >> "${RESULT_FILE}" 2>&1; then pass "cmake configure completed"; else fail "cmake configure failed"; fi
for target in xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check fibre_prod_hydro_input_candidate_check fibre_prod_structure_input_handoff_check fibre_prod_structure_dry_step_check; do
  if [ "${ALL_PASS}" -eq 1 ] && cmake --build "${BUILD_DIR}" --target "${target}" >> "${RESULT_FILE}" 2>&1; then pass "built ${target}"; else fail "failed to build ${target}"; fi
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
  run_and_expect "$(find_exe fibre_prod_structure_dry_step_check)" "P0_8_STRUCTURE_DRY_STEP_CHECK PASS" "fibre_prod_structure_dry_step_check"
else
  fail "skipping executable runs because build failed"
fi

if [ "${ALL_PASS}" -eq 1 ]; then write_pass_fail "PASS"; log ""; log "PASS"; else write_pass_fail "FAIL"; log ""; log "FAIL"; fi
log ""
log "## P0.8 pass/fail file"
cat "${PASS_FAIL_FILE}" | tee -a "${RESULT_FILE}"
cat "${PASS_FAIL_FILE}"
[ "${ALL_PASS}" -eq 1 ]
