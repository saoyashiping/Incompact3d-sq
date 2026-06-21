#!/usr/bin/env bash
set -u

if [ ! -f "src/CMakeLists.txt" ] || [ ! -f "src/xcompact3d.f90" ]; then
  echo "FAIL: run this script from the project root containing src/CMakeLists.txt"
  exit 1
fi

DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
EVIDENCE_DIR="production_recovery/P0_5_evidence"
RESULT_FILE="${EVIDENCE_DIR}/P0_5_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="production_recovery/P0_5_PASS_FAIL.md"
mkdir -p "${EVIDENCE_DIR}"
: > "${RESULT_FILE}"

status=0
log() { echo "$*" | tee -a "${RESULT_FILE}"; }
fail() { log "FAIL: $*"; status=1; }
pass() { log "PASS: $*"; }

write_pass_fail() {
  result="$1"
  cat > "${PASS_FAIL_FILE}" <<EOF2
# P0.5 Production Fibre State Velocity Attachment PASS/FAIL

Result: ${result}

Meaning: PASS means a production fibre state can be initialized, one-way sampled fluid velocity can be attached to the fibre state, and velocity/RHS no-contamination is preserved while P0.2/P0.3/P0.4 validations remain intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0.6 sampled-velocity-to-fluid-force-input candidate diagnostics, no RHS feedback.
EOF2
}

log "# P0.5 Validation Result"
log "DECOMP2D_ROOT=${DECOMP2D_ROOT}"
log ""
log "## Static state velocity attachment audits"
[ -f src/fibre_prod_state_velocity_attachment.f90 ] && pass "state velocity attachment source exists" || fail "state velocity attachment source missing"
grep -q 'sampled_u' src/fibre_prod_state.f90 && grep -q 'has_sampled_velocity' src/fibre_prod_state.f90 && pass "fibre_prod_state sampled velocity storage exists" || fail "fibre_prod_state sampled velocity storage missing"
grep -q 'use fibre_prod_state_velocity_attachment' src/xcompact3d.f90 && pass "xcompact3d imports state velocity attachment" || fail "xcompact3d missing state velocity attachment import"
grep -q 'FIBRE_PROD_STATE_VELOCITY_ATTACH_ENABLE' src/fibre_prod_state_velocity_attachment.f90 && pass "state attachment path is env-gated by FIBRE_PROD_STATE_VELOCITY_ATTACH_ENABLE" || fail "state attachment env gate missing"
if grep -q 'fibre_prod_rhs_adapter_apply' src/fibre_prod_state_velocity_attachment.f90; then fail "state attachment calls RHS adapter"; else pass "state attachment does not call RHS adapter"; fi
if grep -q 'fibre_prod_main_hook_apply_force_buffer' src/fibre_prod_state_velocity_attachment.f90; then fail "state attachment calls force-buffer RHS hook"; else pass "state attachment does not call force-buffer RHS hook"; fi
if grep -qi 'structure.*advance\|advance.*structure\|fibre_prod_structure' src/fibre_prod_state_velocity_attachment.f90; then fail "state attachment references structure advance"; else pass "state attachment does not reference structure advance"; fi
if grep -qF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_state_velocity_attachment.f90 || \
   grep -qF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_state_velocity_attachment.f90; then
  fail "uniform contribution old path is present"
else
  pass "uniform contribution old path is absent"
fi
grep -q 'add_executable(fibre_prod_state_velocity_attachment_check' src/CMakeLists.txt && pass "CMake target fibre_prod_state_velocity_attachment_check exists" || fail "state attachment check target missing"

cat > "${EVIDENCE_DIR}/P0_5_STATE_INIT_AUDIT.txt" <<EOF2
# P0.5 State Init Audit

Static audit confirms sampled velocity storage exists in fibre_prod_state. Dynamic state initialization evidence is provided by fibre_prod_state_velocity_attachment_check when it emits P0_5_STATE_VELOCITY_ATTACHMENT_CHECK PASS.
EOF2
cat > "${EVIDENCE_DIR}/P0_5_VELOCITY_ATTACHMENT_AUDIT.txt" <<EOF2
# P0.5 Velocity Attachment Audit

Static audit confirms the attachment module exists and does not call RHS feedback APIs. Dynamic attachment evidence is provided by the P0.5 check executable.
EOF2
cat > "${EVIDENCE_DIR}/P0_5_XCOMPACT3D_STATE_ATTACHMENT_CALL_PATH_AUDIT.txt" <<EOF2
# P0.5 xcompact3d State Attachment Call Path Audit

Static audit confirms xcompact3d imports the state velocity attachment module and the attachment path is guarded by FIBRE_PROD_STATE_VELOCITY_ATTACH_ENABLE.
EOF2
cat > "${EVIDENCE_DIR}/P0_5_NO_CONTAMINATION_AUDIT.txt" <<EOF2
# P0.5 No-Contamination Audit

Dynamic no-contamination evidence is provided by fibre_prod_state_velocity_attachment_check when it emits P0_5_STATE_VELOCITY_ATTACHMENT_CHECK PASS.
EOF2

log ""
log "## CMake build"
if cmake -S . -B build_p0_5 -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" -DDECOMP2D_ROOT="${DECOMP2D_ROOT}" >> "${RESULT_FILE}" 2>&1; then
  pass "cmake configure completed"
else
  fail "cmake configure failed"
fi
for target in xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check fibre_prod_state_velocity_attachment_check; do
  if [ "${status}" -eq 0 ] && cmake --build build_p0_5 --target "${target}" >> "${RESULT_FILE}" 2>&1; then
    pass "built ${target}"
  else
    fail "failed to build ${target}"
  fi
done

find_exe() { find build_p0_5 -type f -name "$1" -perm -111 2>/dev/null | head -n 1; }
run_and_expect() {
  exe="$1"
  token="$2"
  shift 2
  if [ -z "${exe}" ]; then fail "missing executable for ${token}"; return; fi
  if "$@" "${exe}" 2>&1 | tee -a "${RESULT_FILE}" | grep -q "${token}"; then
    pass "${exe} emitted ${token}"
  else
    fail "${exe} did not emit ${token}"
  fi
}

log ""
log "## Run checks"
if [ "${status}" -eq 0 ]; then
  run_and_expect "$(find_exe fibre_prod_main_hook_check)" "P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS" env
  run_and_expect "$(find_exe fibre_prod_force_buffer_rhs_path_check)" "P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK PASS" \
    env FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-3 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0
  run_and_expect "$(find_exe fibre_prod_runtime_bridge_check)" "P0_3_RUNTIME_BRIDGE_CHECK PASS" \
    env FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=0.0 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0
  run_and_expect "$(find_exe fibre_prod_velocity_bridge_check)" "P0_4_VELOCITY_BRIDGE_CHECK PASS" \
    env FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=0.0 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0
  run_and_expect "$(find_exe fibre_prod_state_velocity_attachment_check)" "P0_5_STATE_VELOCITY_ATTACHMENT_CHECK PASS" \
    env FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=0.0 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0
else
  fail "skipping executable runs because build failed"
fi

if [ "${status}" -eq 0 ]; then
  write_pass_fail "PASS"
  log ""
  log "PASS"
else
  write_pass_fail "FAIL"
  log ""
  log "FAIL"
fi

log ""
log "## P0.5 pass/fail file"
cat "${PASS_FAIL_FILE}" | tee -a "${RESULT_FILE}"
exit "${status}"
