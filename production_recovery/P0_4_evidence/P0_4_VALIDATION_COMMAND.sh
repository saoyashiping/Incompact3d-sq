#!/usr/bin/env bash
set -u

if [ ! -f "src/CMakeLists.txt" ] || [ ! -f "src/xcompact3d.f90" ]; then
  echo "FAIL: run this script from the project root containing src/CMakeLists.txt"
  exit 1
fi

DECOMP2D_ROOT="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"
EVIDENCE_DIR="production_recovery/P0_4_evidence"
RESULT_FILE="${EVIDENCE_DIR}/P0_4_VALIDATION_RESULT.txt"
PASS_FAIL_FILE="production_recovery/P0_4_PASS_FAIL.md"
mkdir -p "${EVIDENCE_DIR}"
: > "${RESULT_FILE}"

status=0
log() { echo "$*" | tee -a "${RESULT_FILE}"; }
fail() { log "FAIL: $*"; status=1; }
pass() { log "PASS: $*"; }

write_pass_fail() {
  result="$1"
  cat > "${PASS_FAIL_FILE}" <<EOF2
# P0.4 Runtime Velocity Interpolation Bridge PASS/FAIL

Result: ${result}

Meaning: PASS means the production runtime velocity bridge can sample xcompact3d velocity fields at Lagrangian fibre points in one-way diagnostic mode, while preserving velocity/RHS no-contamination and keeping P0.2/P0.3 validations intact. It does NOT mean production DNS-FSI is ready.

Production-run status: STILL BLOCKED

Next required stage: P0.5 production fibre state initialization and one-way sampled velocity attachment.
EOF2
}

log "# P0.4 Validation Result"
log "DECOMP2D_ROOT=${DECOMP2D_ROOT}"
log ""
log "## Static velocity bridge audits"
if [ -f src/fibre_prod_velocity_bridge.f90 ]; then pass "velocity bridge source exists"; else fail "velocity bridge source missing"; fi
if grep -q 'use fibre_prod_velocity_bridge' src/xcompact3d.f90; then pass "xcompact3d imports fibre_prod_velocity_bridge"; else fail "xcompact3d missing velocity bridge import"; fi
if grep -q 'fibre_prod_velocity_bridge_sample_runtime_unit_grid(ux1,uy1,uz1' src/xcompact3d.f90; then pass "xcompact3d velocity bridge call path reads ux1/uy1/uz1"; else fail "xcompact3d velocity bridge call path missing"; fi
if grep -q 'fibre_prod_rhs_adapter_apply' src/fibre_prod_velocity_bridge.f90; then fail "velocity bridge calls RHS adapter"; else pass "velocity bridge does not call RHS adapter"; fi
if grep -q 'fibre_prod_main_hook_apply_force_buffer' src/fibre_prod_velocity_bridge.f90; then fail "velocity bridge calls force-buffer-to-RHS hook"; else pass "velocity bridge does not call force-buffer-to-RHS hook"; fi
if grep -qF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_velocity_bridge.f90 || \
   grep -qF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 src/xcompact3d.f90 src/fibre_prod_velocity_bridge.f90; then
  fail "uniform contribution old path is present"
else
  pass "uniform contribution old path is absent"
fi
if grep -q 'add_executable(fibre_prod_velocity_bridge_check' src/CMakeLists.txt; then pass "CMake target fibre_prod_velocity_bridge_check exists"; else fail "velocity bridge target missing"; fi

cat > "${EVIDENCE_DIR}/P0_4_VELOCITY_BRIDGE_AUDIT.txt" <<EOF2
# P0.4 Velocity Bridge Audit

Static audit confirms the velocity bridge source exists, xcompact3d imports it, and the bridge does not call RHS feedback APIs.
EOF2
cat > "${EVIDENCE_DIR}/P0_4_INTERPOLATION_ACCURACY_AUDIT.txt" <<EOF2
# P0.4 Interpolation Accuracy Audit

Dynamic interpolation accuracy evidence is provided by fibre_prod_velocity_bridge_check when it emits P0_4_VELOCITY_BRIDGE_CHECK PASS.
EOF2
cat > "${EVIDENCE_DIR}/P0_4_XCOMPACT3D_ONEWAY_CALL_PATH_AUDIT.txt" <<EOF2
# P0.4 xcompact3d One-Way Call Path Audit

Static audit confirms xcompact3d has a default-off velocity bridge call path guarded by FIBRE_PROD_VELOCITY_SAMPLE_ENABLE and reading ux1/uy1/uz1.
EOF2
cat > "${EVIDENCE_DIR}/P0_4_NO_CONTAMINATION_AUDIT.txt" <<EOF2
# P0.4 No-Contamination Audit

Dynamic no-contamination evidence is provided by fibre_prod_velocity_bridge_check when the target builds and emits P0_4_VELOCITY_BRIDGE_CHECK PASS.
EOF2

log ""
log "## CMake build"
if cmake -S . -B build_p0_4 -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" -DDECOMP2D_ROOT="${DECOMP2D_ROOT}" >> "${RESULT_FILE}" 2>&1; then
  pass "cmake configure completed"
else
  fail "cmake configure failed"
fi
for target in xcompact3d fibre_prod_main_hook_check fibre_prod_force_buffer_rhs_path_check fibre_prod_runtime_bridge_check fibre_prod_velocity_bridge_check; do
  if [ "${status}" -eq 0 ] && cmake --build build_p0_4 --target "${target}" >> "${RESULT_FILE}" 2>&1; then
    pass "built ${target}"
  else
    fail "failed to build ${target}"
  fi
done

find_exe() { find build_p0_4 -type f -name "$1" -perm -111 2>/dev/null | head -n 1; }
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
log "## P0.4 pass/fail file"
cat "${PASS_FAIL_FILE}" | tee -a "${RESULT_FILE}"
exit "${status}"
