#!/usr/bin/env bash
set -u

if [ ! -f "src/CMakeLists.txt" ] || [ ! -f "src/fibre_prod_main_hook.f90" ]; then
  echo "FAIL: run this script from the project root containing src/CMakeLists.txt"
  exit 1
fi

EVIDENCE_DIR="production_recovery/P0_2_evidence"
RESULT_FILE="${EVIDENCE_DIR}/P0_2_VALIDATION_RESULT.txt"
mkdir -p "${EVIDENCE_DIR}"
: > "${RESULT_FILE}"

status=0
log() { echo "$*" | tee -a "${RESULT_FILE}"; }
fail() { log "FAIL: $*"; status=1; }
pass() { log "PASS: $*"; }

log "# P0.2 Validation Result"
log ""
log "## Static source audits"
if grep -q 'fibre_prod_main_hook_apply_force_buffer' src/fibre_prod_main_hook.f90; then
  pass "main hook exposes fibre_prod_main_hook_apply_force_buffer"
else
  fail "main hook buffer API missing"
fi
if grep -qF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90; then
  fail "old uniform contribution formula is present"
else
  pass "old uniform contribution formula is absent"
fi
if grep -qF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 || \
   grep -qF 'rhs_x(i,j,k) = rhs_x(i,j,k) + contribution' src/fibre_prod_rhs_adapter.f90; then
  fail "old uniform rhs_x increment is present"
else
  pass "old uniform rhs_x increment is absent"
fi
if awk '/add_executable\(fibre_prod_force_buffer_rhs_path_check/ {found=1} END {exit found ? 0 : 1}' src/CMakeLists.txt; then
  pass "CMake target fibre_prod_force_buffer_rhs_path_check exists"
else
  fail "CMake target fibre_prod_force_buffer_rhs_path_check missing"
fi

log ""
log "## Build and run checks"
find_exe() {
  find "$1" -type f -name "$2" -perm -111 2>/dev/null | head -n 1
}
run_and_expect() {
  exe="$1"
  token="$2"
  shift 2
  if "$@" "${exe}" 2>&1 | tee -a "${RESULT_FILE}" | grep -q "${token}"; then
    pass "${exe} emitted ${token}"
  else
    fail "${exe} did not emit ${token}"
  fi
}

built=0
if command -v cmake >/dev/null 2>&1; then
  log "Trying CMake build: cmake -S . -B build_p0_2"
  if cmake -S . -B build_p0_2 >> "${RESULT_FILE}" 2>&1 && \
     cmake --build build_p0_2 --target fibre_prod_main_hook_check >> "${RESULT_FILE}" 2>&1 && \
     cmake --build build_p0_2 --target fibre_prod_force_buffer_rhs_path_check >> "${RESULT_FILE}" 2>&1; then
    main_hook_exe="$(find_exe build_p0_2 fibre_prod_main_hook_check)"
    rhs_path_exe="$(find_exe build_p0_2 fibre_prod_force_buffer_rhs_path_check)"
    if [ -n "${main_hook_exe}" ] && [ -n "${rhs_path_exe}" ]; then
      built=1
      run_and_expect "./${main_hook_exe}" "P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS" env
      run_and_expect "./${rhs_path_exe}" "P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK PASS" \
        env FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-3 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0
    else
      fail "CMake build completed but expected executables were not found"
    fi
  else
    log "CMake path unavailable or failed; trying fallback Fortran compiler"
  fi
else
  log "cmake not found; trying fallback Fortran compiler"
fi

if [ "${built}" -eq 0 ]; then
  fc=""
  for candidate in gfortran flang ifort ifx nvfortran; do
    if command -v "${candidate}" >/dev/null 2>&1; then fc="${candidate}"; break; fi
  done
  if [ -z "${fc}" ]; then
    fail "no supported Fortran compiler found for fallback build"
  else
    fallback_dir="build_p0_2_fallback"
    mkdir -p "${fallback_dir}"
    log "Trying fallback compiler: ${fc}"
    if "${fc}" -J "${fallback_dir}" -I "${fallback_dir}" src/fibre_prod_grid_adapter.f90 \
               src/fibre_prod_ibm_delta.f90 \
               src/fibre_prod_ibm_force_buffer.f90 \
               src/fibre_prod_ibm_spreading.f90 \
               src/fibre_prod_runtime_config.f90 \
               src/fibre_prod_main_diagnostics.f90 \
               src/fibre_prod_rhs_adapter.f90 \
               src/fibre_prod_main_hook.f90 \
               src/fibre_prod_main_hook_check.f90 \
               -o "${fallback_dir}/fibre_prod_main_hook_check" >> "${RESULT_FILE}" 2>&1 && \
       "${fc}" -J "${fallback_dir}" -I "${fallback_dir}" src/fibre_prod_grid_adapter.f90 \
               src/fibre_prod_ibm_delta.f90 \
               src/fibre_prod_ibm_force_buffer.f90 \
               src/fibre_prod_ibm_spreading.f90 \
               src/fibre_prod_runtime_config.f90 \
               src/fibre_prod_main_diagnostics.f90 \
               src/fibre_prod_rhs_adapter.f90 \
               src/fibre_prod_main_hook.f90 \
               src/fibre_prod_force_buffer_rhs_path_check.f90 \
               -o "${fallback_dir}/fibre_prod_force_buffer_rhs_path_check" >> "${RESULT_FILE}" 2>&1; then
      run_and_expect "./${fallback_dir}/fibre_prod_main_hook_check" "P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS" env
      run_and_expect "./${fallback_dir}/fibre_prod_force_buffer_rhs_path_check" "P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK PASS" \
        env FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-3 FIBRE_PROD_PENALTY_BETA=2.0 FIBRE_PROD_DIAGNOSTICS=0
    else
      fail "fallback compiler build failed"
    fi
  fi
fi

log ""
log "## P0.2 pass/fail file"
cat production_recovery/P0_2_PASS_FAIL.md | tee -a "${RESULT_FILE}"

if [ "${status}" -eq 0 ]; then
  log ""
  log "PASS"
else
  log ""
  log "FAIL"
fi
exit "${status}"
