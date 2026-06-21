#!/usr/bin/env bash
set -u

if [ ! -f "src/CMakeLists.txt" ] || [ ! -f "src/fibre_prod_rhs_adapter.f90" ]; then
  echo "FAIL: run this script from the project root containing src/CMakeLists.txt"
  exit 1
fi

EVIDENCE_DIR="production_recovery/P0_1_evidence"
RESULT_FILE="${EVIDENCE_DIR}/P0_1_VALIDATION_RESULT.txt"
mkdir -p "${EVIDENCE_DIR}"
: > "${RESULT_FILE}"

status=0
log() { echo "$*" | tee -a "${RESULT_FILE}"; }
fail() { log "FAIL: $*"; status=1; }
pass() { log "PASS: $*"; }

log "# P0.1 Validation Result"
log ""
log "## xcompact3d target module audit"
main_audit="$(mktemp)"
trap 'rm -f "${main_audit}"' EXIT
awk '
  /add_executable\(xcompact3d/ {in_target=1}
  in_target {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); print}
  in_target && /xcompact3d\.f90\)/ {in_target=0}
' src/CMakeLists.txt > "${main_audit}"

required_modules="
fibre_prod_state.f90
fibre_prod_grid_adapter.f90
fibre_prod_ibm_delta.f90
fibre_prod_ibm_interpolation.f90
fibre_prod_ibm_force_buffer.f90
fibre_prod_ibm_spreading.f90
fibre_prod_boundary_conditions.f90
fibre_prod_bending_solver.f90
fibre_prod_tension_solver.f90
fibre_prod_structure_solver.f90
fibre_prod_fsi_config.f90
fibre_prod_fluid_surrogate.f90
fibre_prod_fsi_coupling.f90
fibre_prod_fsi_diagnostics.f90
fibre_prod_wall_geometry.f90
fibre_prod_wall_contact.f90
fibre_prod_fibre_collision.f90
fibre_prod_runtime_config.f90
fibre_prod_main_diagnostics.f90
fibre_prod_rhs_adapter.f90
fibre_prod_main_hook.f90
"
for module in ${required_modules}; do
  if grep -qx "${module}" "${main_audit}"; then
    pass "xcompact3d target includes ${module}"
  else
    fail "xcompact3d target missing ${module}"
  fi
done

log ""
log "## RHS adapter safety audit"
if grep -qF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90; then
  fail "old uniform contribution formula is still present"
else
  pass "old uniform contribution formula is absent"
fi
if grep -qF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 || \
   grep -qF 'rhs_x(i,j,k) = rhs_x(i,j,k) + contribution' src/fibre_prod_rhs_adapter.f90; then
  fail "old uniform rhs_x increment is still present"
else
  pass "old uniform rhs_x increment is absent"
fi
if grep -qF 'fibre_prod_rhs_status_missing_force_buffer = 13' src/fibre_prod_rhs_adapter.f90; then
  pass "missing force-density buffer status 13 is defined"
else
  fail "missing force-density buffer status 13 is not defined"
fi
if grep -q 'optional :: force_x' src/fibre_prod_rhs_adapter.f90 && \
   grep -q 'optional :: force_y' src/fibre_prod_rhs_adapter.f90 && \
   grep -q 'optional :: force_z' src/fibre_prod_rhs_adapter.f90; then
  pass "optional force-density buffer interface exists in rhs adapter"
else
  fail "optional force-density buffer interface missing in rhs adapter"
fi

log ""
log "## Build and hook-check run"
run_hook_check() {
  exe="$1"
  if "${exe}" 2>&1 | tee -a "${RESULT_FILE}" | grep -q 'P0_1_FIBRE_PROD_MAIN_HOOK_CHECK PASS'; then
    pass "fibre_prod_main_hook_check emitted P0_1 PASS token"
    return 0
  fi
  fail "fibre_prod_main_hook_check did not emit P0_1 PASS token"
  return 1
}

built=0
if command -v cmake >/dev/null 2>&1; then
  log "Trying CMake build: cmake -S . -B build_p0_1"
  if cmake -S . -B build_p0_1 >> "${RESULT_FILE}" 2>&1 && \
     cmake --build build_p0_1 --target fibre_prod_main_hook_check >> "${RESULT_FILE}" 2>&1; then
    exe=""
    for candidate in build_p0_1/src/fibre_prod_main_hook_check build_p0_1/fibre_prod_main_hook_check; do
      if [ -x "${candidate}" ]; then exe="${candidate}"; break; fi
    done
    if [ -z "${exe}" ]; then
      exe="$(find build_p0_1 -type f -name fibre_prod_main_hook_check -perm -111 2>/dev/null | head -n 1)"
    fi
    if [ -n "${exe}" ]; then
      built=1
      run_hook_check "./${exe}"
    else
      fail "CMake built target but executable was not found"
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
    fallback_dir="build_p0_1_fallback"
    mkdir -p "${fallback_dir}"
    log "Trying fallback compiler: ${fc}"
    if "${fc}" src/fibre_prod_runtime_config.f90 \
               src/fibre_prod_main_diagnostics.f90 \
               src/fibre_prod_rhs_adapter.f90 \
               src/fibre_prod_main_hook.f90 \
               src/fibre_prod_main_hook_check.f90 \
               -o "${fallback_dir}/fibre_prod_main_hook_check" >> "${RESULT_FILE}" 2>&1; then
      run_hook_check "./${fallback_dir}/fibre_prod_main_hook_check"
    else
      fail "fallback compiler build failed"
    fi
  fi
fi

log ""
log "## P0.1 pass/fail file"
cat production_recovery/P0_1_PASS_FAIL.md | tee -a "${RESULT_FILE}"

if [ "${status}" -eq 0 ]; then
  log ""
  log "PASS"
else
  log ""
  log "FAIL"
fi
exit "${status}"
