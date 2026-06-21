#!/usr/bin/env bash
set -u

ROOT="$(git rev-parse --show-toplevel 2>/dev/null)"
if [ -z "${ROOT}" ]; then
  echo "FAIL: not inside a git repository"
  exit 1
fi
cd "${ROOT}" || exit 1

EVIDENCE_DIR="production_recovery/P0_0_evidence"
mkdir -p "${EVIDENCE_DIR}"

status=0
fail() { echo "FAIL: $*"; status=1; }
pass() { echo "PASS: $*"; }

for f in \
  production_recovery/R10_PASS_FAIL.md \
  production_recovery/R11_PASS_FAIL.md \
  production_recovery/R12_PASS_FAIL.md \
  production_recovery/R12_evidence/R12_FINAL_VALIDATION_MATRIX.md \
  production_recovery/R12_evidence/R12_FINAL_CLOSURE_REPORT.md; do
  echo "--- grep ${f} ---"
  if grep -nE 'BLOCKED|FAIL|PASS' "${f}"; then
    pass "grep status evidence found in ${f}"
  else
    fail "missing status evidence in ${f}"
  fi
done

main_audit="$(mktemp)"
trap 'rm -f "${main_audit}"' EXIT
awk '
  /add_executable\(xcompact3d/ {in_target=1}
  in_target && /fibre_prod/ {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); print}
  in_target && /xcompact3d\.f90\)/ {in_target=0}
' src/CMakeLists.txt > "${main_audit}"

echo "--- xcompact3d direct fibre_prod sources ---"
cat "${main_audit}"
for f in fibre_prod_runtime_config.f90 fibre_prod_main_diagnostics.f90 fibre_prod_rhs_adapter.f90 fibre_prod_main_hook.f90; do
  grep -qx "${f}" "${main_audit}" && pass "xcompact3d includes ${f}" || fail "xcompact3d missing ${f}"
done
for f in \
  fibre_prod_state.f90 fibre_prod_grid_adapter.f90 fibre_prod_ibm_delta.f90 \
  fibre_prod_ibm_interpolation.f90 fibre_prod_ibm_force_buffer.f90 \
  fibre_prod_ibm_spreading.f90 fibre_prod_boundary_conditions.f90 \
  fibre_prod_bending_solver.f90 fibre_prod_tension_solver.f90 \
  fibre_prod_structure_solver.f90 fibre_prod_fsi_config.f90 \
  fibre_prod_fsi_coupling.f90 fibre_prod_fsi_diagnostics.f90 \
  fibre_prod_wall_geometry.f90 fibre_prod_wall_contact.f90 \
  fibre_prod_fibre_collision.f90; do
  if grep -qx "${f}" "${main_audit}"; then
    fail "unexpected production FSI module in xcompact3d target: ${f}"
  else
    pass "confirmed not in xcompact3d target: ${f}"
  fi
done

echo "--- RHS adapter audit ---"
grep -nF 'contribution = config%lambda_fsi * config%penalty_beta * config%dt' src/fibre_prod_rhs_adapter.f90 \
  && pass "small-lambda contribution formula found" \
  || fail "small-lambda contribution formula missing"
grep -nF 'rhs_x(i, j, k) = rhs_x(i, j, k) + contribution' src/fibre_prod_rhs_adapter.f90 \
  && pass "uniform rhs_x injection found" \
  || fail "uniform rhs_x injection missing"

if git diff --quiet -- src; then
  pass "src has no uncommitted P0.0 modifications"
else
  fail "src has uncommitted modifications"
  git diff -- src
fi

for f in \
  production_recovery/P0_0_PLAN.md \
  production_recovery/P0_0_SOURCE_AUDIT.md \
  production_recovery/P0_0_PASS_FAIL.md \
  production_recovery/P0_0_evidence/P0_0_R10_R12_STATUS_AUDIT.txt \
  production_recovery/P0_0_evidence/P0_0_MAIN_TARGET_MODULE_AUDIT.txt \
  production_recovery/P0_0_evidence/P0_0_RHS_ADAPTER_AUDIT.txt \
  production_recovery/P0_0_evidence/P0_0_STANDALONE_VS_PRODUCTION_AUDIT.md; do
  [ -f "${f}" ] && pass "P0.0 file exists: ${f}" || fail "missing P0.0 file: ${f}"
done

echo "--- production_recovery/P0_0_PASS_FAIL.md ---"
cat production_recovery/P0_0_PASS_FAIL.md

if [ "${status}" -eq 0 ]; then
  echo "PASS"
else
  echo "FAIL"
fi
exit "${status}"
