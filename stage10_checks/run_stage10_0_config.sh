#!/usr/bin/env bash
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR=${BUILD_DIR:-build_stage9}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}
fails=()
pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1"; fails+=("$1"); }

if [ ! -d "${BUILD_DIR}" ]; then
  if [ -n "${DECOMP2D_ROOT}" ]; then
    cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" || true
  else
    cmake -S . -B "${BUILD_DIR}" || true
  fi
fi

if cmake --build "${BUILD_DIR}" --target xcompact3d -j; then
  pass "build xcompact3d"
else
  fail "build xcompact3d"
fi

if cmake --build "${BUILD_DIR}" --target fibre_stage10_config_check -j; then
  pass "build fibre_stage10_config_check"
else
  fail "build fibre_stage10_config_check"
fi

if [ -f stage9_checks/STAGE9_CLOSED.md ]; then
  pass "stage9 closure marker present"
else
  echo "[INFO] stage9_checks/STAGE9_CLOSED.md not found; continuing"
fi

mkdir -p stage10_outputs
log="stage10_outputs/stage10_0_config.log"
X3D_STAGE10_PRODUCTION_HOOK=1 \
X3D_STAGE10_FORCE_NOOP=1 \
X3D_STAGE10_MAX_STEPS=3 \
"${BUILD_DIR}/bin/fibre_stage10_config_check" > "${log}" 2>&1 \
  < /dev/null \
  || fail "run fibre_stage10_config_check"

if ! grep -q "STAGE 10.0 CONFIG VERDICT: PASS" "${log}"; then
  fail "missing STAGE 10.0 CONFIG VERDICT: PASS"
else
  pass "config verdict pass line"
fi

dat="stage10_outputs/fibre_stage10_0_config.dat"
for k in \
  "stage10_0_requested_flag[[:space:]]\+1" \
  "stage10_0_noop_mode_status[[:space:]]\+1" \
  "stage10_0_no_fibre_state_status[[:space:]]\+1" \
  "stage10_0_no_force_status[[:space:]]\+1" \
  "stage10_0_no_rhs_injection_status[[:space:]]\+1" \
  "stage10_0_config_status[[:space:]]\+1"
do
  if grep -q "$k" "${dat}"; then
    pass "dat ${k}"
  else
    fail "dat missing ${k}"
  fi
done

if [ ${#fails[@]} -eq 0 ]; then
  echo "STAGE 10.0 FINAL VERDICT: PASS"
else
  echo "STAGE 10.0 FINAL VERDICT: FAIL"
  printf '  - %s\n' "${fails[@]}"
  exit 1
fi
