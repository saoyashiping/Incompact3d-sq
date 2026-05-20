#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR=${BUILD_DIR:-build_stage9}
MPIEXEC=${MPIEXEC:-mpirun}
DECOMP2D_ROOT=${DECOMP2D_ROOT:-}

fails=()

pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1"; fails+=("$1"); }

run_step(){
  local name="$1"; shift
  echo "[INFO] Running: $name"
  if "$@"; then
    pass "$name"
    return 0
  else
    fail "$name"
    return 1
  fi
}

echo "[INFO] Stage 9.2 minimal parallel gate started"
echo "[INFO] REPO_ROOT=${REPO_ROOT}"
echo "[INFO] BUILD_DIR=${BUILD_DIR}"
echo "[INFO] MPIEXEC=${MPIEXEC}"

configure_ok=1
build_ok=1
stage9_1_ok=1
mpi_ok=1

# A) Configure stage
if [ ! -d "${BUILD_DIR}" ]; then
  echo "[INFO] ${BUILD_DIR} does not exist; running CMake configure"
else
  echo "[INFO] ${BUILD_DIR} exists; running CMake configure refresh"
fi

if [ -n "${DECOMP2D_ROOT}" ]; then
  run_step "cmake configure (${BUILD_DIR}) with DECOMP2D_ROOT" \
    cmake -S . -B "${BUILD_DIR}" -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT}" || configure_ok=0
else
  run_step "cmake configure (${BUILD_DIR})" \
    cmake -S . -B "${BUILD_DIR}" || configure_ok=0
fi

# B) Build stage (only if configure succeeded)
if [ ${configure_ok} -eq 1 ]; then
  run_step "xcompact3d build" cmake --build "${BUILD_DIR}" --target xcompact3d -j || build_ok=0
  run_step "fibre_stage9_dependency_gate_check build" cmake --build "${BUILD_DIR}" --target fibre_stage9_dependency_gate_check -j || build_ok=0
  run_step "fibre_stage9_2_minimal_parallel_gate build" cmake --build "${BUILD_DIR}" --target fibre_stage9_2_minimal_parallel_gate -j || build_ok=0
else
  fail "configure stage failed; skipping all build steps"
  build_ok=0
fi

# Stage 9.1 gate (only if configure+build succeeded)
if [ ${configure_ok} -eq 1 ] && [ ${build_ok} -eq 1 ]; then
  run_step "Stage 9.1 interface consistency" bash stage9_checks/run_stage9_1_interface_consistency.sh || stage9_1_ok=0
else
  fail "build preconditions failed; skipping Stage 9.1 interface consistency"
  stage9_1_ok=0
fi

# C) MPI run stage (only if configure+build+stage9.1 all succeeded)
EXE="${BUILD_DIR}/bin/fibre_stage9_2_minimal_parallel_gate"
if [ ${configure_ok} -eq 1 ] && [ ${build_ok} -eq 1 ] && [ ${stage9_1_ok} -eq 1 ]; then
  if [ ! -f "${EXE}" ]; then
    fail "Stage 9.2 executable is missing: ${EXE}"
    mpi_ok=0
  else
    if [ ! -x "${EXE}" ]; then
      echo "[INFO] ${EXE} exists but is not executable; trying chmod +x"
      chmod +x "${EXE}" || true
    fi
    if [ ! -x "${EXE}" ]; then
      fail "Stage 9.2 executable is not executable"
      mpi_ok=0
    fi
  fi

  if [ ${mpi_ok} -eq 1 ]; then
    run_step "np=1 minimal parallel gate" "${MPIEXEC}" -np 1 "${EXE}" || mpi_ok=0
    run_step "np=2 minimal parallel gate" "${MPIEXEC}" -np 2 "${EXE}" || mpi_ok=0
    run_step "np=4 minimal parallel gate" "${MPIEXEC}" -np 4 "${EXE}" || mpi_ok=0
  else
    fail "MPI stage skipped due to non-executable/missing Stage 9.2 binary"
  fi
else
  fail "pre-MPI gates failed (configure/build/stage9.1); skipping MPI runs"
  mpi_ok=0
fi

if [ ${#fails[@]} -eq 0 ]; then
  echo "============================================================"
  echo "STAGE 9.2 FINAL VERDICT: PASS"
  echo "Reason: configure, builds, Stage 9.1 interface gate, and np=1/2/4 MPI minimal parallel gates passed."
  echo "============================================================"
  exit 0
else
  echo "============================================================"
  echo "STAGE 9.2 FINAL VERDICT: FAIL"
  echo "Failed checks:"
  for f in "${fails[@]}"; do
    echo "  - $f"
  done
  echo "============================================================"
  exit 1
fi
