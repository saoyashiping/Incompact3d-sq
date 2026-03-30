#!/usr/bin/env bash
set -euo pipefail

# 安全示例脚本：用于演示如何编译并运行单相槽道 Re_tau≈180 基线。
# 设计原则：
# 1) 不删除任何目录；
# 2) 不覆盖已有运行结果；
# 3) 每一步都可见、可审查。

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASELINE_DIR="${ROOT_DIR}/baseline/channel_re180_single_phase"
BUILD_DIR="${ROOT_DIR}/build"
RUN_DIR="${ROOT_DIR}/run_channel_re180_baseline"
NPROC="${NPROC:-8}"
EXECUTABLE="${BUILD_DIR}/bin/xcompact3d"
LOG_FILE="run.log"

printf '[INFO] Root directory: %s\n' "${ROOT_DIR}"
printf '[INFO] Baseline directory: %s\n' "${BASELINE_DIR}"
printf '[INFO] Build directory: %s\n' "${BUILD_DIR}"
printf '[INFO] Run directory: %s\n' "${RUN_DIR}"
printf '[INFO] MPI ranks (NPROC): %s\n' "${NPROC}"

if [[ ! -f "${BASELINE_DIR}/input.i3d" ]]; then
  printf '[ERROR] Baseline input not found: %s\n' "${BASELINE_DIR}/input.i3d" >&2
  exit 1
fi

if [[ -d "${RUN_DIR}" ]]; then
  printf '[ERROR] Run directory already exists: %s\n' "${RUN_DIR}" >&2
  printf '[ERROR] Refusing to continue to avoid overwriting existing results.\n' >&2
  exit 1
fi

printf '[STEP] Configure project with CMake...\n'
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release

printf '[STEP] Build executable...\n'
cmake --build "${BUILD_DIR}" -j

if [[ ! -x "${EXECUTABLE}" ]]; then
  printf '[ERROR] Executable not found or not executable: %s\n' "${EXECUTABLE}" >&2
  exit 1
fi

printf '[STEP] Create run directory (no overwrite)...\n'
mkdir "${RUN_DIR}"

printf '[STEP] Copy frozen baseline input into run directory...\n'
cp "${BASELINE_DIR}/input.i3d" "${RUN_DIR}/input.i3d"

printf '[STEP] Run solver with MPI and write log to %s...\n' "${LOG_FILE}"
(
  cd "${RUN_DIR}"
  mpirun -np "${NPROC}" "${EXECUTABLE}" > "${LOG_FILE}" 2>&1
)

printf '[DONE] Baseline run finished.\n'
printf '[DONE] Log: %s/%s\n' "${RUN_DIR}" "${LOG_FILE}"
printf '[DONE] You can now run check_outputs.sh for lightweight integrity checks.\n'
