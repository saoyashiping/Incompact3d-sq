#!/usr/bin/env bash
set -euo pipefail

# 轻量完整性检查脚本：
# - 不解析复杂二进制
# - 不做数值后处理
# - 仅检查关键输出是否存在

RUN_DIR="${1:-$(pwd)}"

printf '[INFO] Checking run directory: %s\n' "${RUN_DIR}"

if [[ ! -d "${RUN_DIR}" ]]; then
  printf '[FAIL] Run directory does not exist: %s\n' "${RUN_DIR}" >&2
  exit 1
fi

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0

pass() {
  printf '[PASS] %s\n' "$1"
  PASS_COUNT=$((PASS_COUNT + 1))
}

warn() {
  printf '[WARN] %s\n' "$1"
  WARN_COUNT=$((WARN_COUNT + 1))
}

fail() {
  printf '[FAIL] %s\n' "$1"
  FAIL_COUNT=$((FAIL_COUNT + 1))
}

shopt -s nullglob

# 1) 日志检查
logs=("${RUN_DIR}"/*.log)
if (( ${#logs[@]} > 0 )); then
  pass "Log file(s) found: ${#logs[@]}"
else
  fail "No *.log file found in run directory"
fi

# 2) statistics 目录检查
if [[ -d "${RUN_DIR}/statistics" ]]; then
  pass "statistics/ directory exists"
else
  fail "statistics/ directory is missing"
fi

# 3) 关键统计文件检查（通配）
required_patterns=(
  "umean.dat*"
  "vmean.dat*"
  "wmean.dat*"
  "uumean.dat*"
  "vvmean.dat*"
  "wwmean.dat*"
  "uvmean.dat*"
)

for pattern in "${required_patterns[@]}"; do
  matches=("${RUN_DIR}/statistics"/${pattern})
  if (( ${#matches[@]} > 0 )); then
    pass "Found statistics/${pattern}"
  else
    fail "Missing statistics/${pattern}"
  fi
done

# 4) restart 文件检查（如有）
restart_candidates=(
  "${RUN_DIR}"/restart*
  "${RUN_DIR}"/*.restart*
  "${RUN_DIR}"/backup*
)
restart_found=0
for f in "${restart_candidates[@]}"; do
  if [[ -e "${f}" ]]; then
    restart_found=1
    break
  fi
done

if (( restart_found == 1 )); then
  pass "Restart-related file(s) found"
else
  warn "No restart-related file found (acceptable if checkpoint not reached or disabled)"
fi

printf '\n[SUMMARY] PASS=%d WARN=%d FAIL=%d\n' "${PASS_COUNT}" "${WARN_COUNT}" "${FAIL_COUNT}"

if (( FAIL_COUNT > 0 )); then
  exit 2
fi
