#!/usr/bin/env bash
set -u

failures=()

pass() { echo "[PASS] $1"; }
fail() { echo "[FAIL] $1"; failures+=("$1"); }

search_src() {
  local pattern="$1"
  grep -RInE "$pattern" src --include='*.f90' --include='*.F90' 2>/dev/null
}

check_absent_file() {
  local file="$1" desc="$2"
  if [ -e "$file" ]; then fail "$desc"; else pass "$desc"; fi
}

check_absent_pattern_file() {
  local pattern="$1" file="$2" desc="$3"
  if grep -nE "$pattern" "$file" >/dev/null 2>&1; then fail "$desc"; else pass "$desc"; fi
}

check_absent_pattern_src() {
  local pattern="$1" desc="$2"
  if search_src "$pattern" >/dev/null; then fail "$desc"; else pass "$desc"; fi
}

check_present_pattern_src() {
  local pattern="$1" desc="$2"
  if search_src "$pattern" >/dev/null; then pass "$desc"; else fail "$desc"; fi
}

check_absent_file "src/xcompact3d_decomp_io_compat.f90" "compat source file removed"
check_absent_file "src/fibre_stage_decomp2d_constants_stub.f90" "decomp constants stub removed"
check_absent_pattern_file "xcompact3d_decomp_io_compat\.f90" "src/CMakeLists.txt" "CMakeLists has no compat source entries"
check_absent_pattern_src "^[[:space:]]*(use|USE)[[:space:]]+xcompact3d_decomp_io_compat\b" "no active compat module import in src"
check_absent_pattern_src "^[[:space:]]*(use|USE)[[:space:]]+m_halo\b" "no active m_halo import in src"
check_absent_pattern_src "xcompact3d_coarse_bounds_initialized" "no local coarse-bounds init flag"
check_absent_pattern_src "init_xcompact3d_coarse_bounds" "no local coarse-bounds initializer"
check_present_pattern_src "^[[:space:]]*(use|USE)[[:space:]]+decomp_2d_io\b" "active sources import decomp_2d_io"
check_absent_pattern_src "^[[:space:]]*(use|USE)[[:space:]]+xcompact3d_decomp_io_compat\b.*fine_to_coarse(S|V)" "fine_to_coarseS/V are not imported via compat"
check_absent_pattern_file "^[[:space:]]*use[[:space:]]+param,[[:space:]]*only[[:space:]]*:.*\b(nx|ny|nz)\b" "src/visu.f90" "visu.f90 does not import nx/ny/nz from param"

if [ ${#failures[@]} -eq 0 ]; then
  echo "============================================================"
  echo "STAGE 9.1 FINAL VERDICT: PASS"
  echo "Reason: real 2decomp-fft interface cleanup is consistent."
  echo "============================================================"
  exit 0
else
  echo "============================================================"
  echo "STAGE 9.1 FINAL VERDICT: FAIL"
  echo "Failed checks:"
  for item in "${failures[@]}"; do
    echo "  - $item"
  done
  echo "============================================================"
  exit 1
fi
