#!/usr/bin/env bash
set -euo pipefail

fail() {
  echo "[FAIL] $1" >&2
  exit 1
}

pass() {
  echo "[PASS] $1"
}

[ ! -e src/xcompact3d_decomp_io_compat.f90 ] || fail "src/xcompact3d_decomp_io_compat.f90 must be removed"
pass "compat source removed"

[ ! -e src/fibre_stage_decomp2d_constants_stub.f90 ] || fail "src/fibre_stage_decomp2d_constants_stub.f90 must be removed"
pass "decomp constants stub removed"

if rg -n "xcompact3d_decomp_io_compat\.f90" src/CMakeLists.txt >/dev/null; then
  fail "src/CMakeLists.txt still references xcompact3d_decomp_io_compat.f90"
fi
pass "CMakeLists has no compat source entries"

if rg -n "^[[:space:]]*use[[:space:]]+xcompact3d_decomp_io_compat\b|^[[:space:]]*USE[[:space:]]+xcompact3d_decomp_io_compat\b" src/*.f90 >/dev/null; then
  fail "active source still imports xcompact3d_decomp_io_compat"
fi
pass "no active compat module import in src/*.f90"

if rg -n "^[[:space:]]*use[[:space:]]+m_halo\b|^[[:space:]]*USE[[:space:]]+m_halo\b" src/*.f90 >/dev/null; then
  fail "active source still imports m_halo"
fi
pass "no active m_halo import in src/*.f90"

if rg -n "xcompact3d_coarse_bounds_initialized" src/*.f90 >/dev/null; then
  fail "xcompact3d_coarse_bounds_initialized residual found"
fi
pass "no local coarse-bounds init flag"

if rg -n "init_xcompact3d_coarse_bounds" src/*.f90 >/dev/null; then
  fail "init_xcompact3d_coarse_bounds residual found"
fi
pass "no local coarse-bounds initializer call/subroutine"

if ! rg -n "^[[:space:]]*use[[:space:]]+decomp_2d_io\b|^[[:space:]]*USE[[:space:]]+decomp_2d_io\b" src/*.f90 >/dev/null; then
  fail "decomp_2d_io import not found in active sources"
fi
pass "decomp_2d_io is used by active sources"

if ! rg -n "^[[:space:]]*use[[:space:]]+decomp_2d\b.*fine_to_coarseS|^[[:space:]]*use[[:space:]]+decomp_2d\b.*fine_to_coarseV|^[[:space:]]*USE[[:space:]]+decomp_2d\b.*fine_to_coarseS|^[[:space:]]*USE[[:space:]]+decomp_2d\b.*fine_to_coarseV" src/*.f90 >/dev/null; then
  fail "fine_to_coarseS/fine_to_coarseV are not imported from decomp_2d"
fi
pass "coarse conversion symbols are imported from decomp_2d"

if rg -n "target_link_libraries\([^)]*decomp2d" src/CMakeLists.txt >/dev/null; then
  pass "targets link decomp2d"
else
  fail "no decomp2d target links found in src/CMakeLists.txt"
fi

echo "Stage 9.1 interface consistency checks passed."
