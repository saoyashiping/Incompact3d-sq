#!/usr/bin/env bash
set -euo pipefail

: "${DECOMP2D_ROOT:=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p stage9_outputs

EVIDENCE="stage9_outputs/stage9_1_build_evidence.dat"
: > "$EVIDENCE"

emit(){
  echo "$1 $2" >> "$EVIDENCE"
}

emit stage9_build_evidence_dns_executed_flag 0
emit stage9_build_evidence_rhs_modified_flag 0
emit stage9_build_evidence_projection_called_flag 0
emit stage9_build_evidence_fluid_update_called_flag 0
emit stage9_build_evidence_fibre_advance_called_flag 0
emit stage9_build_evidence_stage_checks_only_flag 0

if [[ -d "$DECOMP2D_ROOT" ]]; then
  emit stage9_build_evidence_decomp_root_exists 1
else
  emit stage9_build_evidence_decomp_root_exists 0
  echo "ERROR: DECOMP2D_ROOT=$DECOMP2D_ROOT does not exist."
  exit 1
fi

if [[ -f "$DECOMP2D_ROOT/include/decomp_2d.mod" ]]; then
  emit stage9_build_evidence_decomp_mod_exists 1
else
  emit stage9_build_evidence_decomp_mod_exists 0
  echo "ERROR: decomp_2d.mod not found under $DECOMP2D_ROOT/include"
  exit 1
fi

if [[ -f "$DECOMP2D_ROOT/include/decomp_2d_io.mod" ]]; then
  emit stage9_build_evidence_decomp_io_mod_exists 1
else
  emit stage9_build_evidence_decomp_io_mod_exists 0
  echo "ERROR: decomp_2d_io.mod not found under $DECOMP2D_ROOT/include"
  exit 1
fi

if [[ -f "$DECOMP2D_ROOT/lib/libdecomp2d.a" ]]; then
  emit stage9_build_evidence_decomp_static_lib_exists 1
else
  emit stage9_build_evidence_decomp_static_lib_exists 0
  echo "ERROR: libdecomp2d.a not found under $DECOMP2D_ROOT/lib"
  exit 1
fi

case "$DECOMP2D_ROOT" in
  *2decomp-fft-xcompact3d-v2.0.4*)
    emit stage9_build_evidence_used_source_matched_decomp_flag 1
    ;;
  *)
    emit stage9_build_evidence_used_source_matched_decomp_flag 0
    echo "ERROR: DECOMP2D_ROOT does not point to source-matched 2decomp-fft-xcompact3d-v2.0.4."
    exit 1
    ;;
esac

compat_hits="$(grep -RIn "xcompact3d_decomp_io_compat" src stage9_checks || true)"
if [[ -z "$compat_hits" ]]; then
  emit stage9_build_evidence_no_compat_remnants_flag 1
else
  emit stage9_build_evidence_no_compat_remnants_flag 0
  echo "ERROR: xcompact3d_decomp_io_compat remnants found:"
  echo "$compat_hits"
  exit 1
fi

rm -rf build_stage9
cmake -S . -B build_stage9 -DCMAKE_PREFIX_PATH="$DECOMP2D_ROOT"
emit stage9_build_evidence_configure_status 1

if grep -R "FIBRE_STAGE_CHECKS_ONLY:BOOL=ON" build_stage9/CMakeCache.txt >/dev/null 2>&1; then
  emit stage9_build_evidence_fibre_stage_checks_only_off_flag 0
  echo "ERROR: FIBRE_STAGE_CHECKS_ONLY is ON."
  exit 1
else
  emit stage9_build_evidence_fibre_stage_checks_only_off_flag 1
fi

cmake --build build_stage9 --target xcompact3d
emit stage9_build_evidence_xcompact3d_compile_status 1

if grep -q "stage8_total_default_production_disabled_status 1" stage8_outputs/fibre_stage8_total_smoke_check.dat; then
  emit stage9_build_evidence_default_production_disabled_flag 1
else
  emit stage9_build_evidence_default_production_disabled_flag 0
fi

cmake --build build_stage9 --target fibre_stage9_build_regression_check

EXE="build_stage9/bin/fibre_stage9_build_regression_check"
if [[ ! -x "$EXE" ]]; then
  EXE="$(find build_stage9 -type f -name fibre_stage9_build_regression_check | head -n 1)"
fi

if [[ -z "${EXE:-}" || ! -x "$EXE" ]]; then
  echo "ERROR: executable fibre_stage9_build_regression_check not found"
  exit 1
fi

"$EXE"

if [[ ! -s stage9_outputs/fibre_stage9_build_regression_check.dat ]]; then
  echo "ERROR: stage9_outputs/fibre_stage9_build_regression_check.dat was not generated"
  exit 1
fi

cat stage9_outputs/fibre_stage9_build_regression_check.dat
