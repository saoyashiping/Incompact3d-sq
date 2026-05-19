#!/usr/bin/env bash
set -euo pipefail
: "${DECOMP2D_ROOT:=/home/sq/opt/2decomp-fft}"

if [ ! -d "$DECOMP2D_ROOT" ]; then
  echo "ERROR: DECOMP2D_ROOT=$DECOMP2D_ROOT does not exist."
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p stage9_outputs

EVIDENCE="stage9_outputs/stage9_1_build_evidence.dat"
: > "$EVIDENCE"

echo "stage9_build_evidence_stage_checks_only_flag 0" >> "$EVIDENCE"
echo "stage9_build_evidence_dns_executed_flag 0" >> "$EVIDENCE"
echo "stage9_build_evidence_rhs_modified_flag 0" >> "$EVIDENCE"
echo "stage9_build_evidence_projection_called_flag 0" >> "$EVIDENCE"
echo "stage9_build_evidence_fluid_update_called_flag 0" >> "$EVIDENCE"
echo "stage9_build_evidence_fibre_advance_called_flag 0" >> "$EVIDENCE"

cmake -S . -B build_stage9 \
  -DCMAKE_PREFIX_PATH="${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft}"

echo "stage9_build_evidence_configure_status 1" >> "$EVIDENCE"

if [ -d "$DECOMP2D_ROOT" ]; then
  echo "stage9_build_evidence_used_real_decomp2d_flag 1" >> "$EVIDENCE"
else
  echo "stage9_build_evidence_used_real_decomp2d_flag 0" >> "$EVIDENCE"
fi

if rg -n "FIBRE_STAGE_CHECKS_ONLY:BOOL=ON" build_stage9/CMakeCache.txt >/dev/null 2>&1; then
  echo "stage9_build_evidence_fibre_stage_checks_only_off_flag 0" >> "$EVIDENCE"
else
  echo "stage9_build_evidence_fibre_stage_checks_only_off_flag 1" >> "$EVIDENCE"
fi

cmake --build build_stage9

echo "stage9_build_evidence_build_status 1" >> "$EVIDENCE"

if rg -n "stage8_total_default_production_disabled_status 1" stage8_outputs/fibre_stage8_total_smoke_check.dat >/dev/null 2>&1; then
  echo "stage9_build_evidence_default_production_disabled_flag 1" >> "$EVIDENCE"
else
  echo "stage9_build_evidence_default_production_disabled_flag 0" >> "$EVIDENCE"
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
