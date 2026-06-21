#!/usr/bin/env bash
set -euo pipefail

cd /home/sq/Incompact3d-sq-Xcompact3D

DECOMP_ROOT="/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4"
BUILD_DIR="build_r10_main_hook"

rm -rf "$BUILD_DIR"
mkdir -p production_recovery/R10_evidence

cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DDECOMP2D_ROOT="$DECOMP_ROOT" \
  -DCMAKE_PREFIX_PATH="$DECOMP_ROOT" \
  2>&1 | tee production_recovery/R10_BUILD_LOG.txt

cmake --build "$BUILD_DIR" --target fibre_prod_main_hook_check -j 2 \
  2>&1 | tee -a production_recovery/R10_BUILD_LOG.txt

if [ -x "$BUILD_DIR/bin/fibre_prod_main_hook_check" ]; then
  "$BUILD_DIR/bin/fibre_prod_main_hook_check" \
    2>&1 | tee production_recovery/R10_RUN_LOG_hook_check.txt
else
  "$BUILD_DIR/src/fibre_prod_main_hook_check" \
    2>&1 | tee production_recovery/R10_RUN_LOG_hook_check.txt
fi

cmake --build "$BUILD_DIR" --target xcompact3d -j 2 \
  2>&1 | tee -a production_recovery/R10_BUILD_LOG.txt

X3D_EXE=""
if [ -x "$BUILD_DIR/bin/xcompact3d" ]; then
  X3D_EXE="$(pwd)/$BUILD_DIR/bin/xcompact3d"
elif [ -x "$BUILD_DIR/src/xcompact3d" ]; then
  X3D_EXE="$(pwd)/$BUILD_DIR/src/xcompact3d"
else
  X3D_EXE="$(find "$BUILD_DIR" -type f -perm -111 -name xcompact3d | head -n 1)"
fi

mkdir -p production_recovery/R10_evidence/r10_run_lambda0
mkdir -p production_recovery/R10_evidence/r10_run_smalllambda

if [ -f input.i3d ]; then
  cp input.i3d production_recovery/R10_evidence/r10_run_lambda0/input.i3d
  cp input.i3d production_recovery/R10_evidence/r10_run_smalllambda/input.i3d
elif [ -f tests/Channel/input_test_x.i3d ]; then
  cp tests/Channel/input_test_x.i3d production_recovery/R10_evidence/r10_run_lambda0/input.i3d
  cp tests/Channel/input_test_x.i3d production_recovery/R10_evidence/r10_run_smalllambda/input.i3d
elif [ -f examples/Channel/input.i3d ]; then
  cp examples/Channel/input.i3d production_recovery/R10_evidence/r10_run_lambda0/input.i3d
  cp examples/Channel/input.i3d production_recovery/R10_evidence/r10_run_smalllambda/input.i3d
else
  echo "BLOCKED: no input.i3d candidate found" | tee production_recovery/R10_PASS_FAIL.md
  exit 1
fi

(
  cd production_recovery/R10_evidence/r10_run_lambda0
  FIBRE_PROD_ENABLE=1 \
  FIBRE_PROD_LAMBDA=0 \
  FIBRE_PROD_DIAGNOSTICS=1 \
  FIBRE_PROD_DIAGNOSTICS_DIR=.. \
  "$X3D_EXE"
) 2>&1 | tee production_recovery/R10_RUN_LOG_lambda0_np1.txt

(
  cd production_recovery/R10_evidence/r10_run_smalllambda
  FIBRE_PROD_ENABLE=1 \
  FIBRE_PROD_LAMBDA=1.0e-8 \
  FIBRE_PROD_DIAGNOSTICS=1 \
  FIBRE_PROD_DIAGNOSTICS_DIR=.. \
  "$X3D_EXE"
) 2>&1 | tee production_recovery/R10_RUN_LOG_smalllambda_np1.txt

HOOK_PASS=0
LAMBDA0_PASS=0
SMALL_PASS=0

grep -q "R10_FIBRE_PROD_MAIN_HOOK_CHECK PASS" production_recovery/R10_RUN_LOG_hook_check.txt && HOOK_PASS=1
grep -q "^Result: PASS" production_recovery/R10_evidence/R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt && LAMBDA0_PASS=1
grep -q "^Result: PASS" production_recovery/R10_evidence/R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt && SMALL_PASS=1

if [ "$HOOK_PASS" -eq 1 ] && [ "$LAMBDA0_PASS" -eq 1 ] && [ "$SMALL_PASS" -eq 1 ]; then
  printf 'PASS\n' > production_recovery/R10_PASS_FAIL.md
else
  printf 'FAIL\n' > production_recovery/R10_PASS_FAIL.md
fi

cat production_recovery/R10_PASS_FAIL.md
