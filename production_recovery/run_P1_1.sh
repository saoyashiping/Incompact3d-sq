#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
export DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
export CMAKE_PREFIX_PATH="$DECOMP2D_ROOT:${CMAKE_PREFIX_PATH:-}"

bash production_recovery/P1_1_evidence/P1_1_VALIDATION_COMMAND.sh

cat production_recovery/P1_1_PASS_FAIL.md
cat production_recovery/P1_1_evidence/P1_1_VALIDATION_RESULT.txt
