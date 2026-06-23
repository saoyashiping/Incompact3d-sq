#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
export DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
production_recovery/P1_2_evidence/P1_2_VALIDATION_COMMAND.sh
cat production_recovery/P1_2_PASS_FAIL.md
cat production_recovery/P1_2_evidence/P1_2_VALIDATION_RESULT.txt
cat production_recovery/P1_2_evidence/P1_2_LAMBDA_SCALING_AUDIT.txt
cat production_recovery/P1_2_evidence/P1_2_RHS_INCREMENT_AUDIT.txt
