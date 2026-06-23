#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
export DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
production_recovery/P1_4_evidence/P1_4_VALIDATION_COMMAND.sh
cat production_recovery/P1_4_PASS_FAIL.md
cat production_recovery/P1_4_evidence/P1_4_VALIDATION_RESULT.txt
cat production_recovery/P1_4_evidence/P1_4_SELF_CONTAINED_VALIDATION_AUDIT.txt
cat production_recovery/P1_4_evidence/P1_4_LAMBDA0_NP_CONSISTENCY_AUDIT.txt
cat production_recovery/P1_4_evidence/P1_4_TWOWAY_NP_CONSISTENCY_AUDIT.txt
cat production_recovery/P1_4_evidence/P1_4_P1_CLOSURE_AUDIT.txt
