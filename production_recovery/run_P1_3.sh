#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
export DECOMP2D_ROOT=${DECOMP2D_ROOT:-/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4}
production_recovery/P1_3_evidence/P1_3_VALIDATION_COMMAND.sh
cat production_recovery/P1_3_PASS_FAIL.md
cat production_recovery/P1_3_evidence/P1_3_VALIDATION_RESULT.txt
cat production_recovery/P1_3_evidence/P1_3_GUARDED_STABILITY_AUDIT.txt
cat production_recovery/P1_3_evidence/P1_3_RESTART_COMPATIBILITY_AUDIT.txt
cat production_recovery/P1_3_evidence/P1_3_STATISTICS_COMPATIBILITY_AUDIT.txt
cat production_recovery/P1_3_evidence/P1_3_VISUALIZATION_COMPATIBILITY_AUDIT.txt
