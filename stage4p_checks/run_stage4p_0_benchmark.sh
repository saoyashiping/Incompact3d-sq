#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p stage4p_outputs build_stage4p
export NP=${NP:-1}
export NX=${NX:-128}
export NY=${NY:-65}
export NZ=${NZ:-128}
export NL=${NL:-129}
export NSTEPS=${NSTEPS:-200}
export DT=${DT:-1.0e-5}
export BETA_DRAG=${BETA_DRAG:-10.0}
cmake -S . -B build_stage4p
cmake --build build_stage4p --target fibre_stage4p_benchmark_check -j
exe="build_stage4p/bin/fibre_stage4p_benchmark_check"
if [[ ! -x "$exe" ]]; then
  exe="$(find build_stage4p -type f -name 'fibre_stage4p_benchmark_check' | head -n 1)"
fi
if [[ "$NP" -gt 1 ]] && command -v mpirun >/dev/null 2>&1; then
  mpirun -np "$NP" "$exe"
else
  "$exe"
fi
cat stage4p_outputs/stage4p_benchmark_summary.dat
