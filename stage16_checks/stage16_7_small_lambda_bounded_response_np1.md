# Stage 16.7 small-lambda bounded response at np=1

Stage 16.7 is a diagnostic-only closure check for the already validated one-fibre closed-loop path:

```text
Stage 11 sampling
-> Stage 12 feedback-force candidate
-> Stage 16.4 structure-side fluid-on-fibre force input
-> controlled structure-state diagnostic update
-> Stage 13 force-density candidate signature
-> Stage 14 RHS diagnostic/controlled injection path
```

The check compares a zero-lambda baseline against a very small positive lambda.  It proves that the approved closed-loop path produces a finite, nonzero, bounded RHS increment for the small-lambda case while keeping the fluid signature finite and bounded.  It does not connect Stage 16.7 into `xcompact3d.f90`, does not alter pressure/projection/Poisson/RK3/channel-forcing numerics, and does not activate wall/contact, multi-fibre, or legacy IBM production forcing.

## User command

```bash
bash stage16_checks/run_stage16_7_small_lambda_bounded_response_np1.sh
```

The wrapper defaults to `BUILD_DIR=build_stage9`, configures the build directory with CMake if needed, builds only `fibre_stage16_small_lambda_response_check`, runs the diagnostic under `np=1`, and then invokes the Stage 16.7 Python assertion helper.

## Default controls

The wrapper supports the Stage 16.7 environment variables requested for this stage, including:

- `STAGE16_7_ZERO_LAMBDA=0.0`
- `STAGE16_7_SMALL_LAMBDA=1.0e-8`
- `STAGE16_7_NP=1`
- `STAGE16_7_NPTS=8`
- `STAGE16_7_MAX_ZERO_RHS_INCREMENT=1.0e-14`
- `STAGE16_7_MIN_RHS_INCREMENT=1.0e-20`
- `STAGE16_7_MAX_RHS_INCREMENT=1.0e-8`
- `STAGE16_7_MAX_FLUID_DELTA=1.0e-8`

## Expected output

The check writes:

```text
stage16_outputs/fibre_stage16_7_small_lambda_bounded_response_np1.dat
```

A passing run prints:

```text
STAGE 16.7 SMALL-LAMBDA BOUNDED RESPONSE NP1 VERDICT: PASS
STAGE 16.7 FINAL VERDICT: PASS
```

## False-positive-safe static audit policy

The Stage 16.7 helper intentionally reuses the corrected Stage 16.6 / Stage 16.5 / Stage 16.4 helper pattern:

- documentation is checked for required-file existence, not scanned as executable regression evidence;
- negative-check strings are not treated as production behavior;
- regex literals such as `rg[[:space:]]` are not treated as real `rg` command usage;
- legitimate Stage 13.5 conservation/sign audit files are not classified as old production diagnostic regressions;
- old Stage 13.5 production force-density names are rejected only in real production/check logic;
- if evidence cannot be distinguished safely, the helper fails closed with explicit reasons rather than silently passing.
