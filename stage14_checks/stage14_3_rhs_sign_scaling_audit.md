# Stage 14.3 — RHS Sign and Scaling Audit

## Target

Stage 14.3 audits the standalone RHS increment sign and scaling that will be used by later controlled Stage 14 production injection work.

This stage verifies that the RHS increment uses the fibre-on-fluid force-density sign:

```text
Delta_RHS_14 = lambda_14 * f_i_cand
```

where `f_i_cand` represents the Stage 13 Eulerian force-density candidate formed from `F_fibre_to_fluid_cand`. It also rejects the wrong fluid-on-fibre sign and verifies linear scaling with `lambda_14`.

Stage 14.3 performs these checks only in controlled standalone arrays. It does not add anything to the production RHS.

## Mathematical / physical meaning

Stage 12 established the sign convention:

```text
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
```

Stage 13 provides the Eulerian force-density candidate:

```text
f_i_cand = S_h[F_fibre_to_fluid_cand]
```

Stage 14 will later use the RHS increment:

```text
Delta_RHS_14 = lambda_14 * f_i_cand
```

Stage 14.3 verifies the volume-integrated standalone identity:

```text
integrated Delta_RHS_14 = lambda_14 * sum_k F_fibre_to_fluid_cand(k) * ds_k
```

The wrong-sign increment is rejected:

```text
Delta_RHS_wrong = lambda_14 * S_h[F_fluid_to_fibre_cand]
```

At Stage 14.3, production FSI remains disabled:

```text
f_fsi = 0
RHS_stage14.3 = RHS_stage14.2 = RHS_stage14.1 = RHS_stage14.0 = RHS_stage13 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Controlled checks

The Stage 14.3 standalone check covers:

- action-reaction consistency pointwise;
- fibre-on-fluid RHS increment sign;
- wrong fluid-on-fibre sign rejection;
- `lambda=0` zero increment;
- `lambda=1` scaling;
- `lambda=0.1` scaling;
- `lambda=-0.1` sign reversal;
- component-wise x/y/z scaling;
- volume-integrated RHS increment sign and scaling;
- finite RHS increment values;
- additive formula consistency through the Stage 14.2 standalone formula module.

## What is intentionally not done

Stage 14.3 intentionally performs no production coupling. In particular, it does not:

- add an `xcompact3d` hook call;
- insert code in the production main loop;
- add any increment to production RHS arrays;
- call production IBM forcing;
- modify pressure, projection, Poisson, RK3, or channel forcing code;
- access or modify production fluid fields;
- advance the fibre structure.

## Pass criteria

The check passes only if all required status keys are `1`, all absolute error diagnostics are within `1.0e-12`, and the wrong-sign separation is at least `1.0e-8`.

The shell gate writes:

```text
stage14_outputs/stage14_3_rhs_sign_scaling_audit_gate.dat
```

and the standalone executable writes:

```text
stage14_outputs/fibre_stage14_3_rhs_sign_scaling_audit.dat
```

Both verdicts must report PASS.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE14_3_RUN_STAGE14_2=0 \
bash stage14_checks/run_stage14_3_rhs_sign_scaling_audit.sh
```
