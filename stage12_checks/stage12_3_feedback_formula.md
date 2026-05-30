# Stage 12.3 feedback formula gate

## Target

Stage 12.3 adds controlled standalone feedback force formula checks for the future
Lagrangian/control-point feedback candidate. It verifies the analytic relation

```text
F_fs_cand = alpha * (U_f - V_f)
```

only in isolated check code and does not connect the formula to the production
solver.

## Mathematical and physical meaning

The production no-fibre DNS remains unchanged:

```text
div(u) = 0
du/dt + u · grad(u) = -grad(p) + nu laplacian(u) + f_drive + f_fsi
f_fsi = 0
RHS_stage12.3 = RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

Stage 12.3 computes `F_fs_cand = alpha * (U_f - V_f)` only for controlled
Lagrangian arrays owned by the standalone check. The sampled velocity `U_f` was
closed in Stage 11, and the prescribed/placeholder velocity `V_f` was introduced
in Stage 12.2, but no production feedback path is activated here.

## Controlled tests

The standalone check covers:

- zero slip: `U_f = V_f` gives `F_fs_cand = 0`;
- positive slip: positive `U_f - V_f` keeps the expected sign;
- negative slip: negative `U_f - V_f` keeps the expected sign;
- `force_sign = -1`: the configured sign reverses the candidate force;
- multicomponent slip: x/y/z components are checked independently;
- gain scaling: doubling `alpha` doubles `F_fs_cand`;
- finite candidate force and finite force norm diagnostics.

## Intentionally not done

Stage 12.3 does not add or activate:

- an `xcompact3d` hook call;
- production main-loop insertion;
- Eulerian force density allocation;
- RHS injection;
- IBM spreading;
- feedback force application to the fluid;
- two-way force density;
- fibre structure advance.

## Pass criteria

The gate passes only if `fibre_stage12_feedback_formula_check` builds, the check
prints `STAGE 12.3 FEEDBACK FORMULA VERDICT: PASS`, all required diagnostic
status keys equal `1`, and all scalar formula errors are less than or equal to
`1.0e-12`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE12_3_RUN_STAGE12_2=0 \
bash stage12_checks/run_stage12_3_feedback_formula.sh
```
