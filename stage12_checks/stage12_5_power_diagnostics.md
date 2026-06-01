# Stage 12.5 power diagnostics gate

## Target

Stage 12.5 adds standalone force work / power diagnostic candidate checks on
controlled Lagrangian arrays. It computes diagnostic scalars only and does not
activate production coupling.

## Mathematical and physical meaning

Using the Stage 12.4 sign convention:

```text
slip = U_f - V_f
F_fluid_to_fibre = alpha * slip
F_fibre_to_fluid = -F_fluid_to_fibre
```

The standalone diagnostics are:

```text
force_norm(k) = sqrt(sum(F_fluid_to_fibre(k,1:3)**2))
total_force_l2 = sqrt(sum_k sum_c F_fluid_to_fibre(k,c)**2)
P_slip = sum_k F_fluid_to_fibre(k,:) · slip(k,:)
P_structure = sum_k F_fluid_to_fibre(k,:) · V_f(k,:)
P_fluid = sum_k F_fibre_to_fluid(k,:) · U_f(k,:)
P_pair = P_structure + P_fluid
P_pair + P_slip = 0
```

For `alpha > 0`, `P_slip = alpha * sum_k |U_f - V_f|^2 >= 0`, and zero slip
requires zero force and zero power. Stage 12.5 keeps `f_fsi = 0`, does not spread
or inject force, and preserves:

```text
RHS_stage12.5 = RHS_stage12.4 = RHS_stage12.3 = RHS_stage12.2
RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Controlled tests

The standalone check covers:

- zero-slip power;
- positive-slip power;
- multicomponent finite diagnostics;
- gain scaling, with `P_slip` and `total_force_l2` linear in `alpha` for
  `F = alpha * slip`;
- action-reaction power consistency;
- pair-power consistency.

## Intentionally not done

Stage 12.5 does not add or activate:

- an `xcompact3d` hook call;
- production main-loop insertion;
- Eulerian force density allocation;
- RHS injection;
- IBM spreading;
- feedback force application to the fluid;
- two-way force density;
- fibre structure advance.

## Pass criteria

The gate passes only if `fibre_stage12_power_diagnostics_check` builds, the check
prints `STAGE 12.5 POWER DIAGNOSTICS VERDICT: PASS`, all required diagnostic
status keys equal `1`, and all scalar audit errors are less than or equal to
`1.0e-12`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE12_5_RUN_STAGE12_4=0 \
bash stage12_checks/run_stage12_5_power_diagnostics.sh
```
