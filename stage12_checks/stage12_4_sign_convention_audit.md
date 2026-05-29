# Stage 12.4 sign convention audit

## Target

Stage 12.4 adds a standalone feedback force sign convention audit. It explicitly
distinguishes the candidate force exerted by the fluid on the fibre/control
point from the equal-and-opposite candidate force exerted by the fibre/control
point on the fluid.

## Mathematical and physical meaning

The Stage 12.4 convention is:

```text
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
```

If `U_f - V_f` is positive in a component, the fluid velocity is larger than the
fibre/control-point velocity in that component, so the fluid-on-fibre candidate
force is positive in that component. The fibre-on-fluid candidate is the
action-reaction counterpart.

For structure-side notation:

- if the structure equation uses `+F_fluid_to_fibre`, the structure-side sign is
  direct;
- if the project equation uses `-F_fs`, then `F_fs` must be interpreted as
  `F_fibre_to_fluid`, so `-F_fs = F_fluid_to_fibre`.

Future Eulerian spreading/RHS stages must use `F_fibre_to_fluid`, not
`F_fluid_to_fibre`. Stage 12.4 does not spread or inject any force, so:

```text
f_IBM = 0
RHS_stage12.4 = RHS_stage12.3 = RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Controlled tests

The standalone audit covers:

- zero slip;
- positive slip;
- negative slip;
- mixed multicomponent slip;
- action-reaction consistency;
- structure-side convention encoding;
- future fluid-side convention encoding;
- finite force and norm diagnostics.

## Intentionally not done

Stage 12.4 does not add or activate:

- an `xcompact3d` hook call;
- production main-loop insertion;
- Eulerian force density allocation;
- RHS injection;
- IBM spreading;
- feedback force application to the fluid;
- two-way force density;
- fibre structure advance.

## Pass criteria

The gate passes only if `fibre_stage12_sign_convention_audit_check` builds, the
check prints `STAGE 12.4 SIGN CONVENTION AUDIT VERDICT: PASS`, all required
diagnostic status keys equal `1`, all convention diagnostics equal `1`, and all
scalar audit errors are less than or equal to `1.0e-12`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE12_4_RUN_STAGE12_3=0 \
bash stage12_checks/run_stage12_4_sign_convention_audit.sh
```
