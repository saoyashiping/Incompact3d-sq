# Stage 12.11 Total Smoke / Closure Gate

## Target

Stage 12.11 is the total-smoke and closure gate for Stage 12. It closes the production feedback-force-candidate diagnostic path by rerunning the stable final evidence chain and generating `stage12_checks/STAGE12_CLOSED.md` only when every required status passes.

Stage 12.11 adds no physics, no source modules, no production solver hooks, and no changes to existing closed-stage files.

## Mathematical and physical meaning

Stage 11 closed the one-way sampled-velocity path:

```text
U_f = I_h[u](X_f)
```

Stage 12 closes the Lagrangian feedback-force-candidate diagnostic path:

```text
slip = U_f - V_f
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
P_slip = sum F_fluid_to_fibre_cand · slip
P_structure = sum F_fluid_to_fibre_cand · V_f
P_fluid = sum F_fibre_to_fluid_cand · U_f
P_pair = P_structure + P_fluid
P_pair + P_slip = 0
```

Stage 12 remains diagnostic-only:

```text
f_fsi = 0
f_IBM = 0
RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

There is no Eulerian force density, no spreading, no RHS injection, and no fibre structure advance. Real Eulerian force-density construction and RHS-side coupling begin after Stage 12.

## Closure evidence

The total-smoke gate runs only the stable final Stage 12 closure gates by default:

- Stage 12.6 production feedback candidate hook smoke.
- Stage 12.7 np=1 force-candidate invariance.
- Stage 12.8 np=1/2/4 parallel force-candidate consistency.
- Stage 12.9 restart / stats / visu / coarse I/O compatibility.
- Stage 12.10 RHS / spreading / structure contamination audit.

The gate forces prerequisite-skip variables for those scripts so closed Stage 12.0--12.5 scripts are not rerun by default.

## Closed-stage warning

During Stage 12.11:

- do not modify Stage 10 files;
- do not modify Stage 11 files;
- do not modify Stage 12.0--12.10 files;
- do not modify production solver sources;
- do not run Stage 10.2 or Stage 10.3 audit logic;
- do not classify `no_*` diagnostic fields as activity.

## What is intentionally not done

Stage 12.11 intentionally does not add or enable:

- Eulerian force density;
- RHS injection;
- IBM spreading;
- feedback force application to the fluid;
- two-way force density;
- fibre structure advance;
- new production physics.

## Pass criteria

The gate passes only if all of the following are true:

1. `xcompact3d` and all required Stage 11 / Stage 12 check targets through Stage 12.6 build.
2. Stage 12.6, 12.7, 12.8, 12.9, and 12.10 final verdict lines are present and PASS.
3. The final Stage 12 feedback-candidate diagnostic file reports hook activity, finite force candidate diagnostics, finite power diagnostics, action-reaction consistency, pair-power consistency, no field modification, no RHS modification, no Eulerian force density, no RHS injection, no IBM spreading, no feedback application, no two-way force, and no structure advance.
4. The final Stage 12.10 audit `.dat` file reports all static and runtime contamination-audit statuses as pass.
5. `stage12_checks/STAGE12_CLOSED.md` is generated only after all Stage 12.11 statuses pass.

The gate writes `stage12_outputs/stage12_11_total_smoke.dat` and prints:

```text
STAGE 12.11 TOTAL SMOKE VERDICT: PASS
STAGE 12.11 FINAL VERDICT: PASS
```

Any failure prints an explicit non-empty `Reasons:` line, and `STAGE12_CLOSED.md` is not generated or overwritten on failure.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
bash stage12_checks/run_stage12_11_total_smoke.sh
```
