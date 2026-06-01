# Stage 13.5 conservation sign audit

## Target

Stage 13.5 adds a controlled force-density conservation and sign audit for the
Eulerian force-density candidate.  It verifies that the standalone spreading
path uses the Stage 12 fluid-side force sign correctly before any production
coupling is considered.

## Mathematical / physical meaning

Stage 12 defines the feedback-force candidate sign convention as

```text
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
```

Stage 13 force-density candidates must spread the fibre-to-fluid force, not the
fluid-to-fibre force:

```text
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

The conservative target is the volume-weighted integral

```text
sum_i f_i_cand * DeltaV_i ~= sum_k F_fibre_to_fluid_cand(k) * Delta_s_k
```

Stage 13.5 also confirms that the same Eulerian integral is opposite to the
wrong fluid-to-fibre total force.  The production no-fibre DNS remains
unchanged:

```text
f_fsi = 0
RHS_stage13.5 = RHS_stage13.4 = RHS_stage13.3 = RHS_stage13.2
              = RHS_stage13.1 = RHS_stage13.0 = RHS_stage12
              = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Controlled tests

The standalone check covers:

- action-reaction pointwise consistency;
- spreading input sign uses `F_fibre_to_fluid_cand`;
- integrated Eulerian force matches the fibre-to-fluid weighted sum;
- wrong-sign rejection against the fluid-to-fibre weighted sum;
- component-wise sign conservation;
- multipoint mixed-`ds` conservation;
- nonuniform-volume sign conservation;
- near-boundary sign conservation;
- finite force-density values;
- clear behavior for all candidate force-density arrays and norms.

## Intentionally not done

Stage 13.5 intentionally performs no production coupling:

- no `xcompact3d` hook call;
- no production main-loop insertion;
- no production RHS injection;
- no production IBM spreading;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance.

## Pass criteria

The gate passes only when the standalone executable prints

```text
STAGE 13.5 CONSERVATION SIGN AUDIT VERDICT: PASS
```

and the diagnostic file reports all required status keys as `1`, all required
conservation/action-reaction scalar errors at or below `1.0e-12`, and the
wrong-sign separation at or above `1.0e-8`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE13_5_RUN_STAGE13_4=0 \
bash stage13_checks/run_stage13_5_conservation_sign_audit.sh
```
