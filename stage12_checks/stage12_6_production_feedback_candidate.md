# Stage 12.6 production feedback candidate hook

## Target

Stage 12.6 adds a guarded production feedback-force candidate hook.  The hook reads production velocity in read-only smoke mode, computes Lagrangian candidate diagnostics, and writes diagnostics proving that no Eulerian force density, fluid-force application, or RHS modification is performed.

## Mathematical / physical meaning

Stage 12.6 keeps the production no-fibre DNS unchanged:

```text
f_fsi = 0
RHS_stage12.6 = RHS_stage12.5 = RHS_stage12.4 = RHS_stage12.3 = RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

The production smoke hook samples a conservative read-only velocity candidate from the production fields:

```text
U_f = sampled production velocity
V_f = prescribed / placeholder velocity
slip = U_f - V_f
F_fluid_to_fibre_cand = alpha * slip
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
P_slip = F_fluid_to_fibre_cand · slip
P_structure = F_fluid_to_fibre_cand · V_f
P_fluid = F_fibre_to_fluid_cand · U_f
P_pair = P_structure + P_fluid
P_pair + P_slip = 0
```

The candidate is diagnostic only.  It is not spread to the Eulerian mesh and it is not injected into the fluid equations.

## Implementation policy

- Use conservative production smoke sampling from safe in-bounds velocity-array indices.
- Treat all production velocity arrays as `intent(in)` in the Stage 12.6 hook API.
- Compute finite force-candidate signatures, action-reaction consistency, and power diagnostics.
- Use the Stage 12.4 sign convention directly: fluid-to-fibre is `+alpha * slip` and fibre-to-fluid is its reaction.
- Do not perform physical feedback, Eulerian spreading, RHS injection, or fibre-structure advancement.

## Intentionally not done

- No Eulerian force density is allocated.
- No RHS injection is performed.
- No IBM spreading is called.
- No feedback force is applied to fluid.
- No two-way force is activated.
- No fibre structure equation is advanced.
- No pressure, projection, RK3, Poisson, or channel-forcing code is modified.

## Pass criteria

The Stage 12.6 gate passes only if all of the following hold:

1. `xcompact3d`, all required Stage 11 check targets, all Stage 12 check targets through Stage 12.6, and the new standalone check target build.
2. The standalone synthetic-array check prints `STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE CHECK VERDICT: PASS`.
3. The Stage 9.9 deterministic no-fibre production smoke run still passes while Stage 11 and Stage 12.6 read-only hooks are enabled.
4. `stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat` reports initialized, sampled, finite force candidate, finite power diagnostics, action-reaction, and pair-power consistency statuses as `1`.
5. `stage12_6_field_modified_status` and `stage12_6_rhs_modified_status` remain `0`.
6. All no-coupling statuses remain `1`: no Eulerian force density, no RHS injection, no IBM spreading, no feedback application, no two-way force, and no structure advance.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE12_6_RUN_STAGE12_5=0 \
bash stage12_checks/run_stage12_6_production_feedback_candidate.sh
```
