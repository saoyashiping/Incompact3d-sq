# Stage 9.10 no-coupling contamination / no-op invariance audit

Stage 9.10 audits that production no-fibre channel DNS remains uncontaminated by fibre/IBM/FSI coupling code paths.

## Mathematical/physical meaning
Stage 9.10 introduces **no new physics**. It is an invariance audit for the existing no-fibre incompressible DNS system:

- `div(u)=0`
- `du/dt + u·grad(u) = -grad(p) + nu*laplacian(u) + f_drive`

The objective is to verify coupling hooks are inactive and production no-fibre output remains unchanged.

## Static audit
The script statically scans active production sources (not comments/docs/diagnostic modules) and checks:

- no active Stage 8 velocity-to-fibre production call;
- no active Stage 8 feedback path;
- no active Stage 8 two-way force-density path;
- no active IBM spreading hook in production RHS path;
- no active fibre structure advance in production time loop;
- no active fibre force injection into production RHS;
- no reintroduction of forbidden Stage 9.5 placeholders.

## Runtime no-op invariance audit
Two deterministic no-fibre np=1 Stage 9.9 runs are compared:

1. baseline run;
2. audit run with inert flag `X3D_STAGE9_10_NO_COUPLING_AUDIT=1`.

Both runs must produce matching final Stage 9.9 signatures within mixed tolerance:

- `STAGE9_10_INVARIANCE_ABS_TOL` (default `1.0e-12`)
- `STAGE9_10_INVARIANCE_REL_TOL` (default `1.0e-14`)

## Pass/fail criteria
Stage 9.10 passes only if all are true:

- build of `xcompact3d` succeeds;
- prerequisite Stage 9.1–9.9 checks pass (unless `STAGE9_SKIP_PREREQS=1`);
- static audit passes;
- runtime baseline passes;
- runtime audit passes;
- baseline vs audit signatures satisfy invariance tolerances;
- dat key `stage9_10_no_coupling_contamination_status` is `1`.

## Output
Script writes:

- `stage9_outputs/stage9_10_no_coupling_contamination_audit.dat`

and prints:

- `STAGE 9.10 NO-COUPLING CONTAMINATION AUDIT VERDICT: PASS|FAIL`
- `STAGE 9.10 FINAL VERDICT: PASS|FAIL`

## Manual command
`bash stage9_checks/run_stage9_10_no_coupling_contamination_audit.sh`

## Intentionally not tested yet
- Real production coupling activation (later stages);
- coupled-fibre/IBM physics behavior;
- long-time coupled dynamics.
