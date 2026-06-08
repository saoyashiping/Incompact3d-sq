# Stage 17.11: Stage 17 total contamination audit and closure

Stage 17.11 is the final diagnostic-only Stage 17 audit.  It verifies that Stage 17.0--17.10 evidence is preserved, Stage 16 closed-loop compatibility remains valid, wall/contact/fibre-fibre placeholder force channels remain zero, and no RHS/IBM/DNS-core/structure contamination has been introduced.

## Scope

The Stage 17.11 wrapper creates `stage17_outputs/`, runs no production physics by default, builds nothing, and invokes only the Stage 17.11 helper unless explicitly configured to run the closed Stage 17.10 wrapper first.  The helper writes `stage17_outputs/fibre_stage17_11_total_contamination_audit_closure.dat` and creates `stage17_checks/STAGE17_CLOSED.md` only after every required Stage 17.11 status passes.

## Closure checks

The audit verifies:

- Stage 17.0 preflight through Stage 17.10 parallel/restart/I/O evidence is present or safely accepted from closed structural evidence.
- Stage 11 sampling, Stage 12 feedback, Stage 16.4 force input, Stage 13 force-density candidate, and Stage 14 RHS-path evidence remain preserved.
- Parallel, restart/I/O, and stats/visu/coarse-I/O compatibility evidence remains preserved.
- Contact and fibre-fibre force/RHS/structure-update norms are zero.
- Production remains single-fibre and fibre-fibre placeholder geometry remains standalone-only.
- No wall/contact/collision force, production multi-fibre logic, direct RHS injection, IBM forcing, pressure/projection/Poisson/RK3/channel-forcing contamination, or structure-dynamics activation is introduced.

## False-positive-safe audit policy

Stage 17.11 continues the corrected Stage 16 and Stage 17.0--17.10 policy without editing closed files.  It does not scan Markdown as code-regression evidence, does not treat documentation or negative-check strings as real behavior, does not require `rg`, and treats missing `.git` metadata in source-only archives as non-contamination.  The Stage 13.6 diagnostic check targets the actual production candidate files, and the Stage 14 small-lambda hook check targets `src/fibre_stage14_production_rhs_injection.f90` plus `src/xcompact3d.f90`.

## Expected verdict

With closed Stage 16 and Stage 17.0--17.10 evidence intact and default Stage 17.11 settings, the wrapper should print:

```text
STAGE 17.11 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: PASS
STAGE 17.11 FINAL VERDICT: PASS
```

Only a full PASS path may create `stage17_checks/STAGE17_CLOSED.md`; any failure emits explicit `reason` lines and does not create closure.
