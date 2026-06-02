Stage 14.11 Total Smoke and Closure
Objective
Stage 14.11 is the Stage 14 total-smoke and closure evidence path. It verifies that the controlled RHS-injection chain remains valid and generates `stage14\_checks/STAGE14\_CLOSED.md` only after every required Stage 14.11 check passes.
The chain under closure is:
Stage 11 one-way sampling.
Stage 12 Lagrangian feedback-force candidate.
Stage 13 Eulerian force-density candidate.
Stage 14 controlled RHS injection.
Stage 14.11 introduces no new physics and does not advance the flexible-fibre structure.
Controlled RHS formula
The Stage 14 controlled update remains:
```text
RHS\_new = RHS\_old + lambda\_14 \* f\_i\_cand
```
Here `f\_i\_cand` is the Stage 13 Eulerian force-density candidate using the already-audited fibre-on-fluid sign convention.
Required regression prohibitions
Stage 14.11 fails if any of these regressions are detected:
The forbidden Stage 14 hook-registration gate `stage14\_get\_injection\_gain() == 0.0` is present.
Small-lambda Stage 14 hook registration is blocked.
Stage 11.5, Stage 13, or Stage 14 production diagnostics are missing or not rank0-safe/race-free.
Stage 13 force-density diagnostics revert to local subdomain center sampling.
Stage 14 touches pressure/projection/Poisson, RK3, or channel forcing code.
Production IBM forcing is activated outside the approved Stage 13/14 diagnostic-to-RHS path.
Any structure advance, bending solve, tension solve, fibre position update, structural velocity update, or wall/contact logic is active in the Stage 14 production path.
Required PASS evidence is silently skipped.
Checks performed
Static regression audit
The wrapper checks that:
Stage 14 hook registration is not gated by `lambda\_14 == 0`.
The production Stage 14 RHS hook is connected to the Stage 13 force-density candidate and Stage 14 request/enable controls.
Stage 11.5, Stage 12.6, Stage 13.6, and Stage 14.5 diagnostic write paths exist.
Stage 11/13/14 production diagnostic writes are rank0-safe or otherwise race-free.
Stage 13 force-density diagnostics do not use local subdomain center sampling.
Stage 14 code has no pressure/projection/Poisson, RK3, channel-forcing, legacy IBM, or structure-advance contamination.
Runtime total smoke
The wrapper runs one small-lambda Stage 14 production smoke case with `STAGE14\_11\_SMALL\_LAMBDA` and `STAGE14\_11\_NP` MPI ranks. It requires:
Stage 11 one-way sampling active.
Stage 12 feedback-force candidate active.
Stage 13 Eulerian force-density candidate active.
Stage 14 RHS hook active.
`lambda\_14` recorded and nonzero.
RHS increment nonzero, finite, sign-correct, and bounded by `STAGE14\_11\_MAX\_RHS\_INCREMENT`.
Fluid response bounded by `STAGE14\_11\_MAX\_FLUID\_DELTA`.
No NaN/Inf in runtime logs or diagnostics.
Prior Stage 14 evidence preservation
Stage 14.11 verifies existing PASS datfiles by default and can rerun closed prerequisites only when explicitly requested:
`STAGE14\_11\_RUN\_STAGE14\_8=1` force-reruns the Stage 14.8 np=1/2/4 parallel small-lambda response gate.
`STAGE14\_11\_RUN\_STAGE14\_9=1` force-reruns the Stage 14.9 restart/statistics/visualization/coarse-I/O compatibility gate.
`STAGE14\_11\_RUN\_STAGE14\_10=1` force-reruns the Stage 14.10 RHS/IBM/structure contamination audit.
`STAGE14\_11\_AUTO\_RUN\_MISSING\_PREREQS=1` regenerates missing prerequisite dat evidence automatically when a freshly unpacked repository does not contain generated Stage 14.8/14.9/14.10 output files. This is the default.
By default, the three force-rerun flags are `0`, but missing Stage 14.8/14.9/14.10 dat files are regenerated rather than silently accepted or faked. Set `STAGE14\_11\_AUTO\_RUN\_MISSING\_PREREQS=0` to fail closed immediately if prerequisite dat files are absent.
Required artifacts
A passing run produces:
`stage14\_outputs/stage14\_11\_total\_smoke\_closure.dat`
`stage14\_outputs/stage14\_11\_static\_audit\_report.txt`
`stage14\_outputs/stage14\_11\_total\_smoke\_runtime.log`
`stage14\_outputs/stage14\_11\_stage11\_oneway.dat`
`stage14\_outputs/stage14\_11\_stage12\_feedback\_candidate.dat`
`stage14\_outputs/stage14\_11\_stage13\_force\_density.dat`
`stage14\_outputs/stage14\_11\_stage14\_rhs\_hook.dat`
`stage14\_checks/STAGE14\_CLOSED.md`
The closure file is generated only after the total-smoke verdict is fully PASS.
Pass criteria
The final status requires all of the following:
Static regression audit PASS.
Runtime total-smoke PASS.
Required Stage 11/12/13/14 diagnostics present and internally PASS.
Small-lambda Stage 14 RHS hook active.
Small-lambda RHS increment nonzero, finite, sign-correct, and bounded.
np=1/2/4 Stage 14.8 parallel consistency evidence present and PASS.
Stage 14.9 I/O/restart evidence present and PASS.
Stage 14.10 contamination-audit evidence present and PASS.
No pressure/projection/Poisson/RK3/channel-forcing contamination.
No production IBM forcing outside the approved Stage 13/14 path.
No fibre structure advance or structural solver path.
The wrapper prints:
```text
STAGE 14.11 TOTAL SMOKE VERDICT: PASS
STAGE 14.11 FINAL VERDICT: PASS
STAGE14\_CLOSED=YES
```
on success, or FAIL with explicit reasons.
Manual command
```bash
DECOMP2D\_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \\
BUILD\_DIR=build\_stage9 \\
MPIEXEC=mpirun \\
MPIEXEC\_FLAGS="--mca btl self,vader,tcp" \\
STAGE14\_11\_RUN\_STAGE14\_8=0 \\
STAGE14\_11\_RUN\_STAGE14\_9=0 \\
STAGE14\_11\_RUN\_STAGE14\_10=0 \\
STAGE14\_11\_AUTO\_RUN\_MISSING\_PREREQS=1 \\
STAGE14\_11\_SMALL\_LAMBDA=1.0e-8 \\
STAGE14\_11\_MAX\_RHS\_INCREMENT=1.0e-4 \\
STAGE14\_11\_MAX\_FLUID\_DELTA=1.0e-4 \\
STAGE14\_11\_MAX\_PARALLEL\_RHS\_DIFF=1.0e-12 \\
STAGE14\_11\_MAX\_PARALLEL\_FORCE\_DENSITY\_DIFF=1.0e-10 \\
STAGE14\_11\_MAX\_RESTART\_DELTA=1.0e-8 \\
STAGE14\_11\_MAX\_IO\_SIGNATURE\_DELTA=1.0e-8 \\
STAGE14\_11\_NP=2 \\
STAGE14\_11\_NP\_LIST="1 2 4" \\
bash stage14\_checks/run\_stage14\_11\_total\_smoke\_closure.sh
```
If prior Stage 14.8/14.9/14.10 datfiles are not already present, rerun them explicitly by setting the corresponding `STAGE14\_11\_RUN\_STAGE14\_\*` flag(s) to `1`.
