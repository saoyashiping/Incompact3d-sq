# R1 Baseline No-Fibre Status

## Fibre/FSI disabled?

Runtime confirmation was not possible because the baseline executable was not built and no R1 run was executed.

## `lambda=0` or equivalent disabled setting found?

Not confirmed at runtime. R1 did not modify any original input file. The intended follow-up is to use an R1-only input copy with fibre/FSI disabled or absent once build and MPI prerequisites exist.

## Fibre hook calls observed?

No runtime observation was possible because the run stage was blocked.

## NaN/Inf status

No NaN/Inf runtime status is available because no baseline run was executed.

## Normal termination status

No normal termination status is available because no baseline run was executed.

## Divergence/CFL/pressure/velocity log status

No divergence, CFL, pressure, or velocity runtime logs were produced because no baseline run was executed.

## R1 validation boundary

R1 does not validate FSI. R1 only attempts to validate the baseline no-fibre build/run path. This attempt is blocked by missing build/MPI prerequisites.
