# Production Recovery R2 Plan

## R2 scope

R2 establishes a standalone production fibre state container foundation. The R2 modules provide configuration, allocation, initialization, reset, destruction, finite-state checks, segment-length residual diagnostics, and force-zeroing diagnostics.

## R2 non-goals

R2 does not connect fibre state to `xcompact3d.f90`, does not enter the production time loop, does not implement IBM interpolation or spreading, does not couple to RHS, does not advance structure positions, and does not implement wall contact or fibre-fibre collision.

## Build strategy

Add only one standalone CMake target, `fibre_prod_state_check`, built from the R2 modules and independent check driver. The target must not be added to the `xcompact3d` executable source list.

## Run strategy

If build succeeds, run the standalone `fibre_prod_state_check` executable directly. It must print `R2_FIBRE_PROD_STATE_CHECK PASS` on success and fail with non-zero stop on any check failure.

## Evidence boundary

R2 PASS only means the production fibre state container passes standalone validation. R2 PASS does not validate IBM, FSI, structure advancement, RHS coupling, wall contact, fibre-fibre collision, MPI consistency, DNS execution, or production-loop integration.
