# Production Recovery R10 Main Hook Audit

R10 adds a minimal environment-gated hook to `src/xcompact3d.f90`.  The hook initializes after existing Stage 15 setup, applies after transient RHS construction and existing Stage 14 RHS injection, and finalizes before existing shutdown paths.

The hook is default-off.  With no `FIBRE_PROD_*` environment variables, the runtime configuration remains disabled and the RHS adapter returns without writing.

With `FIBRE_PROD_ENABLE=1` and `FIBRE_PROD_LAMBDA=0`, the adapter records before/after signatures and returns without writing arrays.

With `FIBRE_PROD_ENABLE=1` and a small positive lambda, the adapter applies a small finite diagnostic RHS contribution to `dux1(:,:,:,1)` only through the controlled R10 path.

R10 does not modify RK3 coefficients, pressure/projection routines, channel forcing, restart, statistics, visualization, or MPI decomposition logic.  R10 does not claim production DNS-FSI final closure.
