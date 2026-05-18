# STAGE8_1 Runtime Grid Source Audit

## Search keywords
`yp|yyp|y_face|yface|yc|ycenter|dy|dy_cell|dx|dz|istret|stretch|mesh|grid|yly|xlx|zlz|xnu|nclx|ncly|nclz|p_row|p_col`

## Candidate files inspected
- `src/parameters.f90`
- `src/variables.f90`
- `src/navier.f90`
- `src/stretching.f90`
- `src/Case-Periodic-hill.f90`
- `src/Case-Sandbox.f90`

## Candidate runtime-grid variables found
- Domain and boundary-control parameters in `parameters.f90`: `nx, ny, nz, xlx, yly, zlz, nclx1/nclxn/ncly1/nclyn/nclz1/nclzn, istret`.
- Grid spacing and runtime arrays referenced from `var`/`variables`/cases: `dx, dy, dz`, and stretched-y arrays `yp` used by case routines.
- Stretching generator in `stretching.f90` (`mod_stret::stretching`) computes `yp`/`ypi` and derivative helper arrays.

## Can these be used directly as a stable Stage-8.1 runtime bridge source?
Not yet as a clean, stable check-target API. The runtime data path is spread across simulation modules (`param`, `variables`, `var`, case setup and decomposition state), and direct use would over-couple the Stage-8 check target to full solver runtime initialization.

## Face/center coordinate availability
- `y` center-like array evidence: `yp` is clearly present/used.
- Explicit module-level `y_face` array with simple standalone access was not identified from this audit path.

## Why explicit-array fallback is used in Stage 8.1
Stage 8.1 uses explicit-array fallback to keep the bridge audit deterministic and decoupled from DNS runtime bootstrapping while still reusing the Stage-7 channel-grid adapter + validator path.

## What to add before Stage-9 DNS regression
- A dedicated exported runtime-grid metadata API (single module) that provides:
  - `nx, ny, nz`
  - `xmin/xmax, zmin/zmax`
  - `y_face(:), y_center(:)`
  - periodic flags in x/z
- Clear initialization ordering contract so check targets can call this API safely without entering full DNS advance.
