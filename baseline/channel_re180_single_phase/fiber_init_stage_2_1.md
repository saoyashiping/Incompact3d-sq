# Stage 2.1: Minimal single-fiber static geometry infrastructure

## Scope implemented in this stage
This stage adds the minimum infrastructure for a **single straight fiber** as a static geometry object:

1. Fiber data structure definition;
2. Input parsing for minimal fiber parameters;
3. One-time initialization of Lagrangian points on a straight centerline;
4. One-time text output of fiber coordinates (`fiber_points.dat`).

## Intentionally not implemented in this stage
The following are intentionally excluded:

- Any fiber dynamics (translation/rotation/time integration);
- Any IBM interpolation/spreading/forcing;
- Any coupling between fiber and flow fields;
- Any collision or wall-contact model;
- Any flexible beam/tension/stiffness model;
- Multi-fiber support.

The fiber is a static geometry object and does not participate in governing equation solves.

## Added files
- `src/fiber_types.f90`
- `src/fiber_init.f90`
- `src/fiber_io.f90`
- `baseline/channel_re180_single_phase/input_fiber_init_test.i3d`
- `baseline/channel_re180_single_phase/fiber_init_stage_2_1.md`

## Minimally modified existing files
- `src/parameters.f90` (optional `FiberParam` namelist parsing + default initialization hook)
- `src/xcompact3d.f90` (initialization-stage call to fiber init and one-time output)
- `src/CMakeLists.txt` (compile new fiber modules)

## How to run minimal checks
### 1) Baseline path unchanged (fiber off by default)
Use original baseline input:

```bash
./xcompact3d baseline/channel_re180_single_phase/input.i3d
```

Expected: no fiber output file is generated and baseline run path remains unchanged.

### 2) Fiber initialization test
Use dedicated test input:

```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_init_test.i3d
```

Expected: `fiber_points.dat` is generated once during initialization.

## How to validate `fiber_points.dat`
Check:

1. Number of points equals `fiber_nl`;
2. Coordinates are uniformly spaced along the specified direction;
3. Endpoints are symmetric around `fiber_center`;
4. For the provided test (`fiber_direction = (1,0,0)`), only `x` varies linearly while `y,z` stay at center values.

Example command:

```bash
head -n 12 fiber_points.dat
```
