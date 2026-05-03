#!/usr/bin/env python3
from pathlib import Path

path = Path('stage4_outputs/fibre_stage4_grid_adapter_check.dat')
if not path.exists():
    raise SystemExit(f'Missing output file: {path}')

values = {}
for line in path.read_text(encoding='utf-8').splitlines():
    line = line.strip()
    if not line:
        continue
    parts = line.split()
    if len(parts) < 2:
        continue
    values[parts[0]] = float(parts[1])


def expect_eq(key: str, expected: float, tol: float = 0.0) -> None:
    if key not in values:
        raise AssertionError(f'Missing key: {key}')
    if abs(values[key] - expected) > tol:
        raise AssertionError(f'{key} expected {expected}, got {values[key]}')


def expect_gt(key: str, threshold: float) -> None:
    if key not in values:
        raise AssertionError(f'Missing key: {key}')
    if values[key] <= threshold:
        raise AssertionError(f'{key} expected > {threshold}, got {values[key]}')

expect_eq('uniform_adapter_nx', 16)
expect_eq('uniform_adapter_ny', 12)
expect_eq('uniform_adapter_nz', 10)
expect_eq('uniform_adapter_periodic_x', 1)
expect_eq('uniform_adapter_periodic_y', 0)
expect_eq('uniform_adapter_periodic_z', 1)
expect_eq('uniform_adapter_uniform_x', 1)
expect_eq('uniform_adapter_uniform_y', 1)
expect_eq('uniform_adapter_uniform_z', 1)
expect_eq('uniform_adapter_uniform_ibm_compatible', 1)
expect_eq('uniform_adapter_velocity_layout_mode', 1)

expect_eq('nonuniform_y_adapter_uniform_x', 1)
expect_eq('nonuniform_y_adapter_uniform_y', 0)
expect_eq('nonuniform_y_adapter_uniform_z', 1)
expect_eq('nonuniform_y_adapter_uniform_ibm_compatible', 0)
expect_gt('nonuniform_y_adapter_dy_max', values['nonuniform_y_adapter_dy_min'])

expect_eq('layout_unknown_mode', 0)
expect_eq('layout_collocated_mode', 1)
expect_eq('layout_staggered_mode', 2)
expect_eq('layout_staggered_requires_component_coordinates', 1)

expect_eq('adapter_stage4_apply_ibm_to_fluid_rhs', 0)
expect_eq('adapter_stage4_rhs_disabled_flag', 1)

expect_eq('uniform_adapter_to_ibm_grid_possible', 1)
expect_eq('nonuniform_adapter_to_ibm_grid_possible', 0)

expect_eq('stage4_grid_adapter_status', 1)

print('STAGE 4.1 GRID ADAPTER CHECK PASSED')
