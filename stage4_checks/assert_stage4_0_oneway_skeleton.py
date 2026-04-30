#!/usr/bin/env python3
from pathlib import Path

path = Path('stage4_outputs/fibre_stage4_oneway_skeleton_check.dat')
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
    key, raw = parts[0], parts[1]
    values[key] = float(raw)


def expect_eq(key: str, expected: float, tol: float = 0.0) -> None:
    if key not in values:
        raise AssertionError(f'Missing key: {key}')
    val = values[key]
    if abs(val - expected) > tol:
        raise AssertionError(f'{key} expected {expected} got {val}')


def expect_le(key: str, upper: float) -> None:
    if key not in values:
        raise AssertionError(f'Missing key: {key}')
    if values[key] > upper:
        raise AssertionError(f'{key} expected <= {upper} got {values[key]}')


def expect_gt(key: str, lower: float) -> None:
    if key not in values:
        raise AssertionError(f'Missing key: {key}')
    if values[key] <= lower:
        raise AssertionError(f'{key} expected > {lower} got {values[key]}')

expect_eq('stage4_config_status', 1)
expect_eq('stage4_enable_oneway', 1)
expect_eq('stage4_enable_ibm_diagnostics', 1)
expect_eq('stage4_apply_ibm_to_fluid_rhs', 0)
expect_eq('stage4_rhs_disabled_flag', 1)
expect_eq('stage4_coupling_mode', 1)
expect_gt('stage4_beta_drag', 0)
expect_gt('stage4_structure_dt', 0)

expect_eq('stage4_grid_nx', 16)
expect_eq('stage4_grid_ny', 12)
expect_eq('stage4_grid_nz', 10)
expect_eq('stage4_grid_periodic_x', 1)
expect_eq('stage4_grid_periodic_y', 0)
expect_eq('stage4_grid_periodic_z', 1)

expect_eq('stage4_fibre_nl', 33)
expect_eq('stage4_lag_nl', 33)
expect_eq('stage4_lag_weight_sum', values['stage4_expected_fibre_length'], tol=1e-12)

expect_eq('stage4_force_buffer_allocated', 1)
expect_le('stage4_force_buffer_initial_max_abs', 1e-14)
expect_le('stage4_force_buffer_initial_l2_norm', 1e-14)

expect_eq('stage4_interpolation_called', 0)
expect_eq('stage4_feedback_called', 0)
expect_eq('stage4_structure_advance_called', 0)
expect_eq('stage4_spreading_called', 0)
expect_eq('stage4_fluid_rhs_modified', 0)
expect_eq('stage4_skeleton_status', 1)

print('STAGE 4.0 ONEWAY SKELETON CHECK PASSED')
