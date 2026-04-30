#!/usr/bin/env python3
from pathlib import Path

path = Path('stage4_outputs/fibre_stage4_interpolation_check.dat')
if not path.exists():
    raise SystemExit(f'Missing output file: {path}')

values = {}
for line in path.read_text(encoding='utf-8').splitlines():
    parts = line.strip().split()
    if len(parts) >= 2:
        values[parts[0]] = float(parts[1])


def eq(k, v, tol=0.0):
    if abs(values[k] - v) > tol:
        raise AssertionError(f'{k} expected {v}, got {values[k]}')


def le(k, v):
    if values[k] > v:
        raise AssertionError(f'{k} expected <= {v}, got {values[k]}')


eq('stage4_const_interp_status', 1)
le('stage4_const_interp_max_error', 1e-12)
le('stage4_const_weight_sum_max_error', 1e-12)

eq('stage4_linear_interp_status', 1)
le('stage4_linear_interp_max_error', 1e-11)

eq('stage4_poiseuille_interp_status', 1)
le('stage4_poiseuille_interp_max_error', 5e-2)
le('stage4_poiseuille_center_value_error', 5e-2)
le('stage4_poiseuille_symmetry_error', 5e-2)
eq('stage4_poiseuille_monotonicity_flag', 1)

eq('stage4_periodic_interp_status', 1)
le('stage4_periodic_interp_error_max', 0.15)
le('stage4_periodic_weight_sum_max_error', 1e-12)

eq('stage4_nonuniform_y_uniform_ibm_compatible', 0)
eq('stage4_nonuniform_y_blocked_flag', 1)

eq('stage4_unknown_layout_blocked_flag', 1)
eq('stage4_staggered_layout_blocked_flag', 1)

eq('stage4_interp_fluid_rhs_modified', 0)
eq('stage4_interpolation_status', 1)

print('STAGE 4.2 INTERPOLATION CHECK PASSED')
