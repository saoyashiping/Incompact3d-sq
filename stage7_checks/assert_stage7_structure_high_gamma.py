#!/usr/bin/env python3
from pathlib import Path

vals = {}
for line in Path('stage7_outputs/fibre_stage7_structure_high_gamma_check.dat').read_text().splitlines():
    parts = line.split()
    if len(parts) >= 2:
        vals[parts[0]] = float(parts[1])

expected_one = [
    'stage7_highgamma_real_structure_solver_called_flag',
    'stage7_highgamma_real_tension_solver_called_flag',
    'stage7_highgamma_real_bending_path_called_flag',
    'stage7_highgamma_freefree_boundary_path_called_flag',
    'stage7_highgamma_energy_finite_flag',
    'stage7_highgamma_curvature_finite_flag',
    'stage7_highgamma_momentum_finite_flag',
    'stage7_highgamma_structure_check_status',
]
for key in expected_one:
    assert vals[key] == 1.0, f'{key} expected 1, got {vals[key]}'

expected_zero = ['stage7_highgamma_nan_detected']
for key in expected_zero:
    assert vals[key] == 0.0, f'{key} expected 0, got {vals[key]}'

assert vals['stage7_highgamma_solver_failure_count'] == 0.0
assert vals['stage7_highgamma_length_error_max'] <= 1e-8
assert vals['stage7_highgamma_stretch_error_max'] <= 1e-8

print('STAGE 7 STRUCTURE HIGH-GAMMA CHECK PASSED')
