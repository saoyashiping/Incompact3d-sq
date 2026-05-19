#!/usr/bin/env python3
from pathlib import Path

vals = {}
for line in Path('stage7_outputs/fibre_stage7_feedback_sign_audit.dat').read_text().splitlines():
    parts = line.split()
    if len(parts) >= 2:
        vals[parts[0]] = float(parts[1])

expected_one = [
    'stage7_feedback_real_module_called_flag',
    'stage7_feedback_force_on_structure_from_module_flag',
    'stage7_feedback_force_on_fluid_from_module_flag',
    'stage7_feedback_zero_slip_flag',
    'stage7_feedback_action_reaction_flag',
    'stage7_feedback_structure_force_slip_dot_positive_flag',
    'stage7_feedback_fluid_force_slip_dot_negative_flag',
    'stage7_feedback_total_power_dissipative_flag',
    'stage7_feedback_sign_audit_status',
]
for key in expected_one:
    assert vals[key] == 1.0, f'{key} expected 1, got {vals[key]}'

assert vals['stage7_feedback_zero_slip_force_norm'] <= 1e-14
assert vals['stage7_feedback_action_reaction_error'] <= 1e-14
assert vals['stage7_feedback_structure_force_slip_dot'] > 0.0
assert vals['stage7_feedback_fluid_force_slip_dot'] < 0.0
assert vals['stage7_feedback_total_power'] <= 0.0
assert vals['stage7_feedback_total_power_error'] <= 1e-14

print('STAGE 7 FEEDBACK SIGN AUDIT PASSED')
