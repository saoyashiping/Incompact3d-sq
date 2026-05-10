#!/usr/bin/env python3
from pathlib import Path

v = {}
for ln in Path('stage5_outputs/fibre_stage5_rhs_audit_check.dat').read_text().splitlines():
    s = ln.split()
    if len(s) >= 2:
        v[s[0]] = float(s[1])

assert v['stage5_rhs_convention_status'] == 1
assert v['stage5_rhs_uses_acceleration_form'] == 1
assert v['stage5_rhs_requires_divide_by_rho'] == 1
assert v['stage5_rhs_rho_positive_flag'] == 1

assert v['stage5_rhs_location_identified'] == 1
assert v['stage5_rhs_insert_before_projection_flag'] == 1
assert v['stage5_rhs_insert_after_projection_flag'] == 0
assert v['stage5_rhs_pressure_poisson_modified_flag'] == 0

assert v['stage5_rhs_rk_substep_policy_status'] == 1
assert v['stage5_rhs_current_substep_velocity_required'] == 1
assert v['stage5_rhs_stale_velocity_forbidden_flag'] == 1
assert v['stage5_rhs_substep_force_recompute_required'] == 1

assert v['stage5_rhs_audit_noop_rhs_change_max'] <= 1e-14
assert v['stage5_rhs_audit_rhs_modified_flag'] == 0

assert v['stage5_rhs_twoway_config_rhs_allowed_flag'] == 1
assert v['stage5_rhs_twoway_config_rhs_modified_in_5_1_flag'] == 0

assert v['stage5_rhs_audit_report_exists'] == 1
assert v['stage5_rhs_audit_status'] == 1
assert v['stage5_rhs_audit_check_status'] == 1

report = Path('stage5_outputs/stage5_rhs_audit_report.md')
assert report.exists()
text = report.read_text()
for key in [
    'acceleration form',
    'f_ibm / rho_fluid',
    'before pressure projection',
    'not directly modified',
    'current substep velocity',
    'audit only',
]:
    assert key in text

print('STAGE 5.1 RHS AUDIT CHECK PASSED')
