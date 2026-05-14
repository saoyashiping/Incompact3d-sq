#!/usr/bin/env python3
from pathlib import Path

vals={}
for line in Path('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2: vals[p[0]]=float(p[1])

expected_one=[
'stage7_interp_stage6_dependency_status','stage7_interp_stage7_0_dependency_status','stage7_interp_stage7_1_dependency_status',
'stage7_interp_weight_sum_status','stage7_interp_constant_status','stage7_interp_linear_status','stage7_interp_poiseuille_status',
'stage7_interp_periodic_x_wrap_status','stage7_interp_periodic_z_wrap_status','stage7_interp_nearwall_blocked_flag','stage7_scalar_interpolation_check_status']
for k in expected_one: assert vals[k]==1.0, f'{k} expected 1 got {vals[k]}'

expected_zero=['stage7_interp_noop_rhs_modified_flag','stage7_interp_pressure_poisson_modified_flag','stage7_interp_projection_modified_flag','stage7_interp_production_dns_called_flag','stage7_interp_fluid_update_called_flag','stage7_interp_fibre_advance_called_flag']
for k in expected_zero: assert vals[k]==0.0, f'{k} expected 0 got {vals[k]}'

assert vals['stage7_interp_weight_sum_error_max']<=1e-12
assert vals['stage7_interp_constant_error_max']<=1e-12
assert vals['stage7_interp_linear_error_max']<=1e-11
assert vals['stage7_interp_poiseuille_error_max']<=1e-11
assert vals['stage7_interp_periodic_wrap_error_max']<=1e-12
assert vals['stage7_interp_noop_rhs_change_max']<=1e-14
assert vals['stage7_interp_valid_weight_count']>=5
assert vals['stage7_interp_nearwall_unsafe_count']>0
print('STAGE 7.2 SCALAR INTERPOLATION CHECK PASSED')
