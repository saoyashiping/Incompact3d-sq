#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2: vals[p[0]]=float(p[1])
expected_one=['stage7_interp_robust_stage6_dependency_status','stage7_interp_robust_stage7_0_dependency_status','stage7_interp_robust_stage7_1_dependency_status','stage7_interp_robust_stage7_2_dependency_status','stage7_interp_robust_repeatability_status','stage7_interp_robust_weight_sum_status','stage7_interp_robust_quadratic_y_status','stage7_interp_robust_cubic_y_status','stage7_interp_robust_mixed_field_status','stage7_interp_robust_periodic_x_shift_status','stage7_interp_robust_periodic_z_shift_status','stage7_interp_robust_y_outside_low_blocked_flag','stage7_interp_robust_y_outside_high_blocked_flag','stage7_interp_robust_y_outside_status','stage7_interp_robust_nearwall_blocked_flag','stage7_scalar_interpolation_robustness_check_status']
for k in expected_one: assert vals[k]==1.0,f'{k} expected1 got {vals[k]}'
expected_zero=['stage7_interp_robust_repeat_index_mismatch_count','stage7_interp_robust_noop_rhs_modified_flag','stage7_interp_robust_pressure_poisson_modified_flag','stage7_interp_robust_projection_modified_flag','stage7_interp_robust_production_dns_called_flag','stage7_interp_robust_fluid_update_called_flag','stage7_interp_robust_fibre_advance_called_flag']
for k in expected_zero: assert vals[k]==0.0,f'{k} expected0 got {vals[k]}'
assert vals['stage7_interp_robust_repeat_weight_error_max']<=1e-14
assert vals['stage7_interp_robust_repeat_weight_sum_error_max']<=1e-14
assert vals['stage7_interp_robust_weight_sum_error_max']<=1e-12
assert vals['stage7_interp_robust_quadratic_y_error_max']<=1e-12
assert vals['stage7_interp_robust_cubic_y_error_max']<=1e-11
assert vals['stage7_interp_robust_mixed_field_error_max']<=1e-11
assert vals['stage7_interp_robust_periodic_shift_error_max']<=1e-12
assert vals['stage7_interp_robust_noop_rhs_change_max']<=1e-14
assert vals['stage7_interp_robust_valid_weight_count']>=5
assert vals['stage7_interp_robust_nearwall_unsafe_count']>0
print('STAGE 7.3 SCALAR INTERPOLATION ROBUSTNESS CHECK PASSED')
