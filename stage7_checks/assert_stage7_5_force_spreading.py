#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_force_spreading_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2: vals[p[0]]=float(p[1])
expected_one=['stage7_spread_stage6_dependency_status','stage7_spread_stage7_0_dependency_status','stage7_spread_stage7_1_dependency_status','stage7_spread_stage7_2_dependency_status','stage7_spread_stage7_3_dependency_status','stage7_spread_stage7_4_dependency_status','stage7_spread_single_force_conservation_status','stage7_spread_multipoint_force_conservation_status','stage7_spread_nonuniform_volume_used_flag','stage7_spread_volume_scaling_status','stage7_spread_periodic_x_wrap_status','stage7_spread_periodic_z_wrap_status','stage7_spread_periodic_status','stage7_spread_nearwall_blocked_flag','stage7_spread_y_outside_low_blocked_flag','stage7_spread_y_outside_high_blocked_flag','stage7_spread_y_outside_status','stage7_spread_force_density_convention_flag','stage7_spread_no_rho_division_flag','stage7_force_spreading_check_status']
for k in expected_one: assert vals[k]==1.0, k
for k in ['stage7_spread_multipoint_unsafe_count','stage7_spread_noop_rhs_modified_flag','stage7_spread_pressure_poisson_modified_flag','stage7_spread_projection_modified_flag','stage7_spread_production_dns_called_flag','stage7_spread_fluid_update_called_flag','stage7_spread_fibre_advance_called_flag']: assert vals[k]==0.0,k
assert vals['stage7_spread_single_force_abs_error']<=1e-12
assert vals['stage7_spread_single_force_relative_error']<=1e-12
assert vals['stage7_spread_multipoint_valid_count']>=5
assert vals['stage7_spread_multipoint_force_abs_error']<=1e-12
assert vals['stage7_spread_multipoint_force_relative_error']<=1e-12
assert vals['stage7_spread_volume_scaling_error_max']<=1e-12
assert vals['stage7_spread_periodic_force_error_max']<=1e-12
assert vals['stage7_spread_nearwall_unsafe_count']>0
assert vals['stage7_spread_nearwall_buffer_change_max']<=1e-14
assert vals['stage7_spread_noop_rhs_change_max']<=1e-14
print('STAGE 7.5 FORCE SPREADING CHECK PASSED')
