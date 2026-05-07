#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage5_outputs/fibre_stage5_spreading_candidate_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage5_spread_zero_force_buffer_max_abs']<=1e-14
assert v['stage5_spread_zero_rhs_candidate_max_abs']<=1e-14
assert v['stage5_spread_zero_force_conservation_error']<=1e-14
assert v['stage5_spread_nonzero_lag_force_norm']>0
assert v['stage5_spread_nonzero_buffer_max_abs']>0
assert v['stage5_spread_nonzero_rhs_candidate_max_abs']>0
assert v['stage5_spread_force_conservation_abs_error']<=1e-10
assert v['stage5_spread_force_conservation_relative_error']<=1e-10
assert v['stage5_spread_total_eulerian_force_norm']>0
assert v['stage5_spread_total_lagrangian_force_norm']>0
assert v['stage5_spread_power_abs_error']<=1e-10
assert v['stage5_spread_power_relative_error']<=1e-10
assert v['stage5_spread_power_error_consistency_check']<=1e-12
assert v['stage5_spread_power_nonzero_flag']==1
assert v['stage5_rhs_candidate_rho_fluid']>0
assert v['stage5_rhs_candidate_expected_error']<=1e-14
assert v['stage5_rhs_candidate_component_x_error']<=1e-14
assert v['stage5_rhs_candidate_component_y_error']<=1e-14
assert v['stage5_rhs_candidate_component_z_error']<=1e-14
assert v['stage5_spread_rhs_candidate_nonzero_flag']==1
assert v['stage5_spread_real_rhs_change_max']<=1e-14
assert v['stage5_spread_real_rhs_modified_flag']==0
assert v['stage5_spread_nearwall_unsafe_count']>0
assert v['stage5_spread_nearwall_blocked_flag']==1
assert v['stage5_spread_nearwall_spreading_called']==0
assert v['stage5_spread_nearwall_buffer_max_abs']<=1e-14
assert v['stage5_spread_pressure_poisson_modified_flag']==0
assert v['stage5_spreading_candidate_check_status']==1
print('STAGE 5.3 SPREADING CANDIDATE CHECK PASSED')
