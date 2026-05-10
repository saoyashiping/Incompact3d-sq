#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage5_outputs/fibre_stage5_momentum_exchange_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage5_momentum_zero_fluid_impulse_norm']<=1e-14
assert v['stage5_momentum_zero_structure_impulse_norm']<=1e-14
assert v['stage5_momentum_zero_exchange_error']<=1e-14
assert v['stage5_momentum_nonzero_force_buffer_max_abs']>0
assert v['stage5_momentum_nonzero_rhs_change_max']>0
assert v['stage5_momentum_fluid_velocity_change_norm']>0
assert v['stage5_momentum_rhs_injected_flag']==1
assert v['stage5_momentum_rhs_modified_flag']==1
assert v['stage5_momentum_action_reaction_error']<=1e-12
assert v['stage5_momentum_fluid_impulse_from_velocity_norm']>0
assert v['stage5_momentum_eulerian_force_impulse_norm']>0
assert v['stage5_momentum_fluid_impulse_error']<=1e-10
assert v['stage5_momentum_fluid_impulse_relative_error']<=1e-10
assert v['stage5_momentum_structure_impulse_norm']>0
assert v['stage5_momentum_exchange_abs_error']<=1e-10
assert v['stage5_momentum_exchange_relative_error']<=1e-10
assert v['stage5_momentum_component_x_error']<=1e-10
assert v['stage5_momentum_component_y_error']<=1e-10
assert v['stage5_momentum_component_z_error']<=1e-10
assert v['stage5_momentum_disabled_velocity_change_norm']<=1e-14
assert v['stage5_momentum_disabled_injected_flag']==0
assert v['stage5_momentum_oneway_velocity_change_norm']<=1e-14
assert v['stage5_momentum_oneway_injected_flag']==0
assert v['stage5_momentum_invalid_rejected_flag']==1
assert v['stage5_momentum_invalid_velocity_change_norm']<=1e-14
assert v['stage5_momentum_invalid_injected_flag']==0
assert v['stage5_momentum_pressure_poisson_modified_flag']==0
assert v['stage5_momentum_main_dns_hooked_flag']==0
assert v['stage5_momentum_synthetic_only_flag']==1
assert v['stage5_momentum_exchange_check_status']==1
print('STAGE 5.5 MOMENTUM EXCHANGE CHECK PASSED')
