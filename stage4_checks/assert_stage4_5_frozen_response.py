#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage4_outputs/fibre_stage4_frozen_response_check.dat').read_text().splitlines():
 s=ln.split();
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage4_frozen_velocity_change_max']<=1e-14
assert v['stage4_frozen_centerline_final_center_velocity_norm']>0
assert v['stage4_frozen_centerline_center_displacement_norm']>0
assert v['stage4_frozen_centerline_f_ext_norm']>0
assert v['stage4_frozen_centerline_length_error']<=1e-8
assert v['stage4_frozen_centerline_unsafe_count']==0
assert v['stage4_frozen_centerline_nan_detected']==0
assert v['stage4_frozen_u_lag_change_norm']>0
assert v['stage4_frozen_f_ext_refresh_error_max']<=1e-12
assert v['stage4_frozen_force_conservation_error']<=1e-10
assert v['stage4_frozen_force_conservation_relative_error']<=1e-10
assert v['stage4_frozen_force_buffer_max_abs']>0
assert v['stage4_frozen_power_abs_error']<=1e-10
assert v['stage4_frozen_power_relative_error']<=1e-10
assert v['stage4_frozen_nearwall_unsafe_count']>0
assert v['stage4_frozen_nearwall_blocked_flag']==1
assert v['stage4_frozen_nearwall_structure_advance_called']==0
assert v['stage4_frozen_fluid_rhs_modified']==0
assert v['stage4_frozen_response_status']==1
print('STAGE 4.5 FROZEN RESPONSE CHECK PASSED')
