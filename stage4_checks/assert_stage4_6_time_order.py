#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage4_outputs/fibre_stage4_time_order_check.dat').read_text().splitlines():
 s=ln.split();
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage4_order_single_step_status']==1
assert v['stage4_order_single_step_f_ext_match_error']<=1e-12
assert v['stage4_order_single_step_advance_called']==1
assert v['stage4_order_single_step_rhs_modified']==0
assert v['stage4_order_multistep_f_ext_refresh_error_max']<=1e-12
assert v['stage4_order_multistep_f_ext_change_norm']>0
assert v['stage4_order_multistep_length_error']<=1e-8
assert v['stage4_order_multistep_nan_detected']==0
assert v['stage4_order_preadvance_force_match_error']<=1e-12
assert v['stage4_order_postadvance_force_staleness_norm']>0
assert v['stage4_substep_force_alpha_order_flag']==1
assert v['stage4_substep_force_distinct_flag']==1
assert v['stage4_substep_f_ext_match_error_max']<=1e-12
assert v['stage4_order_unsafe_count']>0
assert v['stage4_order_unsafe_blocked_flag']==1
assert v['stage4_order_unsafe_advance_called']==0
assert v['stage4_order_unsafe_f_ext_norm']<=1e-14
assert v['stage4_order_velocity_change_max']<=1e-14
assert v['stage4_order_fluid_rhs_modified']==0
assert v['stage4_time_order_status']==1
print('STAGE 4.6 TIME ORDER CHECK PASSED')
