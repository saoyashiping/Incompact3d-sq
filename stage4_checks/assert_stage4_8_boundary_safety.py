#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage4_outputs/fibre_stage4_boundary_safety_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage4_boundary_safe_interior_count']==33
assert v['stage4_boundary_safe_periodic_wrap_count']==0
assert v['stage4_boundary_safe_unsafe_count']==0
assert v['stage4_boundary_safe_outside_count']==0
assert v['stage4_boundary_safe_blocked_flag']==0
assert v['stage4_boundary_periodic_wrap_count']>0
assert v['stage4_boundary_periodic_unsafe_count']==0
assert v['stage4_boundary_periodic_outside_count']==0
assert v['stage4_boundary_periodic_blocked_flag']==0
assert v['stage4_boundary_nearwall_unsafe_count']>0
assert v['stage4_boundary_nearwall_blocked_flag']==1
assert v['stage4_boundary_nearwall_interpolation_called']==0
assert v['stage4_boundary_nearwall_feedback_called']==0
assert v['stage4_boundary_nearwall_structure_advance_called']==0
assert v['stage4_boundary_nearwall_f_ext_norm']<=1e-14
assert v['stage4_boundary_outside_count']>0
assert v['stage4_boundary_outside_blocked_flag']==1
assert v['stage4_boundary_outside_structure_advance_called']==0
assert v['stage4_boundary_nonuniform_y_compatible']==0
assert v['stage4_boundary_nonuniform_y_blocked_flag']==1
assert v['stage4_boundary_unknown_layout_blocked_flag']==1
assert v['stage4_boundary_staggered_layout_blocked_flag']==1
assert v['stage4_boundary_fluid_rhs_modified']==0
assert v['stage4_boundary_safety_status']==1
print('STAGE 4.8 BOUNDARY SAFETY CHECK PASSED')
