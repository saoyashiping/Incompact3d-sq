#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage4_outputs/fibre_stage4_poiseuille_response_check.dat').read_text().splitlines():
 p=ln.split();
 if len(p)>=2:v[p[0]]=float(p[1])
assert v['stage4_zero_flow_preservation_error']<=1e-12
# In zero-flow preservation, x and x_old agree up to position roundoff.
# The feedback force uses central-difference velocity history, so O(1e-14)
# position roundoff divided by dt=1e-5 and multiplied by beta_drag=10 can
# produce O(1e-9) residual f_ext. Preservation and length errors remain
# near machine precision; therefore 1e-8 is the appropriate diagnostic
# tolerance for this residual force norm.
assert v['stage4_zero_flow_f_ext_norm']<=1e-8
assert v['stage4_zero_flow_length_error']<=1e-12
assert v['stage4_zero_flow_solver_failure_count']==0
assert v['stage4_zero_flow_nan_detected']==0
assert v['stage4_centerline_final_center_velocity_x']>0
assert v['stage4_centerline_direction_flag']==1
assert v['stage4_centerline_velocity_error']<=1e-4
assert v['stage4_centerline_shape_error_max']<=1e-8
assert v['stage4_centerline_length_error']<=1e-8
assert v['stage4_centerline_bending_energy_final']<=1e-10
assert v['stage4_centerline_solver_failure_count']==0
assert v['stage4_centerline_nan_detected']==0
assert v['stage4_reverse_center_velocity_x']<0
assert v['stage4_reverse_direction_flag']==1
assert v['stage4_reverse_length_error']<=1e-8
assert v['stage4_reverse_nan_detected']==0
assert v['stage4_vertical_force_variation_norm']>0
assert v['stage4_vertical_center_force_greater_flag']==1
assert v['stage4_vertical_shape_response_max']>0
assert v['stage4_vertical_bending_energy_final']>0
assert v['stage4_vertical_length_error']<=1e-8
assert v['stage4_vertical_unsafe_count']==0
assert v['stage4_vertical_solver_failure_count']==0
assert v['stage4_vertical_nan_detected']==0
assert v['stage4_nearwall_unsafe_count']>0
assert v['stage4_nearwall_blocked_flag']==1
assert v['stage4_nearwall_structure_advance_called']==0
assert v['stage4_poiseuille_fluid_rhs_modified']==0
assert v['stage4_poiseuille_response_status']==1
print('STAGE 4.4 POISEUILLE RESPONSE CHECK PASSED')
