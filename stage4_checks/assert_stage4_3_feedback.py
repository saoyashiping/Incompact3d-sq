#!/usr/bin/env python3
from pathlib import Path
p=Path('stage4_outputs/fibre_stage4_feedback_check.dat')
v={}
for ln in p.read_text().splitlines():
 s=ln.split();
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage4_zero_slip_feedback_status']==1
assert v['stage4_zero_slip_force_structure_norm']<=1e-12
assert v['stage4_zero_slip_force_fluid_norm']<=1e-12
assert v['stage4_zero_slip_f_ext_norm']<=1e-12
assert v['stage4_positive_drag_force_error']<=1e-12
assert v['stage4_positive_drag_mean_fx']>0
assert v['stage4_positive_drag_direction_flag']==1
assert v['stage4_reverse_drag_force_error']<=1e-12
assert v['stage4_reverse_drag_mean_fx']<0
assert v['stage4_reverse_drag_direction_flag']==1
assert v['stage4_action_reaction_pointwise_error']<=1e-12
assert v['stage4_action_reaction_total_error']<=1e-12
assert v['stage4_poiseuille_force_status']==1
assert v['stage4_poiseuille_force_center_greater_flag']==1
assert v['stage4_poiseuille_force_symmetry_error']<=5e-2
assert v['stage4_poiseuille_force_variation_norm']>0
assert v['stage4_f_ext_matches_structure_force_error']<=1e-14
assert v['stage4_structure_advance_called']==0
assert v['stage4_nonuniform_feedback_blocked_flag']==1
assert v['stage4_unknown_layout_feedback_blocked_flag']==1
assert v['stage4_staggered_layout_feedback_blocked_flag']==1
assert v['stage4_feedback_fluid_rhs_modified']==0
assert v['stage4_feedback_status']==1
print('STAGE 4.3 FEEDBACK CHECK PASSED')
