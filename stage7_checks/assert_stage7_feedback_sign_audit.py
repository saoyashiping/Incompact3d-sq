#!/usr/bin/env python3
from pathlib import Path
v={}
for l in Path('stage7_outputs/fibre_stage7_feedback_sign_audit.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2:v[p[0]]=float(p[1])
assert v['stage7_feedback_zero_slip_force_norm']<=1e-14
assert v['stage7_feedback_action_reaction_error']<=1e-14
assert v['stage7_feedback_structure_force_slip_dot_positive_flag']==1
assert v['stage7_feedback_fluid_force_slip_dot_negative_flag']==1
assert v['stage7_feedback_total_power_dissipative_flag']==1
assert v['stage7_feedback_total_power_error']<=1e-14
assert v['stage7_feedback_sign_audit_status']==1
print('STAGE 7 FEEDBACK SIGN AUDIT PASSED')
