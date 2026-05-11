#!/usr/bin/env python3
from pathlib import Path
v={}
for l in Path('stage7_outputs/fibre_stage7_rho_convention_audit.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2:v[p[0]]=float(p[1])
assert v['stage7_rho_force_buffer_independent_of_rho_flag']==1
assert v['stage7_rho_rhs_divides_once_flag']==1
assert v['stage7_rho_scaling_error']<=1e-14
assert v['stage7_rho_double_division_detected_flag']==0
assert v['stage7_rho_invalid_rho_rejected_flag']==1
assert v['stage7_rho_convention_audit_status']==1
print('STAGE 7 RHO CONVENTION AUDIT PASSED')
