#!/usr/bin/env python3
from pathlib import Path
v={}
for l in Path('stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2:v[p[0]]=float(p[1])
assert v['stage5_closed_loop_actual_force_conservation_relative_error']<=1e-12
assert v['stage5_closed_loop_actual_power_relative_error']<=1e-12
assert v['stage5_closed_loop_no_tautological_force_flag']==1
assert v['stage5_closed_loop_no_tautological_power_flag']==1
assert v['stage5_closed_loop_diagnostics_audit_status']==1
print('STAGE 5 CLOSED LOOP DIAGNOSTICS AUDIT PASSED')
