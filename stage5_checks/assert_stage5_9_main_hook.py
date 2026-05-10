#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage5_outputs/fibre_stage5_main_hook_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert vals['stage5_hook_default_rhs_change_max']<=1e-14
assert int(vals['stage5_hook_default_injected_flag'])==0
assert int(vals['stage5_hook_default_modified_flag'])==0
assert int(vals['stage5_hook_default_safe_flag'])==1
assert vals['stage5_hook_oneway_rhs_change_max']<=1e-14
assert int(vals['stage5_hook_oneway_injected_flag'])==0
assert int(vals['stage5_hook_oneway_modified_flag'])==0
assert vals['stage5_hook_twoway_buffer_max_abs']>0
assert vals['stage5_hook_twoway_rhs_change_max']>0
assert vals['stage5_hook_twoway_expected_error']<=1e-14
assert int(vals['stage5_hook_twoway_injected_flag'])==1
assert int(vals['stage5_hook_twoway_modified_flag'])==1
assert vals['stage5_hook_component_x_error']<=1e-14
assert vals['stage5_hook_component_y_error']<=1e-14
assert vals['stage5_hook_component_z_error']<=1e-14
assert int(vals['stage5_main_hook_check_status'])==1
print('STAGE 5.9 MAIN HOOK CHECK PASSED')
