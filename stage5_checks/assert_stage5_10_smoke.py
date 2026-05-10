#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage5_outputs/fibre_stage5_smoke_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert int(vals['stage5_smoke_config_status'])==1
assert int(vals['stage5_smoke_two_way_enabled_flag'])==1
assert vals['stage5_smoke_completed_steps']==vals['stage5_smoke_step_count']
assert vals['stage5_smoke_hook_expected_error_max']<=1e-14
assert int(vals['stage5_smoke_check_status'])==1
print('STAGE 5.10 SMOKE CHECK PASSED')
