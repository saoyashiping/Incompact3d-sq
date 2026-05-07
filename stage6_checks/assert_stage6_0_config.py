#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_config_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert Path('stage5_checks/STAGE5_CLOSED.md').exists()
assert int(vals['stage6_stage5_closed_marker_status'])==1
assert int(vals['stage6_stage5_dependency_status'])==1
assert int(vals['stage6_default_valid_flag'])==1 and int(vals['stage6_default_rhs_allowed_flag'])==0
assert int(vals['stage6_controlled_valid_flag'])==1 and int(vals['stage6_controlled_rhs_allowed_flag'])==1
assert int(vals['stage6_invalid_production_enabled_rejected_flag'])==1
assert int(vals['stage6_invalid_hook_without_controlled_test_rejected_flag'])==1
assert int(vals['stage6_invalid_controlled_without_allow_hook_rejected_flag'])==1
assert int(vals['stage6_invalid_rho_rejected_flag'])==1
assert vals['stage6_config_noop_rhs_change_max']<=1e-14
assert int(vals['stage6_config_noop_rhs_modified_flag'])==0
assert int(vals['stage6_config_check_status'])==1
print('STAGE 6.0 CONFIG CHECK PASSED')
