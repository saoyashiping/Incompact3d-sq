#!/usr/bin/env python3
from pathlib import Path

v = {}
for ln in Path('stage5_outputs/fibre_stage5_config_check.dat').read_text().splitlines():
    s = ln.split()
    if len(s) >= 2:
        v[s[0]] = float(s[1])

assert v['stage5_default_config_status'] == 1
assert v['stage5_default_enable_stage5'] == 0
assert v['stage5_default_coupling_mode'] == 0
assert v['stage5_default_apply_ibm_to_fluid_rhs'] == 0
assert v['stage5_default_two_way_enabled_flag'] == 0
assert v['stage5_default_rhs_allowed_flag'] == 0
assert v['stage5_default_rejected_flag'] == 0
assert v['stage5_default_rhs_modified_flag'] == 0

assert v['stage5_oneway_config_status'] == 1
assert v['stage5_oneway_enable_stage5'] == 1
assert v['stage5_oneway_coupling_mode'] == 1
assert v['stage5_oneway_apply_ibm_to_fluid_rhs'] == 0
assert v['stage5_oneway_two_way_enabled_flag'] == 0
assert v['stage5_oneway_rhs_allowed_flag'] == 0
assert v['stage5_oneway_rejected_flag'] == 0
assert v['stage5_oneway_rhs_modified_flag'] == 0

assert v['stage5_twoway_config_status'] == 1
assert v['stage5_twoway_enable_stage5'] == 1
assert v['stage5_twoway_coupling_mode'] == 2
assert v['stage5_twoway_apply_ibm_to_fluid_rhs'] == 1
assert v['stage5_twoway_two_way_enabled_flag'] == 1
assert v['stage5_twoway_rhs_allowed_flag'] == 1
assert v['stage5_twoway_rejected_flag'] == 0
assert v['stage5_twoway_rhs_modified_flag'] == 0

assert v['stage5_invalid_oneway_rhs_rejected_flag'] == 1
assert v['stage5_invalid_oneway_rhs_valid_flag'] == 0
assert v['stage5_invalid_oneway_rhs_allowed_flag'] == 0

assert v['stage5_invalid_twoway_no_rhs_rejected_flag'] == 1
assert v['stage5_invalid_twoway_no_rhs_valid_flag'] == 0
assert v['stage5_invalid_twoway_no_rhs_allowed_flag'] == 0

assert v['stage5_invalid_twoway_not_allowed_rejected_flag'] == 1
assert v['stage5_invalid_twoway_not_allowed_valid_flag'] == 0

assert v['stage5_invalid_rho_rejected_flag'] == 1
assert v['stage5_invalid_rho_valid_flag'] == 0

assert v['stage5_stage4p_frozen_marker_status'] == 1
assert v['stage5_config_check_status'] == 1
assert Path('stage4p_checks/STAGE4P_FROZEN.md').exists()

print('STAGE 5.0 CONFIG CHECK PASSED')
