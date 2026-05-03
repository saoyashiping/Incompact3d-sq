#!/usr/bin/env python3
from pathlib import Path
v = {}
for ln in Path('stage4_outputs/fibre_stage4_main_contamination_check.dat').read_text().splitlines():
    s = ln.split()
    if len(s) >= 2:
        v[s[0]] = float(s[1])

assert v['stage4_main_config_status'] == 1
assert v['stage4_main_apply_ibm_to_fluid_rhs'] == 0
assert v['stage4_main_rhs_disabled_flag'] == 1
assert v['stage4_main_coupling_mode'] == 1

assert v['stage4_main_zero_buffer_rhs_change_max'] <= 1e-14

assert v['stage4_main_nonzero_buffer_max_abs'] > 0
assert v['stage4_main_nonzero_buffer_rhs_change_max'] <= 1e-14

assert v['stage4_main_disabled_rhs_change_max'] <= 1e-14
assert v['stage4_main_disabled_flag'] == 1

assert v['stage4_main_oneway_rhs_change_max'] <= 1e-14
assert v['stage4_main_oneway_rhs_blocked_flag'] == 1

assert v['stage4_main_twoway_rejected_flag'] == 1
assert v['stage4_main_twoway_rhs_change_max'] <= 1e-14

assert v['stage4_main_default_dns_stage4_autocall_flag'] == 0
assert v['stage4_main_default_dns_safe_flag'] == 1

assert v['stage4_main_contamination_status'] == 1
print('STAGE 4.7 MAIN CONTAMINATION CHECK PASSED')
