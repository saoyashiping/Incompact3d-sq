#!/usr/bin/env python3
from pathlib import Path

v = {}
for ln in Path('stage5_outputs/fibre_stage5_rhs_wrapper_check.dat').read_text().splitlines():
    s = ln.split()
    if len(s) >= 2:
        v[s[0]] = float(s[1])

assert v['stage5_rhs_disabled_change_max'] <= 1e-14
assert v['stage5_rhs_disabled_modified_flag'] == 0
assert v['stage5_rhs_disabled_injected_flag'] == 0

assert v['stage5_rhs_oneway_change_max'] <= 1e-14
assert v['stage5_rhs_oneway_modified_flag'] == 0
assert v['stage5_rhs_oneway_injected_flag'] == 0

assert v['stage5_rhs_twoway_zero_buffer_change_max'] <= 1e-14
assert v['stage5_rhs_twoway_zero_buffer_modified_flag'] == 0
assert v['stage5_rhs_twoway_zero_buffer_injected_flag'] == 1

assert v['stage5_rhs_twoway_nonzero_buffer_max_abs'] > 0
assert v['stage5_rhs_twoway_nonzero_change_max'] > 0
assert v['stage5_rhs_twoway_nonzero_expected_error'] <= 1e-14
assert v['stage5_rhs_twoway_nonzero_modified_flag'] == 1
assert v['stage5_rhs_twoway_nonzero_injected_flag'] == 1

assert v['stage5_rhs_invalid_config_rejected_flag'] == 1
assert v['stage5_rhs_invalid_config_change_max'] <= 1e-14
assert v['stage5_rhs_invalid_config_injected_flag'] == 0

assert v['stage5_rhs_invalid_rho_rejected_flag'] == 1
assert v['stage5_rhs_invalid_rho_change_max'] <= 1e-14
assert v['stage5_rhs_invalid_rho_injected_flag'] == 0

assert v['stage5_rhs_component_x_error'] <= 1e-14
assert v['stage5_rhs_component_y_error'] <= 1e-14
assert v['stage5_rhs_component_z_error'] <= 1e-14

assert v['stage5_rhs_pressure_poisson_modified_flag'] == 0
assert v['stage5_rhs_wrapper_check_status'] == 1

print('STAGE 5.2 RHS WRAPPER CHECK PASSED')
