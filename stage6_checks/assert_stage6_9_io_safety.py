#!/usr/bin/env python3
from pathlib import Path
vals={}
for line in Path('stage6_outputs/fibre_stage6_io_safety_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2:
        vals[p[0]]=float(p[1])

def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]

for k in [
'stage6_io_production_two_way_enabled_by_default','stage6_io_default_rhs_allowed_flag','stage6_io_default_main_hook_enabled_flag','stage6_io_default_config_safe_flag',
'stage6_io_controlled_test_only_flag','stage6_io_controlled_config_valid_flag','stage6_io_controlled_rhs_allowed_flag','stage6_io_controlled_production_enabled_flag',
'stage6_io_invalid_config_rejected_flag','stage6_io_invalid_config_rhs_allowed_flag','stage6_io_production_enabled_rejected_flag','stage6_io_production_enabled_rhs_allowed_flag',
'stage6_io_stage5_closed_marker_exists','stage6_io_stage5_dependency_status','stage6_io_stage6_0_output_exists','stage6_io_stage6_1_output_exists','stage6_io_stage6_2_output_exists','stage6_io_stage6_3_output_exists','stage6_io_stage6_4_output_exists','stage6_io_stage6_5_output_exists','stage6_io_stage6_6_output_exists','stage6_io_stage6_7_output_exists','stage6_io_stage6_8_output_exists','stage6_io_all_prior_outputs_exist','stage6_io_summary_file_exists','stage6_io_summary_file_status','stage6_io_required_keys_present','stage6_io_stage5_outputs_preserved_flag','stage6_io_stage5_closed_marker_preserved_flag','stage6_io_safety_check_status']:
    assert need(k)==1, k

for k in ['stage6_io_pressure_poisson_modified_flag','stage6_io_projection_modified_flag','stage6_io_real_projection_called_flag','stage6_io_production_dns_called_flag']:
    assert need(k)==0, k

summary=Path('stage6_outputs/STAGE6_PROGRESS_SUMMARY.md')
assert summary.exists()
text=summary.read_text()
for key in [
'Stage 6 Progress Summary','Stage 6.0 CONFIG CHECK PASSED','Stage 6.1 RHS AUDIT CHECK PASSED','Stage 6.2 NOOP HOOK CHECK PASSED','Stage 6.3 CONTROLLED RHS CHECK PASSED','Stage 6.4 SINGLEPHASE NOOP CHECK PASSED','Stage 6.5 PROJECTION ORDER CHECK PASSED','Stage 6.6 RK RHS SYNC CHECK PASSED','Stage 6.7 LAYOUT GUARD CHECK PASSED','Stage 6.8 MICRO SMOKE CHECK PASSED','Production two-way enabled by default: false','Controlled test only: true','Pressure Poisson modified: false','Projection modified: false','Real production DNS called: false']:
    assert key in text, key
print('STAGE 6.9 IO SAFETY CHECK PASSED')
