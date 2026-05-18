#!/usr/bin/env python3
from pathlib import Path
vals = {}
for line in Path('stage7_outputs/fibre_stage7_total_smoke_check.dat').read_text().splitlines():
    p = line.split()
    if len(p) >= 2:
        vals[p[0]] = float(p[1])

def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]

expected_one = [
  'stage7_total_stage6_closed_marker_exists','stage7_total_stage6_total_smoke_output_exists','stage7_total_stage6_total_smoke_status','stage7_total_stage6_all_prior_outputs_status','stage7_total_stage6_dependency_status',
  'stage7_total_stage7_0_output_exists','stage7_total_stage7_1_output_exists','stage7_total_stage7_2_output_exists','stage7_total_stage7_3_output_exists','stage7_total_stage7_4_output_exists','stage7_total_stage7_5_output_exists','stage7_total_stage7_6_output_exists','stage7_total_stage7_7_output_exists','stage7_total_stage7_8_output_exists','stage7_total_stage7_9_output_exists','stage7_total_all_stage7_outputs_exist',
  'stage7_total_stage7_0_status','stage7_total_stage7_1_status','stage7_total_stage7_2_status','stage7_total_stage7_3_status','stage7_total_stage7_4_status','stage7_total_stage7_5_status','stage7_total_stage7_6_status','stage7_total_stage7_7_status','stage7_total_stage7_8_status','stage7_total_stage7_9_status','stage7_total_all_stage7_status',
  'stage7_total_default_production_disabled_status','stage7_total_controlled_rhs_path_status','stage7_total_force_density_convention_status','stage7_total_no_double_rho_division_status','stage7_total_production_rejection_status','stage7_total_safety_summary_status',
  'stage7_total_pressure_poisson_untouched_status','stage7_total_projection_untouched_status','stage7_total_production_dns_not_called_status','stage7_total_fluid_update_not_called_status','stage7_total_fibre_advance_not_called_status','stage7_total_noop_safety_summary_status',
  'stage7_total_grid_numeric_status','stage7_total_interpolation_numeric_status','stage7_total_velocity_numeric_status','stage7_total_spreading_numeric_status','stage7_total_power_numeric_status','stage7_total_adapter_numeric_status','stage7_total_rhs_numeric_status','stage7_total_numeric_summary_status',
  'stage7_total_closed_marker_written_flag','stage7_total_closed_marker_status','stage7_total_smoke_check_status'
]
for k in expected_one:
    assert need(k) == 1, k

assert Path('stage7_outputs/STAGE7_CLOSED.md').exists(), 'missing stage7_outputs/STAGE7_CLOSED.md'
print('STAGE 7.10 TOTAL SMOKE CHECK PASSED')
