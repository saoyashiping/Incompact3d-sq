#!/usr/bin/env python3
from pathlib import Path
vals={}
for line in Path('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2:
        vals[p[0]]=float(p[1])

def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]

one=[
'stage8_bridge_stage7_closed_marker_exists','stage8_bridge_stage7_total_smoke_output_exists','stage8_bridge_stage7_total_smoke_status','stage8_bridge_stage7_closed_marker_status','stage8_bridge_stage8_0_output_exists','stage8_bridge_stage8_0_status','stage8_bridge_dependency_status',
'stage8_bridge_explicit_array_fallback_called_flag','stage8_bridge_validate_stage7_grid_called_flag','stage8_bridge_explicit_grid_valid_flag','stage8_bridge_explicit_bridge_status',
'stage8_bridge_metadata_match_status','stage8_bridge_real_source_audit_report_exists','stage8_bridge_explicit_array_fallback_available_flag','stage8_bridge_real_source_evidence_status',
'stage8_bridge_x_uniform_flag','stage8_bridge_z_uniform_flag','stage8_bridge_y_nonuniform_flag','stage8_bridge_y_monotonic_flag','stage8_bridge_dy_positive_flag','stage8_bridge_volume_positive_flag','stage8_bridge_wall_bounds_status','stage8_bridge_grid_validation_detail_status',
'stage8_bridge_scalar_smoke_status','stage8_bridge_velocity_smoke_status','stage8_bridge_spreading_smoke_status','stage8_bridge_boundary_smoke_status','stage8_bridge_stage7_interface_smoke_status',
'stage8_bridge_invalid_yface_rejected_flag','stage8_bridge_invalid_dy_rejected_flag','stage8_bridge_invalid_periodic_x_rejected_flag','stage8_bridge_invalid_periodic_z_rejected_flag','stage8_bridge_invalid_input_status',
'stage8_bridge_noop_safety_status','stage8_runtime_grid_bridge_check_status'
]
zero=['stage8_bridge_explicit_grid_rejected_flag','stage8_bridge_noop_rhs_modified_flag','stage8_bridge_pressure_poisson_modified_flag','stage8_bridge_projection_modified_flag','stage8_bridge_real_projection_called_flag','stage8_bridge_production_dns_called_flag','stage8_bridge_fluid_update_called_flag','stage8_bridge_fibre_advance_called_flag']
for k in one: assert need(k)==1.0,k
for k in zero: assert need(k)==0.0,k
rf=need('stage8_bridge_real_xcompact_grid_source_found_flag'); ra=need('stage8_bridge_real_xcompact_grid_adapter_called_flag')
assert (rf==1.0 and ra==1.0) or (rf==0.0 and ra==0.0)
assert abs(need('stage8_bridge_metadata_match_error_max'))<=1e-12
assert abs(need('stage8_bridge_dy_face_consistency_error_max'))<=1e-12
assert abs(need('stage8_bridge_volume_formula_error_max'))<=1e-12
assert abs(need('stage8_bridge_noop_rhs_change_max'))<=1e-14
print('STAGE 8.1 RUNTIME GRID BRIDGE CHECK PASSED')
