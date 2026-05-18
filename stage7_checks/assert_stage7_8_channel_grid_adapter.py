#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_channel_grid_adapter_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2: vals[p[0]]=float(p[1])
for k in ['stage7_adapter_stage6_dependency_status','stage7_adapter_stage7_0_dependency_status','stage7_adapter_stage7_1_dependency_status','stage7_adapter_stage7_2_dependency_status','stage7_adapter_stage7_3_dependency_status','stage7_adapter_stage7_4_dependency_status','stage7_adapter_stage7_5_dependency_status','stage7_adapter_stage7_6_dependency_status','stage7_adapter_stage7_7_dependency_status','stage7_adapter_from_arrays_valid_flag','stage7_adapter_validate_status','stage7_adapter_metadata_match_status','stage7_adapter_dy_face_consistency_status','stage7_adapter_volume_formula_status','stage7_adapter_total_volume_status','stage7_adapter_periodic_x_flag','stage7_adapter_periodic_z_flag','stage7_adapter_wall_bounds_status','stage7_adapter_invalid_yface_rejected_flag','stage7_adapter_invalid_dy_rejected_flag','stage7_adapter_invalid_periodic_rejected_flag','stage7_adapter_invalid_input_status','stage7_adapter_scalar_smoke_status','stage7_adapter_velocity_smoke_status','stage7_adapter_spreading_smoke_status','stage7_adapter_boundary_smoke_status','stage7_adapter_explicit_array_adapter_called_flag','stage7_channel_grid_adapter_check_status']:
 assert vals[k]==1.0,k
for k in ['stage7_adapter_from_arrays_rejected_flag','stage7_adapter_periodic_y_flag','stage7_adapter_noop_rhs_modified_flag','stage7_adapter_pressure_poisson_modified_flag','stage7_adapter_projection_modified_flag','stage7_adapter_production_dns_called_flag','stage7_adapter_fluid_update_called_flag','stage7_adapter_fibre_advance_called_flag']:
 assert vals[k]==0.0,k
assert vals['stage7_adapter_metadata_match_error_max']<=1e-12
assert vals['stage7_adapter_dy_face_consistency_error_max']<=1e-12
assert vals['stage7_adapter_volume_formula_error_max']<=1e-12
assert vals['stage7_adapter_total_volume_error']<=1e-12
assert vals['stage7_adapter_noop_rhs_change_max']<=1e-14
print('STAGE 7.8 CHANNEL GRID ADAPTER CHECK PASSED')
