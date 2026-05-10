#!/usr/bin/env python3
from pathlib import Path

v = {}
for ln in Path('stage4p_outputs/stage4p_benchmark_summary.dat').read_text().splitlines():
    s = ln.split()
    if len(s) >= 2:
        v[s[0]] = float(s[1])

assert v['stage4p_config_status'] == 1
assert v['stage4p_grid_nx'] > 0 and v['stage4p_grid_ny'] > 0 and v['stage4p_grid_nz'] > 0
assert v['stage4p_fibre_nl'] > 0 and v['stage4p_nsteps'] > 0 and v['stage4p_dt'] > 0 and v['stage4p_beta_drag'] > 0

assert v['stage4p_rank0_output_flag'] == 1
assert v['stage4p_mpi_size'] >= 1
assert v['stage4p_mpi_rank'] == 0

assert v['stage4p_final_center_displacement_norm'] > 0
assert v['stage4p_final_center_velocity_norm'] > 0
assert v['stage4p_max_f_ext_norm'] > 0
assert v['stage4p_max_slip_norm'] > 0
assert v['stage4p_max_length_error'] <= 1e-8
assert v['stage4p_max_unsafe_count'] == 0
assert v['stage4p_nan_detected'] == 0
assert v['stage4p_solver_failure_count'] == 0

assert v['stage4p_force_conservation_relative_error_max'] <= 1e-10
assert v['stage4p_power_relative_error_max'] <= 1e-10
assert v['stage4p_power_error_consistency_check_max'] <= 1e-12
assert v['stage4p_power_nonzero_flag'] == 1

assert v['stage4p_velocity_change_max'] <= 1e-14
assert v['stage4p_fluid_rhs_modified'] == 0
assert v['stage4p_apply_ibm_to_fluid_rhs'] == 0
assert v['stage4p_rhs_disabled_flag'] == 1

assert v['stage4p_time_history_file_exists'] == 1
assert v['stage4p_time_history_line_count'] >= v['stage4p_nsteps'] + 1

assert v['stage4p_benchmark_status'] == 1
print('STAGE 4P.0 BENCHMARK CHECK PASSED')
