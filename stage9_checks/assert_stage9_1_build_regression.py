#!/usr/bin/env python3
from pathlib import Path

p = Path('stage9_outputs/fibre_stage9_build_regression_check.dat')
if not p.is_file():
    raise SystemExit('missing stage9 build regression dat')

vals = {}
for ln in p.read_text().splitlines():
    if not ln.strip():
        continue
    parts = ln.split()
    if len(parts) >= 2:
        vals[parts[0].strip()] = parts[1].strip()

expected_one = [
  'stage9_build_stage9_0_output_exists','stage9_build_stage9_0_status','stage9_build_stage9_entry_allowed_flag','stage9_build_stage9_1_allowed_flag','stage9_build_stage9_0_dependency_status',
  'stage9_build_stage7_closed_marker_exists','stage9_build_stage8_closed_marker_exists','stage9_build_stage8_total_smoke_output_exists','stage9_build_stage8_total_smoke_status','stage9_build_stage8_closed_marker_status','stage9_build_stage8_no_production_side_effect_status','stage9_build_prior_stage_closure_status',
  'stage9_build_evidence_output_exists','stage9_build_configure_status','stage9_build_xcompact3d_compile_status','stage9_build_used_source_matched_decomp_flag','stage9_build_decomp_root_exists','stage9_build_decomp_mod_exists','stage9_build_decomp_io_mod_exists','stage9_build_decomp_static_lib_exists','stage9_build_no_compat_remnants_flag','stage9_build_fibre_stage_checks_only_off_flag','stage9_build_default_production_disabled_status','stage9_build_evidence_status',
  'stage9_build_default_production_disabled_status','stage9_build_production_dns_not_called_status','stage9_build_production_coupling_disabled_status','stage9_build_rhs_untouched_status','stage9_build_projection_untouched_status','stage9_build_fluid_update_not_called_status','stage9_build_fibre_advance_not_called_status','stage9_build_production_disabled_policy_status',
  'stage9_build_noop_safety_status','stage9_build_regression_check_status'
]

expected_zero = [
  'stage9_build_production_coupling_allowed_flag','stage9_build_stage_checks_only_flag','stage9_build_dns_executed_flag','stage9_build_rhs_modified_flag','stage9_build_projection_called_flag','stage9_build_fluid_update_called_flag','stage9_build_fibre_advance_called_flag','stage9_build_rhs_hook_called_flag','stage9_build_pressure_poisson_modified_flag','stage9_build_real_projection_called_flag','stage9_build_production_dns_called_flag'
]

bad = []
for k in expected_one:
    if vals.get(k) != '1':
        bad.append((k, vals.get(k)))
for k in expected_zero:
    if vals.get(k) != '0':
        bad.append((k, vals.get(k)))

if bad:
    raise SystemExit('assert failed: ' + str(bad))

print('STAGE 9.1 BUILD REGRESSION CHECK PASSED')
