#!/usr/bin/env python3
from pathlib import Path

p=Path('stage9_outputs/fibre_stage9_1a_decomp_io_api_audit.dat')
if not p.is_file():
    raise SystemExit('missing audit dat')
vals={}
for ln in p.read_text().splitlines():
    if '=' in ln:
        k,v=ln.split('=',1)
        vals[k.strip()]=v.strip()

expected_one = [
  'stage9_1a_stage9_0_output_exists','stage9_1a_stage9_0_status','stage9_1a_stage9_1_allowed_flag','stage9_1a_stage9_0_dependency_status',
  'stage9_1a_decomp2d_root_exists','stage9_1a_decomp2d_root_nonempty_flag','stage9_1a_real_decomp2d_required_flag','stage9_1a_no_decomp2d_stub_policy_flag','stage9_1a_decomp2d_presence_status',
  'stage9_1a_compat_module_exists','stage9_1a_compat_in_cmake_flag','stage9_1a_compat_before_forces_flag','stage9_1a_compat_before_statistics_flag','stage9_1a_compat_before_visu_flag','stage9_1a_compat_before_tools_flag','stage9_1a_compat_before_les_models_flag','stage9_1a_compat_before_case_files_flag','stage9_1a_compat_cmake_order_status',
  'stage9_1a_forbidden_direct_old_io_imports_absent_flag','stage9_1a_xcompact3d_real_io_init_finalise_allowed_flag','stage9_1a_direct_import_policy_status',
  'stage9_1a_compat_metadata_wrappers_status','stage9_1a_compat_write_one_wrapper_status','stage9_1a_compat_read_one_wrapper_status','stage9_1a_compat_write_plane_wrapper_status','stage9_1a_compat_fine_to_coarse_wrapper_status','stage9_1a_compat_gen_iodir_name_status','stage9_1a_compat_coverage_status',
  'stage9_1a_legacy_io_call_coverage_flag','stage9_1a_legacy_call_coverage_status',
  'stage9_1a_build_only_wrapper_marker_flag','stage9_1a_stage9_10_output_semantics_marker_flag','stage9_1a_wrapper_semantics_warning_status',
  'stage9_1a_stage9_1b_allowed_flag','stage9_1a_next_step_policy_status','stage9_1a_decomp_io_api_audit_status'
]
expected_zero = ['stage9_1a_production_coupling_allowed_flag','stage9_1a_forbidden_direct_old_io_import_count','stage9_1a_legacy_io_uncovered_file_count','stage9_1a_stage9_2_allowed_flag']

bad=[]
for k in expected_one:
    if vals.get(k)!='1': bad.append((k,vals.get(k)))
for k in expected_zero:
    if vals.get(k)!='0': bad.append((k,vals.get(k)))
if not Path('stage9_outputs/STAGE9_1A_DECOMP_IO_API_AUDIT.md').is_file():
    bad.append(('stage9_1a_report_exists','0'))
if bad:
    raise SystemExit('assert failed: '+str(bad))
print('STAGE 9.1A DECOMP IO API AUDIT PASSED')
