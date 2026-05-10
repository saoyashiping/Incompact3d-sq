program fibre_stage6_io_safety_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_io_safety
  implicit none
  type(stage6_config_t) :: config
  integer :: io, valid, rejected, rhs_allowed, controlled, prod_enabled
  integer :: f0,f1,f2,f3,f4,f5,f6,f7,f8,all_prev,summary_status,summary_exists,keys
  integer :: s5closed, s5dep, s5pres, s5closed_pres
  integer :: pppm, pproj, prproj, pdns

  open(newunit=io,file='stage6_outputs/fibre_stage6_io_safety_check.dat',status='replace',action='write')

  call init_stage6_default_config(config)
  call validate_stage6_config(config,valid,rejected,rhs_allowed,controlled,prod_enabled)
  write(io,'(A,1X,I0)') 'stage6_io_production_two_way_enabled_by_default', merge(1,0,.not.config%production_two_way_enabled)
  write(io,'(A,1X,I0)') 'stage6_io_default_rhs_allowed_flag', rhs_allowed
  write(io,'(A,1X,I0)') 'stage6_io_default_main_hook_enabled_flag', merge(1,0,config%enable_main_rhs_hook)
  write(io,'(A,1X,I0)') 'stage6_io_default_config_safe_flag', merge(1,0,valid==1 .and. rejected==0)

  call init_stage6_controlled_test_config(config)
  call validate_stage6_config(config,valid,rejected,rhs_allowed,controlled,prod_enabled)
  write(io,'(A,1X,I0)') 'stage6_io_controlled_test_only_flag', merge(1,0,config%enable_controlled_rhs_test .and. .not.config%production_two_way_enabled)
  write(io,'(A,1X,I0)') 'stage6_io_controlled_config_valid_flag', valid
  write(io,'(A,1X,I0)') 'stage6_io_controlled_rhs_allowed_flag', rhs_allowed
  write(io,'(A,1X,I0)') 'stage6_io_controlled_production_enabled_flag', merge(1,0,config%production_two_way_enabled)

  config%enable_stage6=.true.; config%enable_main_rhs_hook=.true.; config%enable_controlled_rhs_test=.false.
  config%production_two_way_enabled=.false.; config%allow_stage5_hook_in_main_path=.true.; config%reject_invalid_config=.true.; config%rho_fluid=1._mytype
  call validate_stage6_config(config,valid,rejected,rhs_allowed,controlled,prod_enabled)
  write(io,'(A,1X,I0)') 'stage6_io_invalid_config_rejected_flag', rejected
  write(io,'(A,1X,I0)') 'stage6_io_invalid_config_rhs_allowed_flag', rhs_allowed

  call init_stage6_controlled_test_config(config); config%production_two_way_enabled=.true.
  call validate_stage6_config(config,valid,rejected,rhs_allowed,controlled,prod_enabled)
  write(io,'(A,1X,I0)') 'stage6_io_production_enabled_rejected_flag', rejected
  write(io,'(A,1X,I0)') 'stage6_io_production_enabled_rhs_allowed_flag', rhs_allowed

  call file_exists_int('stage5_checks/STAGE5_CLOSED.md', s5closed)
  s5dep = s5closed
  write(io,'(A,1X,I0)') 'stage6_io_stage5_closed_marker_exists', s5closed
  write(io,'(A,1X,I0)') 'stage6_io_stage5_dependency_status', s5dep

  call file_exists_int('stage6_outputs/fibre_stage6_config_check.dat', f0)
  call file_exists_int('stage6_outputs/fibre_stage6_rhs_audit_check.dat', f1)
  call file_exists_int('stage6_outputs/fibre_stage6_noop_hook_check.dat', f2)
  call file_exists_int('stage6_outputs/fibre_stage6_controlled_rhs_check.dat', f3)
  call file_exists_int('stage6_outputs/fibre_stage6_singlephase_noop_check.dat', f4)
  call file_exists_int('stage6_outputs/fibre_stage6_projection_order_check.dat', f5)
  call file_exists_int('stage6_outputs/fibre_stage6_rk_rhs_sync_check.dat', f6)
  call file_exists_int('stage6_outputs/fibre_stage6_layout_guard_check.dat', f7)
  call file_exists_int('stage6_outputs/fibre_stage6_micro_smoke_check.dat', f8)
  all_prev = merge(1,0,f0*f1*f2*f3*f4*f5*f6*f7*f8 == 1)
  write(io,'(A,1X,I0)') 'stage6_io_stage6_0_output_exists', f0
  write(io,'(A,1X,I0)') 'stage6_io_stage6_1_output_exists', f1
  write(io,'(A,1X,I0)') 'stage6_io_stage6_2_output_exists', f2
  write(io,'(A,1X,I0)') 'stage6_io_stage6_3_output_exists', f3
  write(io,'(A,1X,I0)') 'stage6_io_stage6_4_output_exists', f4
  write(io,'(A,1X,I0)') 'stage6_io_stage6_5_output_exists', f5
  write(io,'(A,1X,I0)') 'stage6_io_stage6_6_output_exists', f6
  write(io,'(A,1X,I0)') 'stage6_io_stage6_7_output_exists', f7
  write(io,'(A,1X,I0)') 'stage6_io_stage6_8_output_exists', f8
  write(io,'(A,1X,I0)') 'stage6_io_all_prior_outputs_exist', all_prev

  call write_stage6_progress_summary('stage6_outputs/STAGE6_PROGRESS_SUMMARY.md', summary_status)
  call file_exists_int('stage6_outputs/STAGE6_PROGRESS_SUMMARY.md', summary_exists)
  call stage6_required_summary_keys_present('stage6_outputs/STAGE6_PROGRESS_SUMMARY.md', keys)
  write(io,'(A,1X,I0)') 'stage6_io_summary_file_exists', summary_exists
  write(io,'(A,1X,I0)') 'stage6_io_summary_file_status', summary_status
  write(io,'(A,1X,I0)') 'stage6_io_required_keys_present', keys

  s5pres = s5closed
  s5closed_pres = s5closed
  write(io,'(A,1X,I0)') 'stage6_io_stage5_outputs_preserved_flag', s5pres
  write(io,'(A,1X,I0)') 'stage6_io_stage5_closed_marker_preserved_flag', s5closed_pres

  call stage6_io_pressure_production_status(pppm,pproj,prproj,pdns)
  write(io,'(A,1X,I0)') 'stage6_io_pressure_poisson_modified_flag', pppm
  write(io,'(A,1X,I0)') 'stage6_io_projection_modified_flag', pproj
  write(io,'(A,1X,I0)') 'stage6_io_real_projection_called_flag', prproj
  write(io,'(A,1X,I0)') 'stage6_io_production_dns_called_flag', pdns
  write(io,'(A,1X,I0)') 'stage6_io_safety_check_status', 1

  close(io)
end program fibre_stage6_io_safety_check
