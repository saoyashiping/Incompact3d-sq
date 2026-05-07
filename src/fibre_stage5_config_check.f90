program fibre_stage5_config_check
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_DISABLED, STAGE5_COUPLING_ONE_WAY, STAGE5_COUPLING_TWO_WAY, &
                                   init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config, &
                                   validate_stage5_config, stage5_noop_rhs_guard
  implicit none

  type(stage5_config_t) :: cfg
  integer :: valid, rejected, two_way, rhs_allowed, rhs_modified
  integer :: marker_status, check_status, ios

  marker_status = 0
  check_status = 1

  open(11,file='stage5_outputs/fibre_stage5_config_check.dat',status='replace',iostat=ios)
  if (ios /= 0) stop 1

  call init_stage5_default_config(cfg)
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  call stage5_noop_rhs_guard(cfg, rhs_modified)
  write(11,'(A,1X,I0)') 'stage5_default_config_status', cfg%config_status
  write(11,'(A,1X,I0)') 'stage5_default_enable_stage5', merge(1,0,cfg%enable_stage5)
  write(11,'(A,1X,I0)') 'stage5_default_coupling_mode', cfg%coupling_mode
  write(11,'(A,1X,I0)') 'stage5_default_apply_ibm_to_fluid_rhs', merge(1,0,cfg%apply_ibm_to_fluid_rhs)
  write(11,'(A,1X,I0)') 'stage5_default_two_way_enabled_flag', two_way
  write(11,'(A,1X,I0)') 'stage5_default_rhs_allowed_flag', rhs_allowed
  write(11,'(A,1X,I0)') 'stage5_default_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_default_rhs_modified_flag', rhs_modified

  call init_stage5_oneway_config(cfg)
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  call stage5_noop_rhs_guard(cfg, rhs_modified)
  write(11,'(A,1X,I0)') 'stage5_oneway_config_status', cfg%config_status
  write(11,'(A,1X,I0)') 'stage5_oneway_enable_stage5', merge(1,0,cfg%enable_stage5)
  write(11,'(A,1X,I0)') 'stage5_oneway_coupling_mode', cfg%coupling_mode
  write(11,'(A,1X,I0)') 'stage5_oneway_apply_ibm_to_fluid_rhs', merge(1,0,cfg%apply_ibm_to_fluid_rhs)
  write(11,'(A,1X,I0)') 'stage5_oneway_two_way_enabled_flag', two_way
  write(11,'(A,1X,I0)') 'stage5_oneway_rhs_allowed_flag', rhs_allowed
  write(11,'(A,1X,I0)') 'stage5_oneway_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_oneway_rhs_modified_flag', rhs_modified

  call init_stage5_twoway_config(cfg)
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  call stage5_noop_rhs_guard(cfg, rhs_modified)
  write(11,'(A,1X,I0)') 'stage5_twoway_config_status', cfg%config_status
  write(11,'(A,1X,I0)') 'stage5_twoway_enable_stage5', merge(1,0,cfg%enable_stage5)
  write(11,'(A,1X,I0)') 'stage5_twoway_coupling_mode', cfg%coupling_mode
  write(11,'(A,1X,I0)') 'stage5_twoway_apply_ibm_to_fluid_rhs', merge(1,0,cfg%apply_ibm_to_fluid_rhs)
  write(11,'(A,1X,I0)') 'stage5_twoway_two_way_enabled_flag', two_way
  write(11,'(A,1X,I0)') 'stage5_twoway_rhs_allowed_flag', rhs_allowed
  write(11,'(A,1X,I0)') 'stage5_twoway_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_twoway_rhs_modified_flag', rhs_modified

  call init_stage5_default_config(cfg)
  cfg%enable_stage5 = .true.; cfg%coupling_mode = STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs = .true.
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  write(11,'(A,1X,I0)') 'stage5_invalid_oneway_rhs_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_invalid_oneway_rhs_valid_flag', valid
  write(11,'(A,1X,I0)') 'stage5_invalid_oneway_rhs_allowed_flag', rhs_allowed

  call init_stage5_default_config(cfg)
  cfg%enable_stage5 = .true.; cfg%coupling_mode = STAGE5_COUPLING_TWO_WAY; cfg%apply_ibm_to_fluid_rhs = .false.; cfg%allow_two_way = .true.
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  write(11,'(A,1X,I0)') 'stage5_invalid_twoway_no_rhs_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_invalid_twoway_no_rhs_valid_flag', valid
  write(11,'(A,1X,I0)') 'stage5_invalid_twoway_no_rhs_allowed_flag', rhs_allowed

  call init_stage5_default_config(cfg)
  cfg%enable_stage5 = .true.; cfg%coupling_mode = STAGE5_COUPLING_TWO_WAY; cfg%apply_ibm_to_fluid_rhs = .true.; cfg%allow_two_way = .false.
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  write(11,'(A,1X,I0)') 'stage5_invalid_twoway_not_allowed_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_invalid_twoway_not_allowed_valid_flag', valid

  call init_stage5_default_config(cfg)
  cfg%rho_fluid = 0._mytype
  call validate_stage5_config(cfg, valid, rejected, two_way, rhs_allowed)
  write(11,'(A,1X,I0)') 'stage5_invalid_rho_rejected_flag', rejected
  write(11,'(A,1X,I0)') 'stage5_invalid_rho_valid_flag', valid

  inquire(file='stage4p_checks/STAGE4P_FROZEN.md', exist=cfg%allow_two_way)
  marker_status = merge(1,0,cfg%allow_two_way)
  write(11,'(A,1X,I0)') 'stage5_stage4p_frozen_marker_status', marker_status

  if (marker_status /= 1) check_status = 0
  write(11,'(A,1X,I0)') 'stage5_config_check_status', check_status
  close(11)
end program fibre_stage5_config_check
