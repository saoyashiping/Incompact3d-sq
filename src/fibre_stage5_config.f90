module fibre_stage5_config
  use fibre_parameters, only : mytype
  implicit none

  integer, parameter :: STAGE5_COUPLING_DISABLED = 0
  integer, parameter :: STAGE5_COUPLING_ONE_WAY  = 1
  integer, parameter :: STAGE5_COUPLING_TWO_WAY  = 2

  type stage5_config_t
    logical :: enable_stage5
    integer :: coupling_mode
    logical :: apply_ibm_to_fluid_rhs
    logical :: allow_two_way
    logical :: reject_invalid_config
    real(mytype) :: rho_fluid
    integer :: config_status
  end type stage5_config_t

contains

  subroutine init_stage5_default_config(config)
    type(stage5_config_t), intent(out) :: config
    config%enable_stage5 = .false.
    config%coupling_mode = STAGE5_COUPLING_DISABLED
    config%apply_ibm_to_fluid_rhs = .false.
    config%allow_two_way = .false.
    config%reject_invalid_config = .true.
    config%rho_fluid = 1.0_mytype
    config%config_status = 1
  end subroutine init_stage5_default_config

  subroutine init_stage5_oneway_config(config)
    type(stage5_config_t), intent(out) :: config
    config%enable_stage5 = .true.
    config%coupling_mode = STAGE5_COUPLING_ONE_WAY
    config%apply_ibm_to_fluid_rhs = .false.
    config%allow_two_way = .false.
    config%reject_invalid_config = .true.
    config%rho_fluid = 1.0_mytype
    config%config_status = 1
  end subroutine init_stage5_oneway_config

  subroutine init_stage5_twoway_config(config)
    type(stage5_config_t), intent(out) :: config
    config%enable_stage5 = .true.
    config%coupling_mode = STAGE5_COUPLING_TWO_WAY
    config%apply_ibm_to_fluid_rhs = .true.
    config%allow_two_way = .true.
    config%reject_invalid_config = .true.
    config%rho_fluid = 1.0_mytype
    config%config_status = 1
  end subroutine init_stage5_twoway_config

  subroutine validate_stage5_config(config, valid, rejected_flag, two_way_enabled_flag, rhs_allowed_flag)
    type(stage5_config_t), intent(in) :: config
    integer, intent(out) :: valid, rejected_flag, two_way_enabled_flag, rhs_allowed_flag

    valid = 1
    rejected_flag = 0
    two_way_enabled_flag = 0
    rhs_allowed_flag = 0

    if (config%rho_fluid <= 0._mytype) then
      valid = 0
    else if ((.not. config%enable_stage5) .and. config%coupling_mode==STAGE5_COUPLING_DISABLED .and. (.not. config%apply_ibm_to_fluid_rhs)) then
      valid = 1
    else if (config%enable_stage5 .and. config%coupling_mode==STAGE5_COUPLING_ONE_WAY .and. (.not. config%apply_ibm_to_fluid_rhs)) then
      valid = 1
    else if (config%enable_stage5 .and. config%coupling_mode==STAGE5_COUPLING_TWO_WAY .and. config%apply_ibm_to_fluid_rhs .and. config%allow_two_way) then
      valid = 1
      two_way_enabled_flag = 1
      rhs_allowed_flag = 1
    else
      valid = 0
    end if

    if (valid /= 1) then
      rejected_flag = 1
      two_way_enabled_flag = 0
      rhs_allowed_flag = 0
    end if
  end subroutine validate_stage5_config

  subroutine stage5_noop_rhs_guard(config, rhs_modified_flag)
    type(stage5_config_t), intent(in) :: config
    integer, intent(out) :: rhs_modified_flag
    integer :: valid, rejected_flag, two_way_enabled_flag, rhs_allowed_flag
    call validate_stage5_config(config, valid, rejected_flag, two_way_enabled_flag, rhs_allowed_flag)
    rhs_modified_flag = 0
  end subroutine stage5_noop_rhs_guard

end module fibre_stage5_config
