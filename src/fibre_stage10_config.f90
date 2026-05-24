module fibre_stage10_config
  implicit none
  private

  logical, save :: loaded=.false.
  logical, save :: requested_flag=.false.
  logical, save :: noop_mode_flag=.true.
  integer, save :: max_steps_value=3

  public :: stage10_config_load
  public :: stage10_requested
  public :: stage10_noop_mode
  public :: stage10_get_max_steps
  public :: stage10_get_status_values

contains

  subroutine stage10_config_load()
    character(len=64) :: value
    integer :: status, ios, parsed_steps

    requested_flag=.false.
    noop_mode_flag=.true.
    max_steps_value=3

    value=''
    call get_environment_variable('X3D_STAGE10_PRODUCTION_HOOK', value=value, status=status)
    if (status==0) then
      if (trim(adjustl(value))=='1') requested_flag=.true.
    endif

    value=''
    call get_environment_variable('X3D_STAGE10_FORCE_NOOP', value=value, status=status)
    if (status==0) then
      if (trim(adjustl(value))=='1') noop_mode_flag=.true.
    endif

    value=''
    call get_environment_variable('X3D_STAGE10_MAX_STEPS', value=value, status=status)
    if (status==0) then
      read(value,*,iostat=ios) parsed_steps
      if (ios==0 .and. parsed_steps>0) max_steps_value=parsed_steps
    endif

    if (requested_flag) then
      noop_mode_flag=.true.
    endif

    loaded=.true.
  end subroutine stage10_config_load

  logical function stage10_requested()
    if (.not.loaded) call stage10_config_load()
    stage10_requested=requested_flag
  end function stage10_requested

  logical function stage10_noop_mode()
    if (.not.loaded) call stage10_config_load()
    stage10_noop_mode=noop_mode_flag
  end function stage10_noop_mode

  integer function stage10_get_max_steps()
    if (.not.loaded) call stage10_config_load()
    stage10_get_max_steps=max_steps_value
  end function stage10_get_max_steps

  subroutine stage10_get_status_values(requested, noop_mode_status, disabled_by_default_status, no_fibre_state_status, no_force_status, no_rhs_injection_status)
    integer, intent(out) :: requested
    integer, intent(out) :: noop_mode_status
    integer, intent(out) :: disabled_by_default_status
    integer, intent(out) :: no_fibre_state_status
    integer, intent(out) :: no_force_status
    integer, intent(out) :: no_rhs_injection_status

    if (.not.loaded) call stage10_config_load()

    if (requested_flag) then
      requested=1
    else
      requested=0
    endif

    if (noop_mode_flag) then
      noop_mode_status=1
    else
      noop_mode_status=0
    endif

    if ((.not.requested_flag) .and. noop_mode_flag) then
      disabled_by_default_status=1
    else
      if (requested_flag .and. noop_mode_flag) then
        disabled_by_default_status=1
      else
        disabled_by_default_status=0
      endif
    endif

    no_fibre_state_status=1
    no_force_status=1
    no_rhs_injection_status=1
  end subroutine stage10_get_status_values

end module fibre_stage10_config
