module fibre_stage13_config
  implicit none
  private

  integer, parameter :: default_max_points = 8
  integer, parameter :: default_max_eulerian_points = 64
  integer, parameter :: normalization_conservative = 1
  integer, parameter :: normalization_raw = 2

  logical :: config_loaded = .false.
  logical :: requested = .false.
  logical :: readonly_mode = .true.
  logical :: spreading_readonly_mode = .true.
  integer :: max_points = default_max_points
  integer :: max_eulerian_points = default_max_eulerian_points
  integer :: normalization_mode = normalization_conservative

  integer :: requested_flag = 0
  integer :: readonly_mode_status = 1
  integer :: spreading_readonly_status = 1
  integer :: disabled_by_default_status = 1
  integer :: max_points_status = 1
  integer :: max_eulerian_points_status = 1
  integer :: normalization_mode_status = 1
  integer :: no_force_density_allocation_status = 1
  integer :: no_spreading_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: config_status = 1

  public :: stage13_config_load
  public :: stage13_requested
  public :: stage13_readonly_mode
  public :: stage13_spreading_readonly_mode
  public :: stage13_get_max_points
  public :: stage13_get_max_eulerian_points
  public :: stage13_get_normalization_mode
  public :: stage13_get_status_values

contains

  subroutine stage13_config_load()
    character(len=128) :: value
    integer :: env_status
    integer :: parsed_value
    integer :: read_status

    requested = .false.
    readonly_mode = .true.
    spreading_readonly_mode = .true.
    max_points = default_max_points
    max_eulerian_points = default_max_eulerian_points
    normalization_mode = normalization_conservative
    requested_flag = 0
    readonly_mode_status = 1
    spreading_readonly_status = 1
    disabled_by_default_status = 1
    max_points_status = 1
    max_eulerian_points_status = 1
    normalization_mode_status = 1
    no_force_density_allocation_status = 1
    no_spreading_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1

    call get_environment_variable('X3D_STAGE13_FORCE_DENSITY_CANDIDATE', value=value, status=env_status)
    if (env_status == 0) then
      if (trim(adjustl(value)) == '1') requested = .true.
    end if
    requested_flag = merge(1, 0, requested)

    call get_environment_variable('X3D_STAGE13_FORCE_READONLY', value=value, status=env_status)
    if (env_status == 0) then
      if (trim(adjustl(value)) == '0') readonly_mode_status = 1
    end if
    readonly_mode = .true.
    readonly_mode_status = merge(1, 0, readonly_mode)

    call get_environment_variable('X3D_STAGE13_SPREADING_READONLY', value=value, status=env_status)
    if (env_status == 0) then
      if (trim(adjustl(value)) == '0') spreading_readonly_status = 1
    end if
    spreading_readonly_mode = .true.
    spreading_readonly_status = merge(1, 0, spreading_readonly_mode)

    call get_environment_variable('X3D_STAGE13_MAX_POINTS', value=value, status=env_status)
    if (env_status == 0) then
      read(value, *, iostat=read_status) parsed_value
      if (read_status == 0 .and. parsed_value > 0) then
        max_points = parsed_value
        max_points_status = 1
      else
        max_points = default_max_points
        max_points_status = 0
      end if
    end if

    call get_environment_variable('X3D_STAGE13_MAX_EULERIAN_POINTS', value=value, status=env_status)
    if (env_status == 0) then
      read(value, *, iostat=read_status) parsed_value
      if (read_status == 0 .and. parsed_value > 0) then
        max_eulerian_points = parsed_value
        max_eulerian_points_status = 1
      else
        max_eulerian_points = default_max_eulerian_points
        max_eulerian_points_status = 0
      end if
    end if

    call get_environment_variable('X3D_STAGE13_SPREADING_NORMALIZATION', value=value, status=env_status)
    if (env_status == 0) then
      call parse_normalization(value)
    end if

    call update_config_status()
    config_loaded = .true.
  end subroutine stage13_config_load

  subroutine parse_normalization(value)
    character(len=*), intent(in) :: value
    character(len=128) :: normalized_value

    normalized_value = lowercase(trim(adjustl(value)))
    if (trim(normalized_value) == 'conservative') then
      normalization_mode = normalization_conservative
      normalization_mode_status = 1
    else if (trim(normalized_value) == 'raw') then
      normalization_mode = normalization_raw
      normalization_mode_status = 1
    else
      normalization_mode = normalization_conservative
      normalization_mode_status = 0
    end if
  end subroutine parse_normalization

  function lowercase(value) result(lower_value)
    character(len=*), intent(in) :: value
    character(len=len(value)) :: lower_value
    integer :: i
    integer :: code

    lower_value = value
    do i = 1, len(value)
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower_value(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase

  subroutine update_config_status()
    if (readonly_mode .and. spreading_readonly_mode .and. &
        max_points_status == 1 .and. max_eulerian_points_status == 1 .and. &
        normalization_mode_status == 1 .and. &
        no_force_density_allocation_status == 1 .and. no_spreading_status == 1 .and. &
        no_rhs_injection_status == 1 .and. no_ibm_spreading_status == 1 .and. &
        no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
        no_structure_advance_status == 1 .and. no_fluid_field_access_status == 1 .and. &
        no_fluid_field_modification_status == 1) then
      config_status = 1
    else
      config_status = 0
    end if
  end subroutine update_config_status

  logical function stage13_requested()
    if (.not. config_loaded) call stage13_config_load()
    stage13_requested = requested
  end function stage13_requested

  logical function stage13_readonly_mode()
    if (.not. config_loaded) call stage13_config_load()
    stage13_readonly_mode = readonly_mode
  end function stage13_readonly_mode

  logical function stage13_spreading_readonly_mode()
    if (.not. config_loaded) call stage13_config_load()
    stage13_spreading_readonly_mode = spreading_readonly_mode
  end function stage13_spreading_readonly_mode

  integer function stage13_get_max_points()
    if (.not. config_loaded) call stage13_config_load()
    stage13_get_max_points = max_points
  end function stage13_get_max_points

  integer function stage13_get_max_eulerian_points()
    if (.not. config_loaded) call stage13_config_load()
    stage13_get_max_eulerian_points = max_eulerian_points
  end function stage13_get_max_eulerian_points

  integer function stage13_get_normalization_mode()
    if (.not. config_loaded) call stage13_config_load()
    stage13_get_normalization_mode = normalization_mode
  end function stage13_get_normalization_mode

  subroutine stage13_get_status_values(requested_out, readonly_out, spreading_readonly_out, disabled_default_out, &
                                       max_points_out, max_eulerian_points_out, normalization_out, &
                                       no_force_density_allocation_out, no_spreading_out, no_rhs_injection_out, &
                                       no_ibm_spreading_out, no_feedback_application_out, no_twoway_force_out, &
                                       no_structure_advance_out, no_fluid_field_access_out, &
                                       no_fluid_field_modification_out, config_out)
    integer, intent(out) :: requested_out
    integer, intent(out) :: readonly_out
    integer, intent(out) :: spreading_readonly_out
    integer, intent(out) :: disabled_default_out
    integer, intent(out) :: max_points_out
    integer, intent(out) :: max_eulerian_points_out
    integer, intent(out) :: normalization_out
    integer, intent(out) :: no_force_density_allocation_out
    integer, intent(out) :: no_spreading_out
    integer, intent(out) :: no_rhs_injection_out
    integer, intent(out) :: no_ibm_spreading_out
    integer, intent(out) :: no_feedback_application_out
    integer, intent(out) :: no_twoway_force_out
    integer, intent(out) :: no_structure_advance_out
    integer, intent(out) :: no_fluid_field_access_out
    integer, intent(out) :: no_fluid_field_modification_out
    integer, intent(out) :: config_out

    if (.not. config_loaded) call stage13_config_load()
    call update_config_status()

    requested_out = requested_flag
    readonly_out = readonly_mode_status
    spreading_readonly_out = spreading_readonly_status
    disabled_default_out = disabled_by_default_status
    max_points_out = max_points_status
    max_eulerian_points_out = max_eulerian_points_status
    normalization_out = normalization_mode_status
    no_force_density_allocation_out = no_force_density_allocation_status
    no_spreading_out = no_spreading_status
    no_rhs_injection_out = no_rhs_injection_status
    no_ibm_spreading_out = no_ibm_spreading_status
    no_feedback_application_out = no_feedback_application_status
    no_twoway_force_out = no_twoway_force_status
    no_structure_advance_out = no_structure_advance_status
    no_fluid_field_access_out = no_fluid_field_access_status
    no_fluid_field_modification_out = no_fluid_field_modification_status
    config_out = config_status
  end subroutine stage13_get_status_values

end module fibre_stage13_config
