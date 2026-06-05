module fibre_stage16_one_fibre_case
  use decomp_2d_constants, only : mytype
  use fibre_stage16_config, only : stage16_config_load_from_environment, stage16_diagnostic_only_requested, &
       stage16_config_get_status_values
  implicit none

  private

  public :: stage16_one_fibre_case_reset
  public :: stage16_one_fibre_case_define_default
  public :: stage16_one_fibre_case_load_from_environment
  public :: stage16_one_fibre_case_validate
  public :: stage16_one_fibre_case_get_npts
  public :: stage16_one_fibre_case_get_position
  public :: stage16_one_fibre_case_get_velocity
  public :: stage16_one_fibre_case_get_acceleration
  public :: stage16_one_fibre_case_write_diagnostics
  public :: stage16_one_fibre_case_set_for_test
  public :: stage16_one_fibre_case_get_status_values

  integer, parameter :: default_npts = 8
  integer, parameter :: ndim = 3

  logical :: case_loaded = .false.
  logical :: stage16_2_enabled = .true.
  logical :: one_fibre_case_requested = .true.
  logical :: diagnostic_only_value = .true.
  logical :: wall_contact_requested = .false.
  integer :: npts_value = default_npts
  integer :: fibre_count_value = 1

  real(mytype), allocatable :: x_f(:,:)
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: a_f(:,:)
  real(mytype) :: min_wall_clearance_limit = 1.0e-2_mytype
  real(mytype) :: min_point_spacing_limit = 1.0e-6_mytype
  real(mytype) :: max_initial_velocity_limit = 1.0e-8_mytype
  real(mytype) :: max_initial_acceleration_limit = 1.0e-8_mytype
  real(mytype) :: max_structure_update_limit = 1.0e-12_mytype
  real(mytype) :: max_rhs_increment_limit = 1.0e-8_mytype
  real(mytype) :: fibre_length_value = 0.0_mytype
  real(mytype) :: min_point_spacing_value = 0.0_mytype
  real(mytype) :: min_wall_clearance_value = 0.0_mytype
  real(mytype) :: max_initial_velocity_value = 0.0_mytype
  real(mytype) :: max_initial_acceleration_value = 0.0_mytype

  integer :: stage16_2_requested_status = 1
  integer :: one_fibre_case_definition_status = 0
  integer :: one_fibre_count_status = 1
  integer :: npts_valid_status = 0
  integer :: position_finite_status = 0
  integer :: velocity_finite_status = 0
  integer :: acceleration_finite_status = 0
  integer :: initial_velocity_bound_status = 0
  integer :: initial_acceleration_bound_status = 0
  integer :: point_spacing_status = 0
  integer :: wall_clearance_status = 0
  integer :: domain_containment_status = 0
  integer :: invalid_npts_rejection_status = 1
  integer :: invalid_geometry_rejection_status = 1
  integer :: invalid_velocity_rejection_status = 1
  integer :: invalid_acceleration_rejection_status = 1
  integer :: wall_contact_rejection_status = 1
  integer :: multifibre_rejection_status = 1
  integer :: diagnostic_only_status = 1
  integer :: numeric_parse_status = 1
  integer :: numeric_bounds_status = 1
  integer :: no_production_hook_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_rhs_modification_status = 1
  integer :: no_bending_solve_status = 1
  integer :: no_tension_solve_status = 1
  integer :: no_wall_contact_status = 1
  integer :: no_multifibre_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: no_legacy_ibm_forcing_status = 1
  integer :: stage16_1_config_status = 1
  integer :: case_status = 0

contains

  subroutine stage16_one_fibre_case_reset()
    case_loaded = .false.
    stage16_2_enabled = .true.
    one_fibre_case_requested = .true.
    diagnostic_only_value = .true.
    wall_contact_requested = .false.
    npts_value = default_npts
    fibre_count_value = 1
    if (allocated(x_f)) deallocate(x_f)
    if (allocated(v_f)) deallocate(v_f)
    if (allocated(a_f)) deallocate(a_f)
    call reset_limits_and_statuses()
  end subroutine stage16_one_fibre_case_reset

  subroutine stage16_one_fibre_case_define_default()
    integer :: i
    real(mytype) :: denom

    if (allocated(x_f)) deallocate(x_f)
    if (allocated(v_f)) deallocate(v_f)
    if (allocated(a_f)) deallocate(a_f)
    if (npts_value >= 1) then
      allocate(x_f(npts_value, ndim), v_f(npts_value, ndim), a_f(npts_value, ndim))
      x_f(:,:) = 0.0_mytype
      v_f(:,:) = 0.0_mytype
      a_f(:,:) = 0.0_mytype
      denom = max(real(npts_value - 1, mytype), 1.0_mytype)
      do i = 1, npts_value
        x_f(i,1) = 0.25_mytype + 0.50_mytype * real(i - 1, mytype) / denom
        x_f(i,2) = 0.50_mytype
        x_f(i,3) = 0.50_mytype
      end do
    end if
    case_loaded = .true.
    call stage16_one_fibre_case_validate()
  end subroutine stage16_one_fibre_case_define_default

  subroutine stage16_one_fibre_case_load_from_environment()
    character(len=256) :: env_value
    integer :: env_status
    integer :: parse_status
    integer :: parsed_int
    real(mytype) :: parsed_real
    logical :: parsed_bool
    integer :: default_off_dummy, env_dummy, invalid_flag_dummy
    integer :: numeric_parse_dummy, numeric_bounds_dummy, structure_dummy, rhs_dummy

    call stage16_one_fibre_case_reset()
    call load_bool_env('STAGE16_2_ENABLE', stage16_2_enabled, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_2_ONE_FIBRE_FSI_ENABLE', one_fibre_case_requested, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_2_DIAGNOSTIC_ONLY', diagnostic_only_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_2_WALL_CONTACT_REQUEST', wall_contact_requested, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)

    call get_environment_variable('STAGE16_2_NPTS', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parse_status) parsed_int
      if (parse_status == 0) then
        npts_value = parsed_int
      else
        numeric_parse_status = 0
      end if
    end if

    call load_real_env('STAGE16_2_MIN_WALL_CLEARANCE', min_wall_clearance_limit)
    call load_real_env('STAGE16_2_MIN_POINT_SPACING', min_point_spacing_limit)
    call load_real_env('STAGE16_2_MAX_INITIAL_VELOCITY', max_initial_velocity_limit)
    call load_real_env('STAGE16_2_MAX_INITIAL_ACCELERATION', max_initial_acceleration_limit)
    call load_real_env('STAGE16_2_MAX_STRUCTURE_UPDATE', max_structure_update_limit)
    call load_real_env('STAGE16_2_MAX_RHS_INCREMENT', max_rhs_increment_limit)

    call get_environment_variable('STAGE16_2_FIBRE_COUNT', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parse_status) parsed_int
      if (parse_status == 0) then
        fibre_count_value = parsed_int
      else
        numeric_parse_status = 0
      end if
    end if

    call stage16_config_load_from_environment()
    call stage16_config_get_status_values(default_off_dummy, env_dummy, stage16_1_config_status, &
         invalid_flag_dummy, numeric_parse_dummy, numeric_bounds_dummy, structure_dummy, rhs_dummy)
    if (.not. stage16_diagnostic_only_requested()) diagnostic_only_value = .false.
    call stage16_one_fibre_case_define_default()

  contains
    subroutine load_bool_env(name, value, parse_status_out)
      character(len=*), intent(in) :: name
      logical, intent(inout) :: value
      integer, intent(out) :: parse_status_out
      character(len=256) :: local_value
      integer :: local_status

      parse_status_out = 1
      call get_environment_variable(name, value=local_value, status=local_status)
      if (local_status == 0) then
        call parse_bool(local_value, parsed_bool, parse_status_out)
        if (parse_status_out == 1) value = parsed_bool
      end if
    end subroutine load_bool_env

    subroutine load_real_env(name, value)
      character(len=*), intent(in) :: name
      real(mytype), intent(inout) :: value
      character(len=256) :: local_value
      integer :: local_status

      call get_environment_variable(name, value=local_value, status=local_status)
      if (local_status == 0) then
        read(local_value, *, iostat=parse_status) parsed_real
        if (parse_status == 0 .and. finite_real(parsed_real)) then
          value = parsed_real
        else
          numeric_parse_status = 0
        end if
      end if
    end subroutine load_real_env
  end subroutine stage16_one_fibre_case_load_from_environment

  subroutine stage16_one_fibre_case_set_for_test(npts_in, fibre_count_in, make_bad_geometry, &
       make_bad_velocity, make_bad_acceleration, wall_contact_in)
    integer, intent(in) :: npts_in
    integer, intent(in) :: fibre_count_in
    logical, intent(in) :: make_bad_geometry
    logical, intent(in) :: make_bad_velocity
    logical, intent(in) :: make_bad_acceleration
    logical, intent(in) :: wall_contact_in

    call stage16_one_fibre_case_reset()
    npts_value = npts_in
    fibre_count_value = fibre_count_in
    wall_contact_requested = wall_contact_in
    call stage16_one_fibre_case_define_default()
    if (allocated(x_f) .and. make_bad_geometry) x_f(:,1) = 0.0_mytype
    if (allocated(v_f) .and. make_bad_velocity) v_f(1,1) = 2.0_mytype * max_initial_velocity_limit
    if (allocated(a_f) .and. make_bad_acceleration) a_f(1,1) = 2.0_mytype * max_initial_acceleration_limit
    call stage16_one_fibre_case_validate()
  end subroutine stage16_one_fibre_case_set_for_test

  subroutine stage16_one_fibre_case_validate()
    integer :: i
    real(mytype) :: dx(ndim)
    real(mytype) :: spacing

    call reset_status_values_only()
    stage16_2_requested_status = logical_to_int(stage16_2_enabled .and. one_fibre_case_requested)
    diagnostic_only_status = logical_to_int(diagnostic_only_value)
    one_fibre_count_status = logical_to_int(fibre_count_value == 1)
    multifibre_rejection_status = logical_to_int(fibre_count_value == 1)
    wall_contact_rejection_status = logical_to_int(.not. wall_contact_requested)
    npts_valid_status = logical_to_int(npts_value >= 2 .and. allocated(x_f) .and. size(x_f,1) == npts_value)

    if (allocated(x_f) .and. allocated(v_f) .and. allocated(a_f)) then
      position_finite_status = logical_to_int(all_finite_rank2(x_f))
      velocity_finite_status = logical_to_int(all_finite_rank2(v_f))
      acceleration_finite_status = logical_to_int(all_finite_rank2(a_f))
      max_initial_velocity_value = maxval(abs(v_f))
      max_initial_acceleration_value = maxval(abs(a_f))
      initial_velocity_bound_status = logical_to_int(max_initial_velocity_value <= max_initial_velocity_limit)
      initial_acceleration_bound_status = logical_to_int(max_initial_acceleration_value <= max_initial_acceleration_limit)
      domain_containment_status = logical_to_int(minval(x_f) > 0.0_mytype .and. maxval(x_f) < 1.0_mytype)
      min_wall_clearance_value = min(minval(x_f), 1.0_mytype - maxval(x_f))
      wall_clearance_status = logical_to_int(min_wall_clearance_value > min_wall_clearance_limit)
      if (npts_value >= 2) then
        min_point_spacing_value = huge(min_point_spacing_value)
        fibre_length_value = 0.0_mytype
        do i = 1, npts_value - 1
          dx(:) = x_f(i+1,:) - x_f(i,:)
          spacing = sqrt(sum(dx * dx))
          min_point_spacing_value = min(min_point_spacing_value, spacing)
          fibre_length_value = fibre_length_value + spacing
        end do
        point_spacing_status = logical_to_int(min_point_spacing_value > min_point_spacing_limit)
      end if
    end if

    numeric_bounds_status = logical_to_int(finite_real(min_wall_clearance_limit) .and. &
         finite_real(min_point_spacing_limit) .and. finite_real(max_initial_velocity_limit) .and. &
         finite_real(max_initial_acceleration_limit) .and. finite_real(max_structure_update_limit) .and. &
         finite_real(max_rhs_increment_limit) .and. min_wall_clearance_limit >= 0.0_mytype .and. &
         min_point_spacing_limit >= 0.0_mytype .and. max_initial_velocity_limit >= 0.0_mytype .and. &
         max_initial_acceleration_limit >= 0.0_mytype .and. max_structure_update_limit >= 0.0_mytype .and. &
         max_rhs_increment_limit >= 0.0_mytype)

    invalid_npts_rejection_status = logical_to_int(npts_value >= 2)
    invalid_geometry_rejection_status = logical_to_int(position_finite_status == 1 .and. point_spacing_status == 1 .and. &
         wall_clearance_status == 1 .and. domain_containment_status == 1)
    invalid_velocity_rejection_status = logical_to_int(velocity_finite_status == 1 .and. initial_velocity_bound_status == 1)
    invalid_acceleration_rejection_status = logical_to_int(acceleration_finite_status == 1 .and. &
         initial_acceleration_bound_status == 1)

    one_fibre_case_definition_status = logical_to_int(stage16_2_requested_status == 1 .and. &
         stage16_1_config_status == 1 .and. one_fibre_count_status == 1 .and. &
         npts_valid_status == 1 .and. position_finite_status == 1 .and. &
         velocity_finite_status == 1 .and. acceleration_finite_status == 1 .and. &
         point_spacing_status == 1 .and. &
         wall_clearance_status == 1 .and. domain_containment_status == 1 .and. diagnostic_only_status == 1)

    case_status = logical_to_int(one_fibre_case_definition_status == 1 .and. numeric_parse_status == 1 .and. &
         numeric_bounds_status == 1 .and. initial_velocity_bound_status == 1 .and. &
         initial_acceleration_bound_status == 1 .and. wall_contact_rejection_status == 1 .and. &
         multifibre_rejection_status == 1)
  end subroutine stage16_one_fibre_case_validate

  integer function stage16_one_fibre_case_get_npts()
    if (.not. case_loaded) call stage16_one_fibre_case_load_from_environment()
    stage16_one_fibre_case_get_npts = npts_value
  end function stage16_one_fibre_case_get_npts

  subroutine stage16_one_fibre_case_get_position(index, position)
    integer, intent(in) :: index
    real(mytype), intent(out) :: position(ndim)
    position(:) = 0.0_mytype
    if (.not. case_loaded) call stage16_one_fibre_case_load_from_environment()
    if (allocated(x_f) .and. index >= 1 .and. index <= npts_value) position(:) = x_f(index,:)
  end subroutine stage16_one_fibre_case_get_position

  subroutine stage16_one_fibre_case_get_velocity(index, velocity)
    integer, intent(in) :: index
    real(mytype), intent(out) :: velocity(ndim)
    velocity(:) = 0.0_mytype
    if (.not. case_loaded) call stage16_one_fibre_case_load_from_environment()
    if (allocated(v_f) .and. index >= 1 .and. index <= npts_value) velocity(:) = v_f(index,:)
  end subroutine stage16_one_fibre_case_get_velocity

  subroutine stage16_one_fibre_case_get_acceleration(index, acceleration)
    integer, intent(in) :: index
    real(mytype), intent(out) :: acceleration(ndim)
    acceleration(:) = 0.0_mytype
    if (.not. case_loaded) call stage16_one_fibre_case_load_from_environment()
    if (allocated(a_f) .and. index >= 1 .and. index <= npts_value) acceleration(:) = a_f(index,:)
  end subroutine stage16_one_fibre_case_get_acceleration

  subroutine stage16_one_fibre_case_get_status_values(case_out, count_out, npts_out, spacing_out, &
       geometry_out, velocity_out, acceleration_out, contact_out, multifibre_out, final_out)
    integer, intent(out) :: case_out
    integer, intent(out) :: count_out
    integer, intent(out) :: npts_out
    integer, intent(out) :: spacing_out
    integer, intent(out) :: geometry_out
    integer, intent(out) :: velocity_out
    integer, intent(out) :: acceleration_out
    integer, intent(out) :: contact_out
    integer, intent(out) :: multifibre_out
    integer, intent(out) :: final_out

    case_out = one_fibre_case_definition_status
    count_out = one_fibre_count_status
    npts_out = npts_valid_status
    spacing_out = point_spacing_status
    geometry_out = invalid_geometry_rejection_status
    velocity_out = invalid_velocity_rejection_status
    acceleration_out = invalid_acceleration_rejection_status
    contact_out = wall_contact_rejection_status
    multifibre_out = multifibre_rejection_status
    final_out = case_status
  end subroutine stage16_one_fibre_case_get_status_values

  subroutine stage16_one_fibre_case_write_diagnostics(unit)
    integer, intent(in) :: unit

    write(unit,'(A,1X,I0)') 'stage16_2_requested_status', stage16_2_requested_status
    write(unit,'(A,1X,I0)') 'one_fibre_case_definition_status', one_fibre_case_definition_status
    write(unit,'(A,1X,I0)') 'one_fibre_count_status', one_fibre_count_status
    write(unit,'(A,1X,I0)') 'npts', npts_value
    write(unit,'(A,1X,I0)') 'npts_valid_status', npts_valid_status
    write(unit,'(A,1X,I0)') 'position_finite_status', position_finite_status
    write(unit,'(A,1X,I0)') 'velocity_finite_status', velocity_finite_status
    write(unit,'(A,1X,I0)') 'acceleration_finite_status', acceleration_finite_status
    write(unit,'(A,1X,I0)') 'initial_velocity_bound_status', initial_velocity_bound_status
    write(unit,'(A,1X,I0)') 'initial_acceleration_bound_status', initial_acceleration_bound_status
    write(unit,'(A,1X,ES24.16E3)') 'min_point_spacing', min_point_spacing_value
    write(unit,'(A,1X,I0)') 'point_spacing_status', point_spacing_status
    write(unit,'(A,1X,ES24.16E3)') 'min_wall_clearance', min_wall_clearance_value
    write(unit,'(A,1X,I0)') 'wall_clearance_status', wall_clearance_status
    write(unit,'(A,1X,I0)') 'domain_containment_status', domain_containment_status
    write(unit,'(A,1X,I0)') 'invalid_npts_rejection_status', invalid_npts_rejection_status
    write(unit,'(A,1X,I0)') 'invalid_geometry_rejection_status', invalid_geometry_rejection_status
    write(unit,'(A,1X,I0)') 'invalid_velocity_rejection_status', invalid_velocity_rejection_status
    write(unit,'(A,1X,I0)') 'invalid_acceleration_rejection_status', invalid_acceleration_rejection_status
    write(unit,'(A,1X,I0)') 'wall_contact_rejection_status', wall_contact_rejection_status
    write(unit,'(A,1X,I0)') 'multifibre_rejection_status', multifibre_rejection_status
    write(unit,'(A,1X,I0)') 'diagnostic_only_status', diagnostic_only_status
    write(unit,'(A,1X,I0)') 'stage16_1_config_status', stage16_1_config_status
    write(unit,'(A,1X,ES24.16E3)') 'fibre_length_estimate', fibre_length_value
    write(unit,'(A,1X,ES24.16E3)') 'max_initial_velocity', max_initial_velocity_value
    write(unit,'(A,1X,ES24.16E3)') 'max_initial_acceleration', max_initial_acceleration_value
    write(unit,'(A,1X,ES24.16E3)') 'max_structure_update', max_structure_update_limit
    write(unit,'(A,1X,ES24.16E3)') 'max_rhs_increment', max_rhs_increment_limit
    write(unit,'(A,1X,I0)') 'no_production_hook_status', no_production_hook_status
    write(unit,'(A,1X,I0)') 'no_structure_advance_status', no_structure_advance_status
    write(unit,'(A,1X,I0)') 'no_rhs_modification_status', no_rhs_modification_status
    write(unit,'(A,1X,I0)') 'no_bending_solve_status', no_bending_solve_status
    write(unit,'(A,1X,I0)') 'no_tension_solve_status', no_tension_solve_status
    write(unit,'(A,1X,I0)') 'no_wall_contact_status', no_wall_contact_status
    write(unit,'(A,1X,I0)') 'no_multifibre_status', no_multifibre_status
    write(unit,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit,'(A,1X,I0)') 'no_legacy_ibm_forcing_status', no_legacy_ibm_forcing_status
    write(unit,'(A,1X,I0)') 'case_status', case_status
  end subroutine stage16_one_fibre_case_write_diagnostics

  subroutine reset_limits_and_statuses()
    min_wall_clearance_limit = 1.0e-2_mytype
    min_point_spacing_limit = 1.0e-6_mytype
    max_initial_velocity_limit = 1.0e-8_mytype
    max_initial_acceleration_limit = 1.0e-8_mytype
    max_structure_update_limit = 1.0e-12_mytype
    max_rhs_increment_limit = 1.0e-8_mytype
    numeric_parse_status = 1
    numeric_bounds_status = 1
    no_production_hook_status = 1
    no_structure_advance_status = 1
    no_rhs_modification_status = 1
    no_bending_solve_status = 1
    no_tension_solve_status = 1
    no_wall_contact_status = 1
    no_multifibre_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    no_legacy_ibm_forcing_status = 1
    stage16_1_config_status = 1
    call reset_status_values_only()
  end subroutine reset_limits_and_statuses

  subroutine reset_status_values_only()
    fibre_length_value = 0.0_mytype
    min_point_spacing_value = 0.0_mytype
    min_wall_clearance_value = 0.0_mytype
    max_initial_velocity_value = 0.0_mytype
    max_initial_acceleration_value = 0.0_mytype
    stage16_2_requested_status = 1
    one_fibre_case_definition_status = 0
    one_fibre_count_status = 0
    npts_valid_status = 0
    position_finite_status = 0
    velocity_finite_status = 0
    acceleration_finite_status = 0
    initial_velocity_bound_status = 0
    initial_acceleration_bound_status = 0
    point_spacing_status = 0
    wall_clearance_status = 0
    domain_containment_status = 0
    invalid_npts_rejection_status = 0
    invalid_geometry_rejection_status = 0
    invalid_velocity_rejection_status = 0
    invalid_acceleration_rejection_status = 0
    wall_contact_rejection_status = 0
    multifibre_rejection_status = 0
    diagnostic_only_status = 0
    case_status = 0
  end subroutine reset_status_values_only

  subroutine parse_bool(value, parsed_bool, parse_status)
    character(len=*), intent(in) :: value
    logical, intent(out) :: parsed_bool
    integer, intent(out) :: parse_status
    character(len=256) :: lowered

    lowered = lower_string(trim(adjustl(value)))
    select case (trim(lowered))
    case ('1', 'true', 't', 'yes', 'y', 'on')
      parsed_bool = .true.
      parse_status = 1
    case ('0', 'false', 'f', 'no', 'n', 'off')
      parsed_bool = .false.
      parse_status = 1
    case default
      parsed_bool = .false.
      parse_status = 0
    end select
  end subroutine parse_bool

  character(len=256) function lower_string(value)
    character(len=*), intent(in) :: value
    integer :: i
    integer :: code

    lower_string = ' '
    do i = 1, min(len_trim(value), len(lower_string))
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower_string(i:i) = achar(code + 32)
      else
        lower_string(i:i) = value(i:i)
      end if
    end do
  end function lower_string

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

  logical function finite_real(value)
    real(mytype), intent(in) :: value
    finite_real = (value == value) .and. (abs(value) < huge(value))
  end function finite_real

  logical function all_finite_rank2(values)
    real(mytype), intent(in) :: values(:,:)
    integer :: i, j

    all_finite_rank2 = .true.
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        if (.not. finite_real(values(i,j))) all_finite_rank2 = .false.
      end do
    end do
  end function all_finite_rank2

end module fibre_stage16_one_fibre_case
