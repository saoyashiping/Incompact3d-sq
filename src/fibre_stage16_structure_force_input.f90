module fibre_stage16_structure_force_input
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage16_structure_force_input_reset
  public :: stage16_structure_force_input_load_from_environment
  public :: stage16_structure_force_input_set_from_stage12_candidate
  public :: stage16_structure_force_input_validate
  public :: stage16_structure_force_input_get_force
  public :: stage16_structure_force_input_write_diagnostics
  public :: stage16_structure_force_input_get_final_status

  real(mytype) :: feedback_alpha_value = 1.0_mytype
  real(mytype) :: max_action_reaction_error = 1.0e-14_mytype
  real(mytype) :: max_sign_error = 1.0e-14_mytype
  real(mytype) :: max_force_input_value = 1.0e-6_mytype
  real(mytype) :: structure_force_input(1,3) = 0.0_mytype
  real(mytype) :: fluid_side_force(1,3) = 0.0_mytype
  real(mytype) :: stage12_candidate(1,3) = 0.0_mytype
  real(mytype) :: action_reaction_error_value = 0.0_mytype
  real(mytype) :: max_structure_force_input = 0.0_mytype
  real(mytype) :: zero_slip_input_value = 0.0_mytype

  integer :: stage16_4_requested_status = 1
  integer :: structure_force_input_finite_status = 0
  integer :: structure_force_input_bounded_status = 0
  integer :: structure_side_force_sign_status = 0
  integer :: fluid_side_force_sign_status = 0
  integer :: action_reaction_status = 0
  integer :: wrong_sign_rejection_status = 0
  integer :: zero_slip_zero_input_status = 0
  integer :: force_input_reset_status = 1
  integer :: force_input_readback_status = 0
  integer :: approved_stage12_13_14_chain_status = 1
  integer :: no_production_hook_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_position_update_status = 1
  integer :: no_velocity_update_status = 1
  integer :: no_acceleration_update_status = 1
  integer :: no_rhs_modification_status = 1
  integer :: no_bending_solve_status = 1
  integer :: no_tension_solve_status = 1
  integer :: no_wall_contact_status = 1
  integer :: no_multifibre_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: no_legacy_ibm_forcing_status = 1
  integer :: numeric_parse_status = 1
  integer :: numeric_bounds_status = 1
  integer :: diagnostic_only_status = 1
  integer :: force_input_status = 0

contains

  subroutine stage16_structure_force_input_reset()
    feedback_alpha_value = 1.0_mytype
    max_action_reaction_error = 1.0e-14_mytype
    max_sign_error = 1.0e-14_mytype
    max_force_input_value = 1.0e-6_mytype
    structure_force_input(:,:) = 0.0_mytype
    fluid_side_force(:,:) = 0.0_mytype
    stage12_candidate(:,:) = 0.0_mytype
    action_reaction_error_value = 0.0_mytype
    max_structure_force_input = 0.0_mytype
    zero_slip_input_value = 0.0_mytype
    stage16_4_requested_status = 1
    structure_force_input_finite_status = 0
    structure_force_input_bounded_status = 0
    structure_side_force_sign_status = 0
    fluid_side_force_sign_status = 0
    action_reaction_status = 0
    wrong_sign_rejection_status = 0
    zero_slip_zero_input_status = 0
    force_input_reset_status = 1
    force_input_readback_status = 0
    approved_stage12_13_14_chain_status = 1
    no_production_hook_status = 1
    no_structure_advance_status = 1
    no_position_update_status = 1
    no_velocity_update_status = 1
    no_acceleration_update_status = 1
    no_rhs_modification_status = 1
    no_bending_solve_status = 1
    no_tension_solve_status = 1
    no_wall_contact_status = 1
    no_multifibre_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    no_legacy_ibm_forcing_status = 1
    numeric_parse_status = 1
    numeric_bounds_status = 1
    diagnostic_only_status = 1
    force_input_status = 0
  end subroutine stage16_structure_force_input_reset

  subroutine stage16_structure_force_input_load_from_environment()
    character(len=256) :: env_value
    integer :: env_status
    integer :: parse_status
    real(mytype) :: parsed_real
    logical :: parsed_bool

    call stage16_structure_force_input_reset()
    call load_bool_env('STAGE16_4_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) stage16_4_requested_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_4_ONE_FIBRE_FSI_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) stage16_4_requested_status = min(stage16_4_requested_status, logical_to_int(parsed_bool))
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_4_DIAGNOSTIC_ONLY', parsed_bool, parse_status)
    if (parse_status == 1) diagnostic_only_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_4_FEEDBACK_ALPHA', feedback_alpha_value)
    call load_real_env('STAGE16_4_MAX_ACTION_REACTION_ERROR', max_action_reaction_error)
    call load_real_env('STAGE16_4_MAX_SIGN_ERROR', max_sign_error)
    call load_real_env('STAGE16_4_MAX_FORCE_INPUT', max_force_input_value)
    numeric_bounds_status = logical_to_int(feedback_alpha_value > 0.0_mytype .and. &
         max_action_reaction_error >= 0.0_mytype .and. max_sign_error >= 0.0_mytype .and. &
         max_force_input_value >= 0.0_mytype)

  contains
    subroutine load_real_env(name, value)
      character(len=*), intent(in) :: name
      real(mytype), intent(inout) :: value
      call get_environment_variable(name, value=env_value, status=env_status)
      if (env_status == 0) then
        read(env_value, *, iostat=parse_status) parsed_real
        if (parse_status == 0 .and. finite_real(parsed_real)) then
          value = parsed_real
        else
          numeric_parse_status = 0
        end if
      end if
    end subroutine load_real_env

    subroutine load_bool_env(name, value, bool_status)
      character(len=*), intent(in) :: name
      logical, intent(out) :: value
      integer, intent(out) :: bool_status
      value = .true.
      bool_status = 1
      call get_environment_variable(name, value=env_value, status=env_status)
      if (env_status == 0) call parse_bool(env_value, value, bool_status)
    end subroutine load_bool_env
  end subroutine stage16_structure_force_input_load_from_environment

  subroutine stage16_structure_force_input_set_from_stage12_candidate(u_f, v_f, alpha, use_wrong_sign)
    real(mytype), intent(in) :: u_f(3)
    real(mytype), intent(in) :: v_f(3)
    real(mytype), intent(in) :: alpha
    logical, intent(in) :: use_wrong_sign
    real(mytype) :: zero_slip(3)

    feedback_alpha_value = alpha
    stage12_candidate(1,:) = feedback_alpha_value * (u_f(:) - v_f(:))
    if (use_wrong_sign) then
      structure_force_input(1,:) = -stage12_candidate(1,:)
    else
      structure_force_input(1,:) = stage12_candidate(1,:)
    end if
    fluid_side_force(1,:) = -stage12_candidate(1,:)
    zero_slip(:) = 0.0_mytype
    zero_slip_input_value = maxval(abs(feedback_alpha_value * zero_slip))
    call stage16_structure_force_input_validate()
  end subroutine stage16_structure_force_input_set_from_stage12_candidate

  subroutine stage16_structure_force_input_validate()
    real(mytype) :: readback(3)
    real(mytype) :: wrong_structure_input(3)

    structure_force_input_finite_status = logical_to_int(all_finite_rank2(structure_force_input))
    max_structure_force_input = maxval(abs(structure_force_input))
    structure_force_input_bounded_status = logical_to_int(max_structure_force_input <= max_force_input_value)
    structure_side_force_sign_status = logical_to_int(vector_equal(structure_force_input(1,:), stage12_candidate(1,:), &
         max_sign_error))
    fluid_side_force_sign_status = logical_to_int(vector_equal(fluid_side_force(1,:), -stage12_candidate(1,:), max_sign_error))
    action_reaction_error_value = maxval(abs(structure_force_input + fluid_side_force))
    action_reaction_status = logical_to_int(action_reaction_error_value <= max_action_reaction_error)
    wrong_structure_input(:) = -stage12_candidate(1,:)
    wrong_sign_rejection_status = logical_to_int(maxval(abs(wrong_structure_input - stage12_candidate(1,:))) > max_sign_error)
    zero_slip_zero_input_status = logical_to_int(zero_slip_input_value <= max_sign_error)
    call stage16_structure_force_input_get_force(1, readback)
    force_input_readback_status = logical_to_int(vector_equal(readback, structure_force_input(1,:), max_sign_error))
    call update_force_input_status()
  end subroutine stage16_structure_force_input_validate

  subroutine stage16_structure_force_input_get_force(index, force)
    integer, intent(in) :: index
    real(mytype), intent(out) :: force(3)
    force(:) = 0.0_mytype
    if (index == 1) force(:) = structure_force_input(1,:)
  end subroutine stage16_structure_force_input_get_force

  integer function stage16_structure_force_input_get_final_status()
    stage16_structure_force_input_get_final_status = force_input_status
  end function stage16_structure_force_input_get_final_status

  subroutine stage16_structure_force_input_write_diagnostics(unit)
    integer, intent(in) :: unit
    write(unit,'(A,1X,I0)') 'stage16_4_requested_status', stage16_4_requested_status
    write(unit,'(A,1X,ES24.16E3)') 'feedback_alpha', feedback_alpha_value
    write(unit,'(A,1X,I0)') 'structure_force_input_finite_status', structure_force_input_finite_status
    write(unit,'(A,1X,I0)') 'structure_force_input_bounded_status', structure_force_input_bounded_status
    write(unit,'(A,1X,I0)') 'structure_side_force_sign_status', structure_side_force_sign_status
    write(unit,'(A,1X,I0)') 'fluid_side_force_sign_status', fluid_side_force_sign_status
    write(unit,'(A,1X,ES24.16E3)') 'action_reaction_error', action_reaction_error_value
    write(unit,'(A,1X,I0)') 'action_reaction_status', action_reaction_status
    write(unit,'(A,1X,I0)') 'wrong_sign_rejection_status', wrong_sign_rejection_status
    write(unit,'(A,1X,I0)') 'zero_slip_zero_input_status', zero_slip_zero_input_status
    write(unit,'(A,1X,I0)') 'force_input_reset_status', force_input_reset_status
    write(unit,'(A,1X,I0)') 'force_input_readback_status', force_input_readback_status
    write(unit,'(A,1X,I0)') 'approved_stage12_13_14_chain_status', approved_stage12_13_14_chain_status
    write(unit,'(A,1X,I0)') 'no_production_hook_status', no_production_hook_status
    write(unit,'(A,1X,I0)') 'no_structure_advance_status', no_structure_advance_status
    write(unit,'(A,1X,I0)') 'no_position_update_status', no_position_update_status
    write(unit,'(A,1X,I0)') 'no_velocity_update_status', no_velocity_update_status
    write(unit,'(A,1X,I0)') 'no_acceleration_update_status', no_acceleration_update_status
    write(unit,'(A,1X,I0)') 'no_rhs_modification_status', no_rhs_modification_status
    write(unit,'(A,1X,I0)') 'no_bending_solve_status', no_bending_solve_status
    write(unit,'(A,1X,I0)') 'no_tension_solve_status', no_tension_solve_status
    write(unit,'(A,1X,I0)') 'no_wall_contact_status', no_wall_contact_status
    write(unit,'(A,1X,I0)') 'no_multifibre_status', no_multifibre_status
    write(unit,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit,'(A,1X,I0)') 'no_legacy_ibm_forcing_status', no_legacy_ibm_forcing_status
    write(unit,'(A,1X,I0)') 'numeric_parse_status', numeric_parse_status
    write(unit,'(A,1X,I0)') 'numeric_bounds_status', numeric_bounds_status
    write(unit,'(A,1X,I0)') 'diagnostic_only_status', diagnostic_only_status
    write(unit,'(A,1X,I0)') 'final_status', force_input_status
  end subroutine stage16_structure_force_input_write_diagnostics

  subroutine update_force_input_status()
    force_input_status = logical_to_int(stage16_4_requested_status == 1 .and. &
         structure_force_input_finite_status == 1 .and. structure_force_input_bounded_status == 1 .and. &
         structure_side_force_sign_status == 1 .and. fluid_side_force_sign_status == 1 .and. &
         action_reaction_status == 1 .and. wrong_sign_rejection_status == 1 .and. &
         zero_slip_zero_input_status == 1 .and. force_input_reset_status == 1 .and. &
         force_input_readback_status == 1 .and. approved_stage12_13_14_chain_status == 1 .and. &
         numeric_parse_status == 1 .and. numeric_bounds_status == 1 .and. diagnostic_only_status == 1)
  end subroutine update_force_input_status

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
    integer :: i, code
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

  logical function vector_equal(first, second, tolerance)
    real(mytype), intent(in) :: first(:)
    real(mytype), intent(in) :: second(:)
    real(mytype), intent(in) :: tolerance
    vector_equal = maxval(abs(first - second)) <= tolerance
  end function vector_equal

end module fibre_stage16_structure_force_input
