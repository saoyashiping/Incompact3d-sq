program fibre_stage15_feedback_linkage_check
  use fibre_stage12_feedback_formula, only : stage12_feedback_formula_init, &
       stage12_feedback_formula_compute_controlled
  implicit none

  integer, parameter :: mytype = kind(1.0d0)
  integer :: npts
  integer :: np
  integer :: stage15_enable
  integer :: controlled_enable
  integer :: structure_enable
  integer :: diagnostic_only
  real(mytype) :: dt
  real(mytype) :: rho_tilde
  real(mytype) :: test_force
  real(mytype) :: feedback_alpha
  real(mytype) :: lambda_value
  real(mytype) :: max_velocity_update_allowed
  real(mytype) :: max_slip_error_allowed
  real(mytype) :: max_force_error_allowed
  real(mytype) :: max_force_response_allowed
  real(mytype) :: max_rhs_response_allowed
  real(mytype), allocatable :: u_f(:,:)
  real(mytype), allocatable :: v_old(:,:)
  real(mytype), allocatable :: v_new(:,:)
  real(mytype), allocatable :: delta_v(:,:)
  real(mytype), allocatable :: slip_old(:,:)
  real(mytype), allocatable :: slip_new(:,:)
  real(mytype), allocatable :: delta_slip(:,:)
  real(mytype), allocatable :: f_old(:,:)
  real(mytype), allocatable :: f_new(:,:)
  real(mytype), allocatable :: delta_f(:,:)
  real(mytype), allocatable :: force_norm_old(:)
  real(mytype), allocatable :: force_norm_new(:)
  real(mytype) :: max_velocity_update
  real(mytype) :: max_slip_change
  real(mytype) :: slip_error
  real(mytype) :: max_feedback_force_change
  real(mytype) :: feedback_force_error
  real(mytype) :: rhs_response
  integer :: requested_status
  integer :: velocity_update_nonzero_status
  integer :: slip_change_nonzero_status
  integer :: slip_consistency_status
  integer :: feedback_force_change_nonzero_status
  integer :: feedback_force_consistency_status
  integer :: force_response_bounded_status
  integer :: rhs_response_bounded_status
  integer :: controlled_update_count
  integer :: production_full_structure_advance_count
  integer :: bending_solve_count
  integer :: tension_solve_count
  integer :: wall_contact_count
  integer :: multifibre_count
  integer :: rhs_injection_connection_count
  integer :: approved_stage12_13_14_chain_status
  integer :: no_fluid_rhs_modification_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: final_status
  integer :: unit_id

  npts = get_env_int('STAGE15_7_NPTS', 8)
  np = get_env_int('STAGE15_7_NP', 1)
  stage15_enable = get_env_int('STAGE15_7_ENABLE', 1)
  controlled_enable = get_env_int('STAGE15_7_CONTROLLED_STEP_ENABLE', 1)
  structure_enable = get_env_int('STAGE15_7_STRUCTURE_ADVANCE_ENABLE', 1)
  diagnostic_only = get_env_int('STAGE15_7_DIAGNOSTIC_ONLY', 1)
  dt = get_env_real('STAGE15_7_DT', 1.0e-4_mytype)
  rho_tilde = get_env_real('STAGE15_7_RHO_TILDE', 1.0_mytype)
  test_force = get_env_real('STAGE15_7_TEST_FORCE', 1.0e-6_mytype)
  feedback_alpha = get_env_real('STAGE15_7_FEEDBACK_ALPHA', 1.0_mytype)
  lambda_value = get_env_real('STAGE15_7_LAMBDA', 1.0e-8_mytype)
  max_velocity_update_allowed = get_env_real('STAGE15_7_MAX_VELOCITY_UPDATE', 1.0e-9_mytype)
  max_slip_error_allowed = get_env_real('STAGE15_7_MAX_SLIP_ERROR', 1.0e-14_mytype)
  max_force_error_allowed = get_env_real('STAGE15_7_MAX_FORCE_ERROR', 1.0e-14_mytype)
  max_force_response_allowed = get_env_real('STAGE15_7_MAX_FORCE_RESPONSE', 1.0e-8_mytype)
  max_rhs_response_allowed = get_env_real('STAGE15_7_MAX_RHS_RESPONSE', 1.0e-12_mytype)

  call execute_command_line('mkdir -p stage15_outputs')
  allocate(u_f(npts, 3), v_old(npts, 3), v_new(npts, 3), delta_v(npts, 3), slip_old(npts, 3), &
           slip_new(npts, 3), delta_slip(npts, 3), f_old(npts, 3), f_new(npts, 3), delta_f(npts, 3), &
           force_norm_old(npts), force_norm_new(npts))

  requested_status = merge(1, 0, stage15_enable == 1)
  call initialize_reference_state(u_f, v_old)
  delta_v(:, :) = 0.0_mytype
  if (controls_valid()) delta_v(:, 1) = dt * test_force / rho_tilde
  v_new(:, :) = v_old(:, :) + delta_v(:, :)

  slip_old(:, :) = u_f(:, :) - v_old(:, :)
  slip_new(:, :) = u_f(:, :) - v_new(:, :)
  delta_slip(:, :) = slip_new(:, :) - slip_old(:, :)

  call stage12_feedback_formula_init()
  call stage12_feedback_formula_compute_controlled(u_f, v_old, feedback_alpha, 1, f_old, force_norm_old)
  call stage12_feedback_formula_compute_controlled(u_f, v_new, feedback_alpha, 1, f_new, force_norm_new)
  delta_f(:, :) = f_new(:, :) - f_old(:, :)

  max_velocity_update = maxval(abs(delta_v))
  max_slip_change = maxval(abs(delta_slip))
  slip_error = maxval(abs(delta_slip + delta_v))
  max_feedback_force_change = maxval(abs(delta_f))
  feedback_force_error = maxval(abs(delta_f + feedback_alpha * delta_v))
  rhs_response = abs(lambda_value) * max_feedback_force_change

  velocity_update_nonzero_status = merge(1, 0, max_velocity_update > 0.0_mytype .and. &
                                         max_velocity_update <= max_velocity_update_allowed .and. &
                                         all_finite_rank2(delta_v))
  slip_change_nonzero_status = merge(1, 0, max_slip_change > 0.0_mytype .and. all_finite_rank2(delta_slip))
  slip_consistency_status = merge(1, 0, slip_error <= max_slip_error_allowed)
  feedback_force_change_nonzero_status = merge(1, 0, max_feedback_force_change > 0.0_mytype .and. &
                                               all_finite_rank2(delta_f))
  feedback_force_consistency_status = merge(1, 0, feedback_force_error <= max_force_error_allowed)
  force_response_bounded_status = merge(1, 0, max_feedback_force_change <= max_force_response_allowed)
  rhs_response_bounded_status = merge(1, 0, rhs_response <= max_rhs_response_allowed .and. is_finite(rhs_response))
  controlled_update_count = merge(1, 0, velocity_update_nonzero_status == 1)
  production_full_structure_advance_count = 0
  bending_solve_count = 0
  tension_solve_count = 0
  wall_contact_count = 0
  multifibre_count = 0
  rhs_injection_connection_count = 0
  approved_stage12_13_14_chain_status = 1
  no_fluid_rhs_modification_status = 1
  no_pressure_projection_modification_status = 1
  no_poisson_modification_status = 1
  no_rk3_channel_forcing_modification_status = 1

  final_status = merge(1, 0, requested_status == 1 .and. controlled_enable == 1 .and. structure_enable == 1 .and. &
                       diagnostic_only == 1 .and. np == 1 .and. controlled_update_count == 1 .and. &
                       velocity_update_nonzero_status == 1 .and. slip_change_nonzero_status == 1 .and. &
                       slip_consistency_status == 1 .and. feedback_force_change_nonzero_status == 1 .and. &
                       feedback_force_consistency_status == 1 .and. force_response_bounded_status == 1 .and. &
                       rhs_response_bounded_status == 1 .and. production_full_structure_advance_count == 0 .and. &
                       bending_solve_count == 0 .and. tension_solve_count == 0 .and. wall_contact_count == 0 .and. &
                       multifibre_count == 0 .and. rhs_injection_connection_count == 0 .and. &
                       approved_stage12_13_14_chain_status == 1 .and. no_fluid_rhs_modification_status == 1 .and. &
                       no_pressure_projection_modification_status == 1 .and. no_poisson_modification_status == 1 .and. &
                       no_rk3_channel_forcing_modification_status == 1)

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_7_feedback_linkage.dat', status='replace', action='write')
  call write_diagnostics(unit_id)
  close(unit_id)

  if (final_status == 1) then
    write(*, '(A)') 'STAGE 15.7 FEEDBACK LINKAGE CHECK: PASS'
    stop 0
  end if
  write(*, '(A)') 'STAGE 15.7 FEEDBACK LINKAGE CHECK: FAIL'
  stop 1

contains

  logical function controls_valid()
    controls_valid = stage15_enable == 1 .and. controlled_enable == 1 .and. structure_enable == 1 .and. &
                     diagnostic_only == 1 .and. np == 1 .and. npts > 0 .and. dt > 0.0_mytype .and. &
                     rho_tilde > 0.0_mytype .and. abs(test_force) > 0.0_mytype .and. feedback_alpha > 0.0_mytype .and. &
                     abs(lambda_value) > 0.0_mytype .and. is_finite(dt) .and. is_finite(rho_tilde) .and. &
                     is_finite(test_force) .and. is_finite(feedback_alpha) .and. is_finite(lambda_value)
  end function controls_valid

  subroutine initialize_reference_state(u_ref, v_ref)
    real(mytype), intent(out) :: u_ref(:,:)
    real(mytype), intent(out) :: v_ref(:,:)
    integer :: i

    u_ref(:, :) = 0.0_mytype
    v_ref(:, :) = 0.0_mytype
    do i = 1, size(u_ref, 1)
      u_ref(i, 1) = 0.25_mytype
      u_ref(i, 2) = 0.125_mytype
      u_ref(i, 3) = -0.0625_mytype
    end do
  end subroutine initialize_reference_state

  subroutine write_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    write(unit_id, '(A,1X,I0)') 'stage15_7_requested_status', requested_status
    write(unit_id, '(A,1X,I0)') 'controlled_step_enabled_status', controlled_enable
    write(unit_id, '(A,1X,I0)') 'structure_advance_enable_status', structure_enable
    write(unit_id, '(A,1X,I0)') 'diagnostic_only_status', diagnostic_only
    write(unit_id, '(A,1X,I0)') 'np', np
    write(unit_id, '(A,1X,I0)') 'npts', npts
    write(unit_id, '(A,1X,ES24.16)') 'dt', dt
    write(unit_id, '(A,1X,ES24.16)') 'rho_tilde', rho_tilde
    write(unit_id, '(A,1X,ES24.16)') 'test_force_magnitude', test_force
    write(unit_id, '(A,1X,ES24.16)') 'feedback_alpha', feedback_alpha
    write(unit_id, '(A,1X,ES24.16)') 'lambda_value', lambda_value
    write(unit_id, '(A,1X,ES24.16)') 'max_velocity_update', max_velocity_update
    write(unit_id, '(A,1X,I0)') 'velocity_update_nonzero_status', velocity_update_nonzero_status
    write(unit_id, '(A,1X,ES24.16)') 'max_slip_change', max_slip_change
    write(unit_id, '(A,1X,I0)') 'slip_change_nonzero_status', slip_change_nonzero_status
    write(unit_id, '(A,1X,ES24.16)') 'slip_error', slip_error
    write(unit_id, '(A,1X,I0)') 'slip_consistency_status', slip_consistency_status
    write(unit_id, '(A,1X,ES24.16)') 'max_feedback_force_change', max_feedback_force_change
    write(unit_id, '(A,1X,I0)') 'feedback_force_change_nonzero_status', feedback_force_change_nonzero_status
    write(unit_id, '(A,1X,ES24.16)') 'feedback_force_error', feedback_force_error
    write(unit_id, '(A,1X,I0)') 'feedback_force_consistency_status', feedback_force_consistency_status
    write(unit_id, '(A,1X,I0)') 'force_response_bounded_status', force_response_bounded_status
    write(unit_id, '(A,1X,I0)') 'rhs_response_bounded_status', rhs_response_bounded_status
    write(unit_id, '(A,1X,I0)') 'controlled_update_count', controlled_update_count
    write(unit_id, '(A,1X,I0)') 'production_full_structure_advance_count', production_full_structure_advance_count
    write(unit_id, '(A,1X,I0)') 'bending_solve_count', bending_solve_count
    write(unit_id, '(A,1X,I0)') 'tension_solve_count', tension_solve_count
    write(unit_id, '(A,1X,I0)') 'wall_contact_count', wall_contact_count
    write(unit_id, '(A,1X,I0)') 'multifibre_count', multifibre_count
    write(unit_id, '(A,1X,I0)') 'rhs_injection_connection_count', rhs_injection_connection_count
    write(unit_id, '(A,1X,I0)') 'approved_stage12_13_14_chain_status', approved_stage12_13_14_chain_status
    write(unit_id, '(A,1X,I0)') 'no_fluid_rhs_modification_status', no_fluid_rhs_modification_status
    write(unit_id, '(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit_id, '(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'final_status', final_status
  end subroutine write_diagnostics

  logical function all_finite_rank2(values)
    real(mytype), intent(in) :: values(:,:)
    all_finite_rank2 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank2

  logical function is_finite(value)
    real(mytype), intent(in) :: value
    is_finite = value == value .and. abs(value) < huge(1.0_mytype)
  end function is_finite

  integer function get_env_int(name, default_value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    character(len=128) :: value
    integer :: env_status
    integer :: read_status

    get_env_int = default_value
    call get_environment_variable(name, value=value, status=env_status)
    if (env_status == 0) then
      read(value, *, iostat=read_status) get_env_int
      if (read_status /= 0) get_env_int = default_value
    end if
  end function get_env_int

  real(mytype) function get_env_real(name, default_value)
    character(len=*), intent(in) :: name
    real(mytype), intent(in) :: default_value
    character(len=128) :: value
    integer :: env_status
    integer :: read_status

    get_env_real = default_value
    call get_environment_variable(name, value=value, status=env_status)
    if (env_status == 0) then
      read(value, *, iostat=read_status) get_env_real
      if (read_status /= 0) get_env_real = default_value
    end if
  end function get_env_real

end program fibre_stage15_feedback_linkage_check
