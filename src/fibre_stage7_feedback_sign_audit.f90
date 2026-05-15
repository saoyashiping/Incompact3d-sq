program fibre_stage7_feedback_sign_audit
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_lagrangian_points_t
  use fibre_ibm_feedback, only : compute_ibm_feedback_forces, compute_feedback_pair_power, compute_feedback_expected_dissipation
  implicit none

  real(mytype), parameter :: tol = 1.0e-14_mytype
  real(mytype), parameter :: beta_drag = 2.0_mytype

  type(ibm_lagrangian_points_t) :: lag
  real(mytype), allocatable :: u_lag(:,:), force_on_structure(:,:), force_on_fluid(:,:)
  real(mytype) :: zero_slip_force_norm, action_reaction_error
  real(mytype) :: structure_force_slip_dot, fluid_force_slip_dot
  real(mytype) :: power_structure, power_fluid, total_power, expected_dissipation, expected_total_power, total_power_error
  real(mytype) :: slip(3)
  integer :: zero_slip_flag, action_reaction_flag
  integer :: structure_force_slip_dot_positive_flag, fluid_force_slip_dot_negative_flag
  integer :: total_power_dissipative_flag
  integer :: real_module_called_flag, force_on_structure_from_module_flag, force_on_fluid_from_module_flag
  integer :: audit_status, io

  call init_single_point_lag(lag)
  allocate(u_lag(3, lag%nl), force_on_structure(3, lag%nl), force_on_fluid(3, lag%nl))

  real_module_called_flag = 0
  force_on_structure_from_module_flag = 0
  force_on_fluid_from_module_flag = 0

  ! Case A: zero slip
  u_lag(:, 1) = (/1.0_mytype, -0.5_mytype, 0.25_mytype/)
  lag%v(:, 1) = u_lag(:, 1)
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
  real_module_called_flag = 1
  force_on_structure_from_module_flag = 1
  force_on_fluid_from_module_flag = 1
  zero_slip_force_norm = sqrt(sum(force_on_structure**2) + sum(force_on_fluid**2))
  zero_slip_flag = merge(1, 0, zero_slip_force_norm <= tol)

  ! Case B/C/D: nonzero slip sign + action-reaction + dissipativity
  u_lag(:, 1) = (/1.2_mytype, -0.4_mytype, 0.8_mytype/)
  lag%v(:, 1) = (/0.3_mytype, -0.9_mytype, 0.1_mytype/)
  slip = u_lag(:, 1) - lag%v(:, 1)

  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)

  structure_force_slip_dot = dot_product(force_on_structure(:, 1), slip)
  fluid_force_slip_dot = dot_product(force_on_fluid(:, 1), slip)
  structure_force_slip_dot_positive_flag = merge(1, 0, structure_force_slip_dot > 0._mytype)
  fluid_force_slip_dot_negative_flag = merge(1, 0, fluid_force_slip_dot < 0._mytype)

  action_reaction_error = sqrt(sum((force_on_structure(:, 1) + force_on_fluid(:, 1))**2))
  action_reaction_flag = merge(1, 0, action_reaction_error <= tol)

  call compute_feedback_pair_power(lag, u_lag, force_on_structure, force_on_fluid, power_structure, power_fluid, total_power)
  call compute_feedback_expected_dissipation(lag, u_lag, beta_drag, expected_dissipation)
  expected_total_power = -expected_dissipation
  total_power_error = abs(total_power - expected_total_power)
  total_power_dissipative_flag = merge(1, 0, total_power <= 0._mytype .and. total_power_error <= tol)

  audit_status = merge(1, 0, &
       real_module_called_flag == 1 .and. &
       force_on_structure_from_module_flag == 1 .and. force_on_fluid_from_module_flag == 1 .and. &
       zero_slip_force_norm <= tol .and. zero_slip_flag == 1 .and. &
       action_reaction_error <= tol .and. action_reaction_flag == 1 .and. &
       structure_force_slip_dot > 0._mytype .and. fluid_force_slip_dot < 0._mytype .and. &
       structure_force_slip_dot_positive_flag == 1 .and. fluid_force_slip_dot_negative_flag == 1 .and. &
       total_power <= 0._mytype .and. total_power_error <= tol .and. total_power_dissipative_flag == 1)

  open(newunit=io, file='stage7_outputs/fibre_stage7_feedback_sign_audit.dat', status='replace', action='write')
  write(io,'(A,1X,I0)') 'stage7_feedback_real_module_called_flag', real_module_called_flag
  write(io,'(A,1X,I0)') 'stage7_feedback_force_on_structure_from_module_flag', force_on_structure_from_module_flag
  write(io,'(A,1X,I0)') 'stage7_feedback_force_on_fluid_from_module_flag', force_on_fluid_from_module_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_zero_slip_force_norm', zero_slip_force_norm
  write(io,'(A,1X,I0)') 'stage7_feedback_zero_slip_flag', zero_slip_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_action_reaction_error', action_reaction_error
  write(io,'(A,1X,I0)') 'stage7_feedback_action_reaction_flag', action_reaction_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_structure_force_slip_dot', structure_force_slip_dot
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_fluid_force_slip_dot', fluid_force_slip_dot
  write(io,'(A,1X,I0)') 'stage7_feedback_structure_force_slip_dot_positive_flag', structure_force_slip_dot_positive_flag
  write(io,'(A,1X,I0)') 'stage7_feedback_fluid_force_slip_dot_negative_flag', fluid_force_slip_dot_negative_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_total_power', total_power
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_expected_total_power', expected_total_power
  write(io,'(A,1X,ES24.16E3)') 'stage7_feedback_total_power_error', total_power_error
  write(io,'(A,1X,I0)') 'stage7_feedback_total_power_dissipative_flag', total_power_dissipative_flag
  write(io,'(A,1X,I0)') 'stage7_feedback_sign_audit_status', audit_status
  close(io)

contains
  subroutine init_single_point_lag(lag)
    type(ibm_lagrangian_points_t), intent(out) :: lag

    lag%nl = 1
    allocate(lag%x(3,1), lag%v(3,1), lag%force(3,1), lag%weight(1))
    lag%x(:,1) = (/0.5_mytype, 0.5_mytype, 0.5_mytype/)
    lag%v(:,1) = 0._mytype
    lag%force(:,1) = 0._mytype
    lag%weight(1) = 1._mytype
  end subroutine init_single_point_lag

end program fibre_stage7_feedback_sign_audit
