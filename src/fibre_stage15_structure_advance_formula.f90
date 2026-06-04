module fibre_stage15_structure_advance_formula
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: n_points = 0
  real(mytype) :: dt_value = 0.0_mytype
  real(mytype) :: rho_tilde_value = 0.0_mytype
  real(mytype) :: test_force_magnitude = 0.0_mytype
  real(mytype) :: zero_force_drift = huge(1.0_mytype)
  integer :: zero_force_status = 0
  real(mytype) :: small_force_max_acceleration = huge(1.0_mytype)
  real(mytype) :: small_force_max_velocity_update = huge(1.0_mytype)
  real(mytype) :: small_force_max_position_update = huge(1.0_mytype)
  integer :: small_force_status = 0
  integer :: sign_consistency_status = 0
  real(mytype) :: force_scaling_error = huge(1.0_mytype)
  integer :: force_scaling_status = 0
  real(mytype) :: dt_scaling_error = huge(1.0_mytype)
  integer :: dt_scaling_status = 0
  integer :: finite_value_status = 0
  integer :: production_time_loop_connection_count = 0
  integer :: production_structure_advance_count = 0
  integer :: bending_solve_count = 0
  integer :: tension_solve_count = 0
  integer :: wall_contact_count = 0
  integer :: multifibre_count = 0
  integer :: rhs_injection_connection_count = 0
  integer :: no_fluid_rhs_modification_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: final_status = 0

  public :: stage15_structure_advance_formula_zero_force_check
  public :: stage15_structure_advance_formula_small_force_check
  public :: stage15_structure_advance_formula_sign_check
  public :: stage15_structure_advance_formula_gain_scaling_check
  public :: stage15_structure_advance_formula_dt_scaling_check
  public :: stage15_structure_advance_formula_get_status_values
  public :: stage15_structure_advance_formula_write_diagnostics
  public :: stage15_structure_advance_formula_finalize

contains

  subroutine stage15_structure_advance_formula_zero_force_check(npts, dt, rho_tilde, tolerance)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: tolerance
    real(mytype), allocatable :: x_old(:,:)
    real(mytype), allocatable :: v_old(:,:)
    real(mytype), allocatable :: force_test(:,:)
    real(mytype), allocatable :: a_cand(:,:)
    real(mytype), allocatable :: v_new(:,:)
    real(mytype), allocatable :: x_new(:,:)

    call initialize_common(npts, dt, rho_tilde, 0.0_mytype)
    if (.not. valid_inputs(npts, dt, rho_tilde)) then
      call update_final_status()
      return
    end if

    allocate(x_old(npts, 3), v_old(npts, 3), force_test(npts, 3), a_cand(npts, 3), &
             v_new(npts, 3), x_new(npts, 3))
    call initialize_reference_state(x_old, v_old)
    force_test(:, :) = 0.0_mytype
    call compute_candidate(x_old, v_old, force_test, dt, rho_tilde, a_cand, v_new, x_new)
    zero_force_drift = max(maxval(abs(a_cand)), max(maxval(abs(v_new - v_old)), maxval(abs(x_new - x_old))))
    zero_force_status = merge(1, 0, zero_force_drift <= tolerance .and. all_candidate_finite(a_cand, v_new, x_new))
    call update_final_status()
  end subroutine stage15_structure_advance_formula_zero_force_check

  subroutine stage15_structure_advance_formula_small_force_check(npts, dt, rho_tilde, force_magnitude, max_update)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: max_update
    real(mytype), allocatable :: x_old(:,:)
    real(mytype), allocatable :: v_old(:,:)
    real(mytype), allocatable :: force_test(:,:)
    real(mytype), allocatable :: a_cand(:,:)
    real(mytype), allocatable :: v_new(:,:)
    real(mytype), allocatable :: x_new(:,:)

    n_points = npts
    dt_value = dt
    rho_tilde_value = rho_tilde
    test_force_magnitude = force_magnitude
    small_force_status = 0
    if (.not. valid_inputs(npts, dt, rho_tilde)) then
      call update_final_status()
      return
    end if

    allocate(x_old(npts, 3), v_old(npts, 3), force_test(npts, 3), a_cand(npts, 3), &
             v_new(npts, 3), x_new(npts, 3))
    call initialize_reference_state(x_old, v_old)
    force_test(:, :) = 0.0_mytype
    force_test(:, 1) = force_magnitude
    call compute_candidate(x_old, v_old, force_test, dt, rho_tilde, a_cand, v_new, x_new)
    small_force_max_acceleration = maxval(abs(a_cand))
    small_force_max_velocity_update = maxval(abs(v_new - v_old))
    small_force_max_position_update = maxval(abs(x_new - x_old))
    small_force_status = merge(1, 0, all_candidate_finite(a_cand, v_new, x_new) .and. &
                               small_force_max_acceleration <= max_update .and. &
                               small_force_max_velocity_update <= max_update .and. &
                               small_force_max_position_update <= max_update)
    call update_final_status()
  end subroutine stage15_structure_advance_formula_small_force_check

  subroutine stage15_structure_advance_formula_sign_check(npts, dt, rho_tilde, force_magnitude, tolerance)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: tolerance
    real(mytype) :: expected_acceleration
    real(mytype) :: expected_velocity_update
    real(mytype) :: expected_position_update

    sign_consistency_status = 0
    if (.not. valid_inputs(npts, dt, rho_tilde)) then
      call update_final_status()
      return
    end if
    expected_acceleration = force_magnitude / rho_tilde
    expected_velocity_update = dt * expected_acceleration
    expected_position_update = dt * expected_velocity_update
    if (force_magnitude > 0.0_mytype) then
      sign_consistency_status = merge(1, 0, expected_acceleration >= -tolerance .and. &
                                      expected_velocity_update >= -tolerance .and. &
                                      expected_position_update >= -tolerance)
    else if (force_magnitude < 0.0_mytype) then
      sign_consistency_status = merge(1, 0, expected_acceleration <= tolerance .and. &
                                      expected_velocity_update <= tolerance .and. &
                                      expected_position_update <= tolerance)
    else
      sign_consistency_status = 1
    end if
    call update_final_status()
  end subroutine stage15_structure_advance_formula_sign_check

  subroutine stage15_structure_advance_formula_gain_scaling_check(npts, dt, rho_tilde, force_magnitude, tolerance)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: tolerance
    real(mytype) :: base_velocity_update
    real(mytype) :: double_velocity_update
    real(mytype) :: base_position_update
    real(mytype) :: double_position_update

    force_scaling_error = huge(1.0_mytype)
    force_scaling_status = 0
    if (.not. valid_inputs(npts, dt, rho_tilde) .or. abs(force_magnitude) <= 0.0_mytype) then
      call update_final_status()
      return
    end if
    base_velocity_update = dt * force_magnitude / rho_tilde
    double_velocity_update = dt * (2.0_mytype * force_magnitude) / rho_tilde
    base_position_update = dt * base_velocity_update
    double_position_update = dt * double_velocity_update
    force_scaling_error = max(abs(double_velocity_update - 2.0_mytype * base_velocity_update), &
                              abs(double_position_update - 2.0_mytype * base_position_update))
    force_scaling_status = merge(1, 0, force_scaling_error <= tolerance)
    call update_final_status()
  end subroutine stage15_structure_advance_formula_gain_scaling_check

  subroutine stage15_structure_advance_formula_dt_scaling_check(npts, dt, rho_tilde, force_magnitude, tolerance)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: tolerance
    real(mytype) :: base_velocity_update
    real(mytype) :: double_dt_velocity_update
    real(mytype) :: base_position_update
    real(mytype) :: double_dt_position_update

    dt_scaling_error = huge(1.0_mytype)
    dt_scaling_status = 0
    if (.not. valid_inputs(npts, dt, rho_tilde)) then
      call update_final_status()
      return
    end if
    base_velocity_update = dt * force_magnitude / rho_tilde
    double_dt_velocity_update = (2.0_mytype * dt) * force_magnitude / rho_tilde
    base_position_update = dt * base_velocity_update
    double_dt_position_update = (2.0_mytype * dt) * double_dt_velocity_update
    dt_scaling_error = max(abs(double_dt_velocity_update - 2.0_mytype * base_velocity_update), &
                           abs(double_dt_position_update - 4.0_mytype * base_position_update))
    dt_scaling_status = merge(1, 0, dt_scaling_error <= tolerance)
    call update_final_status()
  end subroutine stage15_structure_advance_formula_dt_scaling_check

  subroutine stage15_structure_advance_formula_get_status_values(npts_out, dt_out, rho_tilde_out, force_out, &
                                                                 zero_drift_out, zero_status_out, &
                                                                 small_accel_out, small_vel_out, small_pos_out, &
                                                                 small_status_out, sign_status_out, force_error_out, &
                                                                 force_status_out, dt_error_out, dt_status_out, &
                                                                 finite_status_out, time_loop_count_out, &
                                                                 structure_count_out, bending_count_out, tension_count_out, &
                                                                 wall_count_out, multifibre_count_out, rhs_count_out, &
                                                                 no_fluid_rhs_out, no_pressure_projection_out, &
                                                                 no_poisson_out, no_rk3_channel_out, final_out)
    integer, intent(out) :: npts_out
    real(mytype), intent(out) :: dt_out
    real(mytype), intent(out) :: rho_tilde_out
    real(mytype), intent(out) :: force_out
    real(mytype), intent(out) :: zero_drift_out
    integer, intent(out) :: zero_status_out
    real(mytype), intent(out) :: small_accel_out
    real(mytype), intent(out) :: small_vel_out
    real(mytype), intent(out) :: small_pos_out
    integer, intent(out) :: small_status_out
    integer, intent(out) :: sign_status_out
    real(mytype), intent(out) :: force_error_out
    integer, intent(out) :: force_status_out
    real(mytype), intent(out) :: dt_error_out
    integer, intent(out) :: dt_status_out
    integer, intent(out) :: finite_status_out
    integer, intent(out) :: time_loop_count_out
    integer, intent(out) :: structure_count_out
    integer, intent(out) :: bending_count_out
    integer, intent(out) :: tension_count_out
    integer, intent(out) :: wall_count_out
    integer, intent(out) :: multifibre_count_out
    integer, intent(out) :: rhs_count_out
    integer, intent(out) :: no_fluid_rhs_out
    integer, intent(out) :: no_pressure_projection_out
    integer, intent(out) :: no_poisson_out
    integer, intent(out) :: no_rk3_channel_out
    integer, intent(out) :: final_out

    call update_final_status()
    npts_out = n_points
    dt_out = dt_value
    rho_tilde_out = rho_tilde_value
    force_out = test_force_magnitude
    zero_drift_out = zero_force_drift
    zero_status_out = zero_force_status
    small_accel_out = small_force_max_acceleration
    small_vel_out = small_force_max_velocity_update
    small_pos_out = small_force_max_position_update
    small_status_out = small_force_status
    sign_status_out = sign_consistency_status
    force_error_out = force_scaling_error
    force_status_out = force_scaling_status
    dt_error_out = dt_scaling_error
    dt_status_out = dt_scaling_status
    finite_status_out = finite_value_status
    time_loop_count_out = production_time_loop_connection_count
    structure_count_out = production_structure_advance_count
    bending_count_out = bending_solve_count
    tension_count_out = tension_solve_count
    wall_count_out = wall_contact_count
    multifibre_count_out = multifibre_count
    rhs_count_out = rhs_injection_connection_count
    no_fluid_rhs_out = no_fluid_rhs_modification_status
    no_pressure_projection_out = no_pressure_projection_modification_status
    no_poisson_out = no_poisson_modification_status
    no_rk3_channel_out = no_rk3_channel_forcing_modification_status
    final_out = final_status
  end subroutine stage15_structure_advance_formula_get_status_values

  subroutine stage15_structure_advance_formula_write_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    call update_final_status()
    write(unit_id,'(A,1X,I0)') 'npts', n_points
    write(unit_id,'(A,1X,ES24.16)') 'dt', dt_value
    write(unit_id,'(A,1X,ES24.16)') 'rho_tilde', rho_tilde_value
    write(unit_id,'(A,1X,ES24.16)') 'test_force_magnitude', test_force_magnitude
    write(unit_id,'(A,1X,ES24.16)') 'zero_force_drift', zero_force_drift
    write(unit_id,'(A,1X,I0)') 'zero_force_status', zero_force_status
    write(unit_id,'(A,1X,ES24.16)') 'small_force_max_acceleration', small_force_max_acceleration
    write(unit_id,'(A,1X,ES24.16)') 'small_force_max_velocity_update', small_force_max_velocity_update
    write(unit_id,'(A,1X,ES24.16)') 'small_force_max_position_update', small_force_max_position_update
    write(unit_id,'(A,1X,I0)') 'small_force_status', small_force_status
    write(unit_id,'(A,1X,I0)') 'sign_consistency_status', sign_consistency_status
    write(unit_id,'(A,1X,ES24.16)') 'force_scaling_error', force_scaling_error
    write(unit_id,'(A,1X,I0)') 'force_scaling_status', force_scaling_status
    write(unit_id,'(A,1X,ES24.16)') 'dt_scaling_error', dt_scaling_error
    write(unit_id,'(A,1X,I0)') 'dt_scaling_status', dt_scaling_status
    write(unit_id,'(A,1X,I0)') 'finite_value_status', finite_value_status
    write(unit_id,'(A,1X,I0)') 'production_time_loop_connection_count', production_time_loop_connection_count
    write(unit_id,'(A,1X,I0)') 'production_structure_advance_count', production_structure_advance_count
    write(unit_id,'(A,1X,I0)') 'bending_solve_count', bending_solve_count
    write(unit_id,'(A,1X,I0)') 'tension_solve_count', tension_solve_count
    write(unit_id,'(A,1X,I0)') 'wall_contact_count', wall_contact_count
    write(unit_id,'(A,1X,I0)') 'multifibre_count', multifibre_count
    write(unit_id,'(A,1X,I0)') 'rhs_injection_connection_count', rhs_injection_connection_count
    write(unit_id,'(A,1X,I0)') 'no_fluid_rhs_modification_status', no_fluid_rhs_modification_status
    write(unit_id,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit_id,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit_id,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit_id,'(A,1X,I0)') 'final_status', final_status
  end subroutine stage15_structure_advance_formula_write_diagnostics

  subroutine stage15_structure_advance_formula_finalize()
    n_points = 0
    dt_value = 0.0_mytype
    rho_tilde_value = 0.0_mytype
    test_force_magnitude = 0.0_mytype
    zero_force_drift = huge(1.0_mytype)
    zero_force_status = 0
    small_force_max_acceleration = huge(1.0_mytype)
    small_force_max_velocity_update = huge(1.0_mytype)
    small_force_max_position_update = huge(1.0_mytype)
    small_force_status = 0
    sign_consistency_status = 0
    force_scaling_error = huge(1.0_mytype)
    force_scaling_status = 0
    dt_scaling_error = huge(1.0_mytype)
    dt_scaling_status = 0
    finite_value_status = 0
    production_time_loop_connection_count = 0
    production_structure_advance_count = 0
    bending_solve_count = 0
    tension_solve_count = 0
    wall_contact_count = 0
    multifibre_count = 0
    rhs_injection_connection_count = 0
    no_fluid_rhs_modification_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    final_status = 0
  end subroutine stage15_structure_advance_formula_finalize

  subroutine initialize_common(npts, dt, rho_tilde, force_magnitude)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude

    n_points = npts
    dt_value = dt
    rho_tilde_value = rho_tilde
    test_force_magnitude = force_magnitude
  end subroutine initialize_common

  logical function valid_inputs(npts, dt, rho_tilde)
    integer, intent(in) :: npts
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde

    valid_inputs = npts > 0 .and. is_finite(dt) .and. is_finite(rho_tilde) .and. &
                   dt > 0.0_mytype .and. rho_tilde > 0.0_mytype
  end function valid_inputs

  subroutine initialize_reference_state(x_old, v_old)
    real(mytype), intent(out) :: x_old(:,:)
    real(mytype), intent(out) :: v_old(:,:)
    integer :: i

    do i = 1, size(x_old, 1)
      if (size(x_old, 1) > 1) then
        x_old(i, 1) = real(i - 1, mytype) / real(size(x_old, 1) - 1, mytype)
      else
        x_old(i, 1) = 0.0_mytype
      end if
      x_old(i, 2) = 0.0_mytype
      x_old(i, 3) = 0.0_mytype
    end do
    v_old(:, :) = 0.0_mytype
  end subroutine initialize_reference_state

  subroutine compute_candidate(x_old, v_old, force_test, dt, rho_tilde, a_cand, v_new, x_new)
    real(mytype), intent(in) :: x_old(:,:)
    real(mytype), intent(in) :: v_old(:,:)
    real(mytype), intent(in) :: force_test(:,:)
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(out) :: a_cand(:,:)
    real(mytype), intent(out) :: v_new(:,:)
    real(mytype), intent(out) :: x_new(:,:)

    a_cand(:, :) = force_test(:, :) / rho_tilde
    v_new(:, :) = v_old(:, :) + dt * a_cand(:, :)
    x_new(:, :) = x_old(:, :) + dt * v_new(:, :)
  end subroutine compute_candidate

  logical function all_candidate_finite(a_cand, v_new, x_new)
    real(mytype), intent(in) :: a_cand(:,:)
    real(mytype), intent(in) :: v_new(:,:)
    real(mytype), intent(in) :: x_new(:,:)

    all_candidate_finite = all_finite_rank2(a_cand) .and. all_finite_rank2(v_new) .and. all_finite_rank2(x_new)
  end function all_candidate_finite

  subroutine update_final_status()
    finite_value_status = merge(1, 0, is_finite(zero_force_drift) .and. &
                                is_finite(small_force_max_acceleration) .and. &
                                is_finite(small_force_max_velocity_update) .and. &
                                is_finite(small_force_max_position_update) .and. &
                                is_finite(force_scaling_error) .and. is_finite(dt_scaling_error))
    final_status = merge(1, 0, zero_force_status == 1 .and. small_force_status == 1 .and. &
                         sign_consistency_status == 1 .and. force_scaling_status == 1 .and. &
                         dt_scaling_status == 1 .and. finite_value_status == 1 .and. &
                         production_time_loop_connection_count == 0 .and. &
                         production_structure_advance_count == 0 .and. bending_solve_count == 0 .and. &
                         tension_solve_count == 0 .and. wall_contact_count == 0 .and. multifibre_count == 0 .and. &
                         rhs_injection_connection_count == 0 .and. no_fluid_rhs_modification_status == 1 .and. &
                         no_pressure_projection_modification_status == 1 .and. &
                         no_poisson_modification_status == 1 .and. &
                         no_rk3_channel_forcing_modification_status == 1)
  end subroutine update_final_status

  logical function all_finite_rank2(values)
    real(mytype), intent(in) :: values(:,:)

    all_finite_rank2 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank2

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = value == value .and. abs(value) < huge(1.0_mytype)
  end function is_finite

end module fibre_stage15_structure_advance_formula
