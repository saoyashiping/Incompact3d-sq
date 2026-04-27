program fibre_frequency_scaling_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_geometry, only : fibre_center_of_mass
  use fibre_diagnostics, only : compute_bending_energy, compute_kinetic_energy
  use fibre_diagnostics, only : compute_total_length_relative_error, compute_linear_momentum
  use fibre_structure_integrator, only : advance_fibre_structure_freefree

  implicit none

  integer, parameter :: nl = 65
  integer, parameter :: steps_per_period = 500
  integer, parameter :: n_periods = 10

  real(mytype), parameter :: pi = acos(-1.0_mytype)
  real(mytype), parameter :: alpha = 22.4_mytype / pi

  real(mytype) :: base_initial_length_error, base_final_length_error
  real(mytype) :: gamma4_initial_length_error, gamma4_final_length_error
  real(mytype) :: rho4_initial_length_error, rho4_final_length_error
  real(mytype) :: length2_initial_length_error, length2_final_length_error

  integer :: base_nan_detected, base_solver_failure_count
  integer :: gamma4_nan_detected, gamma4_solver_failure_count
  integer :: rho4_nan_detected, rho4_solver_failure_count
  integer :: length2_nan_detected, length2_solver_failure_count

  integer :: io_summary

  call run_case(case_name='base', output_file='stage2_outputs/frequency_base.dat', nl=nl, &
                length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=base_initial_length_error, final_length_error=base_final_length_error, &
                nan_detected=base_nan_detected, solver_failure_count=base_solver_failure_count)

  call run_case(case_name='gamma4', output_file='stage2_outputs/frequency_gamma4.dat', nl=nl, &
                length=1.0_mytype, rho_tilde=1.0_mytype, gamma=4.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=gamma4_initial_length_error, final_length_error=gamma4_final_length_error, &
                nan_detected=gamma4_nan_detected, solver_failure_count=gamma4_solver_failure_count)

  call run_case(case_name='rho4', output_file='stage2_outputs/frequency_rho4.dat', nl=nl, &
                length=1.0_mytype, rho_tilde=4.0_mytype, gamma=1.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=rho4_initial_length_error, final_length_error=rho4_final_length_error, &
                nan_detected=rho4_nan_detected, solver_failure_count=rho4_solver_failure_count)

  call run_case(case_name='length2', output_file='stage2_outputs/frequency_length2.dat', nl=nl, &
                length=2.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=length2_initial_length_error, final_length_error=length2_final_length_error, &
                nan_detected=length2_nan_detected, solver_failure_count=length2_solver_failure_count)

  open(newunit=io_summary, file='stage2_outputs/fibre_frequency_scaling_driver_summary.dat', &
       status='replace', action='write', form='formatted')

  write(io_summary, '(A,ES24.16)') 'base_initial_length_error = ', base_initial_length_error
  write(io_summary, '(A,ES24.16)') 'base_final_length_error = ', base_final_length_error
  write(io_summary, '(A,I0)') 'base_nan_detected = ', base_nan_detected
  write(io_summary, '(A,I0)') 'base_solver_failure_count = ', base_solver_failure_count

  write(io_summary, '(A,ES24.16)') 'gamma4_initial_length_error = ', gamma4_initial_length_error
  write(io_summary, '(A,ES24.16)') 'gamma4_final_length_error = ', gamma4_final_length_error
  write(io_summary, '(A,I0)') 'gamma4_nan_detected = ', gamma4_nan_detected
  write(io_summary, '(A,I0)') 'gamma4_solver_failure_count = ', gamma4_solver_failure_count

  write(io_summary, '(A,ES24.16)') 'rho4_initial_length_error = ', rho4_initial_length_error
  write(io_summary, '(A,ES24.16)') 'rho4_final_length_error = ', rho4_final_length_error
  write(io_summary, '(A,I0)') 'rho4_nan_detected = ', rho4_nan_detected
  write(io_summary, '(A,I0)') 'rho4_solver_failure_count = ', rho4_solver_failure_count

  write(io_summary, '(A,ES24.16)') 'length2_initial_length_error = ', length2_initial_length_error
  write(io_summary, '(A,ES24.16)') 'length2_final_length_error = ', length2_final_length_error
  write(io_summary, '(A,I0)') 'length2_nan_detected = ', length2_nan_detected
  write(io_summary, '(A,I0)') 'length2_solver_failure_count = ', length2_solver_failure_count

  close(io_summary)

contains

  subroutine run_case(case_name, output_file, nl, length, rho_tilde, gamma, steps_per_period, n_periods, &
                      initial_length_error, final_length_error, nan_detected, solver_failure_count)
    character(len=*), intent(in) :: case_name, output_file
    integer, intent(in) :: nl, steps_per_period, n_periods
    real(mytype), intent(in) :: length, rho_tilde, gamma
    real(mytype), intent(out) :: initial_length_error, final_length_error
    integer, intent(out) :: nan_detected, solver_failure_count

    type(fibre_t) :: fibre
    integer :: i, j, imid, nseg, nsteps, step, io_ts, ierr
    real(mytype) :: ds, theta0, s_mid, theta
    real(mytype) :: f_ref, t_ref, dt, time
    real(mytype) :: tension_residual, relative_tension_residual
    real(mytype) :: bending_energy, kinetic_energy, total_energy, length_error
    real(mytype) :: signal, momentum_norm
    real(mytype) :: center(3), momentum(3)

    call fibre_init_straight_free_free(fibre, nl=nl, length=length, rho_tilde=rho_tilde, gamma=gamma)

    nseg = fibre%nl - 1
    ds = length / real(nseg, mytype)
    theta0 = 0.03_mytype

    fibre%x(:, 1) = 0._mytype
    do j = 1, nseg
      s_mid = (real(j, mytype) - 0.5_mytype) * ds
      theta = theta0 * sin(pi * s_mid / length)
      fibre%x(1, j + 1) = fibre%x(1, j) + ds * cos(theta)
      fibre%x(2, j + 1) = fibre%x(2, j) + ds * sin(theta)
      fibre%x(3, j + 1) = fibre%x(3, j) + 0._mytype
    end do

    fibre%x_old = fibre%x
    fibre%v = 0._mytype

    f_ref = alpha * sqrt(gamma / (rho_tilde * length**4))
    t_ref = 1.0_mytype / f_ref
    dt = t_ref / real(steps_per_period, mytype)
    nsteps = n_periods * steps_per_period

    call compute_total_length_relative_error(fibre, initial_length_error)

    open(newunit=io_ts, file=output_file, status='replace', action='write', form='formatted')
    write(io_ts, '(A)') '# case ' // trim(case_name)
    write(io_ts, '(A,I0)') '# nl ', fibre%nl
    write(io_ts, '(A,ES24.16)') '# length ', length
    write(io_ts, '(A,ES24.16)') '# rho_tilde ', rho_tilde
    write(io_ts, '(A,ES24.16)') '# gamma ', gamma
    write(io_ts, '(A,ES24.16)') '# dt ', dt
    write(io_ts, '(A,I0)') '# nsteps ', nsteps
    write(io_ts, '(A,ES24.16)') '# f_ref ', f_ref
    write(io_ts, '(A)') '# time signal bending_energy kinetic_energy total_energy length_error momentum_norm center_x center_y center_z'

    solver_failure_count = 0
    nan_detected = 0
    imid = (fibre%nl + 1) / 2

    do step = 0, nsteps
      time = real(step, mytype) * dt

      call compute_bending_energy(fibre, bending_energy)
      call compute_kinetic_energy(fibre, kinetic_energy)
      total_energy = bending_energy + kinetic_energy
      call compute_total_length_relative_error(fibre, length_error)
      call compute_linear_momentum(fibre, momentum)
      momentum_norm = sqrt(sum(momentum**2))
      center = fibre_center_of_mass(fibre)
      signal = fibre%x(2, imid) - center(2)

      write(io_ts, '(10(ES24.16,1X))') time, signal, bending_energy, kinetic_energy, total_energy, &
                                       length_error, momentum_norm, center(1), center(2), center(3)

      if (any(ieee_is_nan(fibre%x)) .or. any(ieee_is_nan(fibre%v)) .or. any(ieee_is_nan(fibre%tension)) .or. &
          ieee_is_nan(bending_energy) .or. ieee_is_nan(kinetic_energy) .or. ieee_is_nan(length_error)) then
        nan_detected = 1
      end if

      if (step < nsteps) then
        call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
        if (ierr /= 0) solver_failure_count = solver_failure_count + 1
      end if
    end do

    close(io_ts)

    call compute_total_length_relative_error(fibre, final_length_error)
  end subroutine run_case

end program fibre_frequency_scaling_check
