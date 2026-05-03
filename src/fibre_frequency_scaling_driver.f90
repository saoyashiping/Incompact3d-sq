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

  real(mytype) :: base_signal_initial, base_signal_min, base_signal_max, base_signal_rms_after_transient
  real(mytype) :: gamma4_signal_initial, gamma4_signal_min, gamma4_signal_max, gamma4_signal_rms_after_transient
  real(mytype) :: rho4_signal_initial, rho4_signal_min, rho4_signal_max, rho4_signal_rms_after_transient
  real(mytype) :: length2_signal_initial, length2_signal_min, length2_signal_max, length2_signal_rms_after_transient

  integer :: base_nan_detected, base_solver_failure_count
  integer :: gamma4_nan_detected, gamma4_solver_failure_count
  integer :: rho4_nan_detected, rho4_solver_failure_count
  integer :: length2_nan_detected, length2_solver_failure_count

  integer :: io_summary

  call run_case(case_name='base', output_file='stage2_outputs/frequency_base.dat', nl=nl, &
                length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=base_initial_length_error, final_length_error=base_final_length_error, &
                nan_detected=base_nan_detected, solver_failure_count=base_solver_failure_count, &
                signal_initial=base_signal_initial, signal_min=base_signal_min, signal_max=base_signal_max, &
                signal_rms_after_transient=base_signal_rms_after_transient)

  call run_case(case_name='gamma4', output_file='stage2_outputs/frequency_gamma4.dat', nl=nl, &
                length=1.0_mytype, rho_tilde=1.0_mytype, gamma=4.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=gamma4_initial_length_error, final_length_error=gamma4_final_length_error, &
                nan_detected=gamma4_nan_detected, solver_failure_count=gamma4_solver_failure_count, &
                signal_initial=gamma4_signal_initial, signal_min=gamma4_signal_min, signal_max=gamma4_signal_max, &
                signal_rms_after_transient=gamma4_signal_rms_after_transient)

  call run_case(case_name='rho4', output_file='stage2_outputs/frequency_rho4.dat', nl=nl, &
                length=1.0_mytype, rho_tilde=4.0_mytype, gamma=1.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=rho4_initial_length_error, final_length_error=rho4_final_length_error, &
                nan_detected=rho4_nan_detected, solver_failure_count=rho4_solver_failure_count, &
                signal_initial=rho4_signal_initial, signal_min=rho4_signal_min, signal_max=rho4_signal_max, &
                signal_rms_after_transient=rho4_signal_rms_after_transient)

  call run_case(case_name='length2', output_file='stage2_outputs/frequency_length2.dat', nl=nl, &
                length=2.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype, &
                steps_per_period=steps_per_period, n_periods=n_periods, &
                initial_length_error=length2_initial_length_error, final_length_error=length2_final_length_error, &
                nan_detected=length2_nan_detected, solver_failure_count=length2_solver_failure_count, &
                signal_initial=length2_signal_initial, signal_min=length2_signal_min, signal_max=length2_signal_max, &
                signal_rms_after_transient=length2_signal_rms_after_transient)

  open(newunit=io_summary, file='stage2_outputs/fibre_frequency_scaling_driver_summary.dat', &
       status='replace', action='write', form='formatted')

  write(io_summary, '(A,ES24.16)') 'base_initial_length_error = ', base_initial_length_error
  write(io_summary, '(A,ES24.16)') 'base_final_length_error = ', base_final_length_error
  write(io_summary, '(A,I0)') 'base_nan_detected = ', base_nan_detected
  write(io_summary, '(A,I0)') 'base_solver_failure_count = ', base_solver_failure_count
  write(io_summary, '(A,ES24.16)') 'base_signal_initial = ', base_signal_initial
  write(io_summary, '(A,ES24.16)') 'base_signal_min = ', base_signal_min
  write(io_summary, '(A,ES24.16)') 'base_signal_max = ', base_signal_max
  write(io_summary, '(A,ES24.16)') 'base_signal_rms_after_transient = ', base_signal_rms_after_transient

  write(io_summary, '(A,ES24.16)') 'gamma4_initial_length_error = ', gamma4_initial_length_error
  write(io_summary, '(A,ES24.16)') 'gamma4_final_length_error = ', gamma4_final_length_error
  write(io_summary, '(A,I0)') 'gamma4_nan_detected = ', gamma4_nan_detected
  write(io_summary, '(A,I0)') 'gamma4_solver_failure_count = ', gamma4_solver_failure_count
  write(io_summary, '(A,ES24.16)') 'gamma4_signal_initial = ', gamma4_signal_initial
  write(io_summary, '(A,ES24.16)') 'gamma4_signal_min = ', gamma4_signal_min
  write(io_summary, '(A,ES24.16)') 'gamma4_signal_max = ', gamma4_signal_max
  write(io_summary, '(A,ES24.16)') 'gamma4_signal_rms_after_transient = ', gamma4_signal_rms_after_transient

  write(io_summary, '(A,ES24.16)') 'rho4_initial_length_error = ', rho4_initial_length_error
  write(io_summary, '(A,ES24.16)') 'rho4_final_length_error = ', rho4_final_length_error
  write(io_summary, '(A,I0)') 'rho4_nan_detected = ', rho4_nan_detected
  write(io_summary, '(A,I0)') 'rho4_solver_failure_count = ', rho4_solver_failure_count
  write(io_summary, '(A,ES24.16)') 'rho4_signal_initial = ', rho4_signal_initial
  write(io_summary, '(A,ES24.16)') 'rho4_signal_min = ', rho4_signal_min
  write(io_summary, '(A,ES24.16)') 'rho4_signal_max = ', rho4_signal_max
  write(io_summary, '(A,ES24.16)') 'rho4_signal_rms_after_transient = ', rho4_signal_rms_after_transient

  write(io_summary, '(A,ES24.16)') 'length2_initial_length_error = ', length2_initial_length_error
  write(io_summary, '(A,ES24.16)') 'length2_final_length_error = ', length2_final_length_error
  write(io_summary, '(A,I0)') 'length2_nan_detected = ', length2_nan_detected
  write(io_summary, '(A,I0)') 'length2_solver_failure_count = ', length2_solver_failure_count
  write(io_summary, '(A,ES24.16)') 'length2_signal_initial = ', length2_signal_initial
  write(io_summary, '(A,ES24.16)') 'length2_signal_min = ', length2_signal_min
  write(io_summary, '(A,ES24.16)') 'length2_signal_max = ', length2_signal_max
  write(io_summary, '(A,ES24.16)') 'length2_signal_rms_after_transient = ', length2_signal_rms_after_transient

  close(io_summary)

contains

  subroutine run_case(case_name, output_file, nl, length, rho_tilde, gamma, steps_per_period, n_periods, &
                      initial_length_error, final_length_error, nan_detected, solver_failure_count, &
                      signal_initial, signal_min, signal_max, signal_rms_after_transient)
    character(len=*), intent(in) :: case_name, output_file
    integer, intent(in) :: nl, steps_per_period, n_periods
    real(mytype), intent(in) :: length, rho_tilde, gamma
    real(mytype), intent(out) :: initial_length_error, final_length_error
    integer, intent(out) :: nan_detected, solver_failure_count
    real(mytype), intent(out) :: signal_initial, signal_min, signal_max, signal_rms_after_transient

    type(fibre_t) :: fibre
    integer :: i, j, nseg, nsteps, step, io_ts, ierr
    integer :: transient_start, rms_count
    real(mytype) :: ds, theta0, s_mid, theta
    real(mytype) :: f_ref, t_ref, dt, time
    real(mytype) :: tension_residual, relative_tension_residual
    real(mytype) :: bending_energy, kinetic_energy, total_energy, length_error
    real(mytype) :: signal, momentum_norm, sigma, mode_norm, mode_norm_sq
    real(mytype) :: signal_rms_sumsq
    real(mytype) :: x_left0(3), x_right0(3), chord0(3)
    real(mytype) :: x_left(3), x_right(3), chord(3)
    real(mytype) :: center(3), momentum(3)
    real(mytype) :: x_chord0(3, nl), r0(3, nl), x_chord(3, nl), r(3, nl), w(nl)

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

    x_left0 = fibre%x(:, 1)
    x_right0 = fibre%x(:, fibre%nl)
    chord0 = x_right0 - x_left0

    w = ds
    w(1) = 0.5_mytype * ds
    w(fibre%nl) = 0.5_mytype * ds

    do i = 1, fibre%nl
      sigma = real(i - 1, mytype) / real(fibre%nl - 1, mytype)
      x_chord0(:, i) = x_left0 + sigma * chord0
      r0(:, i) = fibre%x(:, i) - x_chord0(:, i)
    end do

    mode_norm_sq = 0._mytype
    do i = 1, fibre%nl
      mode_norm_sq = mode_norm_sq + w(i) * dot_product(r0(:, i), r0(:, i))
    end do
    mode_norm = sqrt(mode_norm_sq)

    if (mode_norm <= 1.0e-10_mytype) then
      error stop 'fibre_frequency_scaling_check: mode_norm too small for robust frequency extraction'
    end if

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
    transient_start = int(0.2_mytype * real(nsteps + 1, mytype))
    signal_rms_sumsq = 0._mytype
    rms_count = 0

    do step = 0, nsteps
      time = real(step, mytype) * dt

      x_left = fibre%x(:, 1)
      x_right = fibre%x(:, fibre%nl)
      chord = x_right - x_left

      do i = 1, fibre%nl
        sigma = real(i - 1, mytype) / real(fibre%nl - 1, mytype)
        x_chord(:, i) = x_left + sigma * chord
        r(:, i) = fibre%x(:, i) - x_chord(:, i)
      end do

      signal = 0._mytype
      do i = 1, fibre%nl
        signal = signal + w(i) * dot_product(r(:, i), r0(:, i))
      end do
      signal = signal / mode_norm

      if (step == 0) then
        signal_initial = signal
        signal_min = signal
        signal_max = signal
      else
        signal_min = min(signal_min, signal)
        signal_max = max(signal_max, signal)
      end if

      if (step >= transient_start) then
        signal_rms_sumsq = signal_rms_sumsq + signal * signal
        rms_count = rms_count + 1
      end if

      call compute_bending_energy(fibre, bending_energy)
      call compute_kinetic_energy(fibre, kinetic_energy)
      total_energy = bending_energy + kinetic_energy
      call compute_total_length_relative_error(fibre, length_error)
      call compute_linear_momentum(fibre, momentum)
      momentum_norm = sqrt(sum(momentum**2))
      center = fibre_center_of_mass(fibre)

      write(io_ts, '(10(ES24.16,1X))') time, signal, bending_energy, kinetic_energy, total_energy, &
                                       length_error, momentum_norm, center(1), center(2), center(3)

      if (any(ieee_is_nan(fibre%x)) .or. any(ieee_is_nan(fibre%v)) .or. any(ieee_is_nan(fibre%tension)) .or. &
          ieee_is_nan(bending_energy) .or. ieee_is_nan(kinetic_energy) .or. ieee_is_nan(length_error) .or. &
          ieee_is_nan(signal)) then
        nan_detected = 1
      end if

      if (step < nsteps) then
        call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
        if (ierr /= 0) solver_failure_count = solver_failure_count + 1
      end if
    end do

    close(io_ts)

    if (rms_count > 0) then
      signal_rms_after_transient = sqrt(signal_rms_sumsq / real(rms_count, mytype))
    else
      signal_rms_after_transient = 0._mytype
    end if

    call compute_total_length_relative_error(fibre, final_length_error)
  end subroutine run_case

end program fibre_frequency_scaling_check
