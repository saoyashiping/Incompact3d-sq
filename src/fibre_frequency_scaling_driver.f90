program fibre_frequency_scaling_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_structure_integrator, only : advance_fibre_structure_freefree

  implicit none

  real(mytype) :: freq_gamma1_dt, freq_gamma1_dt_half, freq_gamma4_dt
  real(mytype) :: dt, dt_half
  real(mytype) :: freq_dt_relative_diff, freq_gamma_scaling_ratio, expected_gamma_scaling_ratio
  integer :: solver_failure_count_total, nan_detected
  integer :: io_unit

  dt = 1.0e-5_mytype
  dt_half = 0.5_mytype * dt
  solver_failure_count_total = 0
  nan_detected = 0

  call run_frequency_case(gamma=1.0_mytype, dt=dt, nsteps=30000, freq=freq_gamma1_dt, &
                          solver_failure_count=solver_failure_count_total, nan_detected=nan_detected)
  call run_frequency_case(gamma=1.0_mytype, dt=dt_half, nsteps=60000, freq=freq_gamma1_dt_half, &
                          solver_failure_count=solver_failure_count_total, nan_detected=nan_detected)
  call run_frequency_case(gamma=4.0_mytype, dt=dt, nsteps=30000, freq=freq_gamma4_dt, &
                          solver_failure_count=solver_failure_count_total, nan_detected=nan_detected)

  freq_dt_relative_diff = abs(freq_gamma1_dt_half - freq_gamma1_dt) / max(abs(freq_gamma1_dt_half), 1.0e-30_mytype)
  freq_gamma_scaling_ratio = freq_gamma4_dt / max(freq_gamma1_dt, 1.0e-30_mytype)
  expected_gamma_scaling_ratio = sqrt(4.0_mytype / 1.0_mytype)

  open(newunit=io_unit, file='stage2_outputs/fibre_frequency_scaling_check.dat', status='replace', action='write', form='formatted')
  write(io_unit, '(A,ES24.16)') 'freq_gamma1_dt = ', freq_gamma1_dt
  write(io_unit, '(A,ES24.16)') 'freq_gamma1_dt_half = ', freq_gamma1_dt_half
  write(io_unit, '(A,ES24.16)') 'freq_gamma4_dt = ', freq_gamma4_dt
  write(io_unit, '(A,ES24.16)') 'freq_dt_relative_diff = ', freq_dt_relative_diff
  write(io_unit, '(A,ES24.16)') 'freq_gamma_scaling_ratio = ', freq_gamma_scaling_ratio
  write(io_unit, '(A,ES24.16)') 'expected_gamma_scaling_ratio = ', expected_gamma_scaling_ratio
  write(io_unit, '(A,I0)') 'solver_failure_count_total = ', solver_failure_count_total
  write(io_unit, '(A,I0)') 'nan_detected = ', nan_detected
  close(io_unit)

contains

  subroutine run_frequency_case(gamma, dt, nsteps, freq, solver_failure_count, nan_detected)
    real(mytype), intent(in) :: gamma, dt
    integer, intent(in) :: nsteps
    real(mytype), intent(out) :: freq
    integer, intent(inout) :: solver_failure_count, nan_detected

    type(fibre_t) :: fibre
    integer :: i, step, ierr, imid
    real(mytype) :: s, pi, amp, tension_residual, relative_tension_residual
    real(mytype) :: y_prev2, y_prev1, y_curr, t_prev1, t_curr
    real(mytype) :: peak_times(32), period_sum
    integer :: npeaks, nperiod

    pi = acos(-1.0_mytype)
    amp = 1.0e-3_mytype

    call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=gamma)

    do i = 1, fibre%nl
      s = real(i - 1, mytype) * fibre%ds
      fibre%x(2, i) = amp * sin(2.0_mytype * pi * s / fibre%length)
    end do
    fibre%x_old = fibre%x
    fibre%v = 0._mytype

    imid = fibre%nl / 2 + 1
    y_prev2 = fibre%x(2, imid)
    y_prev1 = y_prev2
    t_prev1 = 0._mytype
    npeaks = 0

    do step = 1, nsteps
      call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
      if (ierr /= 0) solver_failure_count = solver_failure_count + 1

      if (any(ieee_is_nan(fibre%x)) .or. any(ieee_is_nan(fibre%v)) .or. any(ieee_is_nan(fibre%tension))) then
        nan_detected = 1
      end if

      y_curr = fibre%x(2, imid)
      t_curr = real(step, mytype) * dt

      if (y_prev1 > y_prev2 .and. y_prev1 >= y_curr) then
        if (npeaks < size(peak_times)) then
          npeaks = npeaks + 1
          peak_times(npeaks) = t_prev1
        end if
      end if

      y_prev2 = y_prev1
      y_prev1 = y_curr
      t_prev1 = t_curr
    end do

    if (npeaks >= 3) then
      period_sum = 0._mytype
      nperiod = 0
      do i = 2, npeaks
        if (peak_times(i) > peak_times(i - 1)) then
          period_sum = period_sum + (peak_times(i) - peak_times(i - 1))
          nperiod = nperiod + 1
        end if
      end do
      if (nperiod > 0) then
        freq = real(nperiod, mytype) / period_sum
      else
        freq = 0._mytype
      end if
    else
      freq = 0._mytype
    end if
  end subroutine run_frequency_case

end program fibre_frequency_scaling_check
