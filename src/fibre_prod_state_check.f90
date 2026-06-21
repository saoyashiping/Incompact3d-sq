program fibre_prod_state_check
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_config, only : fibre_prod_config_default, fibre_prod_config_is_disabled, &
                                fibre_prod_config_type
  use fibre_prod_diagnostics, only : fibre_prod_diagnostics_collect, &
                                     fibre_prod_state_check_summary
  use fibre_prod_state, only : fibre_prod_state_type, &
                               fibre_prod_state_allocate, &
                               fibre_prod_state_init_straight, &
                               fibre_prod_state_reset_forces, &
                               fibre_prod_state_destroy, &
                               fibre_prod_state_is_allocated, &
                               fibre_prod_state_all_finite, &
                               fibre_prod_state_segment_length_residual, &
                               fibre_prod_state_total_force_norm
  implicit none

  type(fibre_prod_config_type) :: config
  type(fibre_prod_state_type) :: state
  type(fibre_prod_state_check_summary) :: summary
  integer :: ierr
  real(real64), parameter :: tolerance = 1.0e-12_real64
  real(real64), parameter :: segment_length = 0.125_real64
  real(real64) :: origin(3)
  real(real64) :: direction(3)
  real(real64) :: residual
  real(real64) :: force_norm
  real(real64) :: all_force_norm

  config = fibre_prod_config_default()
  if (.not. fibre_prod_config_is_disabled(config)) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: default config does not disable fibre/FSI'
    error stop 1
  end if

  call fibre_prod_state_allocate(state, nfibre=2, nnode=5, stat=ierr)
  if (ierr /= 0 .or. .not. fibre_prod_state_is_allocated(state)) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: allocation failed'
    error stop 2
  end if

  origin = [0.0_real64, 0.0_real64, 0.0_real64]
  direction = [1.0_real64, 0.0_real64, 0.0_real64]
  call fibre_prod_state_init_straight(state, origin, direction, segment_length, stat=ierr)
  if (ierr /= 0) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: straight fibre initialization failed'
    error stop 3
  end if

  if (.not. fibre_prod_state_all_finite(state)) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: state contains non-finite values'
    error stop 4
  end if

  residual = fibre_prod_state_segment_length_residual(state, segment_length)
  if (residual > tolerance) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: segment length residual too large', residual
    error stop 5
  end if

  state%f_fs = 1.0_real64
  state%f_wall = -2.0_real64
  state%f_coll = 3.0_real64
  state%f_total = state%f_fs + state%f_wall + state%f_coll
  call fibre_prod_state_reset_forces(state)
  force_norm = fibre_prod_state_total_force_norm(state)
  all_force_norm = sqrt(sum(state%f_fs * state%f_fs) + sum(state%f_wall * state%f_wall) + &
                        sum(state%f_coll * state%f_coll) + sum(state%f_total * state%f_total))
  if (force_norm /= 0.0_real64 .or. all_force_norm /= 0.0_real64) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: force reset did not zero force arrays', &
             force_norm, all_force_norm
    error stop 6
  end if

  summary = fibre_prod_diagnostics_collect(state, segment_length)
  if (.not. summary%all_finite .or. summary%segment_length_residual > tolerance .or. &
      summary%total_force_norm /= 0.0_real64) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: diagnostics summary mismatch'
    error stop 7
  end if

  call fibre_prod_state_destroy(state)
  if (fibre_prod_state_is_allocated(state)) then
    print *, 'R2_FIBRE_PROD_STATE_CHECK FAIL: destroy did not release arrays'
    error stop 8
  end if

  print *, 'R2_FIBRE_PROD_STATE_CHECK PASS'
end program fibre_prod_state_check
