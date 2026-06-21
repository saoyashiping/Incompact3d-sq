program fibre_prod_synthetic_closed_loop_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_synthetic_closed_loop, only : fibre_prod_synthetic_closed_loop_type, &
                                               fibre_prod_synthetic_closed_loop_init, &
                                               fibre_prod_synthetic_closed_loop_run, &
                                               fibre_prod_synthetic_closed_loop_get_signature, &
                                               fibre_prod_synthetic_closed_loop_finalize
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_synthetic_closed_loop_type) :: loop0
  type(fibre_prod_synthetic_closed_loop_type) :: loop1
  type(fibre_prod_synthetic_closed_loop_type) :: loop2
  type(fibre_prod_synthetic_closed_loop_type) :: loop_beta1
  type(fibre_prod_synthetic_closed_loop_type) :: loop_zero_force_proxy
  real(dp) :: sig0(16)
  real(dp) :: sig1(16)
  real(dp) :: sig2(16)
  real(dp) :: sig_beta1(16)
  integer :: status

  call run_case(loop0, 0.0_dp, 2.0_dp, sig0)
  call require(loop0%completed, 'lambda0 full path did not complete')
  call require(all(ieee_is_finite(sig0)), 'lambda0 signature is nonfinite')
  call require(loop0%max_abs_force_buffer > 0.0_dp, 'lambda0 full path did not create force buffer candidate')
  call require(loop0%max_abs_rhs_increment == 0.0_dp, 'lambda0 full path changed RHS')

  call run_case(loop1, 1.0e-3_dp, 2.0_dp, sig1)
  call require(loop1%completed, 'small-lambda full path did not complete')
  call require(all(ieee_is_finite(sig1)), 'small-lambda signature is nonfinite')
  call require(loop1%max_abs_sampled_u > 0.0_dp, 'sampled velocity is zero')
  call require(loop1%max_abs_hydro_candidate > 0.0_dp, 'hydro candidate is zero')
  call require(loop1%max_abs_structure_input > 0.0_dp, 'structure input is zero')
  call require(loop1%max_abs_dx_trial > 0.0_dp .and. loop1%max_abs_dx_trial <= loop1%max_allowed_dx, &
               'dry-step diagnostics failed')
  call require(loop1%max_abs_reaction_force > 0.0_dp, 'reaction force is zero')
  call require(loop1%max_abs_force_buffer > 0.0_dp, 'force buffer is zero')
  call require(loop1%max_abs_rhs_increment > 0.0_dp .and. loop1%max_abs_rhs_increment < 1.0e3_dp, &
               'small-lambda RHS increment failed bounded/nonzero check')
  call require(loop1%sum_abs_rhs_increment > 0.0_dp, 'small-lambda RHS increment sum is zero')

  call run_case(loop2, 2.0e-3_dp, 2.0_dp, sig2)
  call require(abs(sig2(12) - 2.0_dp * sig1(12)) <= 1.0e-10_dp, 'lambda signature scaling failed')
  call require(abs(sig2(13) - 2.0_dp * sig1(13)) <= 1.0e-10_dp, 'lambda max increment scaling failed')

  call run_case(loop_beta1, 1.0e-3_dp, 1.0_dp, sig_beta1)
  call require(abs(sig1(12) - 2.0_dp * sig_beta1(12)) <= 1.0e-10_dp, 'penalty-beta signature scaling failed')
  call require(abs(sig1(13) - 2.0_dp * sig_beta1(13)) <= 1.0e-10_dp, 'penalty-beta max increment scaling failed')

  call fibre_prod_synthetic_closed_loop_init(loop_zero_force_proxy, 4, 4, 4, 5, status)
  call require(status == 0, 'zero-force proxy init failed')
  call fibre_prod_synthetic_closed_loop_run(loop_zero_force_proxy, 1.0e-3_dp, 2.0_dp, 0.0_dp, 1.0e-4_dp, 1.0_dp, status)
  call require(status == 0, 'zero-force proxy run failed')
  call require(loop_zero_force_proxy%max_abs_force_buffer == 0.0_dp .and. &
               loop_zero_force_proxy%max_abs_rhs_increment == 0.0_dp, 'zero force buffer generated RHS response')

  call require(abs(sig1(16)) < 1.0e2_dp, 'conservation/scale diagnostic exceeded synthetic bound')
  call require(sig1(12) > 0.0_dp .and. sig1(13) > 0.0_dp, 'RHS increment signature is not positive')
  call require(sig1(11) > 0.0_dp, 'force-buffer signature is not positive')

  call print_signature('lambda0', sig0)
  call print_signature('smalllambda', sig1)
  print '(A,ES24.16)', 'P0_12_SIGNATURE sampled_norm=', sig1(1)
  print '(A,ES24.16)', 'P0_12_SIGNATURE force_buffer_norm=', sig1(11)
  print '(A,ES24.16)', 'P0_12_SIGNATURE rhs_increment_norm=', sig1(12)
  print '(A,ES24.16)', 'P0_12_SIGNATURE lambda0_rhs_increment_norm=', sig0(12)
  print '(A,ES24.16)', 'P0_12_SIGNATURE smalllambda_rhs_increment_norm=', sig1(12)
  print '(A,ES24.16)', 'P0_12_SIGNATURE rhs_increment_max=', sig1(13)
  print '(A,ES24.16)', 'P0_12_SIGNATURE conservation_scale_error=', sig1(16)
  print *, 'P0_12_DIAGNOSTIC: synthetic MPI consistency uses deterministic rank-independent setup'
  print *, 'P0_12_SYNTHETIC_CLOSED_LOOP_CHECK PASS'

  call fibre_prod_synthetic_closed_loop_finalize(loop0)
  call fibre_prod_synthetic_closed_loop_finalize(loop1)
  call fibre_prod_synthetic_closed_loop_finalize(loop2)
  call fibre_prod_synthetic_closed_loop_finalize(loop_beta1)
  call fibre_prod_synthetic_closed_loop_finalize(loop_zero_force_proxy)

contains

  subroutine run_case(loop, lambda_fsi, penalty_beta, signature)
    type(fibre_prod_synthetic_closed_loop_type), intent(inout) :: loop
    real(dp), intent(in) :: lambda_fsi
    real(dp), intent(in) :: penalty_beta
    real(dp), intent(out) :: signature(16)
    integer :: local_status

    call fibre_prod_synthetic_closed_loop_init(loop, 4, 4, 4, 5, local_status)
    call require(local_status == 0, 'closed-loop init failed')
    call fibre_prod_synthetic_closed_loop_run(loop, lambda_fsi, penalty_beta, 1.25_dp, 1.0e-4_dp, 2.0_dp, local_status)
    call require(local_status == 0, 'closed-loop run failed')
    call fibre_prod_synthetic_closed_loop_get_signature(loop, signature, local_status)
    call require(local_status == 0, 'closed-loop signature failed')
  end subroutine run_case

  subroutine print_signature(label, signature)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: signature(16)
    integer :: i

    do i = 1, 16
      print '(A,A,A,I2.2,A,ES24.16)', 'P0_12_SIGNATURE ', trim(label), '_s', i, '=', signature(i)
    end do
  end subroutine print_signature

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      print *, 'P0_12_SYNTHETIC_CLOSED_LOOP_CHECK FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine require

end program fibre_prod_synthetic_closed_loop_check
