module fibre_prod_diagnostics
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_state, only : fibre_prod_state_type, &
                               fibre_prod_state_all_finite, &
                               fibre_prod_state_segment_length_residual, &
                               fibre_prod_state_total_force_norm
  implicit none
  private

  public :: fibre_prod_state_check_summary

  type :: fibre_prod_state_check_summary
    logical :: all_finite = .false.
    real(real64) :: segment_length_residual = huge(1.0_real64)
    real(real64) :: total_force_norm = huge(1.0_real64)
  end type fibre_prod_state_check_summary

  public :: fibre_prod_diagnostics_collect

contains

  pure function fibre_prod_diagnostics_collect(state, expected_segment_length) result(summary)
    type(fibre_prod_state_type), intent(in) :: state
    real(real64), intent(in) :: expected_segment_length
    type(fibre_prod_state_check_summary) :: summary

    summary%all_finite = fibre_prod_state_all_finite(state)
    summary%segment_length_residual = fibre_prod_state_segment_length_residual(state, expected_segment_length)
    summary%total_force_norm = fibre_prod_state_total_force_norm(state)
  end function fibre_prod_diagnostics_collect

end module fibre_prod_diagnostics
