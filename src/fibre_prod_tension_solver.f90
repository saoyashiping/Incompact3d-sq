module fibre_prod_tension_solver
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_tension_segment_length_residual
  public :: fibre_prod_tension_max_stretch_error
  public :: fibre_prod_tension_force_diagnostic_candidate

contains

  real(dp) function fibre_prod_tension_segment_length_residual(x, ds, status) result(residual)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: ds
    integer, intent(out) :: status
    integer :: i
    real(dp) :: segment(3)
    real(dp) :: length

    status = 0
    residual = huge(1.0_dp)
    if (size(x, 2) /= 3 .or. size(x, 1) < 2) then
      status = 1
    else if (.not. ieee_is_finite(ds) .or. ds <= 0.0_dp) then
      status = 2
    else if (.not. all(ieee_is_finite(x))) then
      status = 3
    else
      residual = 0.0_dp
      do i = 1, size(x, 1) - 1
        segment = x(i + 1, :) - x(i, :)
        length = sqrt(sum(segment * segment))
        residual = max(residual, abs(length - ds))
      end do
    end if
  end function fibre_prod_tension_segment_length_residual

  real(dp) function fibre_prod_tension_max_stretch_error(x, ds, status) result(max_error)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: ds
    integer, intent(out) :: status
    real(dp) :: residual

    residual = fibre_prod_tension_segment_length_residual(x, ds, status)
    if (status == 0) then
      max_error = residual / ds
    else
      max_error = huge(1.0_dp)
    end if
  end function fibre_prod_tension_max_stretch_error

  subroutine fibre_prod_tension_force_diagnostic_candidate(x, ds, stiffness, force, status)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: ds
    real(dp), intent(in) :: stiffness
    real(dp), intent(out) :: force(:, :)
    integer, intent(out) :: status
    integer :: i
    real(dp) :: segment(3)
    real(dp) :: length
    real(dp) :: magnitude

    force = 0.0_dp
    status = 0
    if (size(x, 2) /= 3 .or. size(force, 1) /= size(x, 1) .or. size(force, 2) /= 3) then
      status = 1
    else if (size(x, 1) < 2) then
      status = 2
    else if (.not. ieee_is_finite(ds) .or. ds <= 0.0_dp) then
      status = 3
    else if (.not. ieee_is_finite(stiffness) .or. stiffness < 0.0_dp) then
      status = 4
    else if (.not. all(ieee_is_finite(x))) then
      status = 5
    else
      do i = 1, size(x, 1) - 1
        segment = x(i + 1, :) - x(i, :)
        length = sqrt(sum(segment * segment))
        if (length > 0.0_dp) then
          magnitude = stiffness * (length - ds)
          force(i, :) = force(i, :) + magnitude * segment / length
          force(i + 1, :) = force(i + 1, :) - magnitude * segment / length
        end if
      end do
    end if
  end subroutine fibre_prod_tension_force_diagnostic_candidate

end module fibre_prod_tension_solver
