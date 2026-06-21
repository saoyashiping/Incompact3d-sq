module fibre_prod_ibm_delta
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_delta_kernel_1d
  public :: fibre_prod_delta_weight_3d
  public :: fibre_prod_delta_inputs_finite
  public :: fibre_prod_delta_compact_support_radius

contains

  pure real(dp) function fibre_prod_delta_compact_support_radius(h) result(radius)
    real(dp), intent(in) :: h

    if (ieee_is_finite(h) .and. h > 0.0_dp) then
      radius = 2.0_dp * h
    else
      radius = -1.0_dp
    end if
  end function fibre_prod_delta_compact_support_radius

  pure logical function fibre_prod_delta_inputs_finite(distance, h) result(is_valid)
    real(dp), intent(in) :: distance
    real(dp), intent(in) :: h

    is_valid = ieee_is_finite(distance) .and. ieee_is_finite(h) .and. h > 0.0_dp
  end function fibre_prod_delta_inputs_finite

  pure real(dp) function fibre_prod_delta_kernel_1d(distance, h) result(weight)
    real(dp), intent(in) :: distance
    real(dp), intent(in) :: h
    real(dp) :: q
    real(dp), parameter :: pi = acos(-1.0_dp)

    weight = 0.0_dp
    if (.not. fibre_prod_delta_inputs_finite(distance, h)) return

    q = abs(distance) / h
    if (q < 2.0_dp) then
      weight = 0.25_dp * (1.0_dp + cos(0.5_dp * pi * q)) / h
    end if
  end function fibre_prod_delta_kernel_1d

  pure real(dp) function fibre_prod_delta_weight_3d(dx, dy, dz, hx, hy, hz) result(weight)
    real(dp), intent(in) :: dx
    real(dp), intent(in) :: dy
    real(dp), intent(in) :: dz
    real(dp), intent(in) :: hx
    real(dp), intent(in) :: hy
    real(dp), intent(in) :: hz

    weight = fibre_prod_delta_kernel_1d(dx, hx) * &
             fibre_prod_delta_kernel_1d(dy, hy) * &
             fibre_prod_delta_kernel_1d(dz, hz)
  end function fibre_prod_delta_weight_3d

end module fibre_prod_ibm_delta
