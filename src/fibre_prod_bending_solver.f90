module fibre_prod_bending_solver
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_boundary_conditions, only : fibre_prod_bc_apply_free_free_bending, &
                                             fibre_prod_bc_positions_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_bending_compute_force
  public :: fibre_prod_bending_force_norm

contains

  subroutine fibre_prod_bending_compute_force(x, gamma, ds, force, status)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: ds
    real(dp), intent(out) :: force(:, :)
    integer, intent(out) :: status
    integer :: nnode
    integer :: i

    status = 0
    force = 0.0_dp
    nnode = size(x, 1)
    if (size(x, 2) /= 3 .or. size(force, 1) /= nnode .or. size(force, 2) /= 3) then
      status = 1
    else if (nnode < 5) then
      status = 2
    else if (.not. ieee_is_finite(gamma) .or. gamma < 0.0_dp) then
      status = 3
    else if (.not. ieee_is_finite(ds) .or. ds <= 0.0_dp) then
      status = 4
    else if (.not. fibre_prod_bc_positions_finite(x)) then
      status = 5
    else
      do i = 3, nnode - 2
        force(i, :) = -gamma * (x(i - 2, :) - 4.0_dp * x(i - 1, :) + 6.0_dp * x(i, :) - &
                                4.0_dp * x(i + 1, :) + x(i + 2, :)) / ds**4
      end do
      call fibre_prod_bc_apply_free_free_bending(force)
    end if
  end subroutine fibre_prod_bending_compute_force

  pure real(dp) function fibre_prod_bending_force_norm(force) result(norm_value)
    real(dp), intent(in) :: force(:, :)

    if (all(ieee_is_finite(force))) then
      norm_value = sqrt(sum(force * force))
    else
      norm_value = huge(1.0_dp)
    end if
  end function fibre_prod_bending_force_norm

end module fibre_prod_bending_solver
