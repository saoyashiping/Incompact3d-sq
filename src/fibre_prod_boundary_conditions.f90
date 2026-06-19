module fibre_prod_boundary_conditions
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_bc_positions_finite
  public :: fibre_prod_bc_apply_free_free_bending
  public :: fibre_prod_bc_endpoint_tangent

contains

  pure logical function fibre_prod_bc_positions_finite(x) result(is_finite)
    real(dp), intent(in) :: x(:, :)

    is_finite = all(ieee_is_finite(x))
  end function fibre_prod_bc_positions_finite

  pure subroutine fibre_prod_bc_apply_free_free_bending(force)
    real(dp), intent(inout) :: force(:, :)
    integer :: nnode

    nnode = size(force, 1)
    if (nnode >= 1) force(1, :) = 0.0_dp
    if (nnode >= 2) force(2, :) = 0.0_dp
    if (nnode >= 3) force(nnode - 1, :) = 0.0_dp
    if (nnode >= 4) force(nnode, :) = 0.0_dp
  end subroutine fibre_prod_bc_apply_free_free_bending

  pure subroutine fibre_prod_bc_endpoint_tangent(x, tangent_start, tangent_end, status)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: tangent_start(3)
    real(dp), intent(out) :: tangent_end(3)
    integer, intent(out) :: status
    integer :: nnode

    nnode = size(x, 1)
    tangent_start = 0.0_dp
    tangent_end = 0.0_dp
    status = 0
    if (size(x, 2) /= 3 .or. nnode < 2 .or. .not. fibre_prod_bc_positions_finite(x)) then
      status = 1
      return
    end if

    tangent_start = x(2, :) - x(1, :)
    tangent_end = x(nnode, :) - x(nnode - 1, :)
  end subroutine fibre_prod_bc_endpoint_tangent

end module fibre_prod_boundary_conditions
