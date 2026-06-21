module fibre_prod_rhs_adapter
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_validate
  use fibre_prod_main_diagnostics, only : fibre_prod_rhs_signature_type, fibre_prod_main_signature
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_rhs_adapter_apply

contains

  subroutine fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
    type(fibre_prod_runtime_config_type), intent(in) :: config
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    type(fibre_prod_rhs_signature_type), intent(out) :: before
    type(fibre_prod_rhs_signature_type), intent(out) :: after
    integer, intent(out) :: modified_cells
    integer, intent(out) :: status
    integer :: i
    integer :: j
    integer :: k
    real(dp) :: contribution

    modified_cells = 0
    status = fibre_prod_runtime_config_validate(config)
    before = fibre_prod_main_signature(rhs_x, rhs_y, rhs_z)
    after = before
    if (status /= 0) return
    if (size(rhs_y, 1) /= size(rhs_x, 1) .or. size(rhs_y, 2) /= size(rhs_x, 2) .or. &
        size(rhs_y, 3) /= size(rhs_x, 3) .or. size(rhs_z, 1) /= size(rhs_x, 1) .or. &
        size(rhs_z, 2) /= size(rhs_x, 2) .or. size(rhs_z, 3) /= size(rhs_x, 3)) then
      status = 9
      return
    end if
    if (.not. before%finite) then
      status = 10
      return
    end if
    if (.not. config%enabled .or. config%lambda_fsi == 0.0_dp) return
    contribution = config%lambda_fsi * config%penalty_beta * config%dt
    if (.not. ieee_is_finite(contribution)) then
      status = 11
      return
    end if
    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          rhs_x(i, j, k) = rhs_x(i, j, k) + contribution
          modified_cells = modified_cells + 1
        end do
      end do
    end do
    after = fibre_prod_main_signature(rhs_x, rhs_y, rhs_z)
    if (.not. after%finite) status = 12
  end subroutine fibre_prod_rhs_adapter_apply

end module fibre_prod_rhs_adapter
