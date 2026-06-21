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

  subroutine fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status, &
                                          force_x, force_y, force_z)
    type(fibre_prod_runtime_config_type), intent(in) :: config
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    type(fibre_prod_rhs_signature_type), intent(out) :: before
    type(fibre_prod_rhs_signature_type), intent(out) :: after
    integer, intent(out) :: modified_cells
    integer, intent(out) :: status
    real(dp), intent(in), optional :: force_x(:, :, :)
    real(dp), intent(in), optional :: force_y(:, :, :)
    real(dp), intent(in), optional :: force_z(:, :, :)
    integer :: i
    integer :: j
    integer :: k
    real(dp) :: scale
    real(dp) :: add_x
    real(dp) :: add_y
    real(dp) :: add_z
    logical :: any_force_argument

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

    ! Disabled and lambda=0 modes must be strict no-ops. This is the safety gate
    ! needed before any production DNS-FSI path is allowed to run.
    if (.not. config%enabled .or. config%lambda_fsi == 0.0_dp) return

    ! P0.1 removes the previous whole-domain scalar RHS perturbation. A nonzero
    ! lambda is now accepted only when a physical Eulerian force-density buffer is
    ! supplied explicitly. This prevents a smoke-test perturbation from being
    ! mistaken for IBM spreading or Lagrangian fibre reaction forcing.
    any_force_argument = present(force_x) .or. present(force_y) .or. present(force_z)
    if (.not. any_force_argument) then
      status = 13
      return
    end if
    if (.not. (present(force_x) .and. present(force_y) .and. present(force_z))) then
      status = 14
      return
    end if
    if (size(force_x, 1) /= size(rhs_x, 1) .or. size(force_x, 2) /= size(rhs_x, 2) .or. &
        size(force_x, 3) /= size(rhs_x, 3) .or. size(force_y, 1) /= size(rhs_x, 1) .or. &
        size(force_y, 2) /= size(rhs_x, 2) .or. size(force_y, 3) /= size(rhs_x, 3) .or. &
        size(force_z, 1) /= size(rhs_x, 1) .or. size(force_z, 2) /= size(rhs_x, 2) .or. &
        size(force_z, 3) /= size(rhs_x, 3)) then
      status = 15
      return
    end if
    if (.not. all(ieee_is_finite(force_x)) .or. .not. all(ieee_is_finite(force_y)) .or. &
        .not. all(ieee_is_finite(force_z))) then
      status = 16
      return
    end if

    scale = config%lambda_fsi * config%penalty_beta
    if (.not. ieee_is_finite(scale)) then
      status = 11
      return
    end if
    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          add_x = scale * force_x(i, j, k)
          add_y = scale * force_y(i, j, k)
          add_z = scale * force_z(i, j, k)
          if (add_x /= 0.0_dp .or. add_y /= 0.0_dp .or. add_z /= 0.0_dp) then
            rhs_x(i, j, k) = rhs_x(i, j, k) + add_x
            rhs_y(i, j, k) = rhs_y(i, j, k) + add_y
            rhs_z(i, j, k) = rhs_z(i, j, k) + add_z
            modified_cells = modified_cells + 1
          end if
        end do
      end do
    end do
    after = fibre_prod_main_signature(rhs_x, rhs_y, rhs_z)
    if (.not. after%finite) status = 12
  end subroutine fibre_prod_rhs_adapter_apply

end module fibre_prod_rhs_adapter
