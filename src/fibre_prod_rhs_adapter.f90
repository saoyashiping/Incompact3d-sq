module fibre_prod_rhs_adapter
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_validate
  use fibre_prod_main_diagnostics, only : fibre_prod_rhs_signature_type, fibre_prod_main_signature
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: fibre_prod_rhs_status_ok = 0
  integer, parameter, public :: fibre_prod_rhs_status_missing_force_buffer = 13
  integer, parameter, public :: fibre_prod_rhs_status_force_shape_mismatch = 14
  integer, parameter, public :: fibre_prod_rhs_status_nonfinite_force_buffer = 15

  public :: fibre_prod_rhs_adapter_apply

contains

  subroutine fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status, &
                                          force_x, force_y, force_z, max_abs_increment, sum_increment, &
                                          zero_force_buffer, missing_force_buffer)
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
    real(dp), intent(out), optional :: max_abs_increment
    real(dp), intent(out), optional :: sum_increment
    logical, intent(out), optional :: zero_force_buffer
    logical, intent(out), optional :: missing_force_buffer
    integer :: i
    integer :: j
    integer :: k
    real(dp) :: local_max_abs_increment
    real(dp) :: local_sum_increment
    logical :: has_force_buffer
    logical :: local_zero_force_buffer
    logical :: local_missing_force_buffer

    modified_cells = 0
    local_max_abs_increment = 0.0_dp
    local_sum_increment = 0.0_dp
    local_zero_force_buffer = .false.
    local_missing_force_buffer = .false.
    status = fibre_prod_runtime_config_validate(config)
    before = fibre_prod_main_signature(rhs_x, rhs_y, rhs_z)
    after = before
    if (present(max_abs_increment)) max_abs_increment = local_max_abs_increment
    if (present(sum_increment)) sum_increment = local_sum_increment
    if (present(zero_force_buffer)) zero_force_buffer = local_zero_force_buffer
    if (present(missing_force_buffer)) missing_force_buffer = local_missing_force_buffer
    if (status /= fibre_prod_rhs_status_ok) return
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

    has_force_buffer = present(force_x) .and. present(force_y) .and. present(force_z)
    if (.not. has_force_buffer) then
      local_missing_force_buffer = .true.
      status = fibre_prod_rhs_status_missing_force_buffer
      if (present(missing_force_buffer)) missing_force_buffer = local_missing_force_buffer
      return
    end if

    if (.not. same_shape(rhs_x, force_x) .or. .not. same_shape(rhs_y, force_y) .or. &
        .not. same_shape(rhs_z, force_z)) then
      status = fibre_prod_rhs_status_force_shape_mismatch
      return
    end if

    if (.not. all(ieee_is_finite(force_x)) .or. .not. all(ieee_is_finite(force_y)) .or. &
        .not. all(ieee_is_finite(force_z))) then
      status = fibre_prod_rhs_status_nonfinite_force_buffer
      return
    end if

    local_zero_force_buffer = all(force_x == 0.0_dp) .and. all(force_y == 0.0_dp) .and. all(force_z == 0.0_dp)
    if (local_zero_force_buffer) then
      if (present(zero_force_buffer)) zero_force_buffer = local_zero_force_buffer
      return
    end if

    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          if (force_x(i, j, k) /= 0.0_dp .or. force_y(i, j, k) /= 0.0_dp .or. force_z(i, j, k) /= 0.0_dp) then
            modified_cells = modified_cells + 1
          end if
          rhs_x(i, j, k) = rhs_x(i, j, k) + force_x(i, j, k)
          rhs_y(i, j, k) = rhs_y(i, j, k) + force_y(i, j, k)
          rhs_z(i, j, k) = rhs_z(i, j, k) + force_z(i, j, k)
          local_sum_increment = local_sum_increment + force_x(i, j, k) + force_y(i, j, k) + force_z(i, j, k)
          local_max_abs_increment = max(local_max_abs_increment, abs(force_x(i, j, k)), &
                                        abs(force_y(i, j, k)), abs(force_z(i, j, k)))
        end do
      end do
    end do
    after = fibre_prod_main_signature(rhs_x, rhs_y, rhs_z)
    if (.not. after%finite) status = 12
    if (present(max_abs_increment)) max_abs_increment = local_max_abs_increment
    if (present(sum_increment)) sum_increment = local_sum_increment
    if (present(zero_force_buffer)) zero_force_buffer = local_zero_force_buffer
    if (present(missing_force_buffer)) missing_force_buffer = local_missing_force_buffer
  end subroutine fibre_prod_rhs_adapter_apply

  pure logical function same_shape(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = size(lhs, 1) == size(rhs, 1) .and. size(lhs, 2) == size(rhs, 2) .and. &
              size(lhs, 3) == size(rhs, 3)
  end function same_shape

end module fibre_prod_rhs_adapter
