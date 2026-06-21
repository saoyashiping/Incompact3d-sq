program fibre_prod_runtime_bridge_check
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init
  use fibre_prod_runtime_bridge, only : fibre_prod_runtime_bridge_type, &
                                        fibre_prod_runtime_bridge_init_from_rhs, &
                                        fibre_prod_runtime_bridge_reset_force_buffer, &
                                        fibre_prod_runtime_bridge_apply_lambda0_noop, &
                                        fibre_prod_runtime_bridge_finalize
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_runtime_bridge_type) :: bridge
  type(fibre_prod_runtime_config_type) :: config
  real(dp) :: rhs_x(5, 6, 7)
  real(dp) :: rhs_y(5, 6, 7)
  real(dp) :: rhs_z(5, 6, 7)
  real(dp) :: rhs_bad(4, 6, 7)
  real(dp) :: base_x(5, 6, 7)
  real(dp) :: base_y(5, 6, 7)
  real(dp) :: base_z(5, 6, 7)
  integer :: status

  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  base_x = rhs_x
  base_y = rhs_y
  base_z = rhs_z
  call fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 1
  if (.not. bridge%initialized) error stop 2
  if (bridge%nx /= 5 .or. bridge%ny /= 6 .or. bridge%nz /= 7) error stop 3
  if (.not. bridge%force_buffer_allocated) error stop 4
  if (.not. allocated(bridge%force_buffer%fx) .or. .not. allocated(bridge%force_buffer%fy) .or. &
      .not. allocated(bridge%force_buffer%fz)) error stop 5
  if (any(shape(bridge%force_buffer%fx) /= shape(rhs_x))) error stop 6
  if (any(shape(bridge%force_buffer%fy) /= shape(rhs_y))) error stop 7
  if (any(shape(bridge%force_buffer%fz) /= shape(rhs_z))) error stop 8
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 9

  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.
  config%lambda_fsi = 0.0_dp
  config%penalty_beta = 2.0_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 10
  call fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 11
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 12

  config%lambda_fsi = 1.0e-3_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 13
  call fibre_prod_runtime_bridge_reset_force_buffer(bridge, status)
  if (status /= 0) error stop 14
  call fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 15
  if (.not. bridge%last_zero_force_buffer) error stop 16
  if (bridge%last_physical_response) error stop 17
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 18

  rhs_bad = 0.0_dp
  call fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_bad, bridge, status)
  if (status == 0) error stop 19
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 20

  call fibre_prod_runtime_bridge_finalize(bridge)
  print *, 'P0_3_RUNTIME_BRIDGE_CHECK PASS'
contains

  subroutine initialize_rhs(rhs_x, rhs_y, rhs_z)
    real(dp), intent(out) :: rhs_x(:, :, :)
    real(dp), intent(out) :: rhs_y(:, :, :)
    real(dp), intent(out) :: rhs_z(:, :, :)
    integer :: i, j, k

    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          rhs_x(i, j, k) = real(i + 10*j + 100*k, dp)
          rhs_y(i, j, k) = -real(2*i + j + k, dp)
          rhs_z(i, j, k) = real(i - 2*j + 3*k, dp)
        end do
      end do
    end do
  end subroutine initialize_rhs

  logical function same_field(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = all(lhs == rhs)
  end function same_field
end program fibre_prod_runtime_bridge_check
