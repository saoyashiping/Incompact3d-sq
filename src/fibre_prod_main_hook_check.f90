program fibre_prod_main_hook_check
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, &
                                        fibre_prod_runtime_config_default, &
                                        fibre_prod_runtime_config_validate
  use fibre_prod_main_diagnostics, only : fibre_prod_main_diagnostics_type
  use fibre_prod_rhs_adapter, only : fibre_prod_rhs_adapter_apply, &
                                     fibre_prod_rhs_status_missing_force_buffer
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init, &
                                   fibre_prod_main_hook_apply, &
                                   fibre_prod_main_hook_apply_force_buffer, &
                                   fibre_prod_main_hook_get_diagnostics
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: tol = 0.0_dp
  type(fibre_prod_runtime_config_type) :: config
  type(fibre_prod_main_diagnostics_type) :: diag
  real(dp) :: rhs_x(4, 3, 2)
  real(dp) :: rhs_y(4, 3, 2)
  real(dp) :: rhs_z(4, 3, 2)
  real(dp) :: base_x(4, 3, 2)
  real(dp) :: base_y(4, 3, 2)
  real(dp) :: base_z(4, 3, 2)
  real(dp) :: force_x(4, 3, 2)
  real(dp) :: force_y(4, 3, 2)
  real(dp) :: force_z(4, 3, 2)
  type(fibre_prod_force_buffer_type) :: force_buffer
  real(dp) :: force_scale
  real(dp) :: max_abs_increment
  real(dp) :: sum_increment
  logical :: zero_force_buffer
  logical :: missing_force_buffer
  integer :: status
  integer :: modified_cells
  integer :: ierr

  call fibre_prod_runtime_config_default(config)
  if (config%enabled .or. config%lambda_fsi /= 0.0_dp .or. config%diagnostics_enabled) error stop 1
  if (fibre_prod_runtime_config_validate(config) /= 0) error stop 2

  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  base_x = rhs_x
  base_y = rhs_y
  base_z = rhs_z
  config%enabled = .true.
  config%lambda_fsi = 0.0_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 3
  call fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status)
  if (status /= 0) error stop 4
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 5
  call fibre_prod_main_hook_get_diagnostics(diag)
  if (.not. diag%no_contamination .or. diag%modified_cells /= 0) error stop 6

  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  base_x = rhs_x
  base_y = rhs_y
  base_z = rhs_z
  config%lambda_fsi = 1.0e-8_dp
  config%penalty_beta = 2.0_dp
  force_scale = config%lambda_fsi * config%penalty_beta
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 7
  call fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status)
  if (status /= fibre_prod_rhs_status_missing_force_buffer) error stop 8
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 9
  call fibre_prod_main_hook_get_diagnostics(diag)
  if (diag%modified_cells /= 0 .or. diag%last_status /= fibre_prod_rhs_status_missing_force_buffer) error stop 10

  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  base_x = rhs_x
  base_y = rhs_y
  base_z = rhs_z
  call initialize_force(force_x, force_y, force_z)
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 11
  call fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status, force_x, force_y, force_z)
  if (status /= 0) error stop 12
  if (.not. same_field(rhs_x, base_x + force_scale * force_x)) error stop 13
  if (.not. same_field(rhs_y, base_y + force_scale * force_y)) error stop 14
  if (.not. same_field(rhs_z, base_z + force_scale * force_z)) error stop 15
  if (is_uniform_increment(force_x, force_y, force_z)) error stop 16
  call fibre_prod_main_hook_get_diagnostics(diag)
  if (diag%modified_cells <= 0 .or. .not. diag%small_lambda_response) error stop 17

  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  base_x = rhs_x
  base_y = rhs_y
  base_z = rhs_z
  force_x = 0.0_dp
  force_y = 0.0_dp
  force_z = 0.0_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 18
  call fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status, force_x, force_y, force_z)
  if (status /= 0) error stop 19
  if (.not. same_field(rhs_x, base_x) .or. .not. same_field(rhs_y, base_y) .or. .not. same_field(rhs_z, base_z)) error stop 20
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, diag%before, diag%after, modified_cells, status, &
                                    force_x, force_y, force_z, max_abs_increment, sum_increment, &
                                    zero_force_buffer, missing_force_buffer)
  if (status /= 0 .or. modified_cells /= 0) error stop 21
  if (.not. zero_force_buffer .or. missing_force_buffer) error stop 22
  if (max_abs_increment /= 0.0_dp .or. sum_increment /= 0.0_dp) error stop 23

  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  rhs_x = 0.0_dp
  rhs_y = 0.0_dp
  rhs_z = 0.0_dp
  base_x = rhs_x
  base_y = rhs_y
  base_z = rhs_z
  allocate(force_buffer%fx(4, 3, 2), force_buffer%fy(4, 3, 2), force_buffer%fz(4, 3, 2), stat=ierr)
  if (ierr /= 0) error stop 24
  force_buffer%nx_local = 4
  force_buffer%ny_local = 3
  force_buffer%nz_local = 2
  force_buffer%allocated = .true.
  call initialize_force(force_buffer%fx, force_buffer%fy, force_buffer%fz)
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 25
  call fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, force_buffer, status)
  if (status /= 0) error stop 26
  if (.not. same_field(rhs_x, base_x + force_scale * force_buffer%fx)) error stop 27
  if (.not. same_field(rhs_y, base_y + force_scale * force_buffer%fy)) error stop 28
  if (.not. same_field(rhs_z, base_z + force_scale * force_buffer%fz)) error stop 29
  if (is_uniform_increment(rhs_x, rhs_y, rhs_z)) error stop 30
  deallocate(force_buffer%fx, force_buffer%fy, force_buffer%fz)
  force_buffer%allocated = .false.

  print *, 'P0_1_FIBRE_PROD_MAIN_HOOK_CHECK PASS'
  print *, 'P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS'
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
          rhs_y(i, j, k) = -real(i + j + k, dp)
          rhs_z(i, j, k) = real(2*i - j + k, dp)
        end do
      end do
    end do
  end subroutine initialize_rhs

  subroutine initialize_force(force_x, force_y, force_z)
    real(dp), intent(out) :: force_x(:, :, :)
    real(dp), intent(out) :: force_y(:, :, :)
    real(dp), intent(out) :: force_z(:, :, :)
    integer :: i, j, k

    do k = 1, size(force_x, 3)
      do j = 1, size(force_x, 2)
        do i = 1, size(force_x, 1)
          force_x(i, j, k) = 0.01_dp * real(i, dp)
          force_y(i, j, k) = -0.02_dp * real(j, dp)
          force_z(i, j, k) = 0.03_dp * real(k + i, dp)
        end do
      end do
    end do
  end subroutine initialize_force

  logical function same_field(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = all(abs(lhs - rhs) <= tol)
  end function same_field

  logical function is_uniform_increment(force_x, force_y, force_z) result(is_uniform)
    real(dp), intent(in) :: force_x(:, :, :)
    real(dp), intent(in) :: force_y(:, :, :)
    real(dp), intent(in) :: force_z(:, :, :)

    is_uniform = all(force_x == force_x(1, 1, 1)) .and. all(force_y == force_y(1, 1, 1)) .and. &
                 all(force_z == force_z(1, 1, 1)) .and. force_x(1, 1, 1) == force_y(1, 1, 1) .and. &
                 force_y(1, 1, 1) == force_z(1, 1, 1)
  end function is_uniform_increment

end program fibre_prod_main_hook_check
