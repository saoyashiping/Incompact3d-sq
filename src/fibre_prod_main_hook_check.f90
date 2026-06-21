program fibre_prod_main_hook_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, &
                                        fibre_prod_runtime_config_default, &
                                        fibre_prod_runtime_config_validate
  use fibre_prod_main_diagnostics, only : fibre_prod_rhs_signature_type, &
                                          fibre_prod_main_diagnostics_type, &
                                          fibre_prod_main_diagnostics_record, &
                                          fibre_prod_main_diagnostics_finite
  use fibre_prod_rhs_adapter, only : fibre_prod_rhs_adapter_apply
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_runtime_config_type) :: config
  type(fibre_prod_rhs_signature_type) :: before
  type(fibre_prod_rhs_signature_type) :: after
  type(fibre_prod_main_diagnostics_type) :: diag
  real(dp) :: rhs_x(4, 3, 2)
  real(dp) :: rhs_y(4, 3, 2)
  real(dp) :: rhs_z(4, 3, 2)
  real(dp) :: force_x(4, 3, 2)
  real(dp) :: force_y(4, 3, 2)
  real(dp) :: force_z(4, 3, 2)
  real(dp) :: before_sum
  integer :: modified_cells
  integer :: status

  call fibre_prod_runtime_config_default(config)
  if (config%enabled .or. config%lambda_fsi /= 0.0_dp .or. config%diagnostics_enabled) error stop 1
  if (fibre_prod_runtime_config_validate(config) /= 0) error stop 2

  rhs_x = 1.0_dp
  rhs_y = 2.0_dp
  rhs_z = 3.0_dp
  before_sum = sum(rhs_x) + sum(rhs_y) + sum(rhs_z)
  config%enabled = .true.
  config%lambda_fsi = 0.0_dp
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
  if (status /= 0 .or. modified_cells /= 0) error stop 3
  if (sum(rhs_x) + sum(rhs_y) + sum(rhs_z) /= before_sum) error stop 4
  call fibre_prod_main_diagnostics_record(diag, config%enabled, config%lambda_fsi, before, after, modified_cells)
  if (.not. diag%no_contamination) error stop 5

  ! P0.1: nonzero lambda without a physical Eulerian force-density buffer must be blocked.
  config%lambda_fsi = 1.0e-8_dp
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
  if (status /= 13 .or. modified_cells /= 0) error stop 6
  if (sum(rhs_x) + sum(rhs_y) + sum(rhs_z) /= before_sum) error stop 7

  ! P0.1: a finite physical force-density buffer is the only allowed nonzero-lambda RHS path.
  force_x = 0.0_dp
  force_y = 0.0_dp
  force_z = 0.0_dp
  force_x(2, 2, 1) = 1.0_dp
  force_y(2, 2, 1) = -0.5_dp
  force_z(2, 2, 1) = 0.25_dp
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status, &
                                    force_x, force_y, force_z)
  if (status /= 0 .or. modified_cells /= 1) error stop 8
  if (sum(rhs_x) + sum(rhs_y) + sum(rhs_z) <= before_sum) error stop 9
  call fibre_prod_main_diagnostics_record(diag, config%enabled, config%lambda_fsi, before, after, modified_cells)
  if (.not. diag%small_lambda_response) error stop 10

  config%lambda_fsi = -1.0_dp
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
  if (status == 0) error stop 11

  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.
  config%lambda_fsi = 1.0e-8_dp
  rhs_x(1, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status, &
                                    force_x, force_y, force_z)
  if (status == 0) error stop 12
  if (.not. fibre_prod_main_diagnostics_finite(diag)) error stop 13

  print *, 'P0_1_FIBRE_PROD_MAIN_HOOK_CHECK PASS'
end program fibre_prod_main_hook_check
