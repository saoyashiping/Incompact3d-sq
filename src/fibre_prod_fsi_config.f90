module fibre_prod_fsi_config
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_fsi_config_type
  public :: fibre_prod_fsi_config_default
  public :: fibre_prod_fsi_config_validate

  type :: fibre_prod_fsi_config_type
    logical :: fsi_enabled = .false.
    real(dp) :: lambda_fsi = 0.0_dp
    real(dp) :: penalty_beta = 1.0_dp
    real(dp) :: dt = 1.0e-3_dp
    real(dp) :: rho_tilde = 1.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: ds = 1.0_dp
    integer :: nstep = 1
  end type fibre_prod_fsi_config_type

contains

  pure function fibre_prod_fsi_config_default() result(config)
    type(fibre_prod_fsi_config_type) :: config

    config%fsi_enabled = .false.
    config%lambda_fsi = 0.0_dp
    config%penalty_beta = 1.0_dp
    config%dt = 1.0e-3_dp
    config%rho_tilde = 1.0_dp
    config%gamma = 0.0_dp
    config%ds = 1.0_dp
    config%nstep = 1
  end function fibre_prod_fsi_config_default

  pure integer function fibre_prod_fsi_config_validate(config) result(status)
    type(fibre_prod_fsi_config_type), intent(in) :: config

    status = 0
    if (.not. ieee_is_finite(config%lambda_fsi) .or. config%lambda_fsi < 0.0_dp) status = 1
    if (status == 0 .and. (.not. ieee_is_finite(config%penalty_beta) .or. &
        config%penalty_beta < 0.0_dp)) status = 2
    if (status == 0 .and. (.not. ieee_is_finite(config%dt) .or. config%dt <= 0.0_dp)) status = 3
    if (status == 0 .and. (.not. ieee_is_finite(config%rho_tilde) .or. &
        config%rho_tilde <= 0.0_dp)) status = 4
    if (status == 0 .and. (.not. ieee_is_finite(config%gamma) .or. config%gamma < 0.0_dp)) status = 5
    if (status == 0 .and. (.not. ieee_is_finite(config%ds) .or. config%ds <= 0.0_dp)) status = 6
    if (status == 0 .and. config%nstep <= 0) status = 7
  end function fibre_prod_fsi_config_validate

end module fibre_prod_fsi_config
