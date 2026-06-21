module fibre_prod_runtime_config
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_runtime_config_type
  public :: fibre_prod_runtime_config_default
  public :: fibre_prod_runtime_config_validate
  public :: fibre_prod_runtime_config_load_env

  type :: fibre_prod_runtime_config_type
    logical :: enabled = .false.
    real(dp) :: lambda_fsi = 0.0_dp
    real(dp) :: penalty_beta = 1.0_dp
    real(dp) :: dt = 1.0_dp
    logical :: diagnostics_enabled = .false.
    logical :: initialized = .false.
  end type fibre_prod_runtime_config_type

contains

  subroutine fibre_prod_runtime_config_default(config)
    type(fibre_prod_runtime_config_type), intent(out) :: config

    config%enabled = .false.
    config%lambda_fsi = 0.0_dp
    config%penalty_beta = 1.0_dp
    config%dt = 1.0_dp
    config%diagnostics_enabled = .false.
    config%initialized = .true.
  end subroutine fibre_prod_runtime_config_default

  integer function fibre_prod_runtime_config_validate(config) result(status)
    type(fibre_prod_runtime_config_type), intent(in) :: config

    status = 0
    if (.not. config%initialized) status = 1
    if (status == 0 .and. (.not. ieee_is_finite(config%lambda_fsi) .or. &
        config%lambda_fsi < 0.0_dp)) status = 2
    if (status == 0 .and. (.not. ieee_is_finite(config%penalty_beta) .or. &
        config%penalty_beta < 0.0_dp)) status = 3
    if (status == 0 .and. (.not. ieee_is_finite(config%dt) .or. config%dt <= 0.0_dp)) status = 4
  end function fibre_prod_runtime_config_validate

  subroutine fibre_prod_runtime_config_load_env(config, status)
    type(fibre_prod_runtime_config_type), intent(out) :: config
    integer, intent(out) :: status

    call fibre_prod_runtime_config_default(config)
    call read_env_logical('FIBRE_PROD_ENABLE', config%enabled)
    call read_env_real('FIBRE_PROD_LAMBDA', config%lambda_fsi)
    call read_env_real('FIBRE_PROD_PENALTY_BETA', config%penalty_beta)
    call read_env_logical('FIBRE_PROD_DIAGNOSTICS', config%diagnostics_enabled)
    status = fibre_prod_runtime_config_validate(config)
  end subroutine fibre_prod_runtime_config_load_env

  subroutine read_env_logical(name, value)
    character(len=*), intent(in) :: name
    logical, intent(inout) :: value
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    call get_environment_variable(name, raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      value = .true.
    case ('0', 'F', 'f', 'FALSE', 'false', 'False', 'NO', 'no', 'OFF', 'off')
      value = .false.
    end select
  end subroutine read_env_logical

  subroutine read_env_real(name, value)
    character(len=*), intent(in) :: name
    real(dp), intent(inout) :: value
    character(len=128) :: raw
    integer :: length
    integer :: env_status
    integer :: read_status
    real(dp) :: parsed

    call get_environment_variable(name, raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    read(raw(1:min(length, len(raw))), *, iostat=read_status) parsed
    if (read_status == 0) value = parsed
  end subroutine read_env_real

end module fibre_prod_runtime_config
