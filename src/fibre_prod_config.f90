module fibre_prod_config
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  private

  public :: fibre_prod_config_type
  public :: fibre_prod_config_default
  public :: fibre_prod_config_is_disabled

  type :: fibre_prod_config_type
    logical :: fibre_enabled = .false.
    logical :: fsi_enabled = .false.
    integer :: nfibre = 0
    integer :: nnode = 0
    integer :: ndim = 3
    real(real64) :: segment_length = 0.0_real64
  end type fibre_prod_config_type

contains

  pure function fibre_prod_config_default() result(config)
    type(fibre_prod_config_type) :: config

    config%fibre_enabled = .false.
    config%fsi_enabled = .false.
    config%nfibre = 0
    config%nnode = 0
    config%ndim = 3
    config%segment_length = 0.0_real64
  end function fibre_prod_config_default

  pure logical function fibre_prod_config_is_disabled(config) result(is_disabled)
    type(fibre_prod_config_type), intent(in) :: config

    is_disabled = (.not. config%fibre_enabled) .and. (.not. config%fsi_enabled)
  end function fibre_prod_config_is_disabled

end module fibre_prod_config
