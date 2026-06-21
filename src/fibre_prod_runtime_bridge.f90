module fibre_prod_runtime_bridge
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, &
                                          fibre_prod_force_buffer_destroy, &
                                          fibre_prod_force_buffer_reset_to_zero
  use fibre_prod_main_hook, only : fibre_prod_main_hook_apply_force_buffer
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_runtime_bridge_type
  public :: fibre_prod_runtime_bridge_init_from_rhs
  public :: fibre_prod_runtime_bridge_reset_force_buffer
  public :: fibre_prod_runtime_bridge_finalize
  public :: fibre_prod_runtime_bridge_apply_lambda0_noop

  type :: fibre_prod_runtime_bridge_type
    logical :: initialized = .false.
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    logical :: force_buffer_allocated = .false.
    logical :: last_zero_force_buffer = .true.
    logical :: last_physical_response = .false.
    type(fibre_prod_force_buffer_type) :: force_buffer
  end type fibre_prod_runtime_bridge_type

contains

  subroutine fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_z, bridge, status)
    real(dp), intent(in) :: rhs_x(:, :, :)
    real(dp), intent(in) :: rhs_y(:, :, :)
    real(dp), intent(in) :: rhs_z(:, :, :)
    type(fibre_prod_runtime_bridge_type), intent(inout) :: bridge
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    if (.not. same_shape(rhs_x, rhs_y) .or. .not. same_shape(rhs_x, rhs_z)) then
      status = 1
      return
    end if

    call fibre_prod_runtime_bridge_finalize(bridge)
    bridge%nx = size(rhs_x, 1)
    bridge%ny = size(rhs_x, 2)
    bridge%nz = size(rhs_x, 3)
    bridge%force_buffer%nx_local = bridge%nx
    bridge%force_buffer%ny_local = bridge%ny
    bridge%force_buffer%nz_local = bridge%nz
    allocate(bridge%force_buffer%fx(bridge%nx, bridge%ny, bridge%nz), &
             bridge%force_buffer%fy(bridge%nx, bridge%ny, bridge%nz), &
             bridge%force_buffer%fz(bridge%nx, bridge%ny, bridge%nz), stat=ierr)
    if (ierr /= 0) then
      status = 2
      call fibre_prod_runtime_bridge_finalize(bridge)
      return
    end if
    bridge%force_buffer%allocated = .true.
    bridge%force_buffer_allocated = .true.
    bridge%initialized = .true.
    call fibre_prod_runtime_bridge_reset_force_buffer(bridge, status)
    if (status /= 0) call fibre_prod_runtime_bridge_finalize(bridge)
  end subroutine fibre_prod_runtime_bridge_init_from_rhs

  subroutine fibre_prod_runtime_bridge_reset_force_buffer(bridge, status)
    type(fibre_prod_runtime_bridge_type), intent(inout) :: bridge
    integer, intent(out) :: status

    call fibre_prod_force_buffer_reset_to_zero(bridge%force_buffer, status)
    if (status == 0) then
      bridge%last_zero_force_buffer = .true.
      bridge%last_physical_response = .false.
    end if
  end subroutine fibre_prod_runtime_bridge_reset_force_buffer

  subroutine fibre_prod_runtime_bridge_finalize(bridge)
    type(fibre_prod_runtime_bridge_type), intent(inout) :: bridge

    call fibre_prod_force_buffer_destroy(bridge%force_buffer)
    bridge%initialized = .false.
    bridge%nx = 0
    bridge%ny = 0
    bridge%nz = 0
    bridge%force_buffer_allocated = .false.
    bridge%last_zero_force_buffer = .true.
    bridge%last_physical_response = .false.
  end subroutine fibre_prod_runtime_bridge_finalize

  subroutine fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    type(fibre_prod_runtime_bridge_type), intent(inout) :: bridge
    integer, intent(out) :: status

    if (.not. bridge%initialized .or. .not. bridge%force_buffer_allocated) then
      status = 3
      return
    end if
    if (.not. same_shape(rhs_x, bridge%force_buffer%fx) .or. &
        .not. same_shape(rhs_y, bridge%force_buffer%fy) .or. &
        .not. same_shape(rhs_z, bridge%force_buffer%fz)) then
      status = 4
      return
    end if
    bridge%last_zero_force_buffer = all(bridge%force_buffer%fx == 0.0_dp) .and. &
                                    all(bridge%force_buffer%fy == 0.0_dp) .and. &
                                    all(bridge%force_buffer%fz == 0.0_dp)
    bridge%last_physical_response = .not. bridge%last_zero_force_buffer
    call fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, bridge%force_buffer, status)
    if (bridge%last_zero_force_buffer) bridge%last_physical_response = .false.
  end subroutine fibre_prod_runtime_bridge_apply_lambda0_noop

  pure logical function same_shape(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = size(lhs, 1) == size(rhs, 1) .and. size(lhs, 2) == size(rhs, 2) .and. &
              size(lhs, 3) == size(rhs, 3)
  end function same_shape

end module fibre_prod_runtime_bridge
