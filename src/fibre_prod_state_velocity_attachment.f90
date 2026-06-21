module fibre_prod_state_velocity_attachment
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity
  use fibre_prod_velocity_bridge, only : fibre_prod_velocity_bridge_sample_points
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_state_velocity_attachment_type
  public :: fibre_prod_state_velocity_attachment_init
  public :: fibre_prod_state_velocity_attachment_set_points
  public :: fibre_prod_state_velocity_attachment_sample
  public :: fibre_prod_state_velocity_attachment_attach_to_state
  public :: fibre_prod_state_velocity_attachment_finalize
  public :: fibre_prod_state_velocity_attachment_env_enabled
  public :: fibre_prod_state_velocity_attachment_runtime_unit_grid

  type :: fibre_prod_state_velocity_attachment_type
    logical :: initialized = .false.
    integer :: nnode = 0
    real(dp), allocatable :: sampled_u(:, :)
    real(dp), allocatable :: sampled_x(:, :)
    real(dp) :: max_abs_sampled_velocity = 0.0_dp
    real(dp) :: sum_abs_sampled_velocity = 0.0_dp
    logical :: has_sampled_velocity = .false.
  end type fibre_prod_state_velocity_attachment_type

contains

  subroutine fibre_prod_state_velocity_attachment_init(attach, nnode, status)
    type(fibre_prod_state_velocity_attachment_type), intent(inout) :: attach
    integer, intent(in) :: nnode
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    call fibre_prod_state_velocity_attachment_finalize(attach)
    if (nnode <= 0) then
      status = 1
      return
    end if
    allocate(attach%sampled_x(nnode, 3), attach%sampled_u(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 2
      call fibre_prod_state_velocity_attachment_finalize(attach)
      return
    end if
    attach%nnode = nnode
    attach%sampled_x = 0.0_dp
    attach%sampled_u = 0.0_dp
    attach%max_abs_sampled_velocity = 0.0_dp
    attach%sum_abs_sampled_velocity = 0.0_dp
    attach%has_sampled_velocity = .false.
    attach%initialized = .true.
  end subroutine fibre_prod_state_velocity_attachment_init

  subroutine fibre_prod_state_velocity_attachment_set_points(attach, points, status)
    type(fibre_prod_state_velocity_attachment_type), intent(inout) :: attach
    real(dp), intent(in) :: points(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. attach%initialized .or. .not. allocated(attach%sampled_x)) then
      status = 3
    else if (size(points, 1) /= attach%nnode .or. size(points, 2) /= 3) then
      status = 4
    else if (.not. all(ieee_is_finite(points))) then
      status = 5
    else
      attach%sampled_x = points
      attach%sampled_u = 0.0_dp
      attach%has_sampled_velocity = .false.
    end if
  end subroutine fibre_prod_state_velocity_attachment_set_points

  subroutine fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    type(fibre_prod_state_velocity_attachment_type), intent(inout) :: attach
    integer, intent(out) :: status

    status = 0
    if (.not. attach%initialized .or. .not. allocated(attach%sampled_x) .or. .not. allocated(attach%sampled_u)) then
      status = 6
      return
    end if
    call fibre_prod_velocity_bridge_sample_points(grid, ux, uy, uz, attach%sampled_x, attach%sampled_u, status)
    if (status == 0) then
      attach%max_abs_sampled_velocity = maxval(abs(attach%sampled_u))
      attach%sum_abs_sampled_velocity = sum(abs(attach%sampled_u))
      attach%has_sampled_velocity = .true.
    end if
  end subroutine fibre_prod_state_velocity_attachment_sample

  subroutine fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)
    type(fibre_prod_state_velocity_attachment_type), intent(in) :: attach
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(out) :: status

    status = 0
    if (.not. attach%initialized .or. .not. attach%has_sampled_velocity .or. .not. allocated(attach%sampled_u)) then
      status = 7
      return
    end if
    call fibre_prod_state_attach_sampled_velocity(state, attach%sampled_u, status)
  end subroutine fibre_prod_state_velocity_attachment_attach_to_state

  subroutine fibre_prod_state_velocity_attachment_finalize(attach)
    type(fibre_prod_state_velocity_attachment_type), intent(inout) :: attach

    if (allocated(attach%sampled_x)) deallocate(attach%sampled_x)
    if (allocated(attach%sampled_u)) deallocate(attach%sampled_u)
    attach%initialized = .false.
    attach%nnode = 0
    attach%max_abs_sampled_velocity = 0.0_dp
    attach%sum_abs_sampled_velocity = 0.0_dp
    attach%has_sampled_velocity = .false.
  end subroutine fibre_prod_state_velocity_attachment_finalize

  logical function fibre_prod_state_velocity_attachment_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_STATE_VELOCITY_ATTACH_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_state_velocity_attachment_env_enabled

  subroutine fibre_prod_state_velocity_attachment_runtime_unit_grid(ux, uy, uz, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_state_velocity_attachment_type) :: attach
    type(fibre_prod_state_type) :: state
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: points(3, 3)

    status = 0
    allocate(x(size(ux, 1)), y(size(ux, 2)), z(size(ux, 3)), stat=status)
    if (status /= 0) return
    call fill_unit_coordinates(x)
    call fill_unit_coordinates(y)
    call fill_unit_coordinates(z)
    call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, size(x), 1, size(y), 1, size(z), &
                                               .false., .false., .false., status)
    points(1, :) = [0.25_dp, 0.5_dp, 0.5_dp]
    points(2, :) = [0.5_dp, 0.5_dp, 0.5_dp]
    points(3, :) = [0.75_dp, 0.5_dp, 0.5_dp]
    if (status == 0) call fibre_prod_state_allocate(state, 1, 3, status)
    if (status == 0) state%x(1, :, :) = points
    if (status == 0) call fibre_prod_state_velocity_attachment_init(attach, 3, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_set_points(attach, points, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)
    call fibre_prod_state_velocity_attachment_finalize(attach)
    call fibre_prod_state_destroy(state)
    call fibre_prod_grid_destroy(grid)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(z)) deallocate(z)
  contains
    subroutine fill_unit_coordinates(coord)
      real(dp), intent(out) :: coord(:)
      integer :: i

      do i = 1, size(coord)
        coord(i) = real(i - 1, dp) / real(max(1, size(coord) - 1), dp)
      end do
    end subroutine fill_unit_coordinates
  end subroutine fibre_prod_state_velocity_attachment_runtime_unit_grid

end module fibre_prod_state_velocity_attachment
