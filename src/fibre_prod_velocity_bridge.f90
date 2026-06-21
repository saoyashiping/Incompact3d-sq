module fibre_prod_velocity_bridge
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_destroy
  use fibre_prod_ibm_interpolation, only : fibre_prod_ibm_interpolate_velocity
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_velocity_sample_type
  public :: fibre_prod_velocity_sample_allocate
  public :: fibre_prod_velocity_sample_set_points
  public :: fibre_prod_velocity_bridge_sample
  public :: fibre_prod_velocity_bridge_sample_points
  public :: fibre_prod_velocity_sample_finalize
  public :: fibre_prod_velocity_bridge_env_enabled
  public :: fibre_prod_velocity_bridge_sample_runtime_unit_grid

  type :: fibre_prod_velocity_sample_type
    integer :: npoint = 0
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: u(:, :)
    logical :: allocated = .false.
    logical :: initialized = .false.
  end type fibre_prod_velocity_sample_type

contains

  subroutine fibre_prod_velocity_sample_allocate(sample, npoint, status)
    type(fibre_prod_velocity_sample_type), intent(inout) :: sample
    integer, intent(in) :: npoint
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    call fibre_prod_velocity_sample_finalize(sample)
    if (npoint <= 0) then
      status = 1
      return
    end if
    allocate(sample%x(npoint, 3), sample%u(npoint, 3), stat=ierr)
    if (ierr /= 0) then
      status = 2
      call fibre_prod_velocity_sample_finalize(sample)
      return
    end if
    sample%npoint = npoint
    sample%x = 0.0_dp
    sample%u = 0.0_dp
    sample%allocated = .true.
    sample%initialized = .true.
  end subroutine fibre_prod_velocity_sample_allocate

  subroutine fibre_prod_velocity_sample_set_points(sample, points, status)
    type(fibre_prod_velocity_sample_type), intent(inout) :: sample
    real(dp), intent(in) :: points(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. sample%allocated .or. .not. allocated(sample%x) .or. .not. allocated(sample%u)) then
      status = 3
    else if (size(points, 1) /= sample%npoint .or. size(points, 2) /= 3) then
      status = 4
    else if (.not. all(ieee_is_finite(points))) then
      status = 5
    else
      sample%x = points
      sample%u = 0.0_dp
      sample%initialized = .true.
    end if
  end subroutine fibre_prod_velocity_sample_set_points

  subroutine fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    type(fibre_prod_velocity_sample_type), intent(inout) :: sample
    integer, intent(out) :: status
    real(dp), allocatable :: velocity_field(:, :, :, :)
    real(dp) :: velocity(3)
    integer :: ipoint
    integer :: ierr

    status = 0
    if (.not. sample%allocated .or. .not. sample%initialized .or. .not. allocated(sample%x) .or. &
        .not. allocated(sample%u)) then
      status = 6
      return
    end if
    if (.not. same_shape(ux, uy) .or. .not. same_shape(ux, uz)) then
      status = 7
      return
    end if
    if (.not. all(ieee_is_finite(ux)) .or. .not. all(ieee_is_finite(uy)) .or. &
        .not. all(ieee_is_finite(uz))) then
      status = 8
      return
    end if

    allocate(velocity_field(size(ux, 1), size(ux, 2), size(ux, 3), 3), stat=ierr)
    if (ierr /= 0) then
      status = 9
      return
    end if
    velocity_field(:, :, :, 1) = ux
    velocity_field(:, :, :, 2) = uy
    velocity_field(:, :, :, 3) = uz
    do ipoint = 1, sample%npoint
      call fibre_prod_ibm_interpolate_velocity(grid, velocity_field, sample%x(ipoint, :), velocity, status)
      if (status /= 0) exit
      sample%u(ipoint, :) = velocity
    end do
    deallocate(velocity_field)
  end subroutine fibre_prod_velocity_bridge_sample


  subroutine fibre_prod_velocity_bridge_sample_points(grid, ux, uy, uz, points, sampled_u, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    real(dp), intent(in) :: points(:, :)
    real(dp), intent(out) :: sampled_u(:, :)
    integer, intent(out) :: status
    type(fibre_prod_velocity_sample_type) :: sample

    status = 0
    if (size(points, 2) /= 3 .or. size(sampled_u, 2) /= 3 .or. size(points, 1) /= size(sampled_u, 1)) then
      status = 10
      return
    end if
    call fibre_prod_velocity_sample_allocate(sample, size(points, 1), status)
    if (status == 0) call fibre_prod_velocity_sample_set_points(sample, points, status)
    if (status == 0) call fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
    if (status == 0) sampled_u = sample%u
    call fibre_prod_velocity_sample_finalize(sample)
  end subroutine fibre_prod_velocity_bridge_sample_points

  subroutine fibre_prod_velocity_sample_finalize(sample)
    type(fibre_prod_velocity_sample_type), intent(inout) :: sample

    if (allocated(sample%x)) deallocate(sample%x)
    if (allocated(sample%u)) deallocate(sample%u)
    sample%npoint = 0
    sample%allocated = .false.
    sample%initialized = .false.
  end subroutine fibre_prod_velocity_sample_finalize

  logical function fibre_prod_velocity_bridge_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_VELOCITY_SAMPLE_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_velocity_bridge_env_enabled

  subroutine fibre_prod_velocity_bridge_sample_runtime_unit_grid(ux, uy, uz, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_velocity_sample_type) :: sample
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: z(:)
    real(dp) :: point(1, 3)
    integer :: i

    status = 0
    allocate(x(size(ux, 1)), y(size(ux, 2)), z(size(ux, 3)), stat=status)
    if (status /= 0) return
    call fill_unit_coordinates(x)
    call fill_unit_coordinates(y)
    call fill_unit_coordinates(z)
    call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, size(x), 1, size(y), 1, size(z), &
                                               .false., .false., .false., status)
    if (status == 0) call fibre_prod_velocity_sample_allocate(sample, 1, status)
    point(1, :) = [0.5_dp, 0.5_dp, 0.5_dp]
    if (status == 0) call fibre_prod_velocity_sample_set_points(sample, point, status)
    if (status == 0) call fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
    call fibre_prod_velocity_sample_finalize(sample)
    call fibre_prod_grid_destroy(grid)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(z)) deallocate(z)
  contains
    subroutine fill_unit_coordinates(coord)
      real(dp), intent(out) :: coord(:)
      integer :: i

      if (size(coord) == 1) then
        coord = 0.0_dp
      else
        do i = 1, size(coord)
          coord(i) = real(i - 1, dp) / real(size(coord) - 1, dp)
        end do
      end if
    end subroutine fill_unit_coordinates
  end subroutine fibre_prod_velocity_bridge_sample_runtime_unit_grid

  pure logical function same_shape(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = size(lhs, 1) == size(rhs, 1) .and. size(lhs, 2) == size(rhs, 2) .and. &
              size(lhs, 3) == size(rhs, 3)
  end function same_shape

end module fibre_prod_velocity_bridge
