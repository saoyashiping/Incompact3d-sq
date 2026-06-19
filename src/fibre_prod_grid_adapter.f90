module fibre_prod_grid_adapter
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_grid_type
  public :: fibre_prod_grid_init_from_coordinates
  public :: fibre_prod_grid_validate
  public :: fibre_prod_grid_destroy
  public :: fibre_prod_grid_is_initialized
  public :: fibre_prod_grid_compute_min_spacing
  public :: fibre_prod_grid_compute_max_spacing
  public :: fibre_prod_grid_compute_total_local_volume
  public :: fibre_prod_grid_find_cell
  public :: fibre_prod_grid_bounds_string

  type :: fibre_prod_grid_type
    integer :: nx_global = 0
    integer :: ny_global = 0
    integer :: nz_global = 0
    integer :: nx_local = 0
    integer :: ny_local = 0
    integer :: nz_local = 0
    integer :: istart = 0
    integer :: iend = 0
    integer :: jstart = 0
    integer :: jend = 0
    integer :: kstart = 0
    integer :: kend = 0
    logical :: periodic_x = .false.
    logical :: periodic_y = .false.
    logical :: periodic_z = .false.
    logical :: initialized = .false.
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: dx(:)
    real(dp), allocatable :: dy(:)
    real(dp), allocatable :: dz(:)
    real(dp), allocatable :: cell_volume(:, :, :)
  end type fibre_prod_grid_type

contains

  subroutine fibre_prod_grid_init_from_coordinates(grid, x_coord, y_coord, z_coord, &
                                                   istart, iend, jstart, jend, kstart, kend, &
                                                   periodic_x, periodic_y, periodic_z, status)
    type(fibre_prod_grid_type), intent(inout) :: grid
    real(dp), intent(in) :: x_coord(:)
    real(dp), intent(in) :: y_coord(:)
    real(dp), intent(in) :: z_coord(:)
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in) :: jstart
    integer, intent(in) :: jend
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    logical, intent(in) :: periodic_x
    logical, intent(in) :: periodic_y
    logical, intent(in) :: periodic_z
    integer, intent(out), optional :: status
    integer :: ierr
    integer :: i
    integer :: j
    integer :: k

    ierr = 0
    call fibre_prod_grid_destroy(grid)

    if (size(x_coord) < 2 .or. size(y_coord) < 2 .or. size(z_coord) < 2) then
      ierr = 1
    else if (.not. coordinates_valid(x_coord)) then
      ierr = 2
    else if (.not. coordinates_valid(y_coord)) then
      ierr = 3
    else if (.not. coordinates_valid(z_coord)) then
      ierr = 4
    else if (.not. valid_range(istart, iend, size(x_coord))) then
      ierr = 5
    else if (.not. valid_range(jstart, jend, size(y_coord))) then
      ierr = 6
    else if (.not. valid_range(kstart, kend, size(z_coord))) then
      ierr = 7
    else
      grid%nx_global = size(x_coord)
      grid%ny_global = size(y_coord)
      grid%nz_global = size(z_coord)
      grid%istart = istart
      grid%iend = iend
      grid%jstart = jstart
      grid%jend = jend
      grid%kstart = kstart
      grid%kend = kend
      grid%nx_local = iend - istart + 1
      grid%ny_local = jend - jstart + 1
      grid%nz_local = kend - kstart + 1
      grid%periodic_x = periodic_x
      grid%periodic_y = periodic_y
      grid%periodic_z = periodic_z

      allocate(grid%x(grid%nx_local), grid%y(grid%ny_local), grid%z(grid%nz_local), stat=ierr)
      if (ierr == 0) allocate(grid%dx(grid%nx_local), grid%dy(grid%ny_local), grid%dz(grid%nz_local), stat=ierr)
      if (ierr == 0) allocate(grid%cell_volume(grid%nx_local, grid%ny_local, grid%nz_local), stat=ierr)

      if (ierr == 0) then
        grid%x = x_coord(istart:iend)
        grid%y = y_coord(jstart:jend)
        grid%z = z_coord(kstart:kend)
        call fill_local_spacing(x_coord, istart, iend, grid%dx)
        call fill_local_spacing(y_coord, jstart, jend, grid%dy)
        call fill_local_spacing(z_coord, kstart, kend, grid%dz)
        do k = 1, grid%nz_local
          do j = 1, grid%ny_local
            do i = 1, grid%nx_local
              grid%cell_volume(i, j, k) = grid%dx(i) * grid%dy(j) * grid%dz(k)
            end do
          end do
        end do
        grid%initialized = fibre_prod_grid_validate(grid)
        if (.not. grid%initialized) ierr = 8
      end if
    end if

    if (ierr /= 0) call fibre_prod_grid_destroy(grid)
    if (present(status)) status = ierr
  end subroutine fibre_prod_grid_init_from_coordinates

  pure logical function fibre_prod_grid_validate(grid) result(valid)
    type(fibre_prod_grid_type), intent(in) :: grid

    valid = grid%initialized .or. allocated(grid%x)
    valid = valid .and. grid%nx_global >= 2 .and. grid%ny_global >= 2 .and. grid%nz_global >= 2
    valid = valid .and. valid_range(grid%istart, grid%iend, grid%nx_global)
    valid = valid .and. valid_range(grid%jstart, grid%jend, grid%ny_global)
    valid = valid .and. valid_range(grid%kstart, grid%kend, grid%nz_global)
    valid = valid .and. grid%nx_local == grid%iend - grid%istart + 1
    valid = valid .and. grid%ny_local == grid%jend - grid%jstart + 1
    valid = valid .and. grid%nz_local == grid%kend - grid%kstart + 1
    valid = valid .and. allocated(grid%x) .and. allocated(grid%y) .and. allocated(grid%z)
    valid = valid .and. allocated(grid%dx) .and. allocated(grid%dy) .and. allocated(grid%dz)
    valid = valid .and. allocated(grid%cell_volume)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%x)) .and. strictly_increasing(grid%x)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%y)) .and. strictly_increasing(grid%y)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%z)) .and. strictly_increasing(grid%z)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%dx)) .and. all(grid%dx > 0.0_dp)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%dy)) .and. all(grid%dy > 0.0_dp)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%dz)) .and. all(grid%dz > 0.0_dp)
    if (valid) valid = valid .and. all(ieee_is_finite(grid%cell_volume))
    if (valid) valid = valid .and. all(grid%cell_volume > 0.0_dp)
  end function fibre_prod_grid_validate

  subroutine fibre_prod_grid_destroy(grid)
    type(fibre_prod_grid_type), intent(inout) :: grid

    if (allocated(grid%x)) deallocate(grid%x)
    if (allocated(grid%y)) deallocate(grid%y)
    if (allocated(grid%z)) deallocate(grid%z)
    if (allocated(grid%dx)) deallocate(grid%dx)
    if (allocated(grid%dy)) deallocate(grid%dy)
    if (allocated(grid%dz)) deallocate(grid%dz)
    if (allocated(grid%cell_volume)) deallocate(grid%cell_volume)
    grid%nx_global = 0
    grid%ny_global = 0
    grid%nz_global = 0
    grid%nx_local = 0
    grid%ny_local = 0
    grid%nz_local = 0
    grid%istart = 0
    grid%iend = 0
    grid%jstart = 0
    grid%jend = 0
    grid%kstart = 0
    grid%kend = 0
    grid%periodic_x = .false.
    grid%periodic_y = .false.
    grid%periodic_z = .false.
    grid%initialized = .false.
  end subroutine fibre_prod_grid_destroy

  pure logical function fibre_prod_grid_is_initialized(grid) result(is_initialized)
    type(fibre_prod_grid_type), intent(in) :: grid

    is_initialized = grid%initialized .and. fibre_prod_grid_validate(grid)
  end function fibre_prod_grid_is_initialized

  pure real(dp) function fibre_prod_grid_compute_min_spacing(grid) result(min_spacing)
    type(fibre_prod_grid_type), intent(in) :: grid

    min_spacing = huge(1.0_dp)
    if (.not. fibre_prod_grid_is_initialized(grid)) return
    min_spacing = min(minval(grid%dx), min(minval(grid%dy), minval(grid%dz)))
  end function fibre_prod_grid_compute_min_spacing

  pure real(dp) function fibre_prod_grid_compute_max_spacing(grid) result(max_spacing)
    type(fibre_prod_grid_type), intent(in) :: grid

    max_spacing = -huge(1.0_dp)
    if (.not. fibre_prod_grid_is_initialized(grid)) return
    max_spacing = max(maxval(grid%dx), max(maxval(grid%dy), maxval(grid%dz)))
  end function fibre_prod_grid_compute_max_spacing

  pure real(dp) function fibre_prod_grid_compute_total_local_volume(grid) result(total_volume)
    type(fibre_prod_grid_type), intent(in) :: grid

    total_volume = -1.0_dp
    if (.not. fibre_prod_grid_is_initialized(grid)) return
    total_volume = sum(grid%cell_volume)
  end function fibre_prod_grid_compute_total_local_volume

  subroutine fibre_prod_grid_find_cell(grid, point, i_cell, j_cell, k_cell, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: point(3)
    integer, intent(out) :: i_cell
    integer, intent(out) :: j_cell
    integer, intent(out) :: k_cell
    integer, intent(out) :: status

    i_cell = -1
    j_cell = -1
    k_cell = -1
    status = 0
    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 1
    else if (.not. all(ieee_is_finite(point))) then
      status = 2
    else
      call find_axis_cell(grid%x, grid%istart, point(1), i_cell, status)
      if (status == 0) call find_axis_cell(grid%y, grid%jstart, point(2), j_cell, status)
      if (status == 0) call find_axis_cell(grid%z, grid%kstart, point(3), k_cell, status)
    end if
  end subroutine fibre_prod_grid_find_cell

  function fibre_prod_grid_bounds_string(grid) result(summary)
    type(fibre_prod_grid_type), intent(in) :: grid
    character(len=256) :: summary

    write(summary, '(A,3(I0,1X),A,6(I0,1X),A,3(L1,1X))') &
      'global=', grid%nx_global, grid%ny_global, grid%nz_global, &
      ' ranges=', grid%istart, grid%iend, grid%jstart, grid%jend, grid%kstart, grid%kend, &
      ' periodic=', grid%periodic_x, grid%periodic_y, grid%periodic_z
  end function fibre_prod_grid_bounds_string

  pure logical function coordinates_valid(coord) result(valid)
    real(dp), intent(in) :: coord(:)

    valid = size(coord) >= 2
    if (valid) valid = all(ieee_is_finite(coord))
    if (valid) valid = strictly_increasing(coord)
  end function coordinates_valid

  pure logical function strictly_increasing(coord) result(valid)
    real(dp), intent(in) :: coord(:)
    integer :: i

    valid = .true.
    do i = 1, size(coord) - 1
      valid = valid .and. coord(i + 1) > coord(i)
    end do
  end function strictly_increasing

  pure logical function valid_range(first_index, last_index, n_global) result(valid)
    integer, intent(in) :: first_index
    integer, intent(in) :: last_index
    integer, intent(in) :: n_global

    valid = first_index >= 1 .and. first_index <= last_index .and. last_index <= n_global
  end function valid_range

  pure subroutine fill_local_spacing(coord, first_index, last_index, spacing)
    real(dp), intent(in) :: coord(:)
    integer, intent(in) :: first_index
    integer, intent(in) :: last_index
    real(dp), intent(out) :: spacing(:)
    integer :: i
    integer :: global_index

    do i = 1, size(spacing)
      global_index = first_index + i - 1
      if (global_index < size(coord)) then
        spacing(i) = coord(global_index + 1) - coord(global_index)
      else
        spacing(i) = coord(global_index) - coord(global_index - 1)
      end if
    end do
  end subroutine fill_local_spacing

  pure subroutine find_axis_cell(coord, global_offset, value, cell_index, status)
    real(dp), intent(in) :: coord(:)
    integer, intent(in) :: global_offset
    real(dp), intent(in) :: value
    integer, intent(out) :: cell_index
    integer, intent(inout) :: status
    integer :: i

    if (value < coord(1) .or. value > coord(size(coord))) then
      status = 3
      cell_index = -1
      return
    end if

    if (value == coord(size(coord))) then
      cell_index = global_offset + size(coord) - 2
      return
    end if

    do i = 1, size(coord) - 1
      if (value >= coord(i) .and. value < coord(i + 1)) then
        cell_index = global_offset + i - 1
        return
      end if
    end do

    status = 4
    cell_index = -1
  end subroutine find_axis_cell

end module fibre_prod_grid_adapter
