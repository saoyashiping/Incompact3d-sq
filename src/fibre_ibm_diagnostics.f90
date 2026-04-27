module fibre_ibm_diagnostics

  use decomp_2d_constants, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t

  implicit none
  private

  public :: compute_lagrangian_weight_sum
  public :: compute_grid_basic_diagnostics
  public :: check_lagrangian_points_inside_box

contains

  subroutine compute_lagrangian_weight_sum(lag, weight_sum)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(out) :: weight_sum

    if (.not. allocated(lag%weight)) then
      error stop 'compute_lagrangian_weight_sum: lag%weight is not allocated'
    end if

    weight_sum = sum(lag%weight)
  end subroutine compute_lagrangian_weight_sum

  subroutine compute_grid_basic_diagnostics(grid, nx, ny, nz, dx, dy, dz, cell_volume, periodic_x, periodic_y, periodic_z)
    type(ibm_grid_t), intent(in) :: grid
    integer, intent(out) :: nx, ny, nz
    real(mytype), intent(out) :: dx, dy, dz, cell_volume
    integer, intent(out) :: periodic_x, periodic_y, periodic_z

    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    dx = grid%dx
    dy = grid%dy
    dz = grid%dz
    cell_volume = grid%cell_volume
    periodic_x = merge(1, 0, grid%periodic_x)
    periodic_y = merge(1, 0, grid%periodic_y)
    periodic_z = merge(1, 0, grid%periodic_z)
  end subroutine compute_grid_basic_diagnostics

  subroutine check_lagrangian_points_inside_box(grid, lag, inside_count, outside_count)
    type(ibm_grid_t), intent(in) :: grid
    type(ibm_lagrangian_points_t), intent(in) :: lag
    integer, intent(out) :: inside_count, outside_count

    integer :: i
    logical :: inside_x, inside_y, inside_z

    if (.not. allocated(lag%x)) then
      error stop 'check_lagrangian_points_inside_box: lag%x is not allocated'
    end if

    inside_count = 0
    outside_count = 0

    do i = 1, lag%nl
      inside_x = (lag%x(1, i) >= grid%xmin .and. lag%x(1, i) <= grid%xmax)
      inside_y = (lag%x(2, i) >= grid%ymin .and. lag%x(2, i) <= grid%ymax)
      inside_z = (lag%x(3, i) >= grid%zmin .and. lag%x(3, i) <= grid%zmax)

      if (inside_x .and. inside_y .and. inside_z) then
        inside_count = inside_count + 1
      else
        outside_count = outside_count + 1
      end if
    end do
  end subroutine check_lagrangian_points_inside_box

end module fibre_ibm_diagnostics
