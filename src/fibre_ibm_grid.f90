module fibre_ibm_grid

  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t

  implicit none
  private

  public :: init_uniform_ibm_grid
  public :: destroy_ibm_grid
  public :: init_lagrangian_points_from_fibre
  public :: destroy_lagrangian_points

contains

  subroutine init_uniform_ibm_grid(grid, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, periodic_x, periodic_y, periodic_z)
    type(ibm_grid_t), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(mytype), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    integer :: i

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
      error stop 'init_uniform_ibm_grid: nx, ny, nz must be > 0'
    end if
    if (xmax <= xmin .or. ymax <= ymin .or. zmax <= zmin) then
      error stop 'init_uniform_ibm_grid: invalid domain bounds'
    end if

    call destroy_ibm_grid(grid)

    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%xmin = xmin
    grid%xmax = xmax
    grid%ymin = ymin
    grid%ymax = ymax
    grid%zmin = zmin
    grid%zmax = zmax
    grid%periodic_x = periodic_x
    grid%periodic_y = periodic_y
    grid%periodic_z = periodic_z

    grid%dx = (xmax - xmin) / real(nx, mytype)
    grid%dy = (ymax - ymin) / real(ny, mytype)
    grid%dz = (zmax - zmin) / real(nz, mytype)
    grid%cell_volume = grid%dx * grid%dy * grid%dz

    allocate(grid%x(nx), grid%y(ny), grid%z(nz))

    do i = 1, nx
      grid%x(i) = xmin + (real(i, mytype) - 0.5_mytype) * grid%dx
    end do
    do i = 1, ny
      grid%y(i) = ymin + (real(i, mytype) - 0.5_mytype) * grid%dy
    end do
    do i = 1, nz
      grid%z(i) = zmin + (real(i, mytype) - 0.5_mytype) * grid%dz
    end do
  end subroutine init_uniform_ibm_grid

  subroutine destroy_ibm_grid(grid)
    type(ibm_grid_t), intent(inout) :: grid

    if (allocated(grid%x)) deallocate(grid%x)
    if (allocated(grid%y)) deallocate(grid%y)
    if (allocated(grid%z)) deallocate(grid%z)

    grid%nx = 0
    grid%ny = 0
    grid%nz = 0
    grid%dx = 0._mytype
    grid%dy = 0._mytype
    grid%dz = 0._mytype
    grid%cell_volume = 0._mytype
  end subroutine destroy_ibm_grid

  subroutine init_lagrangian_points_from_fibre(lag, fibre)
    type(ibm_lagrangian_points_t), intent(inout) :: lag
    type(fibre_t), intent(in) :: fibre

    integer :: i

    call destroy_lagrangian_points(lag)

    lag%nl = fibre%nl
    allocate(lag%x(3, lag%nl), lag%v(3, lag%nl), lag%force(3, lag%nl), lag%weight(lag%nl))

    lag%x = fibre%x
    lag%v = fibre%v
    lag%force = 0._mytype

    lag%weight = fibre%ds
    if (lag%nl >= 1) lag%weight(1) = 0.5_mytype * fibre%ds
    if (lag%nl >= 2) lag%weight(lag%nl) = 0.5_mytype * fibre%ds

    do i = 1, lag%nl
      if (lag%weight(i) < 0._mytype) then
        error stop 'init_lagrangian_points_from_fibre: negative lagrangian weight'
      end if
    end do
  end subroutine init_lagrangian_points_from_fibre

  subroutine destroy_lagrangian_points(lag)
    type(ibm_lagrangian_points_t), intent(inout) :: lag

    if (allocated(lag%x)) deallocate(lag%x)
    if (allocated(lag%v)) deallocate(lag%v)
    if (allocated(lag%force)) deallocate(lag%force)
    if (allocated(lag%weight)) deallocate(lag%weight)
    lag%nl = 0
  end subroutine destroy_lagrangian_points

end module fibre_ibm_grid
