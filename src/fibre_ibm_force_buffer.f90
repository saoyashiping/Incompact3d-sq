module fibre_ibm_force_buffer

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian

  implicit none
  private

  type, public :: ibm_force_buffer_t
    integer :: nx = 0, ny = 0, nz = 0
    real(mytype), allocatable :: fx(:,:,:)
    real(mytype), allocatable :: fy(:,:,:)
    real(mytype), allocatable :: fz(:,:,:)
    logical :: is_allocated = .false.
  end type ibm_force_buffer_t

  public :: allocate_ibm_force_buffer
  public :: destroy_ibm_force_buffer
  public :: clear_ibm_force_buffer
  public :: accumulate_lag_force_to_buffer
  public :: compute_ibm_force_buffer_total_force
  public :: compute_ibm_force_buffer_norms

contains

  subroutine allocate_ibm_force_buffer(buffer, grid)
    type(ibm_force_buffer_t), intent(inout) :: buffer
    type(ibm_grid_t), intent(in) :: grid

    if (buffer%is_allocated) then
      call destroy_ibm_force_buffer(buffer)
    end if

    allocate(buffer%fx(grid%nx, grid%ny, grid%nz))
    allocate(buffer%fy(grid%nx, grid%ny, grid%nz))
    allocate(buffer%fz(grid%nx, grid%ny, grid%nz))

    buffer%nx = grid%nx
    buffer%ny = grid%ny
    buffer%nz = grid%nz

    buffer%fx = 0._mytype
    buffer%fy = 0._mytype
    buffer%fz = 0._mytype

    buffer%is_allocated = .true.
  end subroutine allocate_ibm_force_buffer

  subroutine destroy_ibm_force_buffer(buffer)
    type(ibm_force_buffer_t), intent(inout) :: buffer

    if (allocated(buffer%fx)) deallocate(buffer%fx)
    if (allocated(buffer%fy)) deallocate(buffer%fy)
    if (allocated(buffer%fz)) deallocate(buffer%fz)

    buffer%nx = 0
    buffer%ny = 0
    buffer%nz = 0
    buffer%is_allocated = .false.
  end subroutine destroy_ibm_force_buffer

  subroutine clear_ibm_force_buffer(buffer)
    type(ibm_force_buffer_t), intent(inout) :: buffer

    if (.not. buffer%is_allocated) then
      error stop 'clear_ibm_force_buffer: buffer is not allocated'
    end if

    buffer%fx = 0._mytype
    buffer%fy = 0._mytype
    buffer%fz = 0._mytype
  end subroutine clear_ibm_force_buffer

  subroutine accumulate_lag_force_to_buffer(grid, lag, buffer)
    type(ibm_grid_t), intent(in) :: grid
    type(ibm_lagrangian_points_t), intent(in) :: lag
    type(ibm_force_buffer_t), intent(inout) :: buffer

    if (.not. buffer%is_allocated) then
      error stop 'accumulate_lag_force_to_buffer: buffer is not allocated'
    end if
    if (buffer%nx /= grid%nx .or. buffer%ny /= grid%ny .or. buffer%nz /= grid%nz) then
      error stop 'accumulate_lag_force_to_buffer: buffer shape mismatch with grid'
    end if

    call spread_lag_force_to_eulerian(grid, lag, buffer%fx, buffer%fy, buffer%fz)
  end subroutine accumulate_lag_force_to_buffer

  subroutine compute_ibm_force_buffer_total_force(grid, buffer, total_force)
    type(ibm_grid_t), intent(in) :: grid
    type(ibm_force_buffer_t), intent(in) :: buffer
    real(mytype), intent(out) :: total_force(3)

    if (.not. buffer%is_allocated) then
      error stop 'compute_ibm_force_buffer_total_force: buffer is not allocated'
    end if

    total_force(1) = sum(buffer%fx) * grid%cell_volume
    total_force(2) = sum(buffer%fy) * grid%cell_volume
    total_force(3) = sum(buffer%fz) * grid%cell_volume
  end subroutine compute_ibm_force_buffer_total_force

  subroutine compute_ibm_force_buffer_norms(buffer, max_abs, l2_norm)
    type(ibm_force_buffer_t), intent(in) :: buffer
    real(mytype), intent(out) :: max_abs, l2_norm

    if (.not. buffer%is_allocated) then
      error stop 'compute_ibm_force_buffer_norms: buffer is not allocated'
    end if

    max_abs = max(maxval(abs(buffer%fx)), max(maxval(abs(buffer%fy)), maxval(abs(buffer%fz))))
    l2_norm = sqrt(sum(buffer%fx * buffer%fx + buffer%fy * buffer%fy + buffer%fz * buffer%fz))
  end subroutine compute_ibm_force_buffer_norms

end module fibre_ibm_force_buffer
