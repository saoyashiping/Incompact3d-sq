module fibre_prod_ibm_force_buffer
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_is_initialized
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_force_buffer_type
  public :: fibre_prod_force_buffer_allocate
  public :: fibre_prod_force_buffer_reset_to_zero
  public :: fibre_prod_force_buffer_destroy
  public :: fibre_prod_force_buffer_is_finite
  public :: fibre_prod_force_buffer_total_force
  public :: fibre_prod_force_buffer_l2_norm

  type :: fibre_prod_force_buffer_type
    integer :: nx_local = 0
    integer :: ny_local = 0
    integer :: nz_local = 0
    logical :: allocated = .false.
    real(dp), allocatable :: fx(:, :, :)
    real(dp), allocatable :: fy(:, :, :)
    real(dp), allocatable :: fz(:, :, :)
  end type fibre_prod_force_buffer_type

contains

  subroutine fibre_prod_force_buffer_allocate(buffer, grid, status)
    type(fibre_prod_force_buffer_type), intent(inout) :: buffer
    type(fibre_prod_grid_type), intent(in) :: grid
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    ierr = 0
    call fibre_prod_force_buffer_destroy(buffer)
    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 1
    else
      buffer%nx_local = grid%nx_local
      buffer%ny_local = grid%ny_local
      buffer%nz_local = grid%nz_local
      allocate(buffer%fx(buffer%nx_local, buffer%ny_local, buffer%nz_local), stat=ierr)
      if (ierr == 0) allocate(buffer%fy(buffer%nx_local, buffer%ny_local, buffer%nz_local), stat=ierr)
      if (ierr == 0) allocate(buffer%fz(buffer%nx_local, buffer%ny_local, buffer%nz_local), stat=ierr)
      if (ierr == 0) then
        buffer%allocated = .true.
        call fibre_prod_force_buffer_reset_to_zero(buffer, status)
      else
        status = 2
        call fibre_prod_force_buffer_destroy(buffer)
      end if
    end if
  end subroutine fibre_prod_force_buffer_allocate

  subroutine fibre_prod_force_buffer_reset_to_zero(buffer, status)
    type(fibre_prod_force_buffer_type), intent(inout) :: buffer
    integer, intent(out) :: status

    status = 0
    if (.not. buffer%allocated .or. .not. allocated(buffer%fx) .or. &
        .not. allocated(buffer%fy) .or. .not. allocated(buffer%fz)) then
      status = 1
    else
      buffer%fx = 0.0_dp
      buffer%fy = 0.0_dp
      buffer%fz = 0.0_dp
    end if
  end subroutine fibre_prod_force_buffer_reset_to_zero

  subroutine fibre_prod_force_buffer_destroy(buffer)
    type(fibre_prod_force_buffer_type), intent(inout) :: buffer

    if (allocated(buffer%fx)) deallocate(buffer%fx)
    if (allocated(buffer%fy)) deallocate(buffer%fy)
    if (allocated(buffer%fz)) deallocate(buffer%fz)
    buffer%nx_local = 0
    buffer%ny_local = 0
    buffer%nz_local = 0
    buffer%allocated = .false.
  end subroutine fibre_prod_force_buffer_destroy

  pure logical function fibre_prod_force_buffer_is_finite(buffer) result(is_finite)
    type(fibre_prod_force_buffer_type), intent(in) :: buffer

    is_finite = buffer%allocated .and. allocated(buffer%fx) .and. allocated(buffer%fy) .and. allocated(buffer%fz)
    if (is_finite) is_finite = all(ieee_is_finite(buffer%fx)) .and. &
                                all(ieee_is_finite(buffer%fy)) .and. &
                                all(ieee_is_finite(buffer%fz))
  end function fibre_prod_force_buffer_is_finite

  subroutine fibre_prod_force_buffer_total_force(buffer, grid, total_force, status)
    type(fibre_prod_force_buffer_type), intent(in) :: buffer
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(out) :: total_force(3)
    integer, intent(out) :: status

    total_force = 0.0_dp
    status = 0
    if (.not. buffer_matches_grid(buffer, grid)) then
      status = 1
    else if (.not. fibre_prod_force_buffer_is_finite(buffer)) then
      status = 2
    else
      total_force(1) = sum(buffer%fx * grid%cell_volume)
      total_force(2) = sum(buffer%fy * grid%cell_volume)
      total_force(3) = sum(buffer%fz * grid%cell_volume)
    end if
  end subroutine fibre_prod_force_buffer_total_force

  pure real(dp) function fibre_prod_force_buffer_l2_norm(buffer) result(norm_value)
    type(fibre_prod_force_buffer_type), intent(in) :: buffer

    if (fibre_prod_force_buffer_is_finite(buffer)) then
      norm_value = sqrt(sum(buffer%fx * buffer%fx + buffer%fy * buffer%fy + buffer%fz * buffer%fz))
    else
      norm_value = huge(1.0_dp)
    end if
  end function fibre_prod_force_buffer_l2_norm

  pure logical function buffer_matches_grid(buffer, grid) result(matches)
    type(fibre_prod_force_buffer_type), intent(in) :: buffer
    type(fibre_prod_grid_type), intent(in) :: grid

    matches = fibre_prod_grid_is_initialized(grid) .and. buffer%allocated
    matches = matches .and. buffer%nx_local == grid%nx_local
    matches = matches .and. buffer%ny_local == grid%ny_local
    matches = matches .and. buffer%nz_local == grid%nz_local
  end function buffer_matches_grid

end module fibre_prod_ibm_force_buffer
