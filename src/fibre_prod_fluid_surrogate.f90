module fibre_prod_fluid_surrogate
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_is_initialized
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_fluid_surrogate_type
  public :: fibre_prod_fluid_surrogate_allocate
  public :: fibre_prod_fluid_surrogate_reset_constant
  public :: fibre_prod_fluid_surrogate_apply_force_density_step
  public :: fibre_prod_fluid_surrogate_kinetic_energy
  public :: fibre_prod_fluid_surrogate_is_finite
  public :: fibre_prod_fluid_surrogate_destroy

  type :: fibre_prod_fluid_surrogate_type
    integer :: nx_local = 0
    integer :: ny_local = 0
    integer :: nz_local = 0
    logical :: allocated = .false.
    real(dp), allocatable :: u(:, :, :)
    real(dp), allocatable :: v(:, :, :)
    real(dp), allocatable :: w(:, :, :)
  end type fibre_prod_fluid_surrogate_type

contains

  subroutine fibre_prod_fluid_surrogate_allocate(fluid, grid, status)
    type(fibre_prod_fluid_surrogate_type), intent(inout) :: fluid
    type(fibre_prod_grid_type), intent(in) :: grid
    integer, intent(out) :: status
    integer :: ierr

    call fibre_prod_fluid_surrogate_destroy(fluid)
    status = 0
    ierr = 0
    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 1
    else
      fluid%nx_local = grid%nx_local
      fluid%ny_local = grid%ny_local
      fluid%nz_local = grid%nz_local
      allocate(fluid%u(fluid%nx_local, fluid%ny_local, fluid%nz_local), stat=ierr)
      if (ierr == 0) allocate(fluid%v(fluid%nx_local, fluid%ny_local, fluid%nz_local), stat=ierr)
      if (ierr == 0) allocate(fluid%w(fluid%nx_local, fluid%ny_local, fluid%nz_local), stat=ierr)
      if (ierr == 0) then
        fluid%allocated = .true.
        call fibre_prod_fluid_surrogate_reset_constant(fluid, [0.0_dp, 0.0_dp, 0.0_dp], status)
      else
        status = 2
        call fibre_prod_fluid_surrogate_destroy(fluid)
      end if
    end if
  end subroutine fibre_prod_fluid_surrogate_allocate

  subroutine fibre_prod_fluid_surrogate_reset_constant(fluid, velocity, status)
    type(fibre_prod_fluid_surrogate_type), intent(inout) :: fluid
    real(dp), intent(in) :: velocity(3)
    integer, intent(out) :: status

    status = 0
    if (.not. fluid%allocated .or. .not. all(ieee_is_finite(velocity))) then
      status = 1
    else
      fluid%u = velocity(1)
      fluid%v = velocity(2)
      fluid%w = velocity(3)
    end if
  end subroutine fibre_prod_fluid_surrogate_reset_constant

  subroutine fibre_prod_fluid_surrogate_apply_force_density_step(fluid, buffer, dt, lambda_fsi, status)
    type(fibre_prod_fluid_surrogate_type), intent(inout) :: fluid
    type(fibre_prod_force_buffer_type), intent(in) :: buffer
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: lambda_fsi
    integer, intent(out) :: status

    status = 0
    if (.not. fibre_prod_fluid_surrogate_is_finite(fluid)) then
      status = 1
    else if (.not. fibre_prod_force_buffer_is_finite(buffer)) then
      status = 2
    else if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      status = 3
    else if (.not. ieee_is_finite(lambda_fsi) .or. lambda_fsi < 0.0_dp) then
      status = 4
    else
      fluid%u = fluid%u + dt * lambda_fsi * buffer%fx
      fluid%v = fluid%v + dt * lambda_fsi * buffer%fy
      fluid%w = fluid%w + dt * lambda_fsi * buffer%fz
    end if
  end subroutine fibre_prod_fluid_surrogate_apply_force_density_step

  pure real(dp) function fibre_prod_fluid_surrogate_kinetic_energy(fluid, grid) result(energy)
    type(fibre_prod_fluid_surrogate_type), intent(in) :: fluid
    type(fibre_prod_grid_type), intent(in) :: grid

    if (fibre_prod_fluid_surrogate_is_finite(fluid) .and. fibre_prod_grid_is_initialized(grid)) then
      energy = 0.5_dp * sum((fluid%u * fluid%u + fluid%v * fluid%v + fluid%w * fluid%w) * grid%cell_volume)
    else
      energy = huge(1.0_dp)
    end if
  end function fibre_prod_fluid_surrogate_kinetic_energy

  pure logical function fibre_prod_fluid_surrogate_is_finite(fluid) result(is_finite)
    type(fibre_prod_fluid_surrogate_type), intent(in) :: fluid

    is_finite = fluid%allocated .and. allocated(fluid%u) .and. allocated(fluid%v) .and. allocated(fluid%w)
    if (is_finite) is_finite = all(ieee_is_finite(fluid%u)) .and. &
                                all(ieee_is_finite(fluid%v)) .and. &
                                all(ieee_is_finite(fluid%w))
  end function fibre_prod_fluid_surrogate_is_finite

  subroutine fibre_prod_fluid_surrogate_destroy(fluid)
    type(fibre_prod_fluid_surrogate_type), intent(inout) :: fluid

    if (allocated(fluid%u)) deallocate(fluid%u)
    if (allocated(fluid%v)) deallocate(fluid%v)
    if (allocated(fluid%w)) deallocate(fluid%w)
    fluid%nx_local = 0
    fluid%ny_local = 0
    fluid%nz_local = 0
    fluid%allocated = .false.
  end subroutine fibre_prod_fluid_surrogate_destroy

end module fibre_prod_fluid_surrogate
