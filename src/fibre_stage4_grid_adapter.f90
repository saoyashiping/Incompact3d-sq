module fibre_stage4_grid_adapter

  use fibre_parameters, only : mytype

  implicit none
  private

  type, public :: stage4_grid_adapter_t
    integer :: nx = 0, ny = 0, nz = 0

    real(mytype) :: xmin = 0._mytype, xmax = 0._mytype
    real(mytype) :: ymin = 0._mytype, ymax = 0._mytype
    real(mytype) :: zmin = 0._mytype, zmax = 0._mytype

    real(mytype) :: dx_min = 0._mytype, dx_max = 0._mytype
    real(mytype) :: dy_min = 0._mytype, dy_max = 0._mytype
    real(mytype) :: dz_min = 0._mytype, dz_max = 0._mytype

    logical :: periodic_x = .false.
    logical :: periodic_y = .false.
    logical :: periodic_z = .false.

    logical :: uniform_x = .false.
    logical :: uniform_y = .false.
    logical :: uniform_z = .false.

    logical :: uniform_ibm_compatible = .false.

    integer :: velocity_layout_mode = 0

    real(mytype), allocatable :: x(:)
    real(mytype), allocatable :: y(:)
    real(mytype), allocatable :: z(:)
  end type stage4_grid_adapter_t

  public :: init_stage4_grid_adapter_from_arrays
  public :: destroy_stage4_grid_adapter
  public :: stage4_adapter_rhs_disabled_flag
  public :: write_stage4_grid_adapter_diagnostics

contains

  subroutine init_stage4_grid_adapter_from_arrays(adapter, x_in, y_in, z_in, periodic_x, periodic_y, periodic_z, velocity_layout_mode)
    type(stage4_grid_adapter_t), intent(inout) :: adapter
    real(mytype), intent(in) :: x_in(:), y_in(:), z_in(:)
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(in) :: velocity_layout_mode

    call destroy_stage4_grid_adapter(adapter)

    if (size(x_in) < 2 .or. size(y_in) < 2 .or. size(z_in) < 2) then
      error stop 'init_stage4_grid_adapter_from_arrays: each axis must have at least 2 points'
    end if

    adapter%nx = size(x_in)
    adapter%ny = size(y_in)
    adapter%nz = size(z_in)

    allocate(adapter%x(adapter%nx), adapter%y(adapter%ny), adapter%z(adapter%nz))
    adapter%x = x_in
    adapter%y = y_in
    adapter%z = z_in

    adapter%xmin = minval(adapter%x)
    adapter%xmax = maxval(adapter%x)
    adapter%ymin = minval(adapter%y)
    adapter%ymax = maxval(adapter%y)
    adapter%zmin = minval(adapter%z)
    adapter%zmax = maxval(adapter%z)

    call axis_spacing_stats(adapter%x, adapter%dx_min, adapter%dx_max, adapter%uniform_x)
    call axis_spacing_stats(adapter%y, adapter%dy_min, adapter%dy_max, adapter%uniform_y)
    call axis_spacing_stats(adapter%z, adapter%dz_min, adapter%dz_max, adapter%uniform_z)

    adapter%uniform_ibm_compatible = adapter%uniform_x .and. adapter%uniform_y .and. adapter%uniform_z

    adapter%periodic_x = periodic_x
    adapter%periodic_y = periodic_y
    adapter%periodic_z = periodic_z
    adapter%velocity_layout_mode = velocity_layout_mode
  end subroutine init_stage4_grid_adapter_from_arrays

  subroutine axis_spacing_stats(axis, dmin, dmax, is_uniform)
    real(mytype), intent(in) :: axis(:)
    real(mytype), intent(out) :: dmin, dmax
    logical, intent(out) :: is_uniform

    real(mytype), allocatable :: delta(:)
    real(mytype) :: tol

    allocate(delta(size(axis) - 1))
    delta = axis(2:) - axis(:size(axis)-1)

    dmin = minval(delta)
    dmax = maxval(delta)

    tol = 1.0e-12_mytype * max(1.0_mytype, abs(dmax))
    is_uniform = (dmax - dmin <= tol)

    deallocate(delta)
  end subroutine axis_spacing_stats

  subroutine destroy_stage4_grid_adapter(adapter)
    type(stage4_grid_adapter_t), intent(inout) :: adapter

    if (allocated(adapter%x)) deallocate(adapter%x)
    if (allocated(adapter%y)) deallocate(adapter%y)
    if (allocated(adapter%z)) deallocate(adapter%z)

    adapter%nx = 0
    adapter%ny = 0
    adapter%nz = 0
  end subroutine destroy_stage4_grid_adapter

  subroutine stage4_adapter_rhs_disabled_flag(apply_ibm_to_fluid_rhs, rhs_disabled_flag)
    logical, intent(in) :: apply_ibm_to_fluid_rhs
    integer, intent(out) :: rhs_disabled_flag

    if (.not. apply_ibm_to_fluid_rhs) then
      rhs_disabled_flag = 1
    else
      rhs_disabled_flag = 0
    end if
  end subroutine stage4_adapter_rhs_disabled_flag

  subroutine write_stage4_grid_adapter_diagnostics(adapter, filename)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    character(len=*), intent(in) :: filename

    integer :: unit_id

    open(newunit=unit_id, file=trim(filename), status='replace', action='write', form='formatted')
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_nx', adapter%nx
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_ny', adapter%ny
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_nz', adapter%nz
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_xmin', adapter%xmin
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_xmax', adapter%xmax
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_ymin', adapter%ymin
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_ymax', adapter%ymax
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_zmin', adapter%zmin
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_zmax', adapter%zmax
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_dx_min', adapter%dx_min
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_dx_max', adapter%dx_max
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_dy_min', adapter%dy_min
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_dy_max', adapter%dy_max
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_dz_min', adapter%dz_min
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_adapter_dz_max', adapter%dz_max
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_periodic_x', merge(1, 0, adapter%periodic_x)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_periodic_y', merge(1, 0, adapter%periodic_y)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_periodic_z', merge(1, 0, adapter%periodic_z)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_uniform_x', merge(1, 0, adapter%uniform_x)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_uniform_y', merge(1, 0, adapter%uniform_y)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_uniform_z', merge(1, 0, adapter%uniform_z)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_uniform_ibm_compatible', merge(1, 0, adapter%uniform_ibm_compatible)
    write(unit_id,'(A,1X,I0)') 'stage4_adapter_velocity_layout_mode', adapter%velocity_layout_mode
    close(unit_id)
  end subroutine write_stage4_grid_adapter_diagnostics

end module fibre_stage4_grid_adapter
