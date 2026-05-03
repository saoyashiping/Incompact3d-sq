module fibre_stage4_boundary_policy

  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_external_force, only : clear_fibre_external_force
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_boundary_safety, only : ibm_point_safety_t, init_point_safety, destroy_point_safety, check_ibm_point_boundary_safety, summarize_ibm_boundary_safety
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre, destroy_lagrangian_points, destroy_ibm_grid

  implicit none
  private
  public :: check_stage4_fibre_boundary_policy
  public :: block_stage4_unsafe_fibre

contains

  subroutine check_stage4_fibre_boundary_policy(adapter, fibre, safe_count, periodic_wrap_count, unsafe_count, outside_count, blocked_flag, status)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    type(fibre_t), intent(in) :: fibre
    integer, intent(out) :: safe_count, periodic_wrap_count, unsafe_count, outside_count, blocked_flag, status

    type(ibm_grid_t) :: grid
    type(ibm_lagrangian_points_t) :: lag
    type(ibm_point_safety_t) :: safety

    safe_count = 0; periodic_wrap_count = 0; unsafe_count = 0; outside_count = 0
    blocked_flag = 0; status = 1

    call build_boundary_safety_ibm_grid_from_adapter(adapter, grid, status)
    if (status /= 1) then
      blocked_flag = 1
      return
    end if

    call init_lagrangian_points_from_fibre(lag, fibre)

    call init_point_safety(safety, lag%nl)
    call check_ibm_point_boundary_safety(grid, lag, safety)
    call summarize_ibm_boundary_safety(safety, safe_count, periodic_wrap_count, unsafe_count, outside_count)

    if (unsafe_count + outside_count > 0) blocked_flag = 1

    call destroy_point_safety(safety)
    call destroy_lagrangian_points(lag)
    call destroy_ibm_grid(grid)
  end subroutine

  subroutine block_stage4_unsafe_fibre(fibre, interpolation_called, feedback_called, advance_called)
    type(fibre_t), intent(inout) :: fibre
    integer, intent(out) :: interpolation_called, feedback_called, advance_called

    call clear_fibre_external_force(fibre)
    interpolation_called = 0
    feedback_called = 0
    advance_called = 0
  end subroutine

  subroutine build_boundary_safety_ibm_grid_from_adapter(adapter, grid, status)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    type(ibm_grid_t), intent(inout) :: grid
    integer, intent(out) :: status

    status = 0
    if (.not. adapter%uniform_ibm_compatible) return
    if (adapter%velocity_layout_mode /= 1) return

    call destroy_ibm_grid(grid)
    grid%nx = adapter%nx
    grid%ny = adapter%ny
    grid%nz = adapter%nz
    grid%xmin = adapter%xmin
    grid%xmax = adapter%xmax
    grid%ymin = adapter%ymin
    grid%ymax = adapter%ymax
    grid%zmin = adapter%zmin
    grid%zmax = adapter%zmax
    grid%dx = adapter%dx_max
    grid%dy = adapter%dy_max
    grid%dz = adapter%dz_max
    grid%cell_volume = grid%dx * grid%dy * grid%dz
    grid%periodic_x = adapter%periodic_x
    grid%periodic_y = adapter%periodic_y
    grid%periodic_z = adapter%periodic_z
    allocate(grid%x(grid%nx), grid%y(grid%ny), grid%z(grid%nz))
    grid%x = adapter%x
    grid%y = adapter%y
    grid%z = adapter%z
    status = 1
  end subroutine

end module fibre_stage4_boundary_policy
