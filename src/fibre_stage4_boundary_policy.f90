module fibre_stage4_boundary_policy

  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_external_force, only : clear_fibre_external_force
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_boundary_safety, only : ibm_point_safety_t, init_point_safety, destroy_point_safety, check_ibm_point_boundary_safety, summarize_ibm_boundary_safety
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre, destroy_lagrangian_points
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  use fibre_stage4_interpolation_adapter, only : stage4_adapter_can_use_uniform_collocated_ibm, convert_stage4_adapter_to_ibm_grid

  implicit none
  private
  public :: check_stage4_fibre_boundary_policy
  public :: block_stage4_unsafe_fibre

contains

  subroutine check_stage4_fibre_boundary_policy(adapter, fibre, safe_count, periodic_wrap_count, unsafe_count, outside_count, blocked_flag, status)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    type(fibre_t), intent(in) :: fibre
    integer, intent(out) :: safe_count, periodic_wrap_count, unsafe_count, outside_count, blocked_flag, status

    integer :: can_use
    type(ibm_grid_t) :: grid
    type(ibm_lagrangian_points_t) :: lag
    type(ibm_point_safety_t) :: safety

    safe_count = 0; periodic_wrap_count = 0; unsafe_count = 0; outside_count = 0
    blocked_flag = 0; status = 1

    call stage4_adapter_can_use_uniform_collocated_ibm(adapter, can_use)
    if (can_use /= 1) then
      blocked_flag = 1
      return
    end if

    call init_lagrangian_points_from_fibre(lag, fibre)
    call convert_stage4_adapter_to_ibm_grid(adapter, grid, status)
    if (status /= 1) then
      blocked_flag = 1
      call destroy_lagrangian_points(lag)
      return
    end if

    call init_point_safety(safety, lag%nl)
    call check_ibm_point_boundary_safety(grid, lag, safety)
    call summarize_ibm_boundary_safety(safety, safe_count, periodic_wrap_count, unsafe_count, outside_count)

    if (unsafe_count + outside_count > 0) blocked_flag = 1

    call destroy_point_safety(safety)
    call destroy_lagrangian_points(lag)
  end subroutine

  subroutine block_stage4_unsafe_fibre(fibre, interpolation_called, feedback_called, advance_called)
    type(fibre_t), intent(inout) :: fibre
    integer, intent(out) :: interpolation_called, feedback_called, advance_called

    call clear_fibre_external_force(fibre)
    interpolation_called = 0
    feedback_called = 0
    advance_called = 0
  end subroutine

end module fibre_stage4_boundary_policy
