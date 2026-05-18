module fibre_stage8_runtime_grid_bridge
  use fibre_parameters, only : mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_channel_grid_adapter
  implicit none

  type stage8_runtime_grid_bridge_status_t
    integer :: real_xcompact_grid_source_found_flag
    integer :: real_xcompact_grid_adapter_called_flag
    integer :: explicit_array_fallback_called_flag
    integer :: bridge_grid_valid_flag
    integer :: bridge_grid_rejected_flag
    integer :: validate_stage7_grid_called_flag
    integer :: bridge_status
  end type stage8_runtime_grid_bridge_status_t
contains
  subroutine init_stage8_runtime_bridge_status(st)
    type(stage8_runtime_grid_bridge_status_t), intent(out) :: st
    st%real_xcompact_grid_source_found_flag=0
    st%real_xcompact_grid_adapter_called_flag=0
    st%explicit_array_fallback_called_flag=0
    st%bridge_grid_valid_flag=0
    st%bridge_grid_rejected_flag=0
    st%validate_stage7_grid_called_flag=0
    st%bridge_status=0
  end subroutine

  subroutine init_stage8_grid_from_explicit_arrays_bridge(grid, nx, ny, nz, xmin, xmax, zmin, zmax, y_face, y_center, periodic_x, periodic_z, st)
    type(stage7_channel_grid_t), intent(out) :: grid
    integer, intent(in) :: nx,ny,nz,periodic_x,periodic_z
    real(mytype), intent(in) :: xmin,xmax,zmin,zmax
    real(mytype), intent(in) :: y_face(ny+1), y_center(ny)
    type(stage8_runtime_grid_bridge_status_t), intent(inout) :: st
    integer :: valid, rejected
    call init_stage7_channel_grid_from_arrays(grid,nx,ny,nz,xmin,xmax,zmin,zmax,y_face,y_center,periodic_x,periodic_z,valid,rejected)
    st%explicit_array_fallback_called_flag=1
    st%validate_stage7_grid_called_flag=1
    st%bridge_grid_valid_flag=valid
    st%bridge_grid_rejected_flag=rejected
    st%bridge_status=merge(1,0,valid==1 .and. rejected==0)
  end subroutine

  subroutine init_stage8_grid_from_xcompact_runtime_bridge(grid, st)
    type(stage7_channel_grid_t), intent(out) :: grid
    type(stage8_runtime_grid_bridge_status_t), intent(inout) :: st
    call init_stage7_nonuniform_channel_grid(grid,16,17,12)
    st%real_xcompact_grid_source_found_flag=0
    st%real_xcompact_grid_adapter_called_flag=0
    st%bridge_grid_valid_flag=0
    st%bridge_grid_rejected_flag=0
    st%bridge_status=0
  end subroutine
end module fibre_stage8_runtime_grid_bridge
