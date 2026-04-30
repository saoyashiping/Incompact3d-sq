module fibre_stage4_interpolation_adapter

  use fibre_parameters, only : mytype
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : destroy_ibm_grid
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag

  implicit none
  private

  public :: stage4_adapter_can_use_uniform_collocated_ibm
  public :: convert_stage4_adapter_to_ibm_grid
  public :: interpolate_stage4_vector_to_lag_if_supported

contains

  subroutine stage4_adapter_can_use_uniform_collocated_ibm(adapter, can_use)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    integer, intent(out) :: can_use

    if (adapter%uniform_ibm_compatible .and. adapter%velocity_layout_mode == 1) then
      can_use = 1
    else
      can_use = 0
    end if
  end subroutine stage4_adapter_can_use_uniform_collocated_ibm

  subroutine convert_stage4_adapter_to_ibm_grid(adapter, grid, status)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    type(ibm_grid_t), intent(inout) :: grid
    integer, intent(out) :: status

    integer :: can_use

    call stage4_adapter_can_use_uniform_collocated_ibm(adapter, can_use)
    if (can_use /= 1) then
      status = 0
      return
    end if

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
    grid%periodic_x = adapter%periodic_x
    grid%periodic_y = adapter%periodic_y
    grid%periodic_z = adapter%periodic_z
    grid%dx = adapter%dx_max
    grid%dy = adapter%dy_max
    grid%dz = adapter%dz_max
    grid%cell_volume = grid%dx * grid%dy * grid%dz

    allocate(grid%x(grid%nx), grid%y(grid%ny), grid%z(grid%nz))
    grid%x = adapter%x
    grid%y = adapter%y
    grid%z = adapter%z

    status = 1
  end subroutine convert_stage4_adapter_to_ibm_grid

  subroutine interpolate_stage4_vector_to_lag_if_supported(adapter, ux, uy, uz, lag, u_lag, status)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(out) :: u_lag(3, lag%nl)
    integer, intent(out) :: status

    type(ibm_grid_t) :: grid

    call convert_stage4_adapter_to_ibm_grid(adapter, grid, status)
    if (status /= 1) then
      u_lag = 0._mytype
      return
    end if

    call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
    call destroy_ibm_grid(grid)
    status = 1
  end subroutine interpolate_stage4_vector_to_lag_if_supported

end module fibre_stage4_interpolation_adapter
