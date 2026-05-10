module fibre_stage6_layout_guard
  use fibre_parameters, only : mytype
  use fibre_stage6_config, only : stage6_config_t
  use fibre_stage6_controlled_rhs_hook, only : apply_stage6_controlled_rhs_hook
  implicit none

  integer, parameter :: STAGE6_LAYOUT_UNKNOWN = 0
  integer, parameter :: STAGE6_LAYOUT_UNIFORM_COLLOCATED = 1
  integer, parameter :: STAGE6_LAYOUT_NONUNIFORM_Y = 2
  integer, parameter :: STAGE6_LAYOUT_STAGGERED = 3

  integer, parameter :: STAGE6_LAYOUT_BLOCK_NONE = 0
  integer, parameter :: STAGE6_LAYOUT_BLOCK_UNKNOWN = 1
  integer, parameter :: STAGE6_LAYOUT_BLOCK_NONUNIFORM_Y = 2
  integer, parameter :: STAGE6_LAYOUT_BLOCK_STAGGERED = 3

  type stage6_layout_t
    integer :: layout_kind
    integer :: nx, ny, nz
    integer :: periodic_x
    integer :: periodic_y
    integer :: periodic_z
    integer :: collocated_velocity
    integer :: uniform_x
    integer :: uniform_y
    integer :: uniform_z
    integer :: supported_flag
    integer :: blocked_flag
    integer :: ordinary_path_allowed_flag
    integer :: block_reason_code
  end type stage6_layout_t
contains
  subroutine init_stage6_uniform_collocated_layout(layout, nx, ny, nz)
    type(stage6_layout_t), intent(out) :: layout
    integer, intent(in) :: nx, ny, nz
    layout%layout_kind = STAGE6_LAYOUT_UNIFORM_COLLOCATED
    layout%nx = nx; layout%ny = ny; layout%nz = nz
    layout%periodic_x = 1; layout%periodic_y = 0; layout%periodic_z = 1
    layout%collocated_velocity = 1
    layout%uniform_x = 1; layout%uniform_y = 1; layout%uniform_z = 1
    call validate_stage6_layout_guard(layout)
  end subroutine

  subroutine init_stage6_nonuniform_y_layout(layout, nx, ny, nz)
    type(stage6_layout_t), intent(out) :: layout
    integer, intent(in) :: nx, ny, nz
    layout%layout_kind = STAGE6_LAYOUT_NONUNIFORM_Y
    layout%nx = nx; layout%ny = ny; layout%nz = nz
    layout%periodic_x = 1; layout%periodic_y = 0; layout%periodic_z = 1
    layout%collocated_velocity = 1
    layout%uniform_x = 1; layout%uniform_y = 0; layout%uniform_z = 1
    call validate_stage6_layout_guard(layout)
  end subroutine

  subroutine init_stage6_staggered_layout(layout, nx, ny, nz)
    type(stage6_layout_t), intent(out) :: layout
    integer, intent(in) :: nx, ny, nz
    layout%layout_kind = STAGE6_LAYOUT_STAGGERED
    layout%nx = nx; layout%ny = ny; layout%nz = nz
    layout%periodic_x = 1; layout%periodic_y = 0; layout%periodic_z = 1
    layout%collocated_velocity = 0
    layout%uniform_x = 1; layout%uniform_y = 1; layout%uniform_z = 1
    call validate_stage6_layout_guard(layout)
  end subroutine

  subroutine init_stage6_unknown_layout(layout, nx, ny, nz)
    type(stage6_layout_t), intent(out) :: layout
    integer, intent(in) :: nx, ny, nz
    layout%layout_kind = STAGE6_LAYOUT_UNKNOWN
    layout%nx = nx; layout%ny = ny; layout%nz = nz
    layout%periodic_x = 1; layout%periodic_y = 0; layout%periodic_z = 1
    layout%collocated_velocity = 0
    layout%uniform_x = 0; layout%uniform_y = 0; layout%uniform_z = 0
    call validate_stage6_layout_guard(layout)
  end subroutine

  subroutine validate_stage6_layout_guard(layout)
    type(stage6_layout_t), intent(inout) :: layout
    if (layout%layout_kind == STAGE6_LAYOUT_UNIFORM_COLLOCATED .and. layout%collocated_velocity == 1 .and. &
        layout%uniform_x == 1 .and. layout%uniform_y == 1 .and. layout%uniform_z == 1) then
      layout%supported_flag = 1; layout%blocked_flag = 0; layout%ordinary_path_allowed_flag = 1
      layout%block_reason_code = STAGE6_LAYOUT_BLOCK_NONE
    else if (layout%layout_kind == STAGE6_LAYOUT_NONUNIFORM_Y .or. layout%uniform_y == 0) then
      layout%supported_flag = 0; layout%blocked_flag = 1; layout%ordinary_path_allowed_flag = 0
      layout%block_reason_code = STAGE6_LAYOUT_BLOCK_NONUNIFORM_Y
    else if (layout%layout_kind == STAGE6_LAYOUT_STAGGERED .or. layout%collocated_velocity == 0) then
      layout%supported_flag = 0; layout%blocked_flag = 1; layout%ordinary_path_allowed_flag = 0
      layout%block_reason_code = STAGE6_LAYOUT_BLOCK_STAGGERED
    else
      layout%supported_flag = 0; layout%blocked_flag = 1; layout%ordinary_path_allowed_flag = 0
      layout%block_reason_code = STAGE6_LAYOUT_BLOCK_UNKNOWN
    end if
  end subroutine

  subroutine apply_stage6_layout_guarded_rhs(config, layout, fx, fy, fz, rhsx, rhsy, rhsz, &
                                             interpolation_called, spreading_called, rhs_injection_called, fluid_update_called, &
                                             injected_flag, modified_flag, rejected_flag, rhs_change_max)
    type(stage6_config_t), intent(in) :: config
    type(stage6_layout_t), intent(inout) :: layout
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
    integer, intent(out) :: interpolation_called, spreading_called, rhs_injection_called, fluid_update_called
    integer, intent(out) :: injected_flag, modified_flag, rejected_flag
    real(mytype), intent(out) :: rhs_change_max
    real(mytype), allocatable :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
    integer :: hook_called_flag

    allocate(rhsx0(size(rhsx,1),size(rhsx,2),size(rhsx,3)))
    allocate(rhsy0(size(rhsy,1),size(rhsy,2),size(rhsy,3)))
    allocate(rhsz0(size(rhsz,1),size(rhsz,2),size(rhsz,3)))
    rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz

    call validate_stage6_layout_guard(layout)
    if (layout%ordinary_path_allowed_flag == 0) then
      interpolation_called = 0; spreading_called = 0; rhs_injection_called = 0; fluid_update_called = 0
      injected_flag = 0; modified_flag = 0; rejected_flag = 1
      rhs_change_max = max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
      deallocate(rhsx0, rhsy0, rhsz0)
      return
    end if

    interpolation_called = 0; spreading_called = 0
    call apply_stage6_controlled_rhs_hook(config, fx, fy, fz, rhsx, rhsy, rhsz, hook_called_flag, modified_flag, injected_flag, rejected_flag)
    rhs_injection_called = injected_flag
    fluid_update_called = 0
    rhs_change_max = max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))

    deallocate(rhsx0, rhsy0, rhsz0)
  end subroutine

  subroutine stage6_layout_pressure_status(pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag
    pressure_poisson_modified_flag = 0
    projection_modified_flag = 0
    real_projection_called_flag = 0
  end subroutine
end module fibre_stage6_layout_guard
