module fibre_stage7_boundary_safety
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_velocity_interpolation
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: STAGE7_POINT_SAFE=1,STAGE7_POINT_NEARWALL_BLOCKED=2,STAGE7_POINT_Y_OUTSIDE_BLOCKED=3,STAGE7_POINT_INVALID_COORD_BLOCKED=4,STAGE7_POINT_INVALID_LAYOUT_BLOCKED=5
  type stage7_boundary_safety_result_t
    integer :: safe_flag=0,blocked_flag=0,unsafe_flag=0,nearwall_blocked_flag=0,y_outside_blocked_flag=0,invalid_coord_blocked_flag=0,invalid_layout_blocked_flag=0,periodic_x_allowed_flag=0,periodic_z_allowed_flag=0,status_code=0
  end type
contains
  subroutine init_stage7_boundary_safety_result(res)
    type(stage7_boundary_safety_result_t),intent(out)::res
  end subroutine
  subroutine classify_stage7_boundary_point(grid,layout,x,y,z,require_velocity_layout,res)
    type(stage7_channel_grid_t),intent(in)::grid; type(stage7_velocity_layout_t),intent(in)::layout
    real(mytype),intent(in)::x,y,z; integer,intent(in)::require_velocity_layout
    type(stage7_boundary_safety_result_t),intent(out)::res
    type(stage7_interp_weight_t)::w; integer::v,r; real(mytype)::lx,lz
    call init_stage7_boundary_safety_result(res)
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(y) .or. .not.ieee_is_finite(z)) then
      res%blocked_flag=1;res%unsafe_flag=1;res%invalid_coord_blocked_flag=1;res%status_code=STAGE7_POINT_INVALID_COORD_BLOCKED;return
    end if
    res%periodic_x_allowed_flag=1
    res%periodic_z_allowed_flag=1
    if (y<grid%ymin .or. y>grid%ymax) then
      res%blocked_flag=1;res%unsafe_flag=1;res%y_outside_blocked_flag=1;res%status_code=STAGE7_POINT_Y_OUTSIDE_BLOCKED;return
    end if
    if (require_velocity_layout==1) then
      call validate_stage7_velocity_layout(layout,v,r)
      if (v/=1) then; res%blocked_flag=1;res%unsafe_flag=1;res%invalid_layout_blocked_flag=1;res%status_code=STAGE7_POINT_INVALID_LAYOUT_BLOCKED; return; endif
    end if
    call init_stage7_interp_weight(w,64); call build_stage7_scalar_interp_weight(grid,x,y,z,w)
    if (w%valid_flag/=1 .or. w%unsafe_flag==1) then
      res%blocked_flag=1;res%unsafe_flag=1;res%nearwall_blocked_flag=1;res%status_code=STAGE7_POINT_NEARWALL_BLOCKED
    else
      res%safe_flag=1;res%status_code=STAGE7_POINT_SAFE
    end if
    lx=grid%xmax-grid%xmin; lz=grid%zmax-grid%zmin
    call free_stage7_interp_weight(w)
  end subroutine
end module
