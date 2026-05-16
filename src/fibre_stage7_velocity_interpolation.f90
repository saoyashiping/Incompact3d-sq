module fibre_stage7_velocity_interpolation
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_scalar_interpolation
  implicit none
  type stage7_velocity_layout_t
    integer :: layout_valid_flag=0,rejected_flag=1,collocated_flag=0,component_specific_flag=0
    integer :: u_layout_valid_flag=0,v_layout_valid_flag=0,w_layout_valid_flag=0
    real(mytype) :: u_x_offset=0,u_y_offset=0,u_z_offset=0,v_x_offset=0,v_y_offset=0,v_z_offset=0,w_x_offset=0,w_y_offset=0,w_z_offset=0
  end type
contains
  subroutine init_stage7_collocated_velocity_layout(layout)
    type(stage7_velocity_layout_t),intent(out)::layout
    layout%collocated_flag=1; layout%component_specific_flag=0
    layout%u_layout_valid_flag=1; layout%v_layout_valid_flag=1; layout%w_layout_valid_flag=1
    layout%layout_valid_flag=1; layout%rejected_flag=0
  end subroutine
  subroutine init_stage7_component_velocity_layout(layout)
    type(stage7_velocity_layout_t),intent(out)::layout
    call init_stage7_collocated_velocity_layout(layout)
    layout%collocated_flag=0; layout%component_specific_flag=1
  end subroutine
  subroutine validate_stage7_velocity_layout(layout,valid,rejected)
    type(stage7_velocity_layout_t),intent(in)::layout
    integer,intent(out)::valid,rejected
    valid=1
    if (layout%u_layout_valid_flag/=1 .or. layout%v_layout_valid_flag/=1 .or. layout%w_layout_valid_flag/=1) valid=0
    if (.not.(layout%collocated_flag==1 .or. layout%component_specific_flag==1)) valid=0
    rejected=merge(0,1,valid==1)
  end subroutine
  subroutine build_stage7_component_weight(grid,layout,component_id,x,y,z,weight)
    type(stage7_channel_grid_t),intent(in)::grid
    type(stage7_velocity_layout_t),intent(in)::layout
    integer,intent(in)::component_id
    real(mytype),intent(in)::x,y,z
    type(stage7_interp_weight_t),intent(inout)::weight
    real(mytype)::xc,yc,zc
    if (layout%layout_valid_flag/=1) then; weight%valid_flag=0; weight%unsafe_flag=1; weight%n=0; return; end if
    xc=x; yc=y; zc=z
    select case(component_id)
    case(1); xc=x-layout%u_x_offset; yc=y-layout%u_y_offset; zc=z-layout%u_z_offset
    case(2); xc=x-layout%v_x_offset; yc=y-layout%v_y_offset; zc=z-layout%v_z_offset
    case(3); xc=x-layout%w_x_offset; yc=y-layout%w_y_offset; zc=z-layout%w_z_offset
    case default
      weight%valid_flag=0; weight%unsafe_flag=1; weight%n=0; return
    end select
    call build_stage7_scalar_interp_weight(grid,xc,yc,zc,weight)
  end subroutine
  subroutine interpolate_stage7_velocity(grid,layout,ux,uy,uz,x,y,z,u_lag,valid,unsafe)
    type(stage7_channel_grid_t),intent(in)::grid
    type(stage7_velocity_layout_t),intent(in)::layout
    real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),x,y,z
    real(mytype),intent(out)::u_lag(3)
    integer,intent(out)::valid,unsafe
    type(stage7_interp_weight_t)::wu,wv,ww
    call init_stage7_interp_weight(wu,64); call init_stage7_interp_weight(wv,64); call init_stage7_interp_weight(ww,64)
    call build_stage7_component_weight(grid,layout,1,x,y,z,wu); call build_stage7_component_weight(grid,layout,2,x,y,z,wv); call build_stage7_component_weight(grid,layout,3,x,y,z,ww)
    if (wu%valid_flag/=1 .or. wv%valid_flag/=1 .or. ww%valid_flag/=1) then
      u_lag=0; valid=0; unsafe=1
    else
      call interpolate_stage7_scalar(ux,wu,u_lag(1)); call interpolate_stage7_scalar(uy,wv,u_lag(2)); call interpolate_stage7_scalar(uz,ww,u_lag(3)); valid=1; unsafe=0
    end if
    call free_stage7_interp_weight(wu); call free_stage7_interp_weight(wv); call free_stage7_interp_weight(ww)
  end subroutine
end module
