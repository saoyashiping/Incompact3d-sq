module fibre_stage8_twoway_force_density
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t
  use fibre_stage7_force_spreading, only: spread_stage7_lagrangian_force
  use fibre_stage7_boundary_safety, only: stage7_boundary_safety_result_t, classify_stage7_boundary_point
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  implicit none
  private
  public :: build_stage8_twoway_force_density_candidate
  public :: compute_stage8_eulerian_total_force
  public :: compute_stage8_lagrangian_fluid_force_total
  public :: compute_stage8_force_density_norm
contains
  subroutine build_stage8_twoway_force_density_candidate(grid, layout, state, fx, fy, fz, valid, rejected, valid_count, blocked_count, unsafe_count)
    type(stage7_channel_grid_t),intent(in)::grid
    type(stage7_velocity_layout_t),intent(in)::layout
    type(stage8_lagrangian_state_t),intent(in)::state
    real(mytype),intent(out)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    integer,intent(out)::valid,rejected,valid_count,blocked_count,unsafe_count
    type(stage7_boundary_safety_result_t)::res
    integer::l,sv,su
    fx=0; fy=0; fz=0; valid=0; rejected=1; valid_count=0; blocked_count=0; unsafe_count=0
    if(state%allocated_flag/=1) return
    do l=1,state%nlag
      call classify_stage7_boundary_point(grid,layout,state%x(1,l),state%x(2,l),state%x(3,l),1,res)
      if(res%safe_flag==1) then
        call spread_stage7_lagrangian_force(grid,state%x(1,l),state%x(2,l),state%x(3,l),state%force_fluid(:,l),state%ds(l),fx,fy,fz,sv,su)
        if(sv==1.and.su==0) then; valid_count=valid_count+1; else; unsafe_count=unsafe_count+1; endif
      else
        if(res%blocked_flag==1) blocked_count=blocked_count+1
        if(res%unsafe_flag==1) unsafe_count=unsafe_count+1
      end if
    end do
    if(valid_count>0.and.blocked_count==0.and.unsafe_count==0) then; valid=1; rejected=0; endif
  end
  subroutine compute_stage8_eulerian_total_force(grid,fx,fy,fz,total)
    type(stage7_channel_grid_t),intent(in)::grid
    real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype),intent(out)::total(3)
    integer::i,j,k
    total=0; do k=1,size(fx,3); do j=1,size(fx,2); do i=1,size(fx,1)
      total(1)=total(1)+fx(i,j,k)*grid%volume_y(j); total(2)=total(2)+fy(i,j,k)*grid%volume_y(j); total(3)=total(3)+fz(i,j,k)*grid%volume_y(j)
    end do; end do; end do
  end
  subroutine compute_stage8_lagrangian_fluid_force_total(state,total)
    type(stage8_lagrangian_state_t),intent(in)::state
    real(mytype),intent(out)::total(3)
    integer::l
    total=0; if(state%allocated_flag/=1) return
    do l=1,state%nlag
      if(state%point_valid_flag(l)==1.and.state%point_blocked_flag(l)==0.and.state%point_unsafe_flag(l)==0) total=total+state%force_fluid(:,l)*state%ds(l)
    end do
  end
  subroutine compute_stage8_force_density_norm(fx,fy,fz,norm_max)
    real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype),intent(out)::norm_max
    norm_max=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  end
end module
