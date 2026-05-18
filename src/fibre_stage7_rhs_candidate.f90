module fibre_stage7_rhs_candidate
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_force_spreading
  use fibre_stage7_boundary_safety
  use fibre_stage6_config
  use fibre_stage6_controlled_rhs_hook
  implicit none
contains
  subroutine build_stage7_force_density_candidate(grid,layout,nlag,xlag,flag,ds,fx,fy,fz,valid_count,blocked_count,unsafe_count)
    type(stage7_channel_grid_t),intent(in)::grid
    type(stage7_velocity_layout_t),intent(in)::layout
    integer,intent(in)::nlag
    real(mytype),intent(in)::xlag(3,nlag),flag(3,nlag),ds(nlag)
    real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    integer,intent(out)::valid_count,blocked_count,unsafe_count
    type(stage7_boundary_safety_result_t)::res
    integer::l,v,u
    call clear_stage7_force_buffer(fx,fy,fz); valid_count=0; blocked_count=0; unsafe_count=0
    do l=1,nlag
      call classify_stage7_boundary_point(grid,layout,xlag(1,l),xlag(2,l),xlag(3,l),1,res)
      if (res%safe_flag==1) then
        call spread_stage7_lagrangian_force(grid,xlag(1,l),xlag(2,l),xlag(3,l),flag(:,l),ds(l),fx,fy,fz,v,u)
        if(v==1) valid_count=valid_count+1
      else
        blocked_count=blocked_count+1
        if(res%unsafe_flag==1) unsafe_count=unsafe_count+1
      end if
    end do
  end subroutine
  subroutine apply_stage7_candidate_to_rhs_controlled(config,rhsx,rhsy,rhsz,fx,fy,fz,injected_flag,modified_flag,hook_called_flag,rejected_flag)
    type(stage6_config_t),intent(in)::config
    real(mytype),intent(inout)::rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    integer,intent(out)::injected_flag,modified_flag,hook_called_flag,rejected_flag
    call apply_stage6_controlled_rhs_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,modified_flag,injected_flag,rejected_flag)
  end subroutine
end module
