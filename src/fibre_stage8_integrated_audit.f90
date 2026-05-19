module fibre_stage8_integrated_audit
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  use fibre_stage8_twoway_force_density, only: compute_stage8_eulerian_total_force, compute_stage8_lagrangian_fluid_force_total
  implicit none
contains
  subroutine compute_stage8_integrated_action_reaction_error(state, err_max)
    type(stage8_lagrangian_state_t),intent(in)::state; real(mytype),intent(out)::err_max; integer::l
    err_max=0._mytype; if(state%allocated_flag/=1)return
    do l=1,state%nlag
      if(state%point_valid_flag(l)==1.and.state%point_blocked_flag(l)==0.and.state%point_unsafe_flag(l)==0) err_max=max(err_max,sqrt(sum((state%force_structure(:,l)+state%force_fluid(:,l))**2)))
    end do
  end
  subroutine compute_stage8_integrated_pair_power(state,beta,p_structure,p_fluid_lag,p_pair,p_expected,p_error)
    type(stage8_lagrangian_state_t),intent(in)::state; real(mytype),intent(in)::beta
    real(mytype),intent(out)::p_structure,p_fluid_lag,p_pair,p_expected,p_error; integer::l
    p_structure=0; p_fluid_lag=0; p_pair=0; p_expected=0; p_error=0; if(state%allocated_flag/=1)return
    do l=1,state%nlag
      if(state%point_valid_flag(l)==1.and.state%point_blocked_flag(l)==0.and.state%point_unsafe_flag(l)==0) then
        p_structure=p_structure+dot_product(state%force_structure(:,l),state%v_fibre(:,l))*state%ds(l)
        p_fluid_lag=p_fluid_lag+dot_product(state%force_fluid(:,l),state%u_fluid_lag(:,l))*state%ds(l)
        p_expected=p_expected-beta*dot_product(state%slip(:,l),state%slip(:,l))*state%ds(l)
      end if
    end do
    p_pair=p_structure+p_fluid_lag; p_error=abs(p_pair-p_expected)
  end
  subroutine compute_stage8_integrated_eulerian_power(grid,ux,uy,uz,fx,fy,fz,p_eulerian)
    type(stage7_channel_grid_t),intent(in)::grid; real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:); real(mytype),intent(out)::p_eulerian
    integer::i,j,k
    p_eulerian=0; do k=1,size(fx,3); do j=1,size(fx,2); do i=1,size(fx,1)
      p_eulerian=p_eulerian+(fx(i,j,k)*ux(i,j,k)+fy(i,j,k)*uy(i,j,k)+fz(i,j,k)*uz(i,j,k))*grid%volume_y(j)
    end do; end do; end do
  end
  subroutine compute_stage8_integrated_force_density_convention(grid,fx,fy,fz,state,total_eulerian,total_lagrangian,abs_err,rel_err)
    type(stage7_channel_grid_t),intent(in)::grid; real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:); type(stage8_lagrangian_state_t),intent(in)::state
    real(mytype),intent(out)::total_eulerian(3),total_lagrangian(3),abs_err,rel_err
    call compute_stage8_eulerian_total_force(grid,fx,fy,fz,total_eulerian); call compute_stage8_lagrangian_fluid_force_total(state,total_lagrangian)
    abs_err=sqrt(sum((total_eulerian-total_lagrangian)**2)); rel_err=abs_err/max(sqrt(sum(total_lagrangian**2)),1e-30_mytype)
  end
  subroutine compute_stage8_integrated_buffer_difference(fx1,fy1,fz1,fx2,fy2,fz2,err_max)
    real(mytype),intent(in)::fx1(:,:,:),fy1(:,:,:),fz1(:,:,:),fx2(:,:,:),fy2(:,:,:),fz2(:,:,:); real(mytype),intent(out)::err_max
    err_max=max(maxval(abs(fx1-fx2)),max(maxval(abs(fy1-fy2)),maxval(abs(fz1-fz2))))
  end
end module
