module fibre_stage6_rk_rhs_sync
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_rhs_audit
  use fibre_stage6_controlled_rhs_hook
  implicit none
contains
  subroutine fill_stage6_rk_rhs(substep,rhsx,rhsy,rhsz)
    integer, intent(in) :: substep
    real(mytype), intent(out) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    integer::i,j,k
    do k=1,size(rhsx,3);do j=1,size(rhsx,2);do i=1,size(rhsx,1)
      rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k+0.02_mytype*substep
      rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i-0.01_mytype*substep
      rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j+0.005_mytype*substep
    end do;end do;end do
  end subroutine
  subroutine fill_stage6_rk_force_buffer(substep,fx,fy,fz)
    integer, intent(in) :: substep
    real(mytype), intent(out) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    integer::i,j,k
    do k=1,size(fx,3);do j=1,size(fx,2);do i=1,size(fx,1)
      fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k+0.17_mytype*substep)
      fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k+0.11_mytype*substep)
      fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j+0.07_mytype*substep)
    end do;end do;end do
  end subroutine
  subroutine apply_stage6_rk_substep_rhs(config,substep,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag)
    type(stage6_config_t), intent(in) :: config
    integer, intent(in) :: substep
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    integer, intent(out) :: hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag
    call apply_stage6_controlled_rhs_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag)
  end subroutine
  subroutine compute_stage6_rk_rhs_increment_difference(rhsx0a,rhsy0a,rhsz0a,rhsx1a,rhsy1a,rhsz1a,rhsx0b,rhsy0b,rhsz0b,rhsx1b,rhsy1b,rhsz1b,diff_norm)
    real(mytype), intent(in)::rhsx0a(:,:,:),rhsy0a(:,:,:),rhsz0a(:,:,:),rhsx1a(:,:,:),rhsy1a(:,:,:),rhsz1a(:,:,:),rhsx0b(:,:,:),rhsy0b(:,:,:),rhsz0b(:,:,:),rhsx1b(:,:,:),rhsy1b(:,:,:),rhsz1b(:,:,:)
    real(mytype), intent(out)::diff_norm
    diff_norm=sqrt(sum(((rhsx1a-rhsx0a)-(rhsx1b-rhsx0b))**2+((rhsy1a-rhsy0a)-(rhsy1b-rhsy0b))**2+((rhsz1a-rhsz0a)-(rhsz1b-rhsz0b))**2))
  end subroutine
  subroutine compute_stage6_rk_buffer_difference(fx1,fy1,fz1,fx2,fy2,fz2,diff_norm)
    real(mytype), intent(in)::fx1(:,:,:),fy1(:,:,:),fz1(:,:,:),fx2(:,:,:),fy2(:,:,:),fz2(:,:,:)
    real(mytype), intent(out)::diff_norm
    diff_norm=sqrt(sum((fx2-fx1)**2+(fy2-fy1)**2+(fz2-fz1)**2))
  end subroutine
  subroutine compute_stage6_rk_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,rho,rhsx1,rhsy1,rhsz1,err,ex,ey,ez)
    real(mytype), intent(in)::rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),rho,rhsx1(:,:,:),rhsy1(:,:,:),rhsz1(:,:,:)
    real(mytype), intent(out)::err,ex,ey,ez
    ex=maxval(abs((rhsx1-rhsx0)-fx/rho)); ey=maxval(abs((rhsy1-rhsy0)-fy/rho)); ez=maxval(abs((rhsz1-rhsz0)-fz/rho)); err=max(ex,max(ey,ez))
  end subroutine
  subroutine compute_stage6_rk_stale_force_error(fx_current,fy_current,fz_current,fx_stale,fy_stale,fz_stale,rho,err)
    real(mytype), intent(in)::fx_current(:,:,:),fy_current(:,:,:),fz_current(:,:,:),fx_stale(:,:,:),fy_stale(:,:,:),fz_stale(:,:,:),rho
    real(mytype), intent(out)::err
    err=max(maxval(abs((fx_current-fx_stale)/rho)),max(maxval(abs((fy_current-fy_stale)/rho)),maxval(abs((fz_current-fz_stale)/rho))))
  end subroutine
  subroutine stage6_rk_pressure_status(pressure_poisson_modified_flag,projection_modified_flag,real_projection_called_flag)
    integer, intent(out)::pressure_poisson_modified_flag,projection_modified_flag,real_projection_called_flag
    pressure_poisson_modified_flag=0; projection_modified_flag=0; real_projection_called_flag=0
  end subroutine
end module
