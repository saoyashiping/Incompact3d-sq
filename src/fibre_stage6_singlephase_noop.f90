module fibre_stage6_singlephase_noop
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_main_noop_hook
  implicit none
contains
  subroutine compute_stage6_array_change_max(a0,a1,change_max)
    real(mytype), intent(in) :: a0(:,:,:),a1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max=maxval(abs(a1-a0))
  end subroutine
  subroutine compute_stage6_vector_array_change_max(ax0,ay0,az0,ax1,ay1,az1,change_max,ex,ey,ez)
    real(mytype), intent(in) :: ax0(:,:,:),ay0(:,:,:),az0(:,:,:),ax1(:,:,:),ay1(:,:,:),az1(:,:,:)
    real(mytype), intent(out) :: change_max,ex,ey,ez
    ex=maxval(abs(ax1-ax0)); ey=maxval(abs(ay1-ay0)); ez=maxval(abs(az1-az0)); change_max=max(ex,max(ey,ez))
  end subroutine
  subroutine compute_stage6_synthetic_divergence(ux,uy,uz,div)
    real(mytype), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(mytype), intent(out) :: div(:,:,:)
    integer :: i,j,k,nx,ny,nz
    nx=size(ux,1); ny=size(ux,2); nz=size(ux,3); div=0._mytype
    do k=2,nz-1; do j=2,ny-1; do i=2,nx-1
      div(i,j,k)=0.5_mytype*((ux(i+1,j,k)-ux(i-1,j,k))+(uy(i,j+1,k)-uy(i,j-1,k))+(uz(i,j,k+1)-uz(i,j,k-1)))
    end do; end do; end do
  end subroutine
  subroutine compute_stage6_singlephase_checksum(rhsx,rhsy,rhsz,ux,uy,uz,p,checksum)
    real(mytype), intent(in) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:),p(:,:,:)
    real(mytype), intent(out) :: checksum
    checksum=sum(rhsx)+sum(rhsy)+sum(rhsz)+sum(ux)+sum(uy)+sum(uz)+sum(p)
  end subroutine
  subroutine apply_stage6_default_singlephase_noop_path(config,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,p,prhs,hook_called_flag,injected_flag,modified_flag,fluid_update_called_flag)
    type(stage6_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:),p(:,:,:),prhs(:,:,:)
    integer, intent(out) :: hook_called_flag,injected_flag,modified_flag,fluid_update_called_flag
    integer :: rejected_flag
    call apply_stage6_main_noop_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,modified_flag,injected_flag,rejected_flag)
    fluid_update_called_flag=0
  end subroutine
end module
