module fibre_stage4_frozen_channel
  use fibre_parameters, only : mytype
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  implicit none
  private
  public :: fill_stage4_frozen_channel_velocity, compute_velocity_change_max
contains
  subroutine fill_stage4_frozen_channel_velocity(adapter, uc, amp_u, amp_v, amp_w, ux, uy, uz)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    real(mytype), intent(in) :: uc, amp_u, amp_v, amp_w
    real(mytype), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer :: i,j,k
    real(mytype) :: lx,lz,yc,h,yhat,xx,zz,pi
    pi=acos(-1._mytype); lx=adapter%xmax-adapter%xmin; lz=adapter%zmax-adapter%zmin
    yc=0.5_mytype*(adapter%ymin+adapter%ymax); h=0.5_mytype*(adapter%ymax-adapter%ymin)
    if (lx<=0._mytype .or. lz<=0._mytype .or. h<=0._mytype) error stop 'fill_stage4_frozen_channel_velocity: invalid dimensions'
    do k=1,adapter%nz; zz=(adapter%z(k)-adapter%zmin)/lz
      do j=1,adapter%ny; yhat=(adapter%y(j)-yc)/h
        do i=1,adapter%nx; xx=(adapter%x(i)-adapter%xmin)/lx
          ux(i,j,k)=uc*(1._mytype-yhat**2)+amp_u*sin(2._mytype*pi*xx)*cos(2._mytype*pi*zz)*(1._mytype-yhat**2)
          uy(i,j,k)=amp_v*sin(2._mytype*pi*zz)*(1._mytype-yhat**2)
          uz(i,j,k)=amp_w*cos(2._mytype*pi*xx)*(1._mytype-yhat**2)
        end do
      end do
    end do
  end subroutine
  subroutine compute_velocity_change_max(ux0,uy0,uz0,ux1,uy1,uz1,change_max)
    real(mytype), intent(in) :: ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),ux1(:,:,:),uy1(:,:,:),uz1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max=max(maxval(abs(ux1-ux0)), max(maxval(abs(uy1-uy0)), maxval(abs(uz1-uz0))))
  end subroutine
end module
