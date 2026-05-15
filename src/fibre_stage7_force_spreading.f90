module fibre_stage7_force_spreading
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_scalar_interpolation
  implicit none
contains
  subroutine clear_stage7_force_buffer(fx,fy,fz)
    real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    fx=0._mytype; fy=0._mytype; fz=0._mytype
  end subroutine
  subroutine spread_stage7_lagrangian_force(grid,x,y,z,force_lag,ds_l,fx,fy,fz,valid,unsafe)
    type(stage7_channel_grid_t),intent(in)::grid
    real(mytype),intent(in)::x,y,z,force_lag(3),ds_l
    real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    integer,intent(out)::valid,unsafe
    type(stage7_interp_weight_t)::w
    integer::q,i,j,k
    real(mytype)::dV
    call init_stage7_interp_weight(w,64)
    call build_stage7_scalar_interp_weight(grid,x,y,z,w)
    if (w%valid_flag/=1 .or. w%unsafe_flag==1) then
      valid=0; unsafe=1; call free_stage7_interp_weight(w); return
    end if
    do q=1,w%n
      i=w%i(q); j=w%j(q); k=w%k(q); dV=grid%volume_y(j)
      fx(i,j,k)=fx(i,j,k)+force_lag(1)*ds_l*w%w(q)/dV
      fy(i,j,k)=fy(i,j,k)+force_lag(2)*ds_l*w%w(q)/dV
      fz(i,j,k)=fz(i,j,k)+force_lag(3)*ds_l*w%w(q)/dV
    end do
    valid=1; unsafe=0
    call free_stage7_interp_weight(w)
  end subroutine
  subroutine spread_stage7_lagrangian_force_set(grid,nlag,xlag,flag,ds,fx,fy,fz,valid_count,unsafe_count)
    type(stage7_channel_grid_t),intent(in)::grid
    integer,intent(in)::nlag
    real(mytype),intent(in)::xlag(3,nlag),flag(3,nlag),ds(nlag)
    real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    integer,intent(out)::valid_count,unsafe_count
    integer::l,valid,unsafe
    valid_count=0; unsafe_count=0
    do l=1,nlag
      call spread_stage7_lagrangian_force(grid,xlag(1,l),xlag(2,l),xlag(3,l),flag(:,l),ds(l),fx,fy,fz,valid,unsafe)
      if (valid==1) valid_count=valid_count+1
      if (unsafe==1) unsafe_count=unsafe_count+1
    end do
  end subroutine
  subroutine compute_stage7_eulerian_force_total(grid,fx,fy,fz,total)
    type(stage7_channel_grid_t),intent(in)::grid
    real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype),intent(out)::total(3)
    integer::i,j,k
    total=0
    do k=1,size(fx,3); do j=1,size(fx,2); do i=1,size(fx,1)
      total(1)=total(1)+fx(i,j,k)*grid%volume_y(j)
      total(2)=total(2)+fy(i,j,k)*grid%volume_y(j)
      total(3)=total(3)+fz(i,j,k)*grid%volume_y(j)
    end do; end do; end do
  end subroutine
  subroutine compute_stage7_lagrangian_force_total(nlag,flag,ds,total)
    integer,intent(in)::nlag
    real(mytype),intent(in)::flag(3,nlag),ds(nlag)
    real(mytype),intent(out)::total(3)
    integer::l
    total=0
    do l=1,nlag; total=total+flag(:,l)*ds(l); end do
  end subroutine
end module
