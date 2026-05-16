module fibre_stage7_scalar_interpolation
  use fibre_parameters, only : mytype
  use fibre_stage7_grid_metadata, only : stage7_channel_grid_t
  implicit none

  type stage7_interp_weight_t
    integer :: n=0, valid_flag=0, unsafe_flag=1
    integer, allocatable :: i(:), j(:), k(:)
    real(mytype), allocatable :: w(:)
    real(mytype) :: weight_sum=0._mytype
  end type
contains
  subroutine init_stage7_interp_weight(weight,nmax)
    type(stage7_interp_weight_t), intent(out) :: weight
    integer, intent(in) :: nmax
    allocate(weight%i(nmax),weight%j(nmax),weight%k(nmax),weight%w(nmax))
    weight%n=0; weight%valid_flag=0; weight%unsafe_flag=1; weight%weight_sum=0._mytype
  end subroutine
  subroutine free_stage7_interp_weight(weight)
    type(stage7_interp_weight_t), intent(inout) :: weight
    if (allocated(weight%i)) deallocate(weight%i)
    if (allocated(weight%j)) deallocate(weight%j)
    if (allocated(weight%k)) deallocate(weight%k)
    if (allocated(weight%w)) deallocate(weight%w)
    weight%n=0; weight%valid_flag=0; weight%unsafe_flag=1; weight%weight_sum=0._mytype
  end subroutine
  integer function wrap_periodic_index(idx,n)
    integer, intent(in)::idx,n
    wrap_periodic_index = modulo(idx-1,n)+1
  end function
  subroutine compute_uniform_1d_weights(x,xmin,dx,n,periodic_flag,idx,wx,valid)
    real(mytype), intent(in)::x,xmin,dx
    integer, intent(in)::n,periodic_flag
    integer, intent(out)::idx(4),valid
    real(mytype), intent(out)::wx(4)
    integer :: i0,a,b,p0,pidx
    real(mytype) :: xc(4),den,ldom,xmod,q
    valid=0; wx=0
    if (dx<=0._mytype .or. n<4) return
    if (periodic_flag==1) then
      ldom=real(n,mytype)*dx
      q=modulo(x-xmin,ldom)/dx
      p0=floor(q)-1
      do a=1,4
        pidx=p0+a-1
        idx(a)=modulo(pidx,n)+1
        xc(a)=real(pidx,mytype)
      end do
    else
      q=(x-xmin)/dx
      i0=floor(q)+1
      idx=(/i0-1,i0,i0+1,i0+2/)
      do a=1,4
        if (idx(a)<1 .or. idx(a)>n) return
        xc(a)=real(idx(a)-1,mytype)
      end do
    end if
    do a=1,4
      wx(a)=1._mytype
      do b=1,4
        if (b/=a) then
          den=xc(a)-xc(b)
          if (abs(den)<=tiny(1._mytype)) return
          wx(a)=wx(a)*(q-xc(b))/den
        end if
      end do
    end do
    valid=1
  end subroutine
  subroutine compute_nonuniform_y_weights(y,grid,jy,wy,valid,unsafe)
    real(mytype), intent(in)::y
    type(stage7_channel_grid_t), intent(in)::grid
    integer, intent(out)::jy(4),valid,unsafe
    real(mytype), intent(out)::wy(4)
    integer :: j0,a,b
    real(mytype)::den
    valid=0; unsafe=1; wy=0
    if (y<grid%y_center(2) .or. y>grid%y_center(grid%ny-1)) return
    j0=2
    do while (j0<grid%ny-1 .and. grid%y_center(j0+1)<=y)
      j0=j0+1
    end do
    jy=(/j0-1,j0,j0+1,j0+2/)
    if (jy(1)<1 .or. jy(4)>grid%ny) return
    do a=1,4
      wy(a)=1._mytype
      do b=1,4
        if (b/=a) then
          den=grid%y_center(jy(a))-grid%y_center(jy(b))
          if (abs(den)<=tiny(1._mytype)) return
          wy(a)=wy(a)*(y-grid%y_center(jy(b)))/den
        end if
      end do
    end do
    valid=1; unsafe=0
  end subroutine
  subroutine build_stage7_scalar_interp_weight(grid,x,y,z,weight)
    type(stage7_channel_grid_t), intent(in)::grid
    real(mytype), intent(in)::x,y,z
    type(stage7_interp_weight_t), intent(inout)::weight
    integer :: ix(4),jy(4),kz(4),vx,vy,vz,uy,a,b,c,q
    real(mytype)::wx(4),wy(4),wz(4)
    call compute_uniform_1d_weights(x,grid%xmin,grid%dx,grid%nx,grid%periodic_x,ix,wx,vx)
    call compute_nonuniform_y_weights(y,grid,jy,wy,vy,uy)
    call compute_uniform_1d_weights(z,grid%zmin,grid%dz,grid%nz,grid%periodic_z,kz,wz,vz)
    if (vx==0 .or. vy==0 .or. vz==0) then
      weight%valid_flag=0; weight%unsafe_flag=1; weight%n=0; weight%weight_sum=0; return
    end if
    q=0
    do c=1,4; do b=1,4; do a=1,4
      q=q+1; weight%i(q)=ix(a); weight%j(q)=jy(b); weight%k(q)=kz(c)
      weight%w(q)=wx(a)*wy(b)*wz(c)
    end do; end do; end do
    weight%n=q; weight%weight_sum=sum(weight%w(1:q)); weight%valid_flag=1; weight%unsafe_flag=0
  end subroutine
  subroutine interpolate_stage7_scalar(phi,weight,phi_l)
    real(mytype), intent(in) :: phi(:,:,:)
    type(stage7_interp_weight_t), intent(in)::weight
    real(mytype), intent(out)::phi_l
    integer::q
    if (weight%valid_flag/=1 .or. weight%n<=0) then; phi_l=0._mytype; return; end if
    phi_l=0._mytype
    do q=1,weight%n
      phi_l=phi_l + weight%w(q)*phi(weight%i(q),weight%j(q),weight%k(q))
    end do
  end subroutine
end module
