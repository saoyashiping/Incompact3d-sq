program fibre_stage7_scalar_interpolation_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  implicit none
  type(stage7_channel_grid_t)::grid
  type(stage7_interp_weight_t)::wt
  integer, parameter :: nx=16,ny=17,nz=12,np=6
  real(mytype), allocatable :: phi(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
  real(mytype) :: xp(np),yp(np),zp(np),pl,expv,err
  real(mytype) :: werr_max,cerr_max,lerr_max,perr_max,poerr_max,noop_change
  integer :: i,j,k,p,valid,rej,unsafe_count,valid_cnt,nearwall_unsafe_count
  integer :: dep6,dep70,dep71,wstat,cstat,lstat,pstat,pxstat,pzstat,nwblock
  integer :: noop_mod,pp,proj,rp,dns,flu,fib,status

  call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz)
  call validate_stage7_channel_grid(grid,valid,rej)
  allocate(phi(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz))
  call init_stage7_interp_weight(wt,64)

  call dependency_status(dep6,dep70,dep71)

  xp=(/grid%xmin+3.4_mytype*grid%dx,grid%xmin+7.2_mytype*grid%dx,grid%xmax-0.1_mytype*grid%dx,grid%xmin+0.1_mytype*grid%dx,grid%xmin+9.1_mytype*grid%dx,grid%xmin+12.3_mytype*grid%dx/)
  yp=(/grid%y_center(6),grid%y_center(7),grid%y_center(8),grid%y_center(9),grid%y_center(10),grid%y_center(11)/)
  zp=(/grid%zmin+2.3_mytype*grid%dz,grid%zmin+5.8_mytype*grid%dz,grid%zmax-0.1_mytype*grid%dz,grid%zmin+0.1_mytype*grid%dz,grid%zmin+7.5_mytype*grid%dz,grid%zmin+9.1_mytype*grid%dz/)

  werr_max=0; cerr_max=0; lerr_max=0; poerr_max=0; perr_max=0; unsafe_count=0; valid_cnt=0
  do p=1,np
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),wt)
    if (wt%valid_flag==1) then
      valid_cnt=valid_cnt+1; werr_max=max(werr_max,abs(wt%weight_sum-1._mytype))
    else
      unsafe_count=unsafe_count+1
    end if
  end do
  wstat=merge(1,0,werr_max<=1e-12_mytype .and. valid_cnt>=5 .and. unsafe_count==0)

  phi=2.5_mytype
  do p=1,np
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),wt); call interpolate_stage7_scalar(phi,wt,pl)
    cerr_max=max(cerr_max,abs(pl-2.5_mytype))
  end do
  cstat=merge(1,0,cerr_max<=1e-12_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx
    phi(i,j,k)=1._mytype+0.2_mytype*(grid%xmin+real(i-1,mytype)*grid%dx)-0.3_mytype*grid%y_center(j)+0.1_mytype*(grid%zmin+real(k-1,mytype)*grid%dz)
  end do; end do; end do
  do p=1,2
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),wt); call interpolate_stage7_scalar(phi,wt,pl)
    expv=1._mytype+0.2_mytype*xp(p)-0.3_mytype*yp(p)+0.1_mytype*zp(p)
    lerr_max=max(lerr_max,abs(pl-expv))
  end do
  lstat=merge(1,0,lerr_max<=1e-11_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx
    phi(i,j,k)=1._mytype-grid%y_center(j)**2
  end do; end do; end do
  do p=1,np
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),wt); call interpolate_stage7_scalar(phi,wt,pl)
    poerr_max=max(poerr_max,abs(pl-(1._mytype-yp(p)**2)))
  end do
  pstat=merge(1,0,poerr_max<=1e-11_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx
    phi(i,j,k)=sin(2._mytype*stage7_pi*(grid%xmin+real(i-1,mytype)*grid%dx-grid%xmin)/(grid%xmax-grid%xmin)) + &
               cos(2._mytype*stage7_pi*(grid%zmin+real(k-1,mytype)*grid%dz-grid%zmin)/(grid%zmax-grid%zmin)) + 0.25_mytype*grid%y_center(j)
  end do; end do; end do
  do p=3,4
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),wt); call interpolate_stage7_scalar(phi,wt,pl)
    expv=sin(2._mytype*stage7_pi*(xp(p)-grid%xmin)/(grid%xmax-grid%xmin)) + cos(2._mytype*stage7_pi*(zp(p)-grid%zmin)/(grid%zmax-grid%zmin)) + 0.25_mytype*yp(p)
    perr_max=max(perr_max,abs(pl-expv))
  end do
  pxstat=merge(1,0,perr_max<=1e-10_mytype); pzstat=pxstat

  nearwall_unsafe_count=0
  call build_stage7_scalar_interp_weight(grid,xp(1),grid%y_center(1),zp(1),wt)
  if (wt%unsafe_flag==1 .and. wt%valid_flag==0) nearwall_unsafe_count=nearwall_unsafe_count+1
  nwblock=merge(1,0,nearwall_unsafe_count>0)

  rhsx=1; rhsy=2; rhsz=3
  call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop_change,noop_mod)
  call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)

  status=merge(1,0,dep6==1 .and. dep70==1 .and. dep71==1 .and. wstat==1 .and. cstat==1 .and. lstat==1 .and. pstat==1 .and. pxstat==1 .and. pzstat==1 .and. nwblock==1 .and. noop_mod==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0 .and. noop_change<=1e-14_mytype)

  open(10,file='stage7_outputs/fibre_stage7_scalar_interpolation_check.dat',status='replace')
  write(10,'(A,1X,I0)') 'stage7_interp_stage6_dependency_status',dep6
  write(10,'(A,1X,I0)') 'stage7_interp_stage7_0_dependency_status',dep70
  write(10,'(A,1X,I0)') 'stage7_interp_stage7_1_dependency_status',dep71
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_weight_sum_error_max',werr_max
  write(10,'(A,1X,I0)') 'stage7_interp_weight_sum_status',wstat
  write(10,'(A,1X,I0)') 'stage7_interp_valid_weight_count',valid_cnt
  write(10,'(A,1X,I0)') 'stage7_interp_unsafe_count',unsafe_count
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_constant_error_max',cerr_max
  write(10,'(A,1X,I0)') 'stage7_interp_constant_status',cstat
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_linear_error_max',lerr_max
  write(10,'(A,1X,I0)') 'stage7_interp_linear_status',lstat
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_poiseuille_error_max',poerr_max
  write(10,'(A,1X,I0)') 'stage7_interp_poiseuille_status',pstat
  write(10,'(A,1X,I0)') 'stage7_interp_periodic_x_wrap_status',pxstat
  write(10,'(A,1X,I0)') 'stage7_interp_periodic_z_wrap_status',pzstat
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_periodic_wrap_error_max',perr_max
  write(10,'(A,1X,I0)') 'stage7_interp_nearwall_unsafe_count',nearwall_unsafe_count
  write(10,'(A,1X,I0)') 'stage7_interp_nearwall_blocked_flag',nwblock
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_noop_rhs_change_max',noop_change
  write(10,'(A,1X,I0)') 'stage7_interp_noop_rhs_modified_flag',noop_mod
  write(10,'(A,1X,I0)') 'stage7_interp_pressure_poisson_modified_flag',pp
  write(10,'(A,1X,I0)') 'stage7_interp_projection_modified_flag',proj
  write(10,'(A,1X,I0)') 'stage7_interp_production_dns_called_flag',dns
  write(10,'(A,1X,I0)') 'stage7_interp_fluid_update_called_flag',flu
  write(10,'(A,1X,I0)') 'stage7_interp_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)') 'stage7_scalar_interpolation_check_status',status
  close(10)
contains
  subroutine dependency_status(s6,s70,s71)
    integer,intent(out)::s6,s70,s71
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)* &
       check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype)
    s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype)
    if (.not.file_exists('stage6_outputs/STAGE6_CLOSED.md')) s6=0
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key
    real(mytype),intent(in)::val
    character(len=256)::line,k
    real(mytype)::v
    integer::u,ios
    check_dat=0
    if (.not.file_exists(path)) return
    open(newunit=u,file=path,status='old',iostat=ios); if (ios/=0) return
    do
      read(u,'(A)',iostat=ios) line; if (ios/=0) exit
      read(line,*,iostat=ios) k,v
      if (ios==0 .and. trim(k)==trim(key) .and. abs(v-val)<=1e-12_mytype) then; check_dat=1; exit; end if
    end do
    close(u)
  end function
  logical function file_exists(path)
    character(*),intent(in)::path
    inquire(file=path,exist=file_exists)
  end function
end program
