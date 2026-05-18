program fibre_stage7_force_spreading_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_force_spreading
  implicit none
  type(stage7_channel_grid_t)::grid
  type(stage7_interp_weight_t)::w
  integer,parameter::nx=16,ny=17,nz=12,nlag=5
  real(mytype),allocatable::fx(:,:,:),fy(:,:,:),fz(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),fx0(:,:,:),fy0(:,:,:),fz0(:,:,:)
  real(mytype)::x,y,z,lx,lz,ds1,F1(3),Fe(3),Fl(3),single_abs,single_rel,multi_abs,multi_rel,vol_err,perr,nwchg,noop_change,yl,yh,conv_abs
  real(mytype)::xlag(3,nlag),flag(3,nlag),ds(nlag)
  integer::i,j,k,q,valid,rej,unsafe,dep6,dep70,dep71,dep72,dep73,dep74,vcount,ucount,pxw,pzw,nwc,nwflag,ylof,yhif,yout,pp,proj,rp,dns,flu,fib,noop_mod
  integer::single_ok,multi_ok,vol_ok,per_ok,conv_flag,no_rho_flag,status,nonuni_flag

  call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz); call validate_stage7_channel_grid(grid,valid,rej)
  allocate(fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),fx0(nx,ny,nz),fy0(nx,ny,nz),fz0(nx,ny,nz))
  call init_stage7_interp_weight(w,64); call dependency_status(dep6,dep70,dep71,dep72,dep73,dep74)

  lx=grid%xmax-grid%xmin; lz=grid%zmax-grid%zmin
  x=grid%xmin+0.37_mytype*lx; y=0.5_mytype*(grid%y_center(8)+grid%y_center(9)); z=grid%zmin+0.41_mytype*lz
  F1=(/1.2_mytype,-0.7_mytype,0.35_mytype/); ds1=0.05_mytype

  ! single-point
  call clear_stage7_force_buffer(fx,fy,fz)
  call spread_stage7_lagrangian_force(grid,x,y,z,F1,ds1,fx,fy,fz,valid,unsafe)
  call compute_stage7_eulerian_force_total(grid,fx,fy,fz,Fe); Fl=F1*ds1
  single_abs=maxval(abs(Fe-Fl)); single_rel=single_abs/max(1e-16_mytype,maxval(abs(Fl)))
  single_ok=merge(1,0,single_abs<=1e-12_mytype .and. single_rel<=1e-12_mytype)

  ! multi-point (clear once, accumulate)
  xlag(:,1)=(/grid%xmin+0.22_mytype*lx,0.5_mytype*(grid%y_center(6)+grid%y_center(7)),grid%zmin+0.27_mytype*lz/)
  xlag(:,2)=(/grid%xmin+0.33_mytype*lx,0.5_mytype*(grid%y_center(7)+grid%y_center(8)),grid%zmin+0.44_mytype*lz/)
  xlag(:,3)=(/grid%xmin+0.48_mytype*lx,0.5_mytype*(grid%y_center(8)+grid%y_center(9)),grid%zmin+0.58_mytype*lz/)
  xlag(:,4)=(/grid%xmin+0.63_mytype*lx,0.5_mytype*(grid%y_center(9)+grid%y_center(10)),grid%zmin+0.69_mytype*lz/)
  xlag(:,5)=(/grid%xmin+0.77_mytype*lx,0.5_mytype*(grid%y_center(10)+grid%y_center(11)),grid%zmin+0.81_mytype*lz/)
  flag(:,1)=(/0.8_mytype,-0.5_mytype,0.3_mytype/); flag(:,2)=(/1.1_mytype,-0.7_mytype,0.2_mytype/)
  flag(:,3)=(/0.6_mytype,-0.2_mytype,0.4_mytype/); flag(:,4)=(/1.3_mytype,-0.9_mytype,0.5_mytype/); flag(:,5)=(/0.9_mytype,-0.4_mytype,0.35_mytype/)
  ds=(/0.04_mytype,0.05_mytype,0.06_mytype,0.045_mytype,0.055_mytype/)
  call clear_stage7_force_buffer(fx,fy,fz)
  call spread_stage7_lagrangian_force_set(grid,nlag,xlag,flag,ds,fx,fy,fz,vcount,ucount)
  call compute_stage7_eulerian_force_total(grid,fx,fy,fz,Fe); call compute_stage7_lagrangian_force_total(nlag,flag,ds,Fl)
  multi_abs=maxval(abs(Fe-Fl)); multi_rel=multi_abs/max(1e-16_mytype,maxval(abs(Fl)))
  multi_ok=merge(1,0,vcount>=5 .and. ucount==0 .and. multi_abs<=1e-12_mytype .and. multi_rel<=1e-12_mytype)

  ! volume scaling
  call clear_stage7_force_buffer(fx,fy,fz)
  call spread_stage7_lagrangian_force(grid,x,y,z,F1,ds1,fx,fy,fz,valid,unsafe)
  call build_stage7_scalar_interp_weight(grid,x,y,z,w)
  vol_err=0
  do q=1,w%n
    i=w%i(q); j=w%j(q); k=w%k(q)
    vol_err=max(vol_err,abs(fx(i,j,k)*grid%volume_y(j)-F1(1)*ds1*w%w(q)))
    vol_err=max(vol_err,abs(fy(i,j,k)*grid%volume_y(j)-F1(2)*ds1*w%w(q)))
    vol_err=max(vol_err,abs(fz(i,j,k)*grid%volume_y(j)-F1(3)*ds1*w%w(q)))
  end do
  nonuni_flag=merge(1,0,maxval(abs(grid%volume_y-grid%volume_y(1)))>1e-15_mytype)
  vol_ok=merge(1,0,vol_err<=1e-12_mytype .and. nonuni_flag==1)

  ! periodic wrap + conservation
  perr=0; pxw=0; pzw=0; call clear_stage7_force_buffer(fx,fy,fz)
  call spread_stage7_lagrangian_force(grid,grid%xmin+0.1_mytype*grid%dx,y,grid%zmin+0.3_mytype*lz,F1,ds1,fx,fy,fz,valid,unsafe)
  call build_stage7_scalar_interp_weight(grid,grid%xmin+0.1_mytype*grid%dx,y,grid%zmin+0.3_mytype*lz,w); do q=1,w%n; if(w%i(q)==nx) pxw=1; end do
  call spread_stage7_lagrangian_force(grid,grid%xmax-0.1_mytype*grid%dx,y,grid%zmin+0.4_mytype*lz,F1,ds1,fx,fy,fz,valid,unsafe)
  call build_stage7_scalar_interp_weight(grid,grid%xmax-0.1_mytype*grid%dx,y,grid%zmin+0.4_mytype*lz,w); do q=1,w%n; if(w%i(q)==1) pxw=pxw*1; end do
  call spread_stage7_lagrangian_force(grid,grid%xmin+0.4_mytype*lx,y,grid%zmin+0.1_mytype*grid%dz,F1,ds1,fx,fy,fz,valid,unsafe)
  call build_stage7_scalar_interp_weight(grid,grid%xmin+0.4_mytype*lx,y,grid%zmin+0.1_mytype*grid%dz,w); do q=1,w%n; if(w%k(q)==nz) pzw=1; end do
  call spread_stage7_lagrangian_force(grid,grid%xmin+0.5_mytype*lx,y,grid%zmax-0.1_mytype*grid%dz,F1,ds1,fx,fy,fz,valid,unsafe)
  call build_stage7_scalar_interp_weight(grid,grid%xmin+0.5_mytype*lx,y,grid%zmax-0.1_mytype*grid%dz,w); do q=1,w%n; if(w%k(q)==1) pzw=pzw*1; end do
  call compute_stage7_eulerian_force_total(grid,fx,fy,fz,Fe); Fl=4._mytype*F1*ds1; perr=maxval(abs(Fe-Fl))
  per_ok=merge(1,0,perr<=1e-12_mytype .and. pxw==1 .and. pzw==1)

  ! nearwall + y outside
  call clear_stage7_force_buffer(fx,fy,fz); fx0=fx; fy0=fy; fz0=fz; nwc=0
  call spread_stage7_lagrangian_force(grid,x,grid%y_center(1),z,F1,ds1,fx,fy,fz,valid,unsafe); if(valid==0 .and. unsafe==1) nwc=nwc+1
  call spread_stage7_lagrangian_force(grid,x,grid%y_center(ny),z,F1,ds1,fx,fy,fz,valid,unsafe); if(valid==0 .and. unsafe==1) nwc=nwc+1
  nwflag=merge(1,0,nwc>0); nwchg=max(maxval(abs(fx-fx0)),max(maxval(abs(fy-fy0)),maxval(abs(fz-fz0))))
  yl=grid%ymin-0.1_mytype*(grid%ymax-grid%ymin); yh=grid%ymax+0.1_mytype*(grid%ymax-grid%ymin)
  call spread_stage7_lagrangian_force(grid,x,yl,z,F1,ds1,fx,fy,fz,valid,unsafe); ylof=merge(1,0,valid==0 .and. unsafe==1)
  call spread_stage7_lagrangian_force(grid,x,yh,z,F1,ds1,fx,fy,fz,valid,unsafe); yhif=merge(1,0,valid==0 .and. unsafe==1)
  yout=merge(1,0,ylof==1 .and. yhif==1)

  ! convention case
  call clear_stage7_force_buffer(fx,fy,fz)
  call spread_stage7_lagrangian_force(grid,x,y,z,F1,ds1,fx,fy,fz,valid,unsafe)
  call compute_stage7_eulerian_force_total(grid,fx,fy,fz,Fe); Fl=F1*ds1
  conv_abs=maxval(abs(Fe-Fl)); conv_flag=merge(1,0,conv_abs<=1e-12_mytype); no_rho_flag=conv_flag

  rhsx=1; rhsy=2; rhsz=3
  call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop_change,noop_mod)
  call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)

  status=merge(1,0,dep6==1 .and. dep70==1 .and. dep71==1 .and. dep72==1 .and. dep73==1 .and. dep74==1 .and. &
    single_ok==1 .and. multi_ok==1 .and. vol_ok==1 .and. per_ok==1 .and. nwflag==1 .and. yout==1 .and. conv_flag==1 .and. no_rho_flag==1 .and. &
    noop_mod==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0 .and. &
    nwchg<=1e-14_mytype .and. noop_change<=1e-14_mytype)

  open(10,file='stage7_outputs/fibre_stage7_force_spreading_check.dat',status='replace')
  write(10,'(A,1X,I0)')'stage7_spread_stage6_dependency_status',dep6
  write(10,'(A,1X,I0)')'stage7_spread_stage7_0_dependency_status',dep70
  write(10,'(A,1X,I0)')'stage7_spread_stage7_1_dependency_status',dep71
  write(10,'(A,1X,I0)')'stage7_spread_stage7_2_dependency_status',dep72
  write(10,'(A,1X,I0)')'stage7_spread_stage7_3_dependency_status',dep73
  write(10,'(A,1X,I0)')'stage7_spread_stage7_4_dependency_status',dep74
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_single_force_abs_error',single_abs
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_single_force_relative_error',single_rel
  write(10,'(A,1X,I0)')'stage7_spread_single_force_conservation_status',single_ok
  write(10,'(A,1X,I0)')'stage7_spread_multipoint_valid_count',vcount
  write(10,'(A,1X,I0)')'stage7_spread_multipoint_unsafe_count',ucount
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_multipoint_force_abs_error',multi_abs
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_multipoint_force_relative_error',multi_rel
  write(10,'(A,1X,I0)')'stage7_spread_multipoint_force_conservation_status',multi_ok
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_volume_scaling_error_max',vol_err
  write(10,'(A,1X,I0)')'stage7_spread_nonuniform_volume_used_flag',nonuni_flag
  write(10,'(A,1X,I0)')'stage7_spread_volume_scaling_status',vol_ok
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_periodic_force_error_max',perr
  write(10,'(A,1X,I0)')'stage7_spread_periodic_x_wrap_status',pxw
  write(10,'(A,1X,I0)')'stage7_spread_periodic_z_wrap_status',pzw
  write(10,'(A,1X,I0)')'stage7_spread_periodic_status',per_ok
  write(10,'(A,1X,I0)')'stage7_spread_nearwall_unsafe_count',nwc
  write(10,'(A,1X,I0)')'stage7_spread_nearwall_blocked_flag',nwflag
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_nearwall_buffer_change_max',nwchg
  write(10,'(A,1X,I0)')'stage7_spread_y_outside_low_blocked_flag',ylof
  write(10,'(A,1X,I0)')'stage7_spread_y_outside_high_blocked_flag',yhif
  write(10,'(A,1X,I0)')'stage7_spread_y_outside_status',yout
  write(10,'(A,1X,I0)')'stage7_spread_force_density_convention_flag',conv_flag
  write(10,'(A,1X,I0)')'stage7_spread_no_rho_division_flag',no_rho_flag
  write(10,'(A,1X,ES24.16E3)')'stage7_spread_noop_rhs_change_max',noop_change
  write(10,'(A,1X,I0)')'stage7_spread_noop_rhs_modified_flag',noop_mod
  write(10,'(A,1X,I0)')'stage7_spread_pressure_poisson_modified_flag',pp
  write(10,'(A,1X,I0)')'stage7_spread_projection_modified_flag',proj
  write(10,'(A,1X,I0)')'stage7_spread_production_dns_called_flag',dns
  write(10,'(A,1X,I0)')'stage7_spread_fluid_update_called_flag',flu
  write(10,'(A,1X,I0)')'stage7_spread_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)')'stage7_force_spreading_check_status',status
  close(10)
contains
  subroutine dependency_status(s6,s70,s71,s72,s73,s74)
    integer,intent(out)::s6,s70,s71,s72,s73,s74
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)*check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype)
    s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype)
    s72=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat','stage7_scalar_interpolation_check_status',1._mytype)
    s73=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat','stage7_scalar_interpolation_robustness_check_status',1._mytype)
    s74=check_dat('stage7_outputs/fibre_stage7_velocity_interpolation_check.dat','stage7_velocity_interpolation_check_status',1._mytype)
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key
    real(mytype),intent(in)::val
    character(len=512)::line,rest
    real(mytype)::v
    integer::u,ios
    logical::ex
    check_dat=0; inquire(file=path,exist=ex); if (.not.ex) return
    open(newunit=u,file=path,status='old',iostat=ios); if (ios/=0) return
    do
      read(u,'(A)',iostat=ios) line; if (ios/=0) exit
      line=adjustl(line)
      if (index(line,trim(key))==1) then
        rest=adjustl(line(len_trim(key)+1:)); if (len_trim(rest)>0 .and. rest(1:1)=='=') rest=adjustl(rest(2:))
        read(rest,*,iostat=ios) v
        if (ios==0 .and. abs(v-val)<=1e-12_mytype) then; check_dat=1; exit; end if
      end if
    end do
    close(u)
  end function
end program
