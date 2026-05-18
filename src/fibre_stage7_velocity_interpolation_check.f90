program fibre_stage7_velocity_interpolation_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_velocity_interpolation
  implicit none
  type(stage7_channel_grid_t)::grid
  type(stage7_velocity_layout_t)::layout,layout_bad
  integer,parameter::nx=16,ny=17,nz=12,np=6
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
  real(mytype)::xp(np),yp(np),zp(np),u_lag(3),u_ref(3),err
  real(mytype)::const_err,lin_err,poi_err,ue,ve,we,px_err,pz_err,noop_change
  real(mytype)::lx,lz,xi,yj,zk,ylow,yhigh
  integer::i,j,k,p,valid,rej,unsafe
  integer::dep6,dep70,dep71,dep72,dep73,col_ok,comp_ok,bad_rej,layout_ok
  integer::const_ok,lin_ok,poi_ok,comp_cons_ok,px_ok,pz_ok,nw_count,nw_ok,ylo_ok,yhi_ok,yout_ok
  integer::noop_mod,pp,proj,rp,dns,flu,fib,status

  call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz)
  call validate_stage7_channel_grid(grid,valid,rej)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz))
  call dependency_status(dep6,dep70,dep71,dep72,dep73)

  lx=grid%xmax-grid%xmin; lz=grid%zmax-grid%zmin
  xp=(/grid%xmin+3.2_mytype*grid%dx,grid%xmin+5.4_mytype*grid%dx,grid%xmin+7.1_mytype*grid%dx,grid%xmin+9.3_mytype*grid%dx,grid%xmin+11.2_mytype*grid%dx,grid%xmin+13.0_mytype*grid%dx/)
  yp=(/grid%y_center(6),grid%y_center(7),grid%y_center(8),grid%y_center(9),grid%y_center(10),grid%y_center(11)/)
  zp=(/grid%zmin+2.1_mytype*grid%dz,grid%zmin+3.2_mytype*grid%dz,grid%zmin+4.4_mytype*grid%dz,grid%zmin+6.1_mytype*grid%dz,grid%zmin+7.6_mytype*grid%dz,grid%zmin+9.1_mytype*grid%dz/)

  call init_stage7_collocated_velocity_layout(layout)
  call validate_stage7_velocity_layout(layout,valid,rej); col_ok=valid
  call init_stage7_component_velocity_layout(layout)
  call validate_stage7_velocity_layout(layout,valid,rej); comp_ok=valid
  call validate_stage7_velocity_layout(layout_bad,valid,rej); bad_rej=rej
  layout_ok=merge(1,0,col_ok==1 .and. comp_ok==1 .and. bad_rej==1)

  ux=1.25_mytype; uy=-0.5_mytype; uz=0.75_mytype; const_err=0
  call init_stage7_collocated_velocity_layout(layout)
  do p=1,np
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p),yp(p),zp(p),u_lag,valid,unsafe)
    err=max(abs(u_lag(1)-1.25_mytype),max(abs(u_lag(2)+0.5_mytype),abs(u_lag(3)-0.75_mytype)))
    if (valid/=1 .or. unsafe/=0) err=1._mytype
    const_err=max(const_err,err)
  end do
  const_ok=merge(1,0,const_err<=1e-12_mytype)

  do k=1,nz; zk=grid%zmin+real(k-1,mytype)*grid%dz
    do j=1,ny; yj=grid%y_center(j)
      do i=1,nx; xi=grid%xmin+real(i-1,mytype)*grid%dx
        ux(i,j,k)=1._mytype+0.2_mytype*xi-0.1_mytype*yj+0.05_mytype*zk
        uy(i,j,k)=-0.5_mytype+0.1_mytype*xi+0.3_mytype*yj-0.02_mytype*zk
        uz(i,j,k)=0.25_mytype-0.05_mytype*xi+0.07_mytype*yj+0.15_mytype*zk
      end do
    end do
  end do
  lin_err=0
  do p=1,4
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p),yp(p),zp(p),u_lag,valid,unsafe)
    u_ref(1)=1._mytype+0.2_mytype*xp(p)-0.1_mytype*yp(p)+0.05_mytype*zp(p)
    u_ref(2)=-0.5_mytype+0.1_mytype*xp(p)+0.3_mytype*yp(p)-0.02_mytype*zp(p)
    u_ref(3)=0.25_mytype-0.05_mytype*xp(p)+0.07_mytype*yp(p)+0.15_mytype*zp(p)
    err=maxval(abs(u_lag-u_ref)); if (valid/=1 .or. unsafe/=0) err=1._mytype; lin_err=max(lin_err,err)
  end do
  lin_ok=merge(1,0,lin_err<=1e-11_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=1._mytype-grid%y_center(j)**2; uy(i,j,k)=0._mytype; uz(i,j,k)=0._mytype
  end do; end do; end do
  poi_err=0
  do p=1,np
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p),yp(p),zp(p),u_lag,valid,unsafe)
    u_ref=(/1._mytype-yp(p)**2,0._mytype,0._mytype/)
    err=maxval(abs(u_lag-u_ref)); if (valid/=1 .or. unsafe/=0) err=1._mytype; poi_err=max(poi_err,err)
  end do
  poi_ok=merge(1,0,poi_err<=1e-11_mytype)

  ue=0;ve=0;we=0
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=0.3_mytype+0.1_mytype*grid%y_center(j)+0.05_mytype*grid%y_center(j)**2
    uy(i,j,k)=-0.2_mytype+0.2_mytype*grid%y_center(j)-0.03_mytype*grid%y_center(j)**2
    uz(i,j,k)=0.1_mytype-0.1_mytype*grid%y_center(j)+0.04_mytype*grid%y_center(j)**2
  end do; end do; end do
  do p=1,np
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p),yp(p),zp(p),u_lag,valid,unsafe)
    ue=max(ue,abs(u_lag(1)-(0.3_mytype+0.1_mytype*yp(p)+0.05_mytype*yp(p)**2)))
    ve=max(ve,abs(u_lag(2)-(-0.2_mytype+0.2_mytype*yp(p)-0.03_mytype*yp(p)**2)))
    we=max(we,abs(u_lag(3)-(0.1_mytype-0.1_mytype*yp(p)+0.04_mytype*yp(p)**2)))
  end do
  comp_cons_ok=merge(1,0,ue<=1e-11_mytype .and. ve<=1e-11_mytype .and. we<=1e-11_mytype)

  do k=1,nz; zk=grid%zmin+real(k-1,mytype)*grid%dz
    do j=1,ny; yj=grid%y_center(j)
      do i=1,nx; xi=grid%xmin+real(i-1,mytype)*grid%dx
        ux(i,j,k)=sin(2._mytype*stage7_pi*(xi-grid%xmin)/lx)+0.2_mytype*yj
        uy(i,j,k)=cos(2._mytype*stage7_pi*(zk-grid%zmin)/lz)-0.1_mytype*yj
        uz(i,j,k)=sin(2._mytype*stage7_pi*(xi-grid%xmin)/lx)+cos(2._mytype*stage7_pi*(zk-grid%zmin)/lz)
      end do
    end do
  end do
  px_err=0; pz_err=0
  do p=1,3
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p),yp(p),zp(p),u_lag,valid,unsafe); u_ref=u_lag
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p)+lx,yp(p),zp(p),u_lag,valid,unsafe); px_err=max(px_err,maxval(abs(u_lag-u_ref)))
    call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(p),yp(p),zp(p)+lz,u_lag,valid,unsafe); pz_err=max(pz_err,maxval(abs(u_lag-u_ref)))
  end do
  px_ok=merge(1,0,px_err<=1e-12_mytype); pz_ok=merge(1,0,pz_err<=1e-12_mytype)

  nw_count=0
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(1),grid%y_center(1),zp(1),u_lag,valid,unsafe); if(valid==0 .and. unsafe==1) nw_count=nw_count+1
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(1),grid%y_center(ny),zp(1),u_lag,valid,unsafe); if(valid==0 .and. unsafe==1) nw_count=nw_count+1
  nw_ok=merge(1,0,nw_count>0)

  ylow=grid%ymin-0.1_mytype*(grid%ymax-grid%ymin); yhigh=grid%ymax+0.1_mytype*(grid%ymax-grid%ymin)
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(1),ylow,zp(1),u_lag,valid,unsafe); ylo_ok=merge(1,0,valid==0 .and. unsafe==1)
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xp(1),yhigh,zp(1),u_lag,valid,unsafe); yhi_ok=merge(1,0,valid==0 .and. unsafe==1)
  yout_ok=merge(1,0,ylo_ok==1 .and. yhi_ok==1)

  rhsx=1; rhsy=2; rhsz=3
  call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop_change,noop_mod)
  call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)

  status=merge(1,0,dep6==1 .and. dep70==1 .and. dep71==1 .and. dep72==1 .and. dep73==1 .and. layout_ok==1 .and. const_ok==1 .and. lin_ok==1 .and. poi_ok==1 .and. comp_cons_ok==1 .and. px_ok==1 .and. pz_ok==1 .and. nw_ok==1 .and. yout_ok==1 .and. noop_change<=1e-14_mytype .and. noop_mod==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0)

  open(10,file='stage7_outputs/fibre_stage7_velocity_interpolation_check.dat',status='replace')
  write(10,'(A,1X,I0)') 'stage7_vel_stage6_dependency_status',dep6
  write(10,'(A,1X,I0)') 'stage7_vel_stage7_0_dependency_status',dep70
  write(10,'(A,1X,I0)') 'stage7_vel_stage7_1_dependency_status',dep71
  write(10,'(A,1X,I0)') 'stage7_vel_stage7_2_dependency_status',dep72
  write(10,'(A,1X,I0)') 'stage7_vel_stage7_3_dependency_status',dep73
  write(10,'(A,1X,I0)') 'stage7_vel_collocated_layout_valid_flag',col_ok
  write(10,'(A,1X,I0)') 'stage7_vel_component_layout_valid_flag',comp_ok
  write(10,'(A,1X,I0)') 'stage7_vel_invalid_layout_rejected_flag',bad_rej
  write(10,'(A,1X,I0)') 'stage7_vel_layout_validation_status',layout_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_constant_error_max',const_err
  write(10,'(A,1X,I0)') 'stage7_vel_constant_status',const_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_linear_error_max',lin_err
  write(10,'(A,1X,I0)') 'stage7_vel_linear_status',lin_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_poiseuille_error_max',poi_err
  write(10,'(A,1X,I0)') 'stage7_vel_poiseuille_status',poi_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_u_component_error_max',ue
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_v_component_error_max',ve
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_w_component_error_max',we
  write(10,'(A,1X,I0)') 'stage7_vel_component_consistency_status',comp_cons_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_periodic_x_shift_error_max',px_err
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_periodic_z_shift_error_max',pz_err
  write(10,'(A,1X,I0)') 'stage7_vel_periodic_x_shift_status',px_ok
  write(10,'(A,1X,I0)') 'stage7_vel_periodic_z_shift_status',pz_ok
  write(10,'(A,1X,I0)') 'stage7_vel_nearwall_unsafe_count',nw_count
  write(10,'(A,1X,I0)') 'stage7_vel_nearwall_blocked_flag',nw_ok
  write(10,'(A,1X,I0)') 'stage7_vel_y_outside_low_blocked_flag',ylo_ok
  write(10,'(A,1X,I0)') 'stage7_vel_y_outside_high_blocked_flag',yhi_ok
  write(10,'(A,1X,I0)') 'stage7_vel_y_outside_status',yout_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_vel_noop_rhs_change_max',noop_change
  write(10,'(A,1X,I0)') 'stage7_vel_noop_rhs_modified_flag',noop_mod
  write(10,'(A,1X,I0)') 'stage7_vel_pressure_poisson_modified_flag',pp
  write(10,'(A,1X,I0)') 'stage7_vel_projection_modified_flag',proj
  write(10,'(A,1X,I0)') 'stage7_vel_production_dns_called_flag',dns
  write(10,'(A,1X,I0)') 'stage7_vel_fluid_update_called_flag',flu
  write(10,'(A,1X,I0)') 'stage7_vel_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)') 'stage7_velocity_interpolation_check_status',status
  close(10)
contains
  subroutine dependency_status(s6,s70,s71,s72,s73)
    integer,intent(out)::s6,s70,s71,s72,s73
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)*check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype)
    s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype)
    s72=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat','stage7_scalar_interpolation_check_status',1._mytype)
    s73=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat','stage7_scalar_interpolation_robustness_check_status',1._mytype)
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key
    real(mytype),intent(in)::val
    character(len=512)::line,rest
    real(mytype)::v
    integer::u,ios
    check_dat=0; if (.not.file_exists(path)) return
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
  logical function file_exists(path)
    character(*),intent(in)::path
    inquire(file=path,exist=file_exists)
  end function
end program
