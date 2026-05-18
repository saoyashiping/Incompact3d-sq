program fibre_stage7_power_consistency_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_velocity_interpolation
  use fibre_stage7_force_spreading
  implicit none
  type(stage7_channel_grid_t)::grid
  type(stage7_velocity_layout_t)::layout
  type(stage7_interp_weight_t)::w
  integer,parameter::nx=16,ny=17,nz=12,nlag=5
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),fx0(:,:,:),fy0(:,:,:),fz0(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
  real(mytype)::x,y,z,lx,lz,u_lag(3),F(3),ds1,pe,pl,se,sr,me,mr,ppe,ppl,pee,vol_err,nwchg,noop_change,yl,yh
  real(mytype)::xlag(3,nlag),flag(3,nlag),ds(nlag),ulag(3,nlag)
  integer::i,j,k,l,q,valid,unsafe,rej,vcount,ucount,dep6,dep70,dep71,dep72,dep73,dep74,dep75,pxw,pzw,nwc,nwflag,ylof,yhif,yout,pp,proj,rp,dns,flu,fib,noop_mod
  integer::single_ok,multi_ok,per_ok,vol_ok,conv_flag,no_rho_flag,rhs_hook,status,nonuni
  call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz); call validate_stage7_channel_grid(grid,valid,rej); call init_stage7_collocated_velocity_layout(layout); call init_stage7_interp_weight(w,64)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),fx0(nx,ny,nz),fy0(nx,ny,nz),fz0(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz))
  call dependency_status(dep6,dep70,dep71,dep72,dep73,dep74,dep75)
  lx=grid%xmax-grid%xmin; lz=grid%zmax-grid%zmin
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=1._mytype+0.2_mytype*(grid%xmin+real(i-1,mytype)*grid%dx)-0.1_mytype*grid%y_center(j)+0.05_mytype*(grid%zmin+real(k-1,mytype)*grid%dz)
    uy(i,j,k)=-0.5_mytype+0.1_mytype*(grid%xmin+real(i-1,mytype)*grid%dx)+0.3_mytype*grid%y_center(j)-0.02_mytype*(grid%zmin+real(k-1,mytype)*grid%dz)
    uz(i,j,k)=0.25_mytype-0.05_mytype*(grid%xmin+real(i-1,mytype)*grid%dx)+0.07_mytype*grid%y_center(j)+0.15_mytype*(grid%zmin+real(k-1,mytype)*grid%dz)
  end do; end do; end do
  x=grid%xmin+0.37_mytype*lx; y=0.5_mytype*(grid%y_center(8)+grid%y_center(9)); z=grid%zmin+0.41_mytype*lz; F=(/1.2_mytype,-0.7_mytype,0.35_mytype/); ds1=0.05_mytype
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,x,y,z,u_lag,valid,unsafe); call clear_stage7_force_buffer(fx,fy,fz); call spread_stage7_lagrangian_force(grid,x,y,z,F,ds1,fx,fy,fz,valid,unsafe)
  pe=compute_power(grid,fx,fy,fz,ux,uy,uz); pl=dot_product(F,u_lag)*ds1; se=abs(pe-pl); sr=se/max(1e-16_mytype,abs(pl)); single_ok=merge(1,0,se<=1e-12_mytype .and. sr<=1e-12_mytype)
  xlag(:,1)=(/grid%xmin+0.22_mytype*lx,0.5_mytype*(grid%y_center(6)+grid%y_center(7)),grid%zmin+0.27_mytype*lz/); xlag(:,2)=(/grid%xmin+0.33_mytype*lx,0.5_mytype*(grid%y_center(7)+grid%y_center(8)),grid%zmin+0.44_mytype*lz/); xlag(:,3)=(/grid%xmin+0.48_mytype*lx,0.5_mytype*(grid%y_center(8)+grid%y_center(9)),grid%zmin+0.58_mytype*lz/); xlag(:,4)=(/grid%xmin+0.63_mytype*lx,0.5_mytype*(grid%y_center(9)+grid%y_center(10)),grid%zmin+0.69_mytype*lz/); xlag(:,5)=(/grid%xmin+0.77_mytype*lx,0.5_mytype*(grid%y_center(10)+grid%y_center(11)),grid%zmin+0.81_mytype*lz/)
  flag(:,1)=(/0.8_mytype,-0.5_mytype,0.3_mytype/); flag(:,2)=(/1.1_mytype,-0.7_mytype,0.2_mytype/); flag(:,3)=(/0.6_mytype,-0.2_mytype,0.4_mytype/); flag(:,4)=(/1.3_mytype,-0.9_mytype,0.5_mytype/); flag(:,5)=(/0.9_mytype,-0.4_mytype,0.35_mytype/); ds=(/0.04_mytype,0.05_mytype,0.06_mytype,0.045_mytype,0.055_mytype/)
  call clear_stage7_force_buffer(fx,fy,fz); vcount=0; ucount=0; ppl=0
  do l=1,nlag; call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xlag(1,l),xlag(2,l),xlag(3,l),ulag(:,l),valid,unsafe); if(valid==1) then; vcount=vcount+1; ppl=ppl+dot_product(flag(:,l),ulag(:,l))*ds(l); else; ucount=ucount+1; endif; call spread_stage7_lagrangian_force(grid,xlag(1,l),xlag(2,l),xlag(3,l),flag(:,l),ds(l),fx,fy,fz,valid,unsafe); end do
  pee=compute_power(grid,fx,fy,fz,ux,uy,uz); me=abs(pee-ppl); mr=me/max(1e-16_mytype,abs(ppl)); multi_ok=merge(1,0,vcount>=5 .and. ucount==0 .and. me<=1e-12_mytype .and. mr<=1e-12_mytype)
  ! periodic
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=sin(2._mytype*stage7_pi*((grid%xmin+real(i-1,mytype)*grid%dx)-grid%xmin)/lx)+0.2_mytype*grid%y_center(j)
    uy(i,j,k)=cos(2._mytype*stage7_pi*((grid%zmin+real(k-1,mytype)*grid%dz)-grid%zmin)/lz)-0.1_mytype*grid%y_center(j)
    uz(i,j,k)=sin(2._mytype*stage7_pi*((grid%xmin+real(i-1,mytype)*grid%dx)-grid%xmin)/lx)+cos(2._mytype*stage7_pi*((grid%zmin+real(k-1,mytype)*grid%dz)-grid%zmin)/lz)
  end do; end do; end do
  call clear_stage7_force_buffer(fx,fy,fz); pzw=0; pxw=0; ppl=0
  xlag(:,1)=(/grid%xmin+0.1_mytype*grid%dx,y,grid%zmin+0.3_mytype*lz/); xlag(:,2)=(/grid%xmax-0.1_mytype*grid%dx,y,grid%zmin+0.4_mytype*lz/); xlag(:,3)=(/grid%xmin+0.4_mytype*lx,y,grid%zmin+0.1_mytype*grid%dz/); xlag(:,4)=(/grid%xmin+0.5_mytype*lx,y,grid%zmax-0.1_mytype*grid%dz/); xlag(:,5)=(/grid%xmin+0.35_mytype*lx,y,grid%zmin+0.6_mytype*lz/)
  do l=1,nlag; call interpolate_stage7_velocity(grid,layout,ux,uy,uz,xlag(1,l),xlag(2,l),xlag(3,l),ulag(:,l),valid,unsafe); ppl=ppl+dot_product(flag(:,l),ulag(:,l))*ds(l); call spread_stage7_lagrangian_force(grid,xlag(1,l),xlag(2,l),xlag(3,l),flag(:,l),ds(l),fx,fy,fz,valid,unsafe); call build_stage7_scalar_interp_weight(grid,xlag(1,l),xlag(2,l),xlag(3,l),w); do q=1,w%n; if(l<=2 .and. (w%i(q)==1 .or. w%i(q)==nx)) pxw=1; if(l>=3 .and. l<=4 .and. (w%k(q)==1 .or. w%k(q)==nz)) pzw=1; end do; end do
  ppe=compute_power(grid,fx,fy,fz,ux,uy,uz); pe=abs(ppe-ppl); pl=pe/max(1e-16_mytype,abs(ppl)); per_ok=merge(1,0,pe<=1e-12_mytype .and. pl<=1e-12_mytype .and. pxw==1 .and. pzw==1)
  ! volume audit
  call clear_stage7_force_buffer(fx,fy,fz); call spread_stage7_lagrangian_force(grid,x,y,z,F,ds1,fx,fy,fz,valid,unsafe); call build_stage7_scalar_interp_weight(grid,x,y,z,w); vol_err=0; do q=1,w%n; i=w%i(q); j=w%j(q); k=w%k(q); vol_err=max(vol_err,abs(fx(i,j,k)*grid%volume_y(j)-F(1)*ds1*w%w(q))); vol_err=max(vol_err,abs(fy(i,j,k)*grid%volume_y(j)-F(2)*ds1*w%w(q))); vol_err=max(vol_err,abs(fz(i,j,k)*grid%volume_y(j)-F(3)*ds1*w%w(q))); end do; nonuni=merge(1,0,maxval(abs(grid%volume_y-grid%volume_y(1)))>1e-15_mytype); vol_ok=merge(1,0,nonuni==1 .and. vol_err<=1e-12_mytype)
  ! nearwall/outside
  call clear_stage7_force_buffer(fx,fy,fz); fx0=fx; fy0=fy; fz0=fz; nwc=0
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,x,grid%y_center(1),z,u_lag,valid,unsafe); if(valid==0 .and. unsafe==1) nwc=nwc+1; call spread_stage7_lagrangian_force(grid,x,grid%y_center(1),z,F,ds1,fx,fy,fz,valid,unsafe); if(valid==0 .and. unsafe==1) nwc=nwc+1
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,x,grid%y_center(ny),z,u_lag,valid,unsafe); if(valid==0 .and. unsafe==1) nwc=nwc+1; call spread_stage7_lagrangian_force(grid,x,grid%y_center(ny),z,F,ds1,fx,fy,fz,valid,unsafe); if(valid==0 .and. unsafe==1) nwc=nwc+1
  nwflag=merge(1,0,nwc>0); nwchg=max(maxval(abs(fx-fx0)),max(maxval(abs(fy-fy0)),maxval(abs(fz-fz0))))
  yl=grid%ymin-0.1_mytype*(grid%ymax-grid%ymin); yh=grid%ymax+0.1_mytype*(grid%ymax-grid%ymin)
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,x,yl,z,u_lag,valid,unsafe); ylof=merge(1,0,valid==0 .and. unsafe==1); call spread_stage7_lagrangian_force(grid,x,yl,z,F,ds1,fx,fy,fz,valid,unsafe); ylof=min(ylof,merge(1,0,valid==0 .and. unsafe==1))
  call interpolate_stage7_velocity(grid,layout,ux,uy,uz,x,yh,z,u_lag,valid,unsafe); yhif=merge(1,0,valid==0 .and. unsafe==1); call spread_stage7_lagrangian_force(grid,x,yh,z,F,ds1,fx,fy,fz,valid,unsafe); yhif=min(yhif,merge(1,0,valid==0 .and. unsafe==1))
  yout=merge(1,0,ylof==1 .and. yhif==1)
  conv_flag=single_ok; no_rho_flag=single_ok; rhs_hook=0
  rhsx=1; rhsy=2; rhsz=3; call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop_change,noop_mod); call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)
  status=merge(1,0,dep6==1 .and. dep70==1 .and. dep71==1 .and. dep72==1 .and. dep73==1 .and. dep74==1 .and. dep75==1 .and. single_ok==1 .and. multi_ok==1 .and. per_ok==1 .and. vol_ok==1 .and. nwflag==1 .and. yout==1 .and. conv_flag==1 .and. no_rho_flag==1 .and. rhs_hook==0 .and. noop_mod==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0 .and. nwchg<=1e-14_mytype .and. noop_change<=1e-14_mytype)
  open(10,file='stage7_outputs/fibre_stage7_power_consistency_check.dat',status='replace')
  write(10,'(A,1X,I0)')'stage7_power_stage6_dependency_status',dep6; write(10,'(A,1X,I0)')'stage7_power_stage7_0_dependency_status',dep70; write(10,'(A,1X,I0)')'stage7_power_stage7_1_dependency_status',dep71; write(10,'(A,1X,I0)')'stage7_power_stage7_2_dependency_status',dep72; write(10,'(A,1X,I0)')'stage7_power_stage7_3_dependency_status',dep73; write(10,'(A,1X,I0)')'stage7_power_stage7_4_dependency_status',dep74; write(10,'(A,1X,I0)')'stage7_power_stage7_5_dependency_status',dep75
  write(10,'(A,1X,ES24.16E3)')'stage7_power_single_eulerian_power',pe; write(10,'(A,1X,ES24.16E3)')'stage7_power_single_lagrangian_power',pl; write(10,'(A,1X,ES24.16E3)')'stage7_power_single_abs_error',se; write(10,'(A,1X,ES24.16E3)')'stage7_power_single_relative_error',sr; write(10,'(A,1X,I0)')'stage7_power_single_consistency_status',single_ok
  write(10,'(A,1X,I0)')'stage7_power_multipoint_valid_count',vcount; write(10,'(A,1X,I0)')'stage7_power_multipoint_unsafe_count',ucount; write(10,'(A,1X,ES24.16E3)')'stage7_power_multipoint_eulerian_power',pee; write(10,'(A,1X,ES24.16E3)')'stage7_power_multipoint_lagrangian_power',ppl; write(10,'(A,1X,ES24.16E3)')'stage7_power_multipoint_abs_error',me; write(10,'(A,1X,ES24.16E3)')'stage7_power_multipoint_relative_error',mr; write(10,'(A,1X,I0)')'stage7_power_multipoint_consistency_status',multi_ok
  write(10,'(A,1X,ES24.16E3)')'stage7_power_periodic_abs_error',pe; write(10,'(A,1X,ES24.16E3)')'stage7_power_periodic_relative_error',pl; write(10,'(A,1X,I0)')'stage7_power_periodic_x_wrap_status',pxw; write(10,'(A,1X,I0)')'stage7_power_periodic_z_wrap_status',pzw; write(10,'(A,1X,I0)')'stage7_power_periodic_consistency_status',per_ok
  write(10,'(A,1X,I0)')'stage7_power_nonuniform_volume_used_flag',nonuni; write(10,'(A,1X,ES24.16E3)')'stage7_power_volume_weighting_error_max',vol_err; write(10,'(A,1X,I0)')'stage7_power_volume_weighting_status',vol_ok
  write(10,'(A,1X,I0)')'stage7_power_nearwall_unsafe_count',nwc; write(10,'(A,1X,I0)')'stage7_power_nearwall_blocked_flag',nwflag; write(10,'(A,1X,ES24.16E3)')'stage7_power_nearwall_buffer_change_max',nwchg
  write(10,'(A,1X,I0)')'stage7_power_y_outside_low_blocked_flag',ylof; write(10,'(A,1X,I0)')'stage7_power_y_outside_high_blocked_flag',yhif; write(10,'(A,1X,I0)')'stage7_power_y_outside_status',yout
  write(10,'(A,1X,I0)')'stage7_power_force_density_convention_flag',conv_flag; write(10,'(A,1X,I0)')'stage7_power_no_rho_division_flag',no_rho_flag; write(10,'(A,1X,I0)')'stage7_power_stage6_rhs_hook_called_flag',rhs_hook
  write(10,'(A,1X,ES24.16E3)')'stage7_power_noop_rhs_change_max',noop_change; write(10,'(A,1X,I0)')'stage7_power_noop_rhs_modified_flag',noop_mod; write(10,'(A,1X,I0)')'stage7_power_pressure_poisson_modified_flag',pp; write(10,'(A,1X,I0)')'stage7_power_projection_modified_flag',proj; write(10,'(A,1X,I0)')'stage7_power_production_dns_called_flag',dns; write(10,'(A,1X,I0)')'stage7_power_fluid_update_called_flag',flu; write(10,'(A,1X,I0)')'stage7_power_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)')'stage7_power_consistency_check_status',status; close(10)
contains
  real(mytype) function compute_power(grid,fx,fy,fz,ux,uy,uz)
    type(stage7_channel_grid_t),intent(in)::grid; real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:)
    integer::i,j,k; compute_power=0
    do k=1,size(fx,3); do j=1,size(fx,2); do i=1,size(fx,1); compute_power=compute_power+(fx(i,j,k)*ux(i,j,k)+fy(i,j,k)*uy(i,j,k)+fz(i,j,k)*uz(i,j,k))*grid%volume_y(j); end do; end do; end do
  end function
  subroutine dependency_status(s6,s70,s71,s72,s73,s74,s75)
    integer,intent(out)::s6,s70,s71,s72,s73,s74,s75
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)*check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype)
    s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype)
    s72=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat','stage7_scalar_interpolation_check_status',1._mytype)
    s73=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat','stage7_scalar_interpolation_robustness_check_status',1._mytype)
    s74=check_dat('stage7_outputs/fibre_stage7_velocity_interpolation_check.dat','stage7_velocity_interpolation_check_status',1._mytype)
    s75=check_dat('stage7_outputs/fibre_stage7_force_spreading_check.dat','stage7_force_spreading_check_status',1._mytype)
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key; real(mytype),intent(in)::val; character(len=512)::line,rest; real(mytype)::v; integer::u,ios; logical::ex
    check_dat=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',iostat=ios); if(ios/=0)return
    do; read(u,'(A)',iostat=ios)line; if(ios/=0)exit; line=adjustl(line); if(index(line,trim(key))==1) then; rest=adjustl(line(len_trim(key)+1:)); if(len_trim(rest)>0 .and. rest(1:1)=='=') rest=adjustl(rest(2:)); read(rest,*,iostat=ios)v; if(ios==0 .and. abs(v-val)<=1e-12_mytype) then; check_dat=1; exit; endif; endif; enddo; close(u)
  end function
end program
