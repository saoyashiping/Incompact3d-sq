program fibre_stage7_boundary_safety_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_velocity_interpolation
  use fibre_stage7_force_spreading
  use fibre_stage7_boundary_safety
  use, intrinsic :: ieee_arithmetic
  implicit none
  type(stage7_channel_grid_t)::grid; type(stage7_velocity_layout_t)::layout,layout_bad; type(stage7_boundary_safety_result_t)::res
  integer,parameter::nx=16,ny=17,nz=12,np=5
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),fx0(:,:,:),fy0(:,:,:),fz0(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
  real(mytype)::x(np),y(np),z(np),u_lag(3),chg1,chg2,noop
  real(mytype)::lx,lz,invalid_y
  integer::i,j,k,p,v,r,u,dep6,dep70,dep71,dep72,dep73,dep74,dep75,dep76,safe_count,safe_ok,nw_count,nw_ok,ylo,yhi,yout,invcoord,invstatus,invlay,invlay_status,px,pz,pstat,block_no_write,noop_mod,pp,proj,rp,dns,flu,fib,status
  call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz); call init_stage7_collocated_velocity_layout(layout)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),fx0(nx,ny,nz),fy0(nx,ny,nz),fz0(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz))
  call dependency_status(dep6,dep70,dep71,dep72,dep73,dep74,dep75,dep76)
  lx=grid%xmax-grid%xmin; lz=grid%zmax-grid%zmin
  do k=1,nz;do j=1,ny;do i=1,nx;ux(i,j,k)=1;uy(i,j,k)=2;uz(i,j,k)=3;enddo;enddo;enddo
  x=(/grid%xmin+3.2_mytype*grid%dx,grid%xmin+5.1_mytype*grid%dx,grid%xmin+7.4_mytype*grid%dx,grid%xmin+9.3_mytype*grid%dx,grid%xmin+12.2_mytype*grid%dx/)
  y=(/grid%y_center(6),grid%y_center(7),grid%y_center(8),grid%y_center(9),grid%y_center(10)/)
  z=(/grid%zmin+2.1_mytype*grid%dz,grid%zmin+3.2_mytype*grid%dz,grid%zmin+4.2_mytype*grid%dz,grid%zmin+6.2_mytype*grid%dz,grid%zmin+7.2_mytype*grid%dz/)
  safe_count=0; call clear_stage7_force_buffer(fx,fy,fz)
  do p=1,np
    call classify_stage7_boundary_point(grid,layout,x(p),y(p),z(p),1,res)
    if(res%safe_flag==1) then
      safe_count=safe_count+1; call interpolate_stage7_velocity(grid,layout,ux,uy,uz,x(p),y(p),z(p),u_lag,v,u); call spread_stage7_lagrangian_force(grid,x(p),y(p),z(p),(/1._mytype,0.5_mytype,0.2_mytype/),0.01_mytype,fx,fy,fz,v,u)
    end if
  end do
  safe_ok=merge(1,0,safe_count>=5)
  call clear_stage7_force_buffer(fx,fy,fz); fx0=fx; fy0=fy; fz0=fz; nw_count=0
  call classify_stage7_boundary_point(grid,layout,x(1),grid%y_center(1),z(1),1,res); if(res%nearwall_blocked_flag==1) nw_count=nw_count+1
  call classify_stage7_boundary_point(grid,layout,x(1),grid%y_center(ny),z(1),1,res); if(res%nearwall_blocked_flag==1) nw_count=nw_count+1
  chg1=max(maxval(abs(fx-fx0)),max(maxval(abs(fy-fy0)),maxval(abs(fz-fz0)))); nw_ok=merge(1,0,nw_count>0 .and. chg1<=1e-14_mytype)
  call classify_stage7_boundary_point(grid,layout,x(1),grid%ymin-0.1_mytype*(grid%ymax-grid%ymin),z(1),1,res); ylo=res%y_outside_blocked_flag
  call classify_stage7_boundary_point(grid,layout,x(1),grid%ymax+0.1_mytype*(grid%ymax-grid%ymin),z(1),1,res); yhi=res%y_outside_blocked_flag; yout=merge(1,0,ylo==1 .and. yhi==1)
  invalid_y=ieee_value(0._mytype,ieee_quiet_nan); call classify_stage7_boundary_point(grid,layout,x(1),invalid_y,z(1),1,res); invcoord=res%invalid_coord_blocked_flag; invstatus=merge(1,0,invcoord==1)
  call init_stage7_collocated_velocity_layout(layout_bad); layout_bad%collocated_flag=0; layout_bad%component_specific_flag=0; layout_bad%u_layout_valid_flag=0; layout_bad%v_layout_valid_flag=0; layout_bad%w_layout_valid_flag=0; layout_bad%layout_valid_flag=0
  call classify_stage7_boundary_point(grid,layout_bad,x(1),y(1),z(1),1,res); invlay=res%invalid_layout_blocked_flag; invlay_status=merge(1,0,invlay==1)
  call classify_stage7_boundary_point(grid,layout,x(1)+lx,y(1),z(1)-lz,1,res); px=res%periodic_x_allowed_flag; pz=res%periodic_z_allowed_flag; pstat=merge(1,0,px==1 .and. pz==1 .and. res%safe_flag==1)
  block_no_write=merge(1,0,chg1<=1e-14_mytype); chg2=chg1
  rhsx=1;rhsy=2;rhsz=3; call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop,noop_mod); call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)
  status=merge(1,0,dep6*dep70*dep71*dep72*dep73*dep74*dep75*dep76==1 .and. safe_ok==1 .and. nw_ok==1 .and. yout==1 .and. invstatus==1 .and. invlay_status==1 .and. pstat==1 .and. block_no_write==1 .and. noop_mod==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0 .and. noop<=1e-14_mytype)
  open(10,file='stage7_outputs/fibre_stage7_boundary_safety_check.dat',status='replace')
  write(10,'(A,1X,I0)')'stage7_boundary_stage6_dependency_status',dep6;write(10,'(A,1X,I0)')'stage7_boundary_stage7_0_dependency_status',dep70;write(10,'(A,1X,I0)')'stage7_boundary_stage7_1_dependency_status',dep71;write(10,'(A,1X,I0)')'stage7_boundary_stage7_2_dependency_status',dep72;write(10,'(A,1X,I0)')'stage7_boundary_stage7_3_dependency_status',dep73;write(10,'(A,1X,I0)')'stage7_boundary_stage7_4_dependency_status',dep74;write(10,'(A,1X,I0)')'stage7_boundary_stage7_5_dependency_status',dep75;write(10,'(A,1X,I0)')'stage7_boundary_stage7_6_dependency_status',dep76
  write(10,'(A,1X,I0)')'stage7_boundary_safe_point_count',safe_count;write(10,'(A,1X,I0)')'stage7_boundary_safe_interior_allowed_flag',safe_ok;write(10,'(A,1X,I0)')'stage7_boundary_safe_interior_status',safe_ok
  write(10,'(A,1X,I0)')'stage7_boundary_nearwall_blocked_count',nw_count;write(10,'(A,1X,I0)')'stage7_boundary_nearwall_blocked_flag',merge(1,0,nw_count>0);write(10,'(A,1X,ES24.16E3)')'stage7_boundary_nearwall_buffer_change_max',chg1;write(10,'(A,1X,I0)')'stage7_boundary_nearwall_status',nw_ok
  write(10,'(A,1X,I0)')'stage7_boundary_y_low_blocked_flag',ylo;write(10,'(A,1X,I0)')'stage7_boundary_y_high_blocked_flag',yhi;write(10,'(A,1X,I0)')'stage7_boundary_y_outside_status',yout
  write(10,'(A,1X,I0)')'stage7_boundary_invalid_coord_blocked_flag',invcoord;write(10,'(A,1X,I0)')'stage7_boundary_invalid_coord_status',invstatus
  write(10,'(A,1X,I0)')'stage7_boundary_invalid_layout_blocked_flag',invlay;write(10,'(A,1X,I0)')'stage7_boundary_invalid_layout_status',invlay_status
  write(10,'(A,1X,I0)')'stage7_boundary_periodic_x_allowed_flag',px;write(10,'(A,1X,I0)')'stage7_boundary_periodic_z_allowed_flag',pz;write(10,'(A,1X,I0)')'stage7_boundary_periodic_allowed_status',pstat
  write(10,'(A,1X,ES24.16E3)')'stage7_boundary_blocked_buffer_change_max',chg2;write(10,'(A,1X,I0)')'stage7_boundary_blocked_no_buffer_write_flag',block_no_write
  write(10,'(A,1X,ES24.16E3)')'stage7_boundary_noop_rhs_change_max',noop;write(10,'(A,1X,I0)')'stage7_boundary_noop_rhs_modified_flag',noop_mod;write(10,'(A,1X,I0)')'stage7_boundary_pressure_poisson_modified_flag',pp;write(10,'(A,1X,I0)')'stage7_boundary_projection_modified_flag',proj;write(10,'(A,1X,I0)')'stage7_boundary_production_dns_called_flag',dns;write(10,'(A,1X,I0)')'stage7_boundary_fluid_update_called_flag',flu;write(10,'(A,1X,I0)')'stage7_boundary_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)')'stage7_boundary_safety_check_status',status; close(10)
contains
  subroutine dependency_status(s6,s70,s71,s72,s73,s74,s75,s76)
    integer,intent(out)::s6,s70,s71,s72,s73,s74,s75,s76
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)*check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype); s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype); s72=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat','stage7_scalar_interpolation_check_status',1._mytype); s73=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat','stage7_scalar_interpolation_robustness_check_status',1._mytype); s74=check_dat('stage7_outputs/fibre_stage7_velocity_interpolation_check.dat','stage7_velocity_interpolation_check_status',1._mytype); s75=check_dat('stage7_outputs/fibre_stage7_force_spreading_check.dat','stage7_force_spreading_check_status',1._mytype); s76=check_dat('stage7_outputs/fibre_stage7_power_consistency_check.dat','stage7_power_consistency_check_status',1._mytype)
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key; real(mytype),intent(in)::val; character(len=512)::line,rest; real(mytype)::vv; integer::u,ios; logical::ex
    check_dat=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',iostat=ios); if(ios/=0)return
    do; read(u,'(A)',iostat=ios)line; if(ios/=0)exit; line=adjustl(line); if(index(line,trim(key))==1) then; rest=adjustl(line(len_trim(key)+1:)); if(len_trim(rest)>0 .and. rest(1:1)=='=') rest=adjustl(rest(2:)); read(rest,*,iostat=ios)vv; if(ios==0 .and. abs(vv-val)<=1e-12_mytype) then; check_dat=1; exit; endif; endif; enddo; close(u)
  end function
end program
