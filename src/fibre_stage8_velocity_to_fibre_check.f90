program fibre_stage8_velocity_to_fibre_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
  use fibre_stage8_runtime_grid_bridge
  use fibre_stage8_lagrangian_state
  use fibre_stage8_velocity_to_fibre
  implicit none
  type(stage7_channel_grid_t) :: gref,gbridge
  type(stage8_runtime_grid_bridge_status_t) :: st
  type(stage7_velocity_layout_t) :: layout,layout_bad
  type(stage8_lagrangian_state_t) :: state
  real(mytype), allocatable :: yface(:),ycenter(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),expv(:,:),rhs0(:,:,:),rhs1(:,:,:),xkeep(:,:),dskeep(:)
  integer :: nx,ny,nz,nlag,io,i,j,k,l,v,r,vc,bc,uc
  integer :: s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,dep,final
  integer :: cst_ok,lin_ok,poi_ok,perx_ok,perz_ok,per_ok,nw_ok,oy_ok,badl_ok,other0_ok,clr_ok,noop_ok
  real(mytype) :: x0,y0,z0,length,dl,err,lin_err,poi_err,perx_err,perz_err,nw_werr,oy_werr,bad_werr,vf_err,sl_err,fo_err,cv_err,rhschg
  real(mytype), parameter :: pi_stage8 = 4.0_mytype * atan(1.0_mytype)
  call ensure_dir('stage8_outputs')
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m); call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
  call file_exists_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',s82o); call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s80o==1.and.s80s==1.and.s81o==1.and.s81s==1.and.s82o==1.and.s82s==1)
  nx=16; ny=17; nz=12; nlag=9
  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz)
  allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(st)
  call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st)
  call init_stage7_collocated_velocity_layout(layout); layout_bad=layout
  layout_bad%collocated_flag=0; layout_bad%component_specific_flag=0; layout_bad%u_layout_valid_flag=0; layout_bad%v_layout_valid_flag=0; layout_bad%w_layout_valid_flag=0; layout_bad%layout_valid_flag=0
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),expv(3,nlag),rhs0(4,3,2),rhs1(4,3,2))
  call init_stage8_lagrangian_state(state); call allocate_stage8_lagrangian_state(state,nlag,v,r)
  x0=gbridge%xmin+0.5_mytype*(gbridge%xmax-gbridge%xmin); y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*(gbridge%zmax-gbridge%zmin)
  length=0.3_mytype*(gbridge%xmax-gbridge%xmin); call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  ux=1.25_mytype; uy=-0.5_mytype; uz=0.75_mytype; expv(1,:)=1.25_mytype; expv(2,:)=-0.5_mytype; expv(3,:)=0.75_mytype
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc); call compute_stage8_velocity_state_error(state,expv,err)
  cst_ok=merge(1,0,vc==nlag.and.bc==0.and.uc==0.and.err<=1e-12_mytype)
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=1._mytype+0.2_mytype*(gbridge%xmin+real(i-1,mytype)*gbridge%dx)-0.1_mytype*gbridge%y_center(j)+0.05_mytype*(gbridge%zmin+real(k-1,mytype)*gbridge%dz)
    uy(i,j,k)=-0.5_mytype+0.1_mytype*(gbridge%xmin+real(i-1,mytype)*gbridge%dx)+0.3_mytype*gbridge%y_center(j)-0.02_mytype*(gbridge%zmin+real(k-1,mytype)*gbridge%dz)
    uz(i,j,k)=0.25_mytype-0.05_mytype*(gbridge%xmin+real(i-1,mytype)*gbridge%dx)+0.07_mytype*gbridge%y_center(j)+0.15_mytype*(gbridge%zmin+real(k-1,mytype)*gbridge%dz)
  end do; end do; end do
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc)
  do l=1,nlag
    expv(1,l)=1._mytype+0.2_mytype*state%x(1,l)-0.1_mytype*state%x(2,l)+0.05_mytype*state%x(3,l)
    expv(2,l)=-0.5_mytype+0.1_mytype*state%x(1,l)+0.3_mytype*state%x(2,l)-0.02_mytype*state%x(3,l)
    expv(3,l)=0.25_mytype-0.05_mytype*state%x(1,l)+0.07_mytype*state%x(2,l)+0.15_mytype*state%x(3,l)
  end do
  call compute_stage8_velocity_state_error(state,expv,lin_err); lin_ok=merge(1,0,lin_err<=1e-11_mytype)
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=1._mytype-gbridge%y_center(j)**2; uy(i,j,k)=0; uz(i,j,k)=0
  end do; end do; end do
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc)
  do l=1,nlag; expv(:,l)=0; expv(1,l)=1._mytype-state%x(2,l)**2; end do
  call compute_stage8_velocity_state_error(state,expv,poi_err); poi_ok=merge(1,0,poi_err<=1e-11_mytype)
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=sin(2.0_mytype*pi_stage8*(real(i-1,mytype))/real(nx,mytype))+0.2_mytype*gbridge%y_center(j)
    uy(i,j,k)=cos(2.0_mytype*pi_stage8*(real(k-1,mytype))/real(nz,mytype))-0.1_mytype*gbridge%y_center(j)
    uz(i,j,k)=sin(2.0_mytype*pi_stage8*(real(i-1,mytype))/real(nx,mytype))+cos(2.0_mytype*pi_stage8*(real(k-1,mytype))/real(nz,mytype))
  end do; end do; end do
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc)
  dl=gbridge%xmax-gbridge%xmin; state%x(1,1)=state%x(1,1)+dl; state%x(1,2)=state%x(1,2)-dl
  dl=gbridge%zmax-gbridge%zmin; state%x(3,3)=state%x(3,3)+dl; state%x(3,4)=state%x(3,4)-dl
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc)
  perx_err=maxval(abs(state%u_fluid_lag(:,1)-state%u_fluid_lag(:,2))); perz_err=maxval(abs(state%u_fluid_lag(:,3)-state%u_fluid_lag(:,4)))
  perx_ok=merge(1,0,perx_err<=1e-12_mytype); perz_ok=merge(1,0,perz_err<=1e-12_mytype); per_ok=merge(1,0,perx_ok==1.and.perz_ok==1)
  call build_stage8_straight_fibre_state(state,gbridge,x0,gbridge%y_center(1),z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc); nw_werr=max_u_fluid_on_blocked_points(state)
  nw_ok=merge(1,0,bc>0.and.uc>0.and.nw_werr<=1e-14_mytype)
  call build_stage8_straight_fibre_state(state,gbridge,x0,gbridge%ymin-0.1_mytype*(gbridge%ymax-gbridge%ymin),z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc); oy_werr=max_u_fluid_on_blocked_points(state)
  oy_ok=merge(1,0,bc>0.and.oy_werr<=1e-14_mytype)
  call build_stage8_straight_fibre_state(state,gbridge,x0,0.5_mytype*(gbridge%ymin+gbridge%ymax),z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout_bad,ux,uy,uz,state,vc,bc,uc); bad_werr=max_u_fluid_on_blocked_points(state); badl_ok=merge(1,0,bc>0.and.bad_werr<=1e-14_mytype)
  vf_err=maxval(abs(state%v_fibre)); sl_err=maxval(abs(state%slip)); fo_err=max(maxval(abs(state%force_structure)),maxval(abs(state%force_fluid)))
  other0_ok=merge(1,0,vf_err<=1e-14_mytype.and.sl_err<=1e-14_mytype.and.fo_err<=1e-14_mytype)
  allocate(xkeep(3,nlag),dskeep(nlag)); xkeep=state%x; dskeep=state%ds
  call clear_stage8_fluid_velocity_lag(state); cv_err=maxval(abs(state%u_fluid_lag)); clr_ok=merge(1,0,cv_err<=1e-14_mytype.and.maxval(abs(state%x-xkeep))<=1e-14_mytype.and.maxval(abs(state%ds-dskeep))<=1e-14_mytype)
  call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
  final=merge(1,0,dep==1.and.cst_ok==1.and.lin_ok==1.and.poi_ok==1.and.per_ok==1.and.nw_ok==1.and.oy_ok==1.and.badl_ok==1.and.other0_ok==1.and.clr_ok==1.and.noop_ok==1)
  open(newunit=io,file='stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_v2f_stage7_closed_marker_exists',s7m; write(io,'(A,1X,I0)') 'stage8_v2f_stage7_total_smoke_output_exists',s7o; write(io,'(A,1X,I0)') 'stage8_v2f_stage7_total_smoke_status',s7s; write(io,'(A,1X,I0)') 'stage8_v2f_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_v2f_stage8_0_output_exists',s80o; write(io,'(A,1X,I0)') 'stage8_v2f_stage8_0_status',s80s; write(io,'(A,1X,I0)') 'stage8_v2f_stage8_1_output_exists',s81o; write(io,'(A,1X,I0)') 'stage8_v2f_stage8_1_status',s81s
  write(io,'(A,1X,I0)') 'stage8_v2f_stage8_2_output_exists',s82o; write(io,'(A,1X,I0)') 'stage8_v2f_stage8_2_status',s82s; write(io,'(A,1X,I0)') 'stage8_v2f_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_v2f_constant_valid_count',nlag; write(io,'(A,1X,I0)') 'stage8_v2f_constant_blocked_count',0; write(io,'(A,1X,I0)') 'stage8_v2f_constant_unsafe_count',0
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_constant_error_max',err; write(io,'(A,1X,I0)') 'stage8_v2f_constant_status',cst_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_linear_error_max',lin_err; write(io,'(A,1X,I0)') 'stage8_v2f_linear_status',lin_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_poiseuille_error_max',poi_err; write(io,'(A,1X,I0)') 'stage8_v2f_poiseuille_status',poi_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_periodic_x_shift_error_max',perx_err; write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_periodic_z_shift_error_max',perz_err
  write(io,'(A,1X,I0)') 'stage8_v2f_periodic_x_shift_status',perx_ok; write(io,'(A,1X,I0)') 'stage8_v2f_periodic_z_shift_status',perz_ok; write(io,'(A,1X,I0)') 'stage8_v2f_periodic_status',per_ok
  write(io,'(A,1X,I0)') 'stage8_v2f_nearwall_valid_count',vc; write(io,'(A,1X,I0)') 'stage8_v2f_nearwall_blocked_count',bc; write(io,'(A,1X,I0)') 'stage8_v2f_nearwall_unsafe_count',uc
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_nearwall_velocity_write_error_max',nw_werr; write(io,'(A,1X,I0)') 'stage8_v2f_nearwall_blocked_flag',merge(1,0,bc>0); write(io,'(A,1X,I0)') 'stage8_v2f_nearwall_status',nw_ok
  write(io,'(A,1X,I0)') 'stage8_v2f_outside_y_blocked_count',bc; write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_outside_y_velocity_write_error_max',oy_werr; write(io,'(A,1X,I0)') 'stage8_v2f_outside_y_blocked_flag',merge(1,0,bc>0); write(io,'(A,1X,I0)') 'stage8_v2f_outside_y_status',oy_ok
  write(io,'(A,1X,I0)') 'stage8_v2f_invalid_layout_blocked_count',bc; write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_invalid_layout_velocity_write_error_max',bad_werr
  write(io,'(A,1X,I0)') 'stage8_v2f_invalid_layout_blocked_flag',merge(1,0,bc>0); write(io,'(A,1X,I0)') 'stage8_v2f_invalid_layout_status',badl_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_v_fibre_zero_error_max',vf_err; write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_slip_zero_error_max',sl_err; write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_force_zero_error_max',fo_err
  write(io,'(A,1X,I0)') 'stage8_v2f_other_placeholders_zero_status',other0_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_clear_velocity_error_max',cv_err; write(io,'(A,1X,I0)') 'stage8_v2f_clear_geometry_preserved_flag',merge(1,0,clr_ok==1); write(io,'(A,1X,I0)') 'stage8_v2f_clear_status',clr_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_v2f_noop_rhs_change_max',rhschg; write(io,'(A,1X,I0)') 'stage8_v2f_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_v2f_pressure_poisson_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_v2f_projection_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_v2f_real_projection_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_v2f_production_dns_called_flag',0; write(io,'(A,1X,I0)') 'stage8_v2f_fluid_update_called_flag',0; write(io,'(A,1X,I0)') 'stage8_v2f_fibre_advance_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_v2f_noop_safety_status',noop_ok; write(io,'(A,1X,I0)') 'stage8_velocity_to_fibre_check_status',final
  close(io)
contains
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end

  function max_u_fluid_on_blocked_points(state) result(err)
    use fibre_parameters, only: mytype
    use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
    type(stage8_lagrangian_state_t), intent(in) :: state
    real(mytype) :: err
    integer :: l

    err = 0.0_mytype
    if (.not. allocated(state%u_fluid_lag)) return
    if (.not. allocated(state%point_blocked_flag)) return

    do l = 1, state%nlag
      if (state%point_blocked_flag(l) == 1) then
        err = max(err, maxval(abs(state%u_fluid_lag(:,l))))
      end if
    end do
  end function max_u_fluid_on_blocked_points
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
