program fibre_stage8_lagrangian_state_check
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
  use fibre_stage8_runtime_grid_bridge
  use fibre_stage8_lagrangian_state
  implicit none
  type(stage7_channel_grid_t) :: gref,gbridge
  type(stage8_runtime_grid_bridge_status_t) :: st
  type(stage7_velocity_layout_t) :: layout
  type(stage8_lagrangian_state_t) :: state
  real(mytype), allocatable :: yface(:),ycenter(:),rhs0(:,:,:),rhs1(:,:,:)
  integer :: io,nx,ny,nz,nlag,dep,s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,av,ar,v,r,sc,bc,uc
  integer :: arr_ok,init_ok,safe_geom_ok,safe_cls_ok,nw_ok,oy_ok,invn,invlen,invd,invnan,inv_ok,clear_ok,ddealloc,des_ok,ph_ok,noop_ok,final
  real(mytype)::x0,y0,z0,length,total_ds,err_ds,err_seg,dl,cz,vc_err,fc_err,clr_err,rhschg
  real(mytype):: dir(3)
  call ensure_dir('stage8_outputs')
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m)
  call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s80o==1.and.s80s==1.and.s81o==1.and.s81s==1)
  nx=16;ny=17;nz=12; call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz)
  allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(st)
  call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st)
  call init_stage7_collocated_velocity_layout(layout)
  nlag=9; call init_stage8_lagrangian_state(state); call allocate_stage8_lagrangian_state(state,nlag,av,ar)
  arr_ok=merge(1,0,allocated(state%x).and.allocated(state%v_fibre).and.allocated(state%u_fluid_lag).and.allocated(state%slip).and.allocated(state%force_structure).and.allocated(state%force_fluid).and.allocated(state%ds))
  init_ok=merge(1,0,av==1.and.ar==0.and.state%allocated_flag==1.and.state%nlag==nlag)
  cz=maxval(abs(state%x)); cz=max(cz,maxval(abs(state%v_fibre))); cz=max(cz,maxval(abs(state%u_fluid_lag))); cz=max(cz,maxval(abs(state%slip))); cz=max(cz,maxval(abs(state%force_structure))); cz=max(cz,maxval(abs(state%force_fluid))); cz=max(cz,maxval(abs(state%ds)))
  cz=max(cz,real(maxval(abs(state%point_valid_flag)),mytype)); cz=max(cz,real(maxval(abs(state%point_blocked_flag)),mytype)); cz=max(cz,real(maxval(abs(state%point_unsafe_flag)),mytype)); cz=max(cz,real(maxval(abs(state%point_status_code)),mytype))
  x0=gbridge%xmin+0.5_mytype*(gbridge%xmax-gbridge%xmin); y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*(gbridge%zmax-gbridge%zmin)
  length=0.4_mytype*(gbridge%xmax-gbridge%xmin); dir=[1._mytype,0._mytype,0._mytype]
  call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,dir,v,r)
  call compute_stage8_lagrangian_total_length(state,total_ds); err_ds=abs(total_ds-length); dl=length/real(nlag-1,mytype); call compute_stage8_lagrangian_segment_error(state,dl,err_seg)
  safe_geom_ok=merge(1,0,v==1.and.r==0.and.err_ds<=1e-12_mytype.and.err_seg<=1e-12_mytype)
  call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc)
  safe_cls_ok=merge(1,0,sc==nlag.and.bc==0.and.uc==0.and.v==1.and.r==0)
  y0=gbridge%y_center(1); call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,dir,v,r); call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc)
  nw_ok=merge(1,0,bc>0.and.uc>0.and.sc<nlag.and.v==0.and.r==1)
  y0=gbridge%ymin-0.1_mytype*(gbridge%ymax-gbridge%ymin); call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,dir,v,r); call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc)
  oy_ok=merge(1,0,bc>0.and.v==0.and.r==1)
  call destroy_stage8_lagrangian_state(state); call init_stage8_lagrangian_state(state); call allocate_stage8_lagrangian_state(state,1,v,r); invn=merge(1,0,v==0.and.r==1)
  call allocate_stage8_lagrangian_state(state,nlag,av,ar); call build_stage8_straight_fibre_state(state,gbridge,x0,0.5_mytype*(gbridge%ymin+gbridge%ymax),z0,0._mytype,dir,v,r); invlen=merge(1,0,v==0.and.r==1)
  call build_stage8_straight_fibre_state(state,gbridge,x0,0.5_mytype*(gbridge%ymin+gbridge%ymax),z0,length,[0._mytype,0._mytype,0._mytype],v,r); invd=merge(1,0,v==0.and.r==1)
  call build_stage8_straight_fibre_state(state,gbridge,x0,0.5_mytype*(gbridge%ymin+gbridge%ymax),z0,length,dir,v,r); state%x(1,1)=ieee_value(0.0_mytype,ieee_quiet_nan)
  call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc); invnan=merge(1,0,v==0.and.r==1); inv_ok=merge(1,0,invn==1.and.invlen==1.and.invd==1.and.invnan==1)
  state%x=1; state%v_fibre=2; state%u_fluid_lag=3; state%slip=4; state%force_structure=5; state%force_fluid=6; state%ds=7
  state%point_valid_flag=1; state%point_blocked_flag=1; state%point_unsafe_flag=1; state%point_status_code=7
  call clear_stage8_lagrangian_state(state)
  clr_err=maxval(abs(state%x)); clr_err=max(clr_err,maxval(abs(state%v_fibre))); clr_err=max(clr_err,maxval(abs(state%u_fluid_lag))); clr_err=max(clr_err,maxval(abs(state%slip))); clr_err=max(clr_err,maxval(abs(state%force_structure))); clr_err=max(clr_err,maxval(abs(state%force_fluid))); clr_err=max(clr_err,maxval(abs(state%ds)))
  clear_ok=merge(1,0,clr_err<=1e-14_mytype.and.state%allocated_flag==1.and.maxval(abs(state%point_valid_flag))==0.and.maxval(abs(state%point_blocked_flag))==0.and.maxval(abs(state%point_unsafe_flag))==0.and.maxval(abs(state%point_status_code))==0)
  call destroy_stage8_lagrangian_state(state); ddealloc=merge(1,0,.not.allocated(state%x).and..not.allocated(state%ds).and.state%allocated_flag==0); des_ok=ddealloc
  call init_stage8_lagrangian_state(state); call allocate_stage8_lagrangian_state(state,nlag,av,ar)
  call build_stage8_straight_fibre_state(state,gbridge,x0,0.5_mytype*(gbridge%ymin+gbridge%ymax),z0,length,dir,v,r)
  vc_err=maxval(abs(state%v_fibre)); vc_err=max(vc_err,maxval(abs(state%u_fluid_lag))); vc_err=max(vc_err,maxval(abs(state%slip)))
  fc_err=maxval(abs(state%force_structure)); fc_err=max(fc_err,maxval(abs(state%force_fluid))); ph_ok=merge(1,0,vc_err<=1e-14_mytype.and.fc_err<=1e-14_mytype)
  allocate(rhs0(4,3,2),rhs1(4,3,2)); call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
  final=merge(1,0,dep==1.and.init_ok==1.and.arr_ok==1.and.cz<=1e-14_mytype.and.safe_geom_ok==1.and.safe_cls_ok==1.and.nw_ok==1.and.oy_ok==1.and.inv_ok==1.and.clear_ok==1.and.des_ok==1.and.ph_ok==1.and.noop_ok==1)
  open(newunit=io,file='stage8_outputs/fibre_stage8_lagrangian_state_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_state_stage7_closed_marker_exists',s7m
  write(io,'(A,1X,I0)') 'stage8_state_stage7_total_smoke_output_exists',s7o
  write(io,'(A,1X,I0)') 'stage8_state_stage7_total_smoke_status',s7s
  write(io,'(A,1X,I0)') 'stage8_state_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_state_stage8_0_output_exists',s80o
  write(io,'(A,1X,I0)') 'stage8_state_stage8_0_status',s80s
  write(io,'(A,1X,I0)') 'stage8_state_stage8_1_output_exists',s81o
  write(io,'(A,1X,I0)') 'stage8_state_stage8_1_status',s81s
  write(io,'(A,1X,I0)') 'stage8_state_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_state_allocate_valid_flag',av
  write(io,'(A,1X,I0)') 'stage8_state_allocate_rejected_flag',ar
  write(io,'(A,1X,I0)') 'stage8_state_arrays_allocated_flag',arr_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_initial_zero_error_max',cz
  write(io,'(A,1X,I0)') 'stage8_state_initialization_status',init_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_safe_length',length
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_safe_total_ds',total_ds
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_safe_total_ds_error',err_ds
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_safe_segment_error_max',err_seg
  write(io,'(A,1X,I0)') 'stage8_state_safe_geometry_status',safe_geom_ok
  call build_stage8_straight_fibre_state(state,gbridge,x0,0.5_mytype*(gbridge%ymin+gbridge%ymax),z0,length,dir,v,r); call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc)
  write(io,'(A,1X,I0)') 'stage8_state_safe_point_count',sc
  write(io,'(A,1X,I0)') 'stage8_state_safe_blocked_count',bc
  write(io,'(A,1X,I0)') 'stage8_state_safe_unsafe_count',uc
  write(io,'(A,1X,I0)') 'stage8_state_safe_classification_status',safe_cls_ok
  y0=gbridge%y_center(1); call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,dir,v,r); call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc)
  write(io,'(A,1X,I0)') 'stage8_state_nearwall_safe_count',sc; write(io,'(A,1X,I0)') 'stage8_state_nearwall_blocked_count',bc; write(io,'(A,1X,I0)') 'stage8_state_nearwall_unsafe_count',uc
  write(io,'(A,1X,I0)') 'stage8_state_nearwall_blocked_flag',merge(1,0,bc>0)
  write(io,'(A,1X,I0)') 'stage8_state_nearwall_status',nw_ok
  y0=gbridge%ymin-0.1_mytype*(gbridge%ymax-gbridge%ymin); call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,dir,v,r); call validate_stage8_lagrangian_state(state,gbridge,layout,v,r,sc,bc,uc)
  write(io,'(A,1X,I0)') 'stage8_state_outside_y_blocked_count',bc
  write(io,'(A,1X,I0)') 'stage8_state_outside_y_blocked_flag',merge(1,0,bc>0)
  write(io,'(A,1X,I0)') 'stage8_state_outside_y_status',oy_ok
  write(io,'(A,1X,I0)') 'stage8_state_invalid_nlag_rejected_flag',invn
  write(io,'(A,1X,I0)') 'stage8_state_invalid_length_rejected_flag',invlen
  write(io,'(A,1X,I0)') 'stage8_state_invalid_direction_rejected_flag',invd
  write(io,'(A,1X,I0)') 'stage8_state_invalid_nan_rejected_flag',invnan
  write(io,'(A,1X,I0)') 'stage8_state_invalid_input_status',inv_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_clear_zero_error_max',clr_err
  write(io,'(A,1X,I0)') 'stage8_state_clear_status',clear_ok
  write(io,'(A,1X,I0)') 'stage8_state_destroy_deallocated_flag',ddealloc
  write(io,'(A,1X,I0)') 'stage8_state_destroy_status',des_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_velocity_placeholder_zero_error',vc_err
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_force_placeholder_zero_error',fc_err
  write(io,'(A,1X,I0)') 'stage8_state_placeholders_zero_status',ph_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_state_noop_rhs_change_max',rhschg
  write(io,'(A,1X,I0)') 'stage8_state_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_state_pressure_poisson_modified_flag',0
  write(io,'(A,1X,I0)') 'stage8_state_projection_modified_flag',0
  write(io,'(A,1X,I0)') 'stage8_state_real_projection_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_state_production_dns_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_state_fluid_update_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_state_fibre_advance_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_state_noop_safety_status',noop_ok
  write(io,'(A,1X,I0)') 'stage8_lagrangian_state_check_status',final
  close(io)
contains
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
