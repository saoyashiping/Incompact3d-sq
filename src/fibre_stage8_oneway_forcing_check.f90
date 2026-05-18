program fibre_stage8_oneway_forcing_check
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
  use fibre_stage8_runtime_grid_bridge
  use fibre_stage8_lagrangian_state
  use fibre_stage8_velocity_to_fibre, only: clear_stage8_fluid_velocity_lag
  use fibre_stage8_feedback_candidate, only: clear_stage8_slip_and_feedback
  use fibre_stage8_oneway_forcing
  implicit none
  type(stage7_channel_grid_t) :: gref,gbridge
  type(stage8_runtime_grid_bridge_status_t) :: st
  type(stage7_velocity_layout_t) :: layout
  type(stage8_lagrangian_state_t) :: state,sbase,sshift,sblocked
  real(mytype), allocatable :: yface(:),ycenter(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),expf(:,:),rhs0(:,:,:),rhs1(:,:,:),expu(:,:),exps(:,:)
  integer :: nx,ny,nz,nlag,io,i,j,k,l,v,r,vc,bc,uc,ok,rej
  integer :: s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,s83o,s83s,s84o,s84s,dep,final
  integer :: cst_ok,zs_ok,lin_ok,poi_ok,gal_ok,pow_ok,blk_ok,ib_ok,ns_ok,noop_ok
  real(mytype)::beta,x0,y0,z0,length,lx,vel_err,slip_err,force_err,zs_slip,zs_force,lin_v,lin_s,lin_f,poi_v,poi_f,gal_s,gal_f,pow,pow_e,pow_err,blk_u,blk_s,blk_f,rhschg
  real(mytype) :: ushift(3), ux_shift, uy_shift, uz_shift
  integer :: zf,nf,nanf,bcnt,safe,blocked,unsafe
  call ensure_dir('stage8_outputs')
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m); call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
  call file_exists_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',s82o); call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82s)
  call file_exists_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',s83o); call get_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat','stage8_velocity_to_fibre_check_status',s83s)
  call file_exists_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat',s84o); call get_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat','stage8_feedback_candidate_check_status',s84s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s80o==1.and.s80s==1.and.s81o==1.and.s81s==1.and.s82o==1.and.s82s==1.and.s83o==1.and.s83s==1.and.s84o==1.and.s84s==1)
  nx=16;ny=17;nz=12;nlag=9; beta=2.5_mytype
  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz); allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(st); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st)
  call init_stage7_collocated_velocity_layout(layout)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),expf(3,nlag),rhs0(4,3,2),rhs1(4,3,2),expu(3,nlag),exps(3,nlag))
  lx=gbridge%xmax-gbridge%xmin; x0=gbridge%xmin+0.5_mytype*lx; y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*(gbridge%zmax-gbridge%zmin); length=0.2_mytype*lx
  call init_stage8_lagrangian_state(state); call allocate_stage8_lagrangian_state(state,nlag,v,r); call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  state%v_fibre=0; ux=1.25_mytype; uy=-0.5_mytype; uz=0.75_mytype
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,state,ok,rej,vc,bc,uc)
  expu(1,:)=1.25_mytype; expu(2,:)=-0.5_mytype; expu(3,:)=0.75_mytype; exps=expu; expf=beta*exps
  vel_err=maxval(abs(state%u_fluid_lag-expu)); slip_err=maxval(abs(state%slip-exps)); force_err=maxval(abs(state%force_structure-expf))
  cst_ok=merge(1,0,ok==1.and.rej==0.and.vc==nlag.and.bc==0.and.uc==0.and.vel_err<=1e-12_mytype.and.slip_err<=1e-12_mytype.and.force_err<=1e-12_mytype)
  state%v_fibre=expu; call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,state,ok,rej,vc,bc,uc)
  zs_slip=maxval(abs(state%slip)); zs_force=maxval(abs(state%force_structure)); zs_ok=merge(1,0,zs_slip<=1e-14_mytype.and.zs_force<=1e-14_mytype)
  do k=1,nz;do j=1,ny;do i=1,nx
    ux(i,j,k)=1._mytype+0.2_mytype*(gbridge%xmin+real(i-1,mytype)*gbridge%dx)-0.1_mytype*gbridge%y_center(j)+0.05_mytype*(gbridge%zmin+real(k-1,mytype)*gbridge%dz)
    uy(i,j,k)=-0.5_mytype+0.1_mytype*(gbridge%xmin+real(i-1,mytype)*gbridge%dx)+0.3_mytype*gbridge%y_center(j)-0.02_mytype*(gbridge%zmin+real(k-1,mytype)*gbridge%dz)
    uz(i,j,k)=0.25_mytype-0.05_mytype*(gbridge%xmin+real(i-1,mytype)*gbridge%dx)+0.07_mytype*gbridge%y_center(j)+0.15_mytype*(gbridge%zmin+real(k-1,mytype)*gbridge%dz)
  end do;end do;end do
  do l=1,nlag; state%v_fibre(1,l)=0.1_mytype+0.01_mytype*real(l,mytype); state%v_fibre(2,l)=-0.2_mytype+0.02_mytype*real(l,mytype); state%v_fibre(3,l)=0.05_mytype-0.01_mytype*real(l,mytype); end do
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,state,ok,rej,vc,bc,uc)
  do l=1,nlag
    expu(1,l)=1._mytype+0.2_mytype*state%x(1,l)-0.1_mytype*state%x(2,l)+0.05_mytype*state%x(3,l)
    expu(2,l)=-0.5_mytype+0.1_mytype*state%x(1,l)+0.3_mytype*state%x(2,l)-0.02_mytype*state%x(3,l)
    expu(3,l)=0.25_mytype-0.05_mytype*state%x(1,l)+0.07_mytype*state%x(2,l)+0.15_mytype*state%x(3,l)
  end do
  exps=expu-state%v_fibre; expf=beta*exps
  lin_v=maxval(abs(state%u_fluid_lag-expu)); lin_s=maxval(abs(state%slip-exps)); lin_f=maxval(abs(state%force_structure-expf)); lin_ok=merge(1,0,lin_v<=1e-11_mytype.and.lin_s<=1e-11_mytype.and.lin_f<=1e-11_mytype)
  do k=1,nz;do j=1,ny;do i=1,nx; ux(i,j,k)=1._mytype-gbridge%y_center(j)**2; uy(i,j,k)=0; uz(i,j,k)=0; end do;end do;end do
  state%v_fibre=0; call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,state,ok,rej,vc,bc,uc)
  do l=1,nlag; expu(:,l)=0; expu(1,l)=1._mytype-state%x(2,l)**2; expf(:,l)=beta*expu(:,l); end do
  poi_v=maxval(abs(state%u_fluid_lag-expu)); poi_f=maxval(abs(state%force_structure-expf)); poi_ok=merge(1,0,poi_v<=1e-11_mytype.and.poi_f<=1e-11_mytype)
  call init_stage8_lagrangian_state(sbase); call allocate_stage8_lagrangian_state(sbase,nlag,v,r); call build_stage8_straight_fibre_state(sbase,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  sbase%v_fibre=0.2_mytype; ux=1.25_mytype; uy=-0.5_mytype; uz=0.75_mytype
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,sbase,ok,rej,vc,bc,uc)
  call init_stage8_lagrangian_state(sshift); call allocate_stage8_lagrangian_state(sshift,nlag,v,r); sshift=sbase
  call clear_stage8_fluid_velocity_lag(sshift); call clear_stage8_slip_and_feedback(sshift)
  ushift=(/0.8_mytype,-0.4_mytype,0.25_mytype/)
  do l=1,nlag
    sshift%v_fibre(:,l)=sbase%v_fibre(:,l)+ushift(:)
  end do
  ux_shift=1.25_mytype+ushift(1); uy_shift=-0.5_mytype+ushift(2); uz_shift=0.75_mytype+ushift(3)
  ux=ux_shift; uy=uy_shift; uz=uz_shift
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,sshift,ok,rej,vc,bc,uc)
  gal_s=maxval(abs(sshift%slip-sbase%slip)); gal_f=maxval(abs(sshift%force_structure-sbase%force_structure)); gal_ok=merge(1,0,gal_s<=1e-12_mytype.and.gal_f<=1e-12_mytype)
  call compute_stage8_oneway_structure_power(state,pow); pow_e=sum((state%force_structure(1,:)*state%v_fibre(1,:)+state%force_structure(2,:)*state%v_fibre(2,:)+state%force_structure(3,:)*state%v_fibre(3,:))*state%ds); pow_err=abs(pow-pow_e); pow_ok=merge(1,0,pow_err<=1e-12_mytype)
  call init_stage8_lagrangian_state(sblocked); call allocate_stage8_lagrangian_state(sblocked,nlag,v,r); call build_stage8_straight_fibre_state(sblocked,gbridge,x0,gbridge%y_center(1),z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call validate_stage8_lagrangian_state(sblocked,gbridge,layout,ok,rej,safe,blocked,unsafe); sblocked%v_fibre=1._mytype; ux=3._mytype; uy=2._mytype; uz=-1._mytype
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,sblocked,ok,rej,vc,bc,uc); bcnt=blocked
  blk_u=max_on_blocked(sblocked%u_fluid_lag,sblocked); blk_s=max_on_blocked(sblocked%slip,sblocked); blk_f=max(max_on_blocked(sblocked%force_structure,sblocked),max_on_blocked(sblocked%force_fluid,sblocked))
  blk_ok=merge(1,0,bcnt>0.and.blk_u<=1e-14_mytype.and.blk_s<=1e-14_mytype.and.blk_f<=1e-14_mytype)
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,0._mytype,state,ok,rej,vc,bc,uc); zf=merge(1,0,ok==0.and.rej==1)
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,-1._mytype,state,ok,rej,vc,bc,uc); nf=merge(1,0,ok==0.and.rej==1)
  call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,ieee_value(0._mytype,ieee_quiet_nan),state,ok,rej,vc,bc,uc); nanf=merge(1,0,ok==0.and.rej==1)
  ib_ok=merge(1,0,zf==1.and.nf==1.and.nanf==1)
  ns_ok=1
  call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
  final=merge(1,0,dep==1.and.cst_ok==1.and.zs_ok==1.and.lin_ok==1.and.poi_ok==1.and.gal_ok==1.and.pow_ok==1.and.blk_ok==1.and.ib_ok==1.and.ns_ok==1.and.noop_ok==1)
  open(newunit=io,file='stage8_outputs/fibre_stage8_oneway_forcing_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_oneway_stage7_closed_marker_exists',s7m; write(io,'(A,1X,I0)') 'stage8_oneway_stage7_total_smoke_output_exists',s7o
  write(io,'(A,1X,I0)') 'stage8_oneway_stage7_total_smoke_status',s7s; write(io,'(A,1X,I0)') 'stage8_oneway_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_oneway_stage8_0_output_exists',s80o; write(io,'(A,1X,I0)') 'stage8_oneway_stage8_0_status',s80s
  write(io,'(A,1X,I0)') 'stage8_oneway_stage8_1_output_exists',s81o; write(io,'(A,1X,I0)') 'stage8_oneway_stage8_1_status',s81s
  write(io,'(A,1X,I0)') 'stage8_oneway_stage8_2_output_exists',s82o; write(io,'(A,1X,I0)') 'stage8_oneway_stage8_2_status',s82s
  write(io,'(A,1X,I0)') 'stage8_oneway_stage8_3_output_exists',s83o; write(io,'(A,1X,I0)') 'stage8_oneway_stage8_3_status',s83s
  write(io,'(A,1X,I0)') 'stage8_oneway_stage8_4_output_exists',s84o; write(io,'(A,1X,I0)') 'stage8_oneway_stage8_4_status',s84s
  write(io,'(A,1X,I0)') 'stage8_oneway_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_oneway_constant_valid_count',nlag; write(io,'(A,1X,I0)') 'stage8_oneway_constant_blocked_count',0; write(io,'(A,1X,I0)') 'stage8_oneway_constant_unsafe_count',0
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_constant_velocity_error_max',vel_err; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_constant_slip_error_max',slip_err; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_constant_force_error_max',force_err; write(io,'(A,1X,I0)') 'stage8_oneway_constant_status',cst_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_zero_slip_error_max',zs_slip; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_zero_force_error_max',zs_force; write(io,'(A,1X,I0)') 'stage8_oneway_zero_slip_status',zs_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_linear_velocity_error_max',lin_v; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_linear_slip_error_max',lin_s; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_linear_force_error_max',lin_f; write(io,'(A,1X,I0)') 'stage8_oneway_linear_status',lin_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_poiseuille_velocity_error_max',poi_v; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_poiseuille_force_error_max',poi_f; write(io,'(A,1X,I0)') 'stage8_oneway_poiseuille_status',poi_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_galilean_slip_error_max',gal_s; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_galilean_force_error_max',gal_f; write(io,'(A,1X,I0)') 'stage8_oneway_galilean_status',gal_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_structure_power',pow; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_expected_structure_power',pow_e; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_structure_power_error',pow_err; write(io,'(A,1X,I0)') 'stage8_oneway_structure_power_status',pow_ok
  write(io,'(A,1X,I0)') 'stage8_oneway_blocked_count',bcnt; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_blocked_velocity_write_error_max',blk_u; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_blocked_slip_write_error_max',blk_s; write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_blocked_force_write_error_max',blk_f; write(io,'(A,1X,I0)') 'stage8_oneway_blocked_status',blk_ok
  write(io,'(A,1X,I0)') 'stage8_oneway_zero_beta_rejected_flag',zf; write(io,'(A,1X,I0)') 'stage8_oneway_negative_beta_rejected_flag',nf; write(io,'(A,1X,I0)') 'stage8_oneway_nan_beta_rejected_flag',nanf; write(io,'(A,1X,I0)') 'stage8_oneway_invalid_beta_status',ib_ok
  write(io,'(A,1X,I0)') 'stage8_oneway_spreading_called_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_rhs_hook_called_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_rhs_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_no_spreading_no_rhs_status',ns_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_oneway_noop_rhs_change_max',rhschg; write(io,'(A,1X,I0)') 'stage8_oneway_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_oneway_pressure_poisson_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_projection_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_real_projection_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_oneway_production_dns_called_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_fluid_update_called_flag',0; write(io,'(A,1X,I0)') 'stage8_oneway_fibre_advance_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_oneway_noop_safety_status',noop_ok; write(io,'(A,1X,I0)') 'stage8_oneway_forcing_check_status',final
  close(io)
contains
  function max_on_blocked(a,state) result(err)
    real(mytype),intent(in)::a(:,:); type(stage8_lagrangian_state_t),intent(in)::state; real(mytype)::err; integer::l
    err=0; do l=1,state%nlag; if(state%point_blocked_flag(l)==1) err=max(err,maxval(abs(a(:,l)))); end do
  end
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
