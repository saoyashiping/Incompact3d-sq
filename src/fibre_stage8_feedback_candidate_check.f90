program fibre_stage8_feedback_candidate_check
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
  use fibre_stage8_runtime_grid_bridge
  use fibre_stage8_lagrangian_state
  use fibre_stage8_velocity_to_fibre
  use fibre_stage8_feedback_candidate
  implicit none
  type(stage7_channel_grid_t) :: gref,gbridge
  type(stage8_runtime_grid_bridge_status_t) :: st
  type(stage7_velocity_layout_t) :: layout
  type(stage8_lagrangian_state_t) :: state,sbase,sshift
  real(mytype), allocatable :: yface(:),ycenter(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),rhs0(:,:,:),rhs1(:,:,:),xkeep(:,:),dskeep(:),vfkeep(:,:),ufkeep(:,:)
  integer :: nx,ny,nz,nlag,io,i,j,k,v,r,vc,bc,uc,fv,fr,fvc,fbc,fuc
  integer :: s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,s83o,s83s,dep,final
  integer :: cst_ok,zero_ok,ar_ok,sgn_ok,pow_ok,gal_ok,blk_ok,ib_ok,clr_ok,noop_ok
  integer :: zbrej,nbrej,nanrej
  real(mytype) :: beta,err_slip,err_force,err_zero_slip,err_zero_force,err_ar,dot_s,dot_f,pair_pow,exp_pow,pow_err,gal_slip_err,gal_force_err
  real(mytype) :: blk_slip_err,blk_force_err,clr_err,rhschg,lx,x0,y0,z0,length
  call ensure_dir('stage8_outputs')
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m); call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
  call file_exists_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',s82o); call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82s)
  call file_exists_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',s83o); call get_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat','stage8_velocity_to_fibre_check_status',s83s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s80o==1.and.s80s==1.and.s81o==1.and.s81s==1.and.s82o==1.and.s82s==1.and.s83o==1.and.s83s==1)
  nx=16;ny=17;nz=12;nlag=9
  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz)
  allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(st)
  call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st)
  call init_stage7_collocated_velocity_layout(layout)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),rhs0(4,3,2),rhs1(4,3,2))
  ux=1.25_mytype; uy=-0.5_mytype; uz=0.75_mytype
  call init_stage8_lagrangian_state(state); call allocate_stage8_lagrangian_state(state,nlag,v,r)
  x0=gbridge%xmin+0.5_mytype*(gbridge%xmax-gbridge%xmin); y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*(gbridge%zmax-gbridge%zmin)
  length=0.2_mytype*(gbridge%xmax-gbridge%xmin)
  call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc)
  beta=2.5_mytype
  call assemble_stage8_slip_velocity(state,vc,bc,uc)
  err_slip=maxval(abs(state%slip-state%u_fluid_lag))
  call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc)
  err_force=maxval(abs(state%force_structure-beta*state%slip)); err_force=max(err_force,maxval(abs(state%force_fluid+beta*state%slip)))
  cst_ok=merge(1,0,vc==nlag.and.bc==0.and.uc==0.and.err_slip<=1e-12_mytype.and.err_force<=1e-12_mytype)
  state%v_fibre=state%u_fluid_lag
  call assemble_stage8_slip_velocity(state,vc,bc,uc); call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc)
  err_zero_slip=maxval(abs(state%slip)); err_zero_force=max(maxval(abs(state%force_structure)),maxval(abs(state%force_fluid)))
  zero_ok=merge(1,0,err_zero_slip<=1e-14_mytype.and.err_zero_force<=1e-14_mytype)
  state%v_fibre=0; call assemble_stage8_slip_velocity(state,vc,bc,uc); call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc)
  call compute_stage8_feedback_action_reaction_error(state,err_ar); call compute_stage8_feedback_sign_metrics(state,dot_s,dot_f,pair_pow)
  ar_ok=merge(1,0,err_ar<=1e-12_mytype); sgn_ok=merge(1,0,dot_s>0._mytype.and.dot_f<0._mytype)
  exp_pow=-beta*sum((state%slip(1,:)**2+state%slip(2,:)**2+state%slip(3,:)**2)*state%ds); pow_err=abs(pair_pow-exp_pow)
  pow_ok=merge(1,0,pair_pow<=0._mytype.and.pow_err<=1e-12_mytype)
  call init_stage8_lagrangian_state(sbase); call allocate_stage8_lagrangian_state(sbase,nlag,v,r)
  call build_stage8_straight_fibre_state(sbase,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,sbase,vc,bc,uc); sbase%v_fibre=0.2_mytype
  call assemble_stage8_slip_velocity(sbase,vc,bc,uc); call assemble_stage8_linear_feedback_candidate(sbase,beta,fv,fr,fvc,fbc,fuc)
  call init_stage8_lagrangian_state(sshift); call allocate_stage8_lagrangian_state(sshift,nlag,v,r)
  sshift=sbase; sshift%u_fluid_lag=sbase%u_fluid_lag+reshape([0.8_mytype,-0.4_mytype,0.25_mytype],[3,1]); sshift%v_fibre=sbase%v_fibre+reshape([0.8_mytype,-0.4_mytype,0.25_mytype],[3,1])
  call assemble_stage8_slip_velocity(sshift,vc,bc,uc); call assemble_stage8_linear_feedback_candidate(sshift,beta,fv,fr,fvc,fbc,fuc)
  gal_slip_err=maxval(abs(sshift%slip-sbase%slip)); gal_force_err=max(maxval(abs(sshift%force_structure-sbase%force_structure)),maxval(abs(sshift%force_fluid-sbase%force_fluid)))
  gal_ok=merge(1,0,gal_slip_err<=1e-12_mytype.and.gal_force_err<=1e-12_mytype)
  call build_stage8_straight_fibre_state(state,gbridge,x0,gbridge%y_center(1),z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc); state%u_fluid_lag=3._mytype; state%v_fibre=1._mytype
  call assemble_stage8_slip_velocity(state,vc,bc,uc); call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc)
  blk_slip_err=max_on_blocked(state%slip,state); blk_force_err=max(max_on_blocked(state%force_structure,state),max_on_blocked(state%force_fluid,state))
  blk_ok=merge(1,0,bc>0.and.blk_slip_err<=1e-14_mytype.and.blk_force_err<=1e-14_mytype)
  zbrej=0; nbrej=0; nanrej=0
  call assemble_stage8_linear_feedback_candidate(state,0._mytype,fv,fr,fvc,fbc,fuc); zbrej=merge(1,0,fv==0.and.fr==1)
  call assemble_stage8_linear_feedback_candidate(state,-1._mytype,fv,fr,fvc,fbc,fuc); nbrej=merge(1,0,fv==0.and.fr==1)
  call assemble_stage8_linear_feedback_candidate(state,ieee_value(0._mytype,ieee_quiet_nan),fv,fr,fvc,fbc,fuc); nanrej=merge(1,0,fv==0.and.fr==1)
  ib_ok=merge(1,0,zbrej==1.and.nbrej==1.and.nanrej==1)
  call build_stage8_straight_fibre_state(state,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call interpolate_stage8_velocity_to_state(gbridge,layout,ux,uy,uz,state,vc,bc,uc); state%v_fibre=0.2_mytype
  call assemble_stage8_slip_velocity(state,vc,bc,uc); call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc)
  allocate(xkeep(3,nlag),dskeep(nlag),vfkeep(3,nlag),ufkeep(3,nlag)); xkeep=state%x; dskeep=state%ds; vfkeep=state%v_fibre; ufkeep=state%u_fluid_lag
  call clear_stage8_slip_and_feedback(state)
  clr_err=max(maxval(abs(state%slip)),max(maxval(abs(state%force_structure)),maxval(abs(state%force_fluid))))
  clr_ok=merge(1,0,clr_err<=1e-14_mytype.and.maxval(abs(state%x-xkeep))<=1e-14_mytype.and.maxval(abs(state%ds-dskeep))<=1e-14_mytype.and.maxval(abs(state%v_fibre-vfkeep))<=1e-14_mytype.and.maxval(abs(state%u_fluid_lag-ufkeep))<=1e-14_mytype)
  call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
  final=merge(1,0,dep==1.and.cst_ok==1.and.zero_ok==1.and.ar_ok==1.and.sgn_ok==1.and.pow_ok==1.and.gal_ok==1.and.blk_ok==1.and.ib_ok==1.and.clr_ok==1.and.noop_ok==1)
  open(newunit=io,file='stage8_outputs/fibre_stage8_feedback_candidate_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_feedback_stage7_closed_marker_exists',s7m
  write(io,'(A,1X,I0)') 'stage8_feedback_stage7_total_smoke_output_exists',s7o
  write(io,'(A,1X,I0)') 'stage8_feedback_stage7_total_smoke_status',s7s
  write(io,'(A,1X,I0)') 'stage8_feedback_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_0_output_exists',s80o
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_0_status',s80s
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_1_output_exists',s81o
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_1_status',s81s
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_2_output_exists',s82o
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_2_status',s82s
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_3_output_exists',s83o
  write(io,'(A,1X,I0)') 'stage8_feedback_stage8_3_status',s83s
  write(io,'(A,1X,I0)') 'stage8_feedback_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_feedback_constant_valid_count',nlag
  write(io,'(A,1X,I0)') 'stage8_feedback_constant_blocked_count',0
  write(io,'(A,1X,I0)') 'stage8_feedback_constant_unsafe_count',0
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_constant_slip_error_max',err_slip
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_constant_force_error_max',err_force
  write(io,'(A,1X,I0)') 'stage8_feedback_constant_status',cst_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_zero_slip_error_max',err_zero_slip
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_zero_force_error_max',err_zero_force
  write(io,'(A,1X,I0)') 'stage8_feedback_zero_slip_status',zero_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_action_reaction_error_max',err_ar
  write(io,'(A,1X,I0)') 'stage8_feedback_action_reaction_status',ar_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_structure_force_slip_dot',dot_s
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_fluid_force_slip_dot',dot_f
  write(io,'(A,1X,I0)') 'stage8_feedback_structure_force_slip_dot_positive_flag',merge(1,0,dot_s>0._mytype)
  write(io,'(A,1X,I0)') 'stage8_feedback_fluid_force_slip_dot_negative_flag',merge(1,0,dot_f<0._mytype)
  write(io,'(A,1X,I0)') 'stage8_feedback_sign_status',sgn_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_total_pair_power',pair_pow
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_expected_pair_power',exp_pow
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_pair_power_error',pow_err
  write(io,'(A,1X,I0)') 'stage8_feedback_pair_power_dissipative_flag',merge(1,0,pair_pow<=0._mytype)
  write(io,'(A,1X,I0)') 'stage8_feedback_pair_power_status',pow_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_galilean_slip_error_max',gal_slip_err
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_galilean_force_error_max',gal_force_err
  write(io,'(A,1X,I0)') 'stage8_feedback_galilean_status',gal_ok
  write(io,'(A,1X,I0)') 'stage8_feedback_blocked_count',bc
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_blocked_slip_write_error_max',blk_slip_err
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_blocked_force_write_error_max',blk_force_err
  write(io,'(A,1X,I0)') 'stage8_feedback_blocked_status',blk_ok
  write(io,'(A,1X,I0)') 'stage8_feedback_zero_beta_rejected_flag',zbrej
  write(io,'(A,1X,I0)') 'stage8_feedback_negative_beta_rejected_flag',nbrej
  write(io,'(A,1X,I0)') 'stage8_feedback_nan_beta_rejected_flag',nanrej
  write(io,'(A,1X,I0)') 'stage8_feedback_invalid_beta_status',ib_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_clear_slip_force_error_max',clr_err
  write(io,'(A,1X,I0)') 'stage8_feedback_clear_geometry_preserved_flag',merge(1,0,maxval(abs(state%x-xkeep))<=1e-14_mytype.and.maxval(abs(state%ds-dskeep))<=1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_feedback_clear_velocity_preserved_flag',merge(1,0,maxval(abs(state%v_fibre-vfkeep))<=1e-14_mytype.and.maxval(abs(state%u_fluid_lag-ufkeep))<=1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_feedback_clear_status',clr_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_feedback_noop_rhs_change_max',rhschg
  write(io,'(A,1X,I0)') 'stage8_feedback_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_feedback_pressure_poisson_modified_flag',0
  write(io,'(A,1X,I0)') 'stage8_feedback_projection_modified_flag',0
  write(io,'(A,1X,I0)') 'stage8_feedback_real_projection_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_feedback_production_dns_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_feedback_fluid_update_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_feedback_fibre_advance_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_feedback_noop_safety_status',noop_ok
  write(io,'(A,1X,I0)') 'stage8_feedback_candidate_check_status',final
  close(io)
contains
  function max_on_blocked(a,state) result(err)
    real(mytype), intent(in) :: a(:,:)
    type(stage8_lagrangian_state_t), intent(in) :: state
    real(mytype) :: err
    integer :: l
    err=0
    do l=1,state%nlag
      if(state%point_blocked_flag(l)==1) err=max(err,maxval(abs(a(:,l))))
    end do
  end function
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
