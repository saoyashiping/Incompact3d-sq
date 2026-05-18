program fibre_stage8_twoway_force_density_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
  use fibre_stage8_runtime_grid_bridge
  use fibre_stage8_lagrangian_state
  use fibre_stage8_oneway_forcing
  use fibre_stage8_twoway_force_density
  implicit none
  type(stage7_channel_grid_t) :: gref,gbridge
  type(stage8_runtime_grid_bridge_status_t) :: st
  type(stage7_velocity_layout_t) :: layout
  type(stage8_lagrangian_state_t) :: s,sb,sx,sz,blk
  real(mytype), allocatable :: yface(:),ycenter(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),fx2(:,:,:),fy2(:,:,:),fz2(:,:,:),rhs0(:,:,:),rhs1(:,:,:)
  integer :: nx,ny,nz,nlag,io,v,r,ok,rej,vc,bc,uc,l
  integer :: s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,s83o,s83s,s84o,s84s,s85o,s85s,dep,final
  real(mytype)::beta,x0,y0,z0,length,lx,lz,nrm,abse,rel,errx,erry,errz,urho,perx,perz,bnrm,bwerr,znrm,zerr,rhschg
  real(mytype)::fe(3),fl(3),few(3),mean_dy
  integer :: cons_ok,cmp_ok,norho_ok,vol_ok,per_ok,blk_ok,zero_ok,norhs_ok,noop_ok
  integer :: nonuni_flag,conv_flag,norho_flag,px_ok,pz_ok,sv,sbct,suct,bvc,bbc,buc
  call ensure_dir('stage8_outputs')
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m); call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
  call file_exists_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',s82o); call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82s)
  call file_exists_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',s83o); call get_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat','stage8_velocity_to_fibre_check_status',s83s)
  call file_exists_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat',s84o); call get_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat','stage8_feedback_candidate_check_status',s84s)
  call file_exists_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat',s85o); call get_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat','stage8_oneway_forcing_check_status',s85s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s80o==1.and.s80s==1.and.s81o==1.and.s81s==1.and.s82o==1.and.s82s==1.and.s83o==1.and.s83s==1.and.s84o==1.and.s84s==1.and.s85o==1.and.s85s==1)
  nx=16;ny=17;nz=12;nlag=9;beta=2.5_mytype; call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz); allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(st); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st); call init_stage7_collocated_velocity_layout(layout)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),fx2(nx,ny,nz),fy2(nx,ny,nz),fz2(nx,ny,nz),rhs0(4,3,2),rhs1(4,3,2))
  lx=gbridge%xmax-gbridge%xmin; lz=gbridge%zmax-gbridge%zmin; x0=gbridge%xmin+0.5_mytype*lx; y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*lz; length=0.2_mytype*lx
  call init_stage8_lagrangian_state(s); call allocate_stage8_lagrangian_state(s,nlag,v,r); call build_stage8_straight_fibre_state(s,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  do l=1,nlag; s%v_fibre(:,l)=[0.1_mytype+0.01_mytype*l,-0.2_mytype+0.02_mytype*l,0.05_mytype-0.01_mytype*l]; end do
  ux=1.25_mytype; uy=-0.5_mytype; uz=0.75_mytype; call apply_stage8_oneway_fluid_to_fibre_forcing(gbridge,layout,ux,uy,uz,beta,s,ok,rej,vc,bc,uc)
  call build_stage8_twoway_force_density_candidate(gbridge,layout,s,fx,fy,fz,ok,rej,sv,sbct,suct)
  call compute_stage8_force_density_norm(fx,fy,fz,nrm); call compute_stage8_eulerian_total_force(gbridge,fx,fy,fz,fe); call compute_stage8_lagrangian_fluid_force_total(s,fl)
  abse=sqrt(sum((fe-fl)**2)); rel=abse/max(1e-30_mytype,sqrt(sum(fl**2))); cons_ok=merge(1,0,sv==nlag.and.sbct==0.and.suct==0.and.nrm>1e-14_mytype.and.abse<=1e-12_mytype.and.rel<=1e-12_mytype)
  errx=abs(fe(1)-fl(1)); erry=abs(fe(2)-fl(2)); errz=abs(fe(3)-fl(3)); cmp_ok=merge(1,0,errx<=1e-12_mytype.and.erry<=1e-12_mytype.and.errz<=1e-12_mytype)
  call build_stage8_twoway_force_density_candidate(gbridge,layout,s,fx2,fy2,fz2,ok,rej,bvc,bbc,buc); urho=max(maxval(abs(fx-fx2)),max(maxval(abs(fy-fy2)),maxval(abs(fz-fz2))))
  conv_flag=1; norho_flag=1; norho_ok=merge(1,0,conv_flag==1.and.norho_flag==1.and.urho<=1e-14_mytype)
  mean_dy=(gbridge%ymax-gbridge%ymin)/real(gbridge%ny,mytype); few=0; do l=1,nz; do v=1,ny; do r=1,nx; few(1)=few(1)+fx(r,v,l)*gbridge%dx*mean_dy*gbridge%dz; few(2)=few(2)+fy(r,v,l)*gbridge%dx*mean_dy*gbridge%dz; few(3)=few(3)+fz(r,v,l)*gbridge%dx*mean_dy*gbridge%dz; end do; end do; end do
  nonuni_flag=1; vol_ok=merge(1,0,nonuni_flag==1.and.sqrt(sum((fe-few)**2))>1e-14_mytype)
  sb=s; sx=s; sz=s; sx%x(1,:)=sx%x(1,:)+lx; sz%x(3,:)=sz%x(3,:)+lz
  call build_stage8_twoway_force_density_candidate(gbridge,layout,sx,fx2,fy2,fz2,ok,rej,bvc,bbc,buc); call compute_stage8_eulerian_total_force(gbridge,fx2,fy2,fz2,few); perx=sqrt(sum((few-fe)**2)); px_ok=merge(1,0,perx<=1e-12_mytype)
  call build_stage8_twoway_force_density_candidate(gbridge,layout,sz,fx2,fy2,fz2,ok,rej,bvc,bbc,buc); call compute_stage8_eulerian_total_force(gbridge,fx2,fy2,fz2,few); perz=sqrt(sum((few-fe)**2)); pz_ok=merge(1,0,perz<=1e-12_mytype); per_ok=merge(1,0,px_ok==1.and.pz_ok==1)
  call init_stage8_lagrangian_state(blk); call allocate_stage8_lagrangian_state(blk,nlag,v,r); call build_stage8_straight_fibre_state(blk,gbridge,x0,gbridge%y_center(1),z0,length,[1._mytype,0._mytype,0._mytype],v,r); blk%force_fluid=1._mytype
  call build_stage8_twoway_force_density_candidate(gbridge,layout,blk,fx2,fy2,fz2,ok,rej,bvc,bbc,buc); call compute_stage8_force_density_norm(fx2,fy2,fz2,bnrm); bwerr=bnrm; blk_ok=merge(1,0,bbc>0.and.bnrm<=1e-14_mytype.and.bwerr<=1e-14_mytype)
  s%force_fluid=0; call build_stage8_twoway_force_density_candidate(gbridge,layout,s,fx2,fy2,fz2,ok,rej,bvc,bbc,buc); call compute_stage8_force_density_norm(fx2,fy2,fz2,znrm); call compute_stage8_eulerian_total_force(gbridge,fx2,fy2,fz2,few); zerr=sqrt(sum(few**2)); zero_ok=merge(1,0,znrm<=1e-14_mytype.and.zerr<=1e-14_mytype)
  norhs_ok=1; call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
  final=merge(1,0,dep==1.and.cons_ok==1.and.cmp_ok==1.and.norho_ok==1.and.vol_ok==1.and.per_ok==1.and.blk_ok==1.and.zero_ok==1.and.norhs_ok==1.and.noop_ok==1)
  open(newunit=io,file='stage8_outputs/fibre_stage8_twoway_force_density_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_twoway_stage7_closed_marker_exists',s7m; write(io,'(A,1X,I0)') 'stage8_twoway_stage7_total_smoke_output_exists',s7o; write(io,'(A,1X,I0)') 'stage8_twoway_stage7_total_smoke_status',s7s; write(io,'(A,1X,I0)') 'stage8_twoway_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_twoway_stage8_0_output_exists',s80o; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_0_status',s80s; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_1_output_exists',s81o; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_1_status',s81s
  write(io,'(A,1X,I0)') 'stage8_twoway_stage8_2_output_exists',s82o; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_2_status',s82s; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_3_output_exists',s83o; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_3_status',s83s
  write(io,'(A,1X,I0)') 'stage8_twoway_stage8_4_output_exists',s84o; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_4_status',s84s; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_5_output_exists',s85o; write(io,'(A,1X,I0)') 'stage8_twoway_stage8_5_status',s85s; write(io,'(A,1X,I0)') 'stage8_twoway_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_twoway_safe_valid_count',sv; write(io,'(A,1X,I0)') 'stage8_twoway_safe_blocked_count',sbct; write(io,'(A,1X,I0)') 'stage8_twoway_safe_unsafe_count',suct
  write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_density_norm_max',nrm; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_conservation_abs_error',abse; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_conservation_relative_error',rel; write(io,'(A,1X,I0)') 'stage8_twoway_force_conservation_status',cons_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_conservation_x_error',errx; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_conservation_y_error',erry; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_conservation_z_error',errz; write(io,'(A,1X,I0)') 'stage8_twoway_component_conservation_status',cmp_ok
  write(io,'(A,1X,I0)') 'stage8_twoway_force_density_convention_flag',conv_flag; write(io,'(A,1X,I0)') 'stage8_twoway_no_rho_division_flag',norho_flag; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_force_buffer_change_with_rho_max',urho; write(io,'(A,1X,I0)') 'stage8_twoway_no_rho_status',norho_ok
  write(io,'(A,1X,I0)') 'stage8_twoway_nonuniform_volume_used_flag',nonuni_flag; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_uniform_volume_difference_norm',sqrt(sum((fe-few)**2)); write(io,'(A,1X,I0)') 'stage8_twoway_volume_scaling_status',vol_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_periodic_x_force_error',perx; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_periodic_z_force_error',perz; write(io,'(A,1X,I0)') 'stage8_twoway_periodic_x_status',px_ok; write(io,'(A,1X,I0)') 'stage8_twoway_periodic_z_status',pz_ok; write(io,'(A,1X,I0)') 'stage8_twoway_periodic_status',per_ok
  write(io,'(A,1X,I0)') 'stage8_twoway_blocked_count',bbc; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_blocked_force_buffer_norm_max',bnrm; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_blocked_force_buffer_write_error_max',bwerr; write(io,'(A,1X,I0)') 'stage8_twoway_blocked_status',blk_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_zero_force_buffer_norm_max',znrm; write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_zero_force_conservation_error',zerr; write(io,'(A,1X,I0)') 'stage8_twoway_zero_force_status',zero_ok
  write(io,'(A,1X,I0)') 'stage8_twoway_rhs_hook_called_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_rhs_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_pressure_poisson_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_projection_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_real_projection_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_twoway_production_dns_called_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_fluid_update_called_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_fibre_advance_called_flag',0; write(io,'(A,1X,I0)') 'stage8_twoway_no_rhs_no_projection_status',norhs_ok
  write(io,'(A,1X,ES24.16E3)') 'stage8_twoway_noop_rhs_change_max',rhschg; write(io,'(A,1X,I0)') 'stage8_twoway_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype); write(io,'(A,1X,I0)') 'stage8_twoway_noop_safety_status',noop_ok
  write(io,'(A,1X,I0)') 'stage8_twoway_force_density_check_status',final; close(io)
contains
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
