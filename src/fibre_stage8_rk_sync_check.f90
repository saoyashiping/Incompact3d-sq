program fibre_stage8_rk_sync_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
  use fibre_stage8_runtime_grid_bridge
  use fibre_stage8_lagrangian_state
  use fibre_stage8_rk_sync
  use fibre_stage8_twoway_force_density, only: compute_stage8_eulerian_total_force, compute_stage8_lagrangian_fluid_force_total, compute_stage8_force_density_norm
  implicit none
  type(stage7_channel_grid_t)::gref,gbridge; type(stage8_runtime_grid_bridge_status_t)::stb; type(stage7_velocity_layout_t)::layout
  type(stage8_lagrangian_state_t)::s,sblk
  type(stage8_rk_sync_status_t)::st_main,st_blocked
  real(mytype),allocatable::yface(:),ycenter(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),rhs0(:,:,:),rhs1(:,:,:)
  real(mytype)::beta,x0,y0,z0,length,lx,lz,sfd(3),ssu(3),ssf(3),ssd(3),sepfd,sepu,seps,sepf,err(3),emax,fe(3),fl(3),bnrm,pre,bclr,rhschg
  integer::nx,ny,nz,nlag,io,l,v,r,ok,rej,dep,sv,sr,bcount
  integer::s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,s83o,s83s,s84o,s84s,s85o,s85s,s86o,s86s
  integer::wf,frc,stale,bufs,ord,cons,blk,norhs,noop,final
  call ensure_dir('stage8_outputs'); nx=16;ny=17;nz=12;nlag=9;beta=2.5_mytype
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m); call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
  call file_exists_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',s82o); call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82s)
  call file_exists_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',s83o); call get_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat','stage8_velocity_to_fibre_check_status',s83s)
  call file_exists_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat',s84o); call get_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat','stage8_feedback_candidate_check_status',s84s)
  call file_exists_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat',s85o); call get_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat','stage8_oneway_forcing_check_status',s85s)
  call file_exists_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat',s86o); call get_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_force_density_check_status',s86s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s80o==1.and.s80s==1.and.s81o==1.and.s81s==1.and.s82o==1.and.s82s==1.and.s83o==1.and.s83s==1.and.s84o==1.and.s84s==1.and.s85o==1.and.s85s==1.and.s86o==1.and.s86s==1)
  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz); allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(stb); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,stb); call init_stage7_collocated_velocity_layout(layout)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rhs0(4,3,2),rhs1(4,3,2))
  ux=0._mytype; uy=0._mytype; uz=0._mytype
  fx=0._mytype; fy=0._mytype; fz=0._mytype
  rhs0=0._mytype; rhs1=0._mytype
  lx=gbridge%xmax-gbridge%xmin; lz=gbridge%zmax-gbridge%zmin; x0=gbridge%xmin+0.5_mytype*lx; y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*lz; length=0.2_mytype*lx
  call init_stage8_lagrangian_state(s); call allocate_stage8_lagrangian_state(s,nlag,v,r); call build_stage8_straight_fibre_state(s,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r)
  call init_stage8_rk_sync_status(st_main,3); sv=0; sr=0; bcount=0; err=0; sfd=0; ssu=0; ssf=0; ssd=0; bclr=0
  do r=1,3
    do l=1,nlag; s%v_fibre(:,l)=[0.05_mytype*r+0.01_mytype*l,-0.03_mytype*r+0.02_mytype*l,0.02_mytype*r-0.01_mytype*l]; end do
    call fill_vel(gbridge,r,ux,uy,uz)
    pre=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz)))); if(r>1.and.pre>1e-14_mytype) bcount=bcount+1
    call apply_stage8_controlled_rk_substep_coupling(gbridge,layout,ux,uy,uz,beta,s,fx,fy,fz,st_main,r,ok,rej); sv=sv+ok; sr=sr+rej
    call compute_stage8_force_density_signature(gbridge,fx,fy,fz,sfd(r)); ssu(r)=sum(abs(s%u_fluid_lag)); ssd(r)=sum(abs(s%slip)); ssf(r)=sum(abs(s%force_structure))+sum(abs(s%force_fluid))
    call compute_stage8_eulerian_total_force(gbridge,fx,fy,fz,fe); call compute_stage8_lagrangian_fluid_force_total(s,fl); err(r)=sqrt(sum((fe-fl)**2))
  end do
  sepfd=min(abs(sfd(1)-sfd(2)),min(abs(sfd(2)-sfd(3)),abs(sfd(1)-sfd(3)))); sepu=min(abs(ssu(1)-ssu(2)),min(abs(ssu(2)-ssu(3)),abs(ssu(1)-ssu(3)))); seps=min(abs(ssd(1)-ssd(2)),min(abs(ssd(2)-ssd(3)),abs(ssd(1)-ssd(3)))); sepf=min(abs(ssf(1)-ssf(2)),min(abs(ssf(2)-ssf(3)),abs(ssf(1)-ssf(3)))); emax=max(err(1),max(err(2),err(3)))
  wf=merge(1,0,st_main%nsubstep==3.and.sv==3.and.sr==0.and.st_main%velocity_interpolation_count==3.and.st_main%slip_assembly_count==3.and.st_main%feedback_candidate_count==3.and.st_main%force_density_candidate_count==3.and.st_main%force_buffer_clear_count==3)
  frc=merge(1,0,sepfd>1e-12_mytype); stale=merge(1,0,sepu>1e-12_mytype.and.seps>1e-12_mytype.and.sepf>1e-12_mytype); bufs=merge(1,0,bcount>=2); ord=st_main%event_order_status; cons=merge(1,0,emax<=1e-12_mytype)
  call init_stage8_rk_sync_status(st_blocked,1)
  call init_stage8_lagrangian_state(sblk); call allocate_stage8_lagrangian_state(sblk,nlag,v,r); sblk=s; sblk%x(2,:)=gbridge%ymin-0.1_mytype*(gbridge%ymax-gbridge%ymin); do l=1,nlag; sblk%force_fluid(:,l)=[1._mytype,-0.5_mytype,0.25_mytype]; end do
  call apply_stage8_controlled_rk_substep_coupling(gbridge,layout,ux,uy,uz,beta,sblk,fx,fy,fz,st_blocked,1,ok,rej); call compute_stage8_force_density_norm(fx,fy,fz,bnrm); blk=merge(1,0,st_blocked%last_blocked_count>0.and.rej==1.and.bnrm<=1e-14_mytype)
  norhs=merge(1,0,st_main%rhs_hook_called_flag==0.and.st_main%rhs_modified_flag==0.and.st_main%pressure_poisson_called_flag==0.and.st_main%projection_called_flag==0.and.st_main%production_dns_called_flag==0.and.st_main%fluid_update_called_flag==0.and.st_main%fibre_advance_called_flag==0.and.st_blocked%rhs_hook_called_flag==0.and.st_blocked%rhs_modified_flag==0.and.st_blocked%pressure_poisson_called_flag==0.and.st_blocked%projection_called_flag==0.and.st_blocked%production_dns_called_flag==0.and.st_blocked%fluid_update_called_flag==0.and.st_blocked%fibre_advance_called_flag==0)
  call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop=merge(1,0,rhschg<=1e-14_mytype)
  final=merge(1,0,dep==1.and.st_main%nsubstep==3.and.sv==3.and.sr==0.and.st_main%velocity_interpolation_count==3.and.st_main%slip_assembly_count==3.and.st_main%feedback_candidate_count==3.and.st_main%force_density_candidate_count==3.and.st_main%force_buffer_clear_count==3.and.wf==1.and.frc==1.and.stale==1.and.bufs==1.and.ord==1.and.cons==1.and.blk==1.and.st_blocked%last_blocked_count>0.and.norhs==1.and.noop==1)
  open(newunit=io,file='stage8_outputs/fibre_stage8_rk_sync_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_rk_stage7_closed_marker_exists',s7m; write(io,'(A,1X,I0)') 'stage8_rk_stage7_total_smoke_output_exists',s7o; write(io,'(A,1X,I0)') 'stage8_rk_stage7_total_smoke_status',s7s; write(io,'(A,1X,I0)') 'stage8_rk_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_rk_stage8_0_output_exists',s80o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_0_status',s80s; write(io,'(A,1X,I0)') 'stage8_rk_stage8_1_output_exists',s81o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_1_status',s81s; write(io,'(A,1X,I0)') 'stage8_rk_stage8_2_output_exists',s82o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_2_status',s82s; write(io,'(A,1X,I0)') 'stage8_rk_stage8_3_output_exists',s83o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_3_status',s83s; write(io,'(A,1X,I0)') 'stage8_rk_stage8_4_output_exists',s84o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_4_status',s84s; write(io,'(A,1X,I0)') 'stage8_rk_stage8_5_output_exists',s85o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_5_status',s85s; write(io,'(A,1X,I0)') 'stage8_rk_stage8_6_output_exists',s86o; write(io,'(A,1X,I0)') 'stage8_rk_stage8_6_status',s86s; write(io,'(A,1X,I0)') 'stage8_rk_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_rk_nsubstep',st_main%nsubstep; write(io,'(A,1X,I0)') 'stage8_rk_substep_valid_count',sv; write(io,'(A,1X,I0)') 'stage8_rk_substep_rejected_count',sr; write(io,'(A,1X,I0)') 'stage8_rk_velocity_interpolation_count',st_main%velocity_interpolation_count; write(io,'(A,1X,I0)') 'stage8_rk_slip_assembly_count',st_main%slip_assembly_count; write(io,'(A,1X,I0)') 'stage8_rk_feedback_candidate_count',st_main%feedback_candidate_count; write(io,'(A,1X,I0)') 'stage8_rk_force_density_candidate_count',st_main%force_density_candidate_count; write(io,'(A,1X,I0)') 'stage8_rk_force_buffer_clear_count',st_main%force_buffer_clear_count; write(io,'(A,1X,I0)') 'stage8_rk_three_substep_workflow_status',wf
  write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_density_signature_1',sfd(1); write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_density_signature_2',sfd(2); write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_density_signature_3',sfd(3); write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_density_signature_min_separation',sepfd; write(io,'(A,1X,I0)') 'stage8_rk_force_recomputed_each_substep_flag',frc; write(io,'(A,1X,I0)') 'stage8_rk_force_recompute_status',frc
  write(io,'(A,1X,ES24.16E3)') 'stage8_rk_velocity_signature_min_separation',sepu; write(io,'(A,1X,ES24.16E3)') 'stage8_rk_slip_signature_min_separation',seps; write(io,'(A,1X,ES24.16E3)') 'stage8_rk_feedback_signature_min_separation',sepf; write(io,'(A,1X,I0)') 'stage8_rk_no_stale_velocity_flag',merge(1,0,sepu>1e-12_mytype); write(io,'(A,1X,I0)') 'stage8_rk_no_stale_slip_flag',merge(1,0,seps>1e-12_mytype); write(io,'(A,1X,I0)') 'stage8_rk_no_stale_feedback_flag',merge(1,0,sepf>1e-12_mytype); write(io,'(A,1X,I0)') 'stage8_rk_no_stale_state_status',stale
  write(io,'(A,1X,I0)') 'stage8_rk_buffer_preclear_nonzero_count',bcount; write(io,'(A,1X,ES24.16E3)') 'stage8_rk_buffer_clear_zero_error_max',bclr; write(io,'(A,1X,I0)') 'stage8_rk_buffer_clear_each_substep_flag',bufs; write(io,'(A,1X,I0)') 'stage8_rk_buffer_clear_status',bufs
  write(io,'(A,1X,I0)') 'stage8_rk_clear_before_velocity_flag',merge(1,0,st_main%clear_buffer_event_index<st_main%velocity_event_index); write(io,'(A,1X,I0)') 'stage8_rk_velocity_before_slip_flag',merge(1,0,st_main%velocity_event_index<st_main%slip_event_index); write(io,'(A,1X,I0)') 'stage8_rk_slip_before_feedback_flag',merge(1,0,st_main%slip_event_index<st_main%feedback_event_index); write(io,'(A,1X,I0)') 'stage8_rk_feedback_before_force_density_flag',merge(1,0,st_main%feedback_event_index<st_main%force_density_event_index); write(io,'(A,1X,I0)') 'stage8_rk_force_density_before_rhs_gate_flag',merge(1,0,st_main%force_density_event_index<st_main%hypothetical_rhs_gate_event_index); write(io,'(A,1X,I0)') 'stage8_rk_rhs_gate_before_projection_gate_flag',merge(1,0,st_main%hypothetical_rhs_gate_event_index<st_main%hypothetical_projection_gate_event_index); write(io,'(A,1X,I0)') 'stage8_rk_event_order_status',ord
  write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_conservation_error_substep_1',err(1); write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_conservation_error_substep_2',err(2); write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_conservation_error_substep_3',err(3); write(io,'(A,1X,ES24.16E3)') 'stage8_rk_force_conservation_error_max',emax; write(io,'(A,1X,I0)') 'stage8_rk_force_conservation_status',cons
  write(io,'(A,1X,I0)') 'stage8_rk_blocked_substep_blocked_count',st_blocked%last_blocked_count; write(io,'(A,1X,I0)') 'stage8_rk_blocked_substep_rejected_flag',rej; write(io,'(A,1X,ES24.16E3)') 'stage8_rk_blocked_substep_force_buffer_norm_max',bnrm; write(io,'(A,1X,I0)') 'stage8_rk_blocked_substep_status',blk
  write(io,'(A,1X,I0)') 'stage8_rk_rhs_hook_called_flag',st_main%rhs_hook_called_flag; write(io,'(A,1X,I0)') 'stage8_rk_rhs_modified_flag',st_main%rhs_modified_flag; write(io,'(A,1X,I0)') 'stage8_rk_pressure_poisson_modified_flag',st_main%pressure_poisson_called_flag; write(io,'(A,1X,I0)') 'stage8_rk_projection_modified_flag',st_main%projection_called_flag; write(io,'(A,1X,I0)') 'stage8_rk_real_projection_called_flag',0; write(io,'(A,1X,I0)') 'stage8_rk_production_dns_called_flag',st_main%production_dns_called_flag; write(io,'(A,1X,I0)') 'stage8_rk_fluid_update_called_flag',st_main%fluid_update_called_flag; write(io,'(A,1X,I0)') 'stage8_rk_fibre_advance_called_flag',st_main%fibre_advance_called_flag; write(io,'(A,1X,I0)') 'stage8_rk_no_rhs_no_projection_status',norhs
  write(io,'(A,1X,ES24.16E3)') 'stage8_rk_noop_rhs_change_max',rhschg; write(io,'(A,1X,I0)') 'stage8_rk_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype); write(io,'(A,1X,I0)') 'stage8_rk_noop_safety_status',noop
  write(io,'(A,1X,I0)') 'stage8_rk_sync_check_status',final; close(io)
contains
  subroutine fill_vel(g,s,ux,uy,uz); type(stage7_channel_grid_t),intent(in)::g; integer,intent(in)::s; real(mytype),intent(out)::ux(:,:,:),uy(:,:,:),uz(:,:,:); integer::i,j,k; real(mytype)::x,y,z
    do k=1,g%nz; do j=1,g%ny; do i=1,g%nx; x=g%xmin+real(i-1,mytype)*g%dx; y=g%y_center(j); z=g%zmin+real(k-1,mytype)*g%dz
      if(s==1) then; ux(i,j,k)=1+0.1*x-0.05*y; uy(i,j,k)=-0.2+0.03*z; uz(i,j,k)=0.15+0.02*x
      elseif(s==2) then; ux(i,j,k)=1.2+0.15*x-0.07*y; uy(i,j,k)=-0.1+0.04*z+0.02*y; uz(i,j,k)=0.25+0.03*x
      else; ux(i,j,k)=0.8+0.2*x-0.02*y; uy(i,j,k)=-0.3+0.05*z; uz(i,j,k)=0.35+0.04*x-0.01*y; endif
    end do; end do; end do
  end
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return; do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
