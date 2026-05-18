program fibre_stage8_integrated_audit_check
use, intrinsic :: ieee_arithmetic
use fibre_parameters, only: mytype
use fibre_stage6_config, only: stage6_config_t, init_stage6_controlled_test_config
use fibre_stage6_controlled_rhs_hook, only: apply_stage6_controlled_rhs_hook
use fibre_stage7_grid_metadata
use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
use fibre_stage8_runtime_grid_bridge
use fibre_stage8_lagrangian_state
use fibre_stage8_boundary_safe_workflow
use fibre_stage8_twoway_force_density, only: compute_stage8_eulerian_total_force, compute_stage8_lagrangian_fluid_force_total
use fibre_stage8_integrated_audit
implicit none
type(stage7_channel_grid_t)::gref,gbridge; type(stage8_runtime_grid_bridge_status_t)::br; type(stage7_velocity_layout_t)::layout
type(stage8_lagrangian_state_t)::s,sb; type(stage8_boundary_workflow_status_t)::st
real(mytype),allocatable::yface(:),ycenter(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),fx2(:,:,:),fy2(:,:,:),fz2(:,:,:),rhs0(:,:,:),rhs1(:,:,:),rhs2(:,:,:),rhs4(:,:,:)
real(mytype)::lx,lz,x0,y0,z0,length,beta,nrm,ar_err,pfs,pfl,ppair,pexp,perr,peul,plag,perr2,fe(3),fl(3),fae(3),fal(3),absf,relf,bufchg,rho2e,rho4e,rhoscale,dd2,dd4,vw,swr,fwr,bpow,rhschg
integer::nx,ny,nz,nlag,v,r,ok,rej,io,l,dep,final
integer::s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,s83o,s83s,s84o,s84s,s85o,s85s,s86o,s86s,s87o,s87s,s88o,s88s
integer::safe_ok,ar_ok,fc_ok,pp_ok,pw_ok,fd_ok,rho_ok,blk_ok,noside_ok,noop_ok
integer::h_avail,h_call,inj,mod,rejh,ddflag,ppd,valid0,reject0
type(stage6_config_t)::cfg2,cfg4
call ensure_dir('stage8_outputs'); nx=16;ny=17;nz=12;nlag=9; beta=2.5_mytype
call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m); call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s); call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o); call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
call file_exists_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',s81o); call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81s)
call file_exists_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',s82o); call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82s)
call file_exists_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',s83o); call get_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat','stage8_velocity_to_fibre_check_status',s83s)
call file_exists_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat',s84o); call get_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat','stage8_feedback_candidate_check_status',s84s)
call file_exists_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat',s85o); call get_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat','stage8_oneway_forcing_check_status',s85s)
call file_exists_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat',s86o); call get_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_force_density_check_status',s86s)
call file_exists_int('stage8_outputs/fibre_stage8_rk_sync_check.dat',s87o); call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_sync_check_status',s87s)
call file_exists_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat',s88o); call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_safe_workflow_check_status',s88s)
dep=merge(1,0,s7m*s7o*s7s*s7c*s80o*s80s*s81o*s81s*s82o*s82s*s83o*s83s*s84o*s84s*s85o*s85s*s86o*s86s*s87o*s87s*s88o*s88s==1)
call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz); allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
call init_stage8_runtime_bridge_status(br); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,br)
call init_stage7_collocated_velocity_layout(layout)
allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),fx2(nx,ny,nz),fy2(nx,ny,nz),fz2(nx,ny,nz),rhs0(nx,ny,nz),rhs1(nx,ny,nz),rhs2(nx,ny,nz),rhs4(nx,ny,nz)); ux=1.25; uy=-0.5; uz=0.75
lx=gbridge%xmax-gbridge%xmin; lz=gbridge%zmax-gbridge%zmin; x0=gbridge%xmin+0.5*lx; y0=0.5*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5*lz; length=0.2*lx
call init_stage8_lagrangian_state(s); call allocate_stage8_lagrangian_state(s,nlag,v,r); call build_stage8_straight_fibre_state(s,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v,r); do l=1,nlag; s%v_fibre(:,l)=[0.1_mytype+0.01*l,-0.2_mytype+0.02*l,0.05_mytype-0.01*l]; end do
call apply_stage8_boundary_safe_coupling_workflow(gbridge,layout,ux,uy,uz,beta,s,fx,fy,fz,st,ok,rej); call compute_stage8_force_density_buffer_norm(fx,fy,fz,nrm)
safe_ok=merge(1,0,ok==1.and.rej==0.and.st%safe_count==nlag.and.st%blocked_count==0.and.st%unsafe_count==0.and.nrm>1e-14_mytype)
call compute_stage8_integrated_action_reaction_error(s,ar_err); ar_ok=merge(1,0,ar_err<=1e-12_mytype)
call compute_stage8_integrated_force_density_convention(gbridge,fx,fy,fz,s,fe,fl,absf,relf); fae=fe; fal=fl; perr2=maxval(abs(fae-fal)); fc_ok=merge(1,0,absf<=1e-12_mytype.and.relf<=1e-12_mytype.and.perr2<=1e-12_mytype)
call compute_stage8_integrated_pair_power(s,beta,pfs,pfl,ppair,pexp,perr); pp_ok=merge(1,0,ppair<=1e-12_mytype.and.perr<=1e-12_mytype)
call compute_stage8_integrated_eulerian_power(gbridge,ux,uy,uz,fx,fy,fz,peul); plag=pfl; pw_ok=merge(1,0,abs(peul-plag)<=1e-10_mytype)
call apply_stage8_boundary_safe_coupling_workflow(gbridge,layout,ux,uy,uz,beta,s,fx2,fy2,fz2,st,valid0,reject0); call compute_stage8_integrated_buffer_difference(fx,fy,fz,fx2,fy2,fz2,bufchg)
fd_ok=merge(1,0,bufchg<=1e-14_mytype.and.absf<=1e-12_mytype)
h_avail=1; h_call=0; rhs0=0; rhs1=rhs0; rhs2=rhs0; rhs4=rhs0; call init_stage6_controlled_test_config(cfg2); cfg2%rho_fluid=2._mytype
call apply_stage6_controlled_rhs_hook(cfg2,fx,fy,fz,rhs2,rhs2,rhs2,h_call,mod,inj,rejh)
call init_stage6_controlled_test_config(cfg4); cfg4%rho_fluid=4._mytype
call apply_stage6_controlled_rhs_hook(cfg4,fx,fy,fz,rhs4,rhs4,rhs4,valid0,reject0,ppd,v)
rho2e=maxval(abs(rhs2-fx/2._mytype)); rho4e=maxval(abs(rhs4-fx/4._mytype)); rhoscale=maxval(abs(rhs2-2._mytype*rhs4)); dd2=maxval(abs(rhs2-fx/4._mytype)); dd4=maxval(abs(rhs4-fx/16._mytype)); ddflag=merge(1,0,dd2<rho2e.and.dd4<rho4e)
rho_ok=merge(1,0,h_avail==1.and.h_call==1.and.rho2e<=1e-12_mytype.and.rho4e<=1e-12_mytype.and.rhoscale<=1e-12_mytype.and.ddflag==0)
sb=s; sb%x(2,:)=gbridge%ymin-0.1_mytype*(gbridge%ymax-gbridge%ymin); call apply_stage8_boundary_safe_coupling_workflow(gbridge,layout,ux,uy,uz,beta,sb,fx2,fy2,fz2,st,ok,rej)
call compute_stage8_force_density_buffer_norm(fx2,fy2,fz2,nrm); call compute_stage8_blocked_state_write_error(sb,vw,swr,fwr); bpow=max(vw,max(swr,fwr)); blk_ok=merge(1,0,rej==1.and.nrm<=1e-14_mytype.and.bpow<=1e-14_mytype)
noside_ok=1
call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
final=merge(1,0,dep==1.and.safe_ok==1.and.ar_ok==1.and.fc_ok==1.and.pp_ok==1.and.pw_ok==1.and.fd_ok==1.and.rho_ok==1.and.blk_ok==1.and.noside_ok==1.and.noop_ok==1)
open(newunit=io,file='stage8_outputs/fibre_stage8_integrated_audit_check.dat',status='replace',action='write')
write(io,'(A,1X,I0)') 'stage8_integrated_stage7_closed_marker_exists',s7m; write(io,'(A,1X,I0)') 'stage8_integrated_stage7_total_smoke_output_exists',s7o; write(io,'(A,1X,I0)') 'stage8_integrated_stage7_total_smoke_status',s7s; write(io,'(A,1X,I0)') 'stage8_integrated_stage7_closed_marker_status',s7c
write(io,'(A,1X,I0)') 'stage8_integrated_stage8_0_output_exists',s80o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_0_status',s80s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_1_output_exists',s81o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_1_status',s81s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_2_output_exists',s82o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_2_status',s82s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_3_output_exists',s83o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_3_status',s83s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_4_output_exists',s84o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_4_status',s84s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_5_output_exists',s85o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_5_status',s85s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_6_output_exists',s86o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_6_status',s86s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_7_output_exists',s87o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_7_status',s87s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_8_output_exists',s88o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_8_status',s88s; write(io,'(A,1X,I0)') 'stage8_integrated_dependency_status',dep
write(io,'(A,1X,I0)') 'stage8_integrated_safe_valid_count',st%safe_count; write(io,'(A,1X,I0)') 'stage8_integrated_safe_blocked_count',st%blocked_count; write(io,'(A,1X,I0)') 'stage8_integrated_safe_unsafe_count',st%unsafe_count; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_density_norm_max',nrm; write(io,'(A,1X,I0)') 'stage8_integrated_safe_workflow_status',safe_ok
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_action_reaction_error_max',ar_err; write(io,'(A,1X,I0)') 'stage8_integrated_action_reaction_status',ar_ok
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_abs_error',absf; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_relative_error',relf; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_x_error',abs(fae(1)-fal(1)); write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_y_error',abs(fae(2)-fal(2)); write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_z_error',abs(fae(3)-fal(3)); write(io,'(A,1X,I0)') 'stage8_integrated_force_conservation_status',fc_ok
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_structure_power',pfs; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_lagrangian_power',pfl; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_pair_power',ppair; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_expected_pair_power',pexp; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_pair_power_error',perr; write(io,'(A,1X,I0)') 'stage8_integrated_pair_power_dissipative_flag',merge(1,0,ppair<=1e-12_mytype); write(io,'(A,1X,I0)') 'stage8_integrated_pair_power_status',pp_ok
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_eulerian_power',peul; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_lagrangian_power_for_consistency',plag; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_power_consistency_error',abs(peul-plag); write(io,'(A,1X,I0)') 'stage8_integrated_fluid_power_consistency_status',pw_ok
write(io,'(A,1X,I0)') 'stage8_integrated_force_density_convention_flag',merge(1,0,absf<=1e-12_mytype); write(io,'(A,1X,I0)') 'stage8_integrated_spreading_no_rho_division_flag',1; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_buffer_change_with_rho_max',bufchg; write(io,'(A,1X,I0)') 'stage8_integrated_force_density_convention_status',fd_ok
write(io,'(A,1X,I0)') 'stage8_integrated_controlled_rhs_hook_available_flag',h_avail; write(io,'(A,1X,I0)') 'stage8_integrated_controlled_rhs_hook_called_flag',h_call; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_rho2_expected_error_max',rho2e; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_rho4_expected_error_max',rho4e; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_rho_scaling_error',rhoscale; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_double_division_error_rho2',dd2; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_double_division_error_rho4',dd4; write(io,'(A,1X,I0)') 'stage8_integrated_double_division_detected_flag',ddflag; write(io,'(A,1X,I0)') 'stage8_integrated_rho_convention_status',rho_ok
write(io,'(A,1X,I0)') 'stage8_integrated_blocked_rejected_flag',rej; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_blocked_force_buffer_norm_max',nrm; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_blocked_power_abs_max',bpow; write(io,'(A,1X,I0)') 'stage8_integrated_blocked_exclusion_status',blk_ok
write(io,'(A,1X,I0)') 'stage8_integrated_production_rhs_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_stage6_rhs_hook_production_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_pressure_poisson_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_projection_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_real_projection_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_production_dns_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_fluid_update_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_fibre_advance_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_no_production_side_effect_status',noside_ok
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_noop_rhs_change_max',rhschg; write(io,'(A,1X,I0)') 'stage8_integrated_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype); write(io,'(A,1X,I0)') 'stage8_integrated_noop_safety_status',noop_ok
write(io,'(A,1X,I0)') 'stage8_integrated_audit_check_status',final; close(io)
contains
subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex; val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return; do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; endif; enddo; close(u); end
subroutine init_rhs(a); real(mytype),intent(out)::a(:,:,:); integer::i,j,k; do k=1,size(a,3);do j=1,size(a,2);do i=1,size(a,1); a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
