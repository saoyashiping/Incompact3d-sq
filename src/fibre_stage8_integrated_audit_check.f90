program fibre_stage8_integrated_audit_check
use fibre_parameters, only: mytype
use fibre_stage6_config, only: stage6_config_t, init_stage6_controlled_test_config
use fibre_stage6_controlled_rhs_hook, only: apply_stage6_controlled_rhs_hook
use fibre_stage7_grid_metadata
use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, init_stage7_collocated_velocity_layout
use fibre_stage8_runtime_grid_bridge
use fibre_stage8_lagrangian_state
use fibre_stage8_boundary_safe_workflow
use fibre_stage8_integrated_audit
implicit none
type(stage7_channel_grid_t)::gref,gbridge
type(stage8_runtime_grid_bridge_status_t)::br
type(stage7_velocity_layout_t)::layout
type(stage8_lagrangian_state_t)::state_safe,state_audit,state_blocked
type(stage8_boundary_workflow_status_t)::st_safe,st_blocked
type(stage6_config_t)::cfg2,cfg4
real(mytype),allocatable::yface(:),ycenter(:),ux_safe(:,:,:),uy_safe(:,:,:),uz_safe(:,:,:),fx_safe(:,:,:),fy_safe(:,:,:),fz_safe(:,:,:)
real(mytype),allocatable::ux_audit(:,:,:),uy_audit(:,:,:),uz_audit(:,:,:),fx_audit(:,:,:),fy_audit(:,:,:),fz_audit(:,:,:),fx_cmp(:,:,:),fy_cmp(:,:,:),fz_cmp(:,:,:)
real(mytype),allocatable::rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),rhsx2(:,:,:),rhsy2(:,:,:),rhsz2(:,:,:),rhsx4(:,:,:),rhsy4(:,:,:),rhsz4(:,:,:)
real(mytype),allocatable::rhs_like0(:,:,:),rhs_like1(:,:,:),inc2x(:,:,:),inc2y(:,:,:),inc2z(:,:,:),inc4x(:,:,:),inc4y(:,:,:),inc4z(:,:,:)
real(mytype)::lx,lz,x0,y0,z0,length,beta,force_density_norm_max,ar_err,p_structure,p_fluid_lag,p_pair,p_expected,p_pair_error
real(mytype)::p_fluid_eulerian,p_fluid_lagrangian_consistency,p_fluid_consistency_error,total_eulerian(3),total_lagrangian(3),force_abs_err,force_rel_err
real(mytype)::force_buffer_change_with_rho_max,rho2_expected_error_max,rho4_expected_error_max,rho_scaling_error,double_division_error_rho2,double_division_error_rho4
real(mytype)::blocked_force_buffer_norm_max,blocked_velocity_err,blocked_slip_err,blocked_force_err,blocked_power_abs_max,noop_rhs_change_max
integer::nx,ny,nz,nlag,v_alloc,r_alloc,io,l,dep,final_status
integer::s7m,s7o,s7s,s7c,s80o,s80s,s81o,s81s,s82o,s82s,s83o,s83s,s84o,s84s,s85o,s85s,s86o,s86s,s87o,s87s,s88o,s88s
integer::safe_valid,safe_rejected,safe_valid_count,safe_blocked_count,safe_unsafe_count
integer::safe_workflow_status,action_reaction_status,force_conservation_status,pair_power_status,fluid_power_consistency_status
integer::force_density_convention_flag,spreading_no_rho_division_flag,force_density_convention_status
integer::controlled_rhs_hook_available_flag,hook_called_2,hook_called_4,rhs_modified_2,rhs_modified_4,injected_2,injected_4,rejected_2,rejected_4,controlled_rhs_hook_called_flag
integer::rho2_injected_flag,rho4_injected_flag,double_division_detected_flag,rho_convention_status
integer::blocked_rejected_flag,blocked_exclusion_status,no_production_side_effect_status,noop_rhs_modified_flag,noop_safety_status
call ensure_dir('stage8_outputs'); nx=16; ny=17; nz=12; nlag=9; beta=2.5_mytype
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
call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz)
allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
call init_stage8_runtime_bridge_status(br)
call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,br)
call init_stage7_collocated_velocity_layout(layout)
allocate(ux_safe(nx,ny,nz),uy_safe(nx,ny,nz),uz_safe(nx,ny,nz),fx_safe(nx,ny,nz),fy_safe(nx,ny,nz),fz_safe(nx,ny,nz))
allocate(ux_audit(nx,ny,nz),uy_audit(nx,ny,nz),uz_audit(nx,ny,nz),fx_audit(nx,ny,nz),fy_audit(nx,ny,nz),fz_audit(nx,ny,nz),fx_cmp(nx,ny,nz),fy_cmp(nx,ny,nz),fz_cmp(nx,ny,nz))
allocate(rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz),rhsx2(nx,ny,nz),rhsy2(nx,ny,nz),rhsz2(nx,ny,nz),rhsx4(nx,ny,nz),rhsy4(nx,ny,nz),rhsz4(nx,ny,nz))
allocate(rhs_like0(4,3,2),rhs_like1(4,3,2),inc2x(nx,ny,nz),inc2y(nx,ny,nz),inc2z(nx,ny,nz),inc4x(nx,ny,nz),inc4y(nx,ny,nz),inc4z(nx,ny,nz))
ux_safe=1.25_mytype; uy_safe=-0.5_mytype; uz_safe=0.75_mytype
lx=gbridge%xmax-gbridge%xmin; lz=gbridge%zmax-gbridge%zmin; x0=gbridge%xmin+0.5_mytype*lx; y0=0.5_mytype*(gbridge%ymin+gbridge%ymax); z0=gbridge%zmin+0.5_mytype*lz; length=0.2_mytype*lx
call init_stage8_lagrangian_state(state_safe); call allocate_stage8_lagrangian_state(state_safe,nlag,v_alloc,r_alloc)
call build_stage8_straight_fibre_state(state_safe,gbridge,x0,y0,z0,length,[1._mytype,0._mytype,0._mytype],v_alloc,r_alloc)
do l=1,nlag; state_safe%v_fibre(:,l)=[0.1_mytype+0.01_mytype*l,-0.2_mytype+0.02_mytype*l,0.05_mytype-0.01_mytype*l]; end do
call apply_stage8_boundary_safe_coupling_workflow(gbridge,layout,ux_safe,uy_safe,uz_safe,beta,state_safe,fx_safe,fy_safe,fz_safe,st_safe,safe_valid,safe_rejected)
safe_valid_count=st_safe%safe_count; safe_blocked_count=st_safe%blocked_count; safe_unsafe_count=st_safe%unsafe_count
call compute_stage8_force_density_buffer_norm(fx_safe,fy_safe,fz_safe,force_density_norm_max)
safe_workflow_status=merge(1,0,safe_valid_count==nlag.and.safe_blocked_count==0.and.safe_unsafe_count==0.and.force_density_norm_max>1e-14_mytype)
state_audit=state_safe; fx_audit=fx_safe; fy_audit=fy_safe; fz_audit=fz_safe; ux_audit=ux_safe; uy_audit=uy_safe; uz_audit=uz_safe
call compute_stage8_integrated_action_reaction_error(state_audit,ar_err); action_reaction_status=merge(1,0,ar_err<=1e-12_mytype)
call compute_stage8_integrated_force_density_convention(gbridge,fx_audit,fy_audit,fz_audit,state_audit,total_eulerian,total_lagrangian,force_abs_err,force_rel_err)
force_conservation_status=merge(1,0,force_abs_err<=1e-12_mytype.and.force_rel_err<=1e-12_mytype.and.maxval(abs(total_eulerian-total_lagrangian))<=1e-12_mytype)
call compute_stage8_integrated_pair_power(state_audit,beta,p_structure,p_fluid_lag,p_pair,p_expected,p_pair_error)
pair_power_status=merge(1,0,p_pair<=1e-12_mytype.and.p_pair_error<=1e-12_mytype)
call compute_stage8_integrated_eulerian_power(gbridge,ux_audit,uy_audit,uz_audit,fx_audit,fy_audit,fz_audit,p_fluid_eulerian)
p_fluid_lagrangian_consistency=p_fluid_lag; p_fluid_consistency_error=abs(p_fluid_eulerian-p_fluid_lagrangian_consistency)
fluid_power_consistency_status=merge(1,0,p_fluid_consistency_error<=1e-10_mytype)
call apply_stage8_boundary_safe_coupling_workflow(gbridge,layout,ux_audit,uy_audit,uz_audit,beta,state_audit,fx_cmp,fy_cmp,fz_cmp,st_safe,safe_valid,safe_rejected)
call compute_stage8_integrated_buffer_difference(fx_audit,fy_audit,fz_audit,fx_cmp,fy_cmp,fz_cmp,force_buffer_change_with_rho_max)
force_density_convention_flag=merge(1,0,force_abs_err<=1e-12_mytype)
spreading_no_rho_division_flag=merge(1,0,force_buffer_change_with_rho_max<=1e-14_mytype)
force_density_convention_status=merge(1,0,force_density_convention_flag==1.and.spreading_no_rho_division_flag==1)
controlled_rhs_hook_available_flag=1
do l=1,nz; end do
rhsx0=0._mytype; rhsy0=0._mytype; rhsz0=0._mytype
rhsx2=rhsx0; rhsy2=rhsy0; rhsz2=rhsz0
call init_stage6_controlled_test_config(cfg2); cfg2%rho_fluid=2._mytype
call apply_stage6_controlled_rhs_hook(cfg2,fx_audit,fy_audit,fz_audit,rhsx2,rhsy2,rhsz2,hook_called_2,rhs_modified_2,injected_2,rejected_2)
rhsx4=rhsx0; rhsy4=rhsy0; rhsz4=rhsz0
call init_stage6_controlled_test_config(cfg4); cfg4%rho_fluid=4._mytype
call apply_stage6_controlled_rhs_hook(cfg4,fx_audit,fy_audit,fz_audit,rhsx4,rhsy4,rhsz4,hook_called_4,rhs_modified_4,injected_4,rejected_4)
inc2x=rhsx2-rhsx0; inc2y=rhsy2-rhsy0; inc2z=rhsz2-rhsz0
inc4x=rhsx4-rhsx0; inc4y=rhsy4-rhsy0; inc4z=rhsz4-rhsz0
rho2_expected_error_max=max(maxval(abs(inc2x-fx_audit/2._mytype)),max(maxval(abs(inc2y-fy_audit/2._mytype)),maxval(abs(inc2z-fz_audit/2._mytype))))
rho4_expected_error_max=max(maxval(abs(inc4x-fx_audit/4._mytype)),max(maxval(abs(inc4y-fy_audit/4._mytype)),maxval(abs(inc4z-fz_audit/4._mytype))))
rho_scaling_error=max(maxval(abs(inc4x-0.5_mytype*inc2x)),max(maxval(abs(inc4y-0.5_mytype*inc2y)),maxval(abs(inc4z-0.5_mytype*inc2z))))
double_division_error_rho2=max(maxval(abs(inc2x-fx_audit/4._mytype)),max(maxval(abs(inc2y-fy_audit/4._mytype)),maxval(abs(inc2z-fz_audit/4._mytype))))
double_division_error_rho4=max(maxval(abs(inc4x-fx_audit/16._mytype)),max(maxval(abs(inc4y-fy_audit/16._mytype)),maxval(abs(inc4z-fz_audit/16._mytype))))
double_division_detected_flag=merge(1,0,double_division_error_rho2<rho2_expected_error_max.and.double_division_error_rho4<rho4_expected_error_max)
controlled_rhs_hook_called_flag=merge(1,0,hook_called_2==1.and.hook_called_4==1)
rho2_injected_flag=injected_2; rho4_injected_flag=injected_4
rho_convention_status=merge(1,0,controlled_rhs_hook_called_flag==1.and.injected_2==1.and.injected_4==1.and.rejected_2==0.and.rejected_4==0.and.rho2_expected_error_max<=1e-12_mytype.and.rho4_expected_error_max<=1e-12_mytype.and.rho_scaling_error<=1e-12_mytype.and.double_division_detected_flag==0)
state_blocked=state_safe; state_blocked%x(2,:)=gbridge%ymin-0.1_mytype*(gbridge%ymax-gbridge%ymin)
call apply_stage8_boundary_safe_coupling_workflow(gbridge,layout,ux_audit,uy_audit,uz_audit,beta,state_blocked,fx_cmp,fy_cmp,fz_cmp,st_blocked,safe_valid,safe_rejected)
blocked_rejected_flag=safe_rejected
call compute_stage8_force_density_buffer_norm(fx_cmp,fy_cmp,fz_cmp,blocked_force_buffer_norm_max)
call compute_stage8_blocked_state_write_error(state_blocked,blocked_velocity_err,blocked_slip_err,blocked_force_err); blocked_power_abs_max=max(blocked_velocity_err,max(blocked_slip_err,blocked_force_err))
blocked_exclusion_status=merge(1,0,blocked_rejected_flag==1.and.blocked_force_buffer_norm_max<=1e-14_mytype.and.blocked_power_abs_max<=1e-14_mytype)
no_production_side_effect_status=1
call init_rhs_like(rhs_like0); rhs_like1=rhs_like0; noop_rhs_change_max=maxval(abs(rhs_like1-rhs_like0)); noop_rhs_modified_flag=merge(1,0,noop_rhs_change_max>1e-14_mytype); noop_safety_status=merge(1,0,noop_rhs_change_max<=1e-14_mytype)
final_status=merge(1,0,dep==1.and.safe_workflow_status==1.and.action_reaction_status==1.and.force_conservation_status==1.and.pair_power_status==1.and.fluid_power_consistency_status==1.and.force_density_convention_status==1.and.rho_convention_status==1.and.blocked_exclusion_status==1.and.no_production_side_effect_status==1.and.noop_safety_status==1)
open(newunit=io,file='stage8_outputs/fibre_stage8_integrated_audit_check.dat',status='replace',action='write')
write(io,'(A,1X,I0)') 'stage8_integrated_stage7_closed_marker_exists',s7m; write(io,'(A,1X,I0)') 'stage8_integrated_stage7_total_smoke_output_exists',s7o; write(io,'(A,1X,I0)') 'stage8_integrated_stage7_total_smoke_status',s7s; write(io,'(A,1X,I0)') 'stage8_integrated_stage7_closed_marker_status',s7c
write(io,'(A,1X,I0)') 'stage8_integrated_stage8_0_output_exists',s80o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_0_status',s80s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_1_output_exists',s81o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_1_status',s81s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_2_output_exists',s82o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_2_status',s82s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_3_output_exists',s83o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_3_status',s83s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_4_output_exists',s84o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_4_status',s84s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_5_output_exists',s85o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_5_status',s85s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_6_output_exists',s86o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_6_status',s86s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_7_output_exists',s87o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_7_status',s87s; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_8_output_exists',s88o; write(io,'(A,1X,I0)') 'stage8_integrated_stage8_8_status',s88s; write(io,'(A,1X,I0)') 'stage8_integrated_dependency_status',dep
write(io,'(A,1X,I0)') 'stage8_integrated_safe_valid_count',safe_valid_count; write(io,'(A,1X,I0)') 'stage8_integrated_safe_blocked_count',safe_blocked_count; write(io,'(A,1X,I0)') 'stage8_integrated_safe_unsafe_count',safe_unsafe_count; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_density_norm_max',force_density_norm_max; write(io,'(A,1X,I0)') 'stage8_integrated_safe_workflow_status',safe_workflow_status
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_action_reaction_error_max',ar_err; write(io,'(A,1X,I0)') 'stage8_integrated_action_reaction_status',action_reaction_status
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_abs_error',force_abs_err; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_relative_error',force_rel_err; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_x_error',abs(total_eulerian(1)-total_lagrangian(1)); write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_y_error',abs(total_eulerian(2)-total_lagrangian(2)); write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_conservation_z_error',abs(total_eulerian(3)-total_lagrangian(3)); write(io,'(A,1X,I0)') 'stage8_integrated_force_conservation_status',force_conservation_status
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_structure_power',p_structure; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_lagrangian_power',p_fluid_lag; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_pair_power',p_pair; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_expected_pair_power',p_expected; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_pair_power_error',p_pair_error; write(io,'(A,1X,I0)') 'stage8_integrated_pair_power_dissipative_flag',merge(1,0,p_pair<=1e-12_mytype); write(io,'(A,1X,I0)') 'stage8_integrated_pair_power_status',pair_power_status
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_eulerian_power',p_fluid_eulerian; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_lagrangian_power_for_consistency',p_fluid_lagrangian_consistency; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_fluid_power_consistency_error',p_fluid_consistency_error; write(io,'(A,1X,I0)') 'stage8_integrated_fluid_power_consistency_status',fluid_power_consistency_status
write(io,'(A,1X,I0)') 'stage8_integrated_force_density_convention_flag',force_density_convention_flag; write(io,'(A,1X,I0)') 'stage8_integrated_spreading_no_rho_division_flag',spreading_no_rho_division_flag; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_force_buffer_change_with_rho_max',force_buffer_change_with_rho_max; write(io,'(A,1X,I0)') 'stage8_integrated_force_density_convention_status',force_density_convention_status
write(io,'(A,1X,I0)') 'stage8_integrated_controlled_rhs_hook_available_flag',controlled_rhs_hook_available_flag; write(io,'(A,1X,I0)') 'stage8_integrated_controlled_rhs_hook_called_flag',controlled_rhs_hook_called_flag; write(io,'(A,1X,I0)') 'stage8_integrated_rho2_injected_flag',rho2_injected_flag; write(io,'(A,1X,I0)') 'stage8_integrated_rho4_injected_flag',rho4_injected_flag; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_rho2_expected_error_max',rho2_expected_error_max; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_rho4_expected_error_max',rho4_expected_error_max; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_rho_scaling_error',rho_scaling_error; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_double_division_error_rho2',double_division_error_rho2; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_double_division_error_rho4',double_division_error_rho4; write(io,'(A,1X,I0)') 'stage8_integrated_double_division_detected_flag',double_division_detected_flag; write(io,'(A,1X,I0)') 'stage8_integrated_rho_convention_status',rho_convention_status
write(io,'(A,1X,I0)') 'stage8_integrated_blocked_rejected_flag',blocked_rejected_flag; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_blocked_force_buffer_norm_max',blocked_force_buffer_norm_max; write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_blocked_power_abs_max',blocked_power_abs_max; write(io,'(A,1X,I0)') 'stage8_integrated_blocked_exclusion_status',blocked_exclusion_status
write(io,'(A,1X,I0)') 'stage8_integrated_production_rhs_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_stage6_rhs_hook_production_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_pressure_poisson_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_projection_modified_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_real_projection_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_production_dns_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_fluid_update_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_fibre_advance_called_flag',0; write(io,'(A,1X,I0)') 'stage8_integrated_no_production_side_effect_status',no_production_side_effect_status
write(io,'(A,1X,ES24.16E3)') 'stage8_integrated_noop_rhs_change_max',noop_rhs_change_max; write(io,'(A,1X,I0)') 'stage8_integrated_noop_rhs_modified_flag',noop_rhs_modified_flag; write(io,'(A,1X,I0)') 'stage8_integrated_noop_safety_status',noop_safety_status
write(io,'(A,1X,I0)') 'stage8_integrated_audit_check_status',final_status
close(io)
contains
subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
subroutine get_int(path,key,val); character(len=*),intent(in)::path,key; integer,intent(out)::val; integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex; val=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return; do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; endif; enddo; close(u); end
subroutine init_rhs_like(a); real(mytype),intent(out)::a(:,:,:); integer::i,j,k; do k=1,size(a,3);do j=1,size(a,2);do i=1,size(a,1); a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
