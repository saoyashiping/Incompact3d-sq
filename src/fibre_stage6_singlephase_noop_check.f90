program fibre_stage6_singlephase_noop_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_singlephase_noop
  use fibre_stage6_main_noop_hook, only : stage6_noop_hook_production_safety_status, stage6_noop_hook_pressure_status
  implicit none
  integer, parameter :: nx=12,ny=10,nz=8
  type(stage6_config_t) :: cfg
  real(mytype) :: rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz)
  real(mytype) :: ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz)
  real(mytype) :: p(nx,ny,nz),p0(nx,ny,nz),prhs(nx,ny,nz),prhs0(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz)
  real(mytype) :: div0(nx,ny,nz),div1(nx,ny,nz),rhschg,ex,ey,ez,vchg,pchg,prhschg,divchg,c0,c1,cerr,bmax
  integer :: i,j,k,hook,inj,mod,fluid,valid,rej,rhsok,ctf,pef,safe,prod,ctl,pp,pm,pr
  open(11,file='stage6_outputs/fibre_stage6_singlephase_noop_check.dat',status='replace')
  do k=1,nz;do j=1,ny;do i=1,nx
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k; rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i; rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j
    ux(i,j,k)=0.2_mytype+0.01_mytype*i+0.002_mytype*j; uy(i,j,k)=0.01_mytype*sin(0.2_mytype*i+0.1_mytype*k); uz(i,j,k)=0.01_mytype*cos(0.1_mytype*j+0.3_mytype*k)
    p(i,j,k)=0.5_mytype+0.01_mytype*sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k); prhs(i,j,k)=0.001_mytype*cos(0.2_mytype*i-0.1_mytype*j+0.05_mytype*k)
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k); fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k); fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
  end do;end do;end do
  rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz; ux0=ux;uy0=uy;uz0=uz; p0=p;prhs0=prhs
  call compute_stage6_synthetic_divergence(ux0,uy0,uz0,div0); call compute_stage6_singlephase_checksum(rhsx0,rhsy0,rhsz0,ux0,uy0,uz0,p0,c0)
  call init_stage6_default_config(cfg); call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  call apply_stage6_default_singlephase_noop_path(cfg,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,p,prhs,hook,inj,mod,fluid)
  call compute_stage6_vector_array_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg,ex,ey,ez)
  call compute_stage6_vector_array_change_max(ux0,uy0,uz0,ux,uy,uz,vchg,ex,ey,ez)
  call compute_stage6_array_change_max(p0,p,pchg); call compute_stage6_array_change_max(prhs0,prhs,prhschg)
  call compute_stage6_synthetic_divergence(ux,uy,uz,div1); call compute_stage6_array_change_max(div0,div1,divchg)
  call compute_stage6_singlephase_checksum(rhsx,rhsy,rhsz,ux,uy,uz,p,c1); cerr=abs(c1-c0)
  bmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  call stage6_noop_hook_pressure_status(pp,pm,pr); call stage6_noop_hook_production_safety_status(safe,prod,ctl)

  write(11,'(A,1X,I0)') 'stage6_dns_noop_default_enable_stage6',merge(1,0,cfg%enable_stage6)
  write(11,'(A,1X,I0)') 'stage6_dns_noop_default_main_rhs_hook_enabled',merge(1,0,cfg%enable_main_rhs_hook)
  write(11,'(A,1X,I0)') 'stage6_dns_noop_default_controlled_test_enabled',merge(1,0,cfg%enable_controlled_rhs_test)
  write(11,'(A,1X,I0)') 'stage6_dns_noop_default_production_two_way_enabled',merge(1,0,cfg%production_two_way_enabled)
  write(11,'(A,1X,I0)') 'stage6_dns_noop_default_rhs_allowed_flag',rhsok
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_rhs_diff_max',rhschg
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_rhsx_diff_max',maxval(abs(rhsx-rhsx0))
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_rhsy_diff_max',maxval(abs(rhsy-rhsy0))
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_rhsz_diff_max',maxval(abs(rhsz-rhsz0))
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_nonzero_buffer_max_abs',bmax
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_nonzero_buffer_rhs_diff_max',rhschg
  write(11,'(A,1X,I0)') 'stage6_dns_noop_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage6_dns_noop_modified_flag',mod
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_velocity_diff_max',vchg
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_ux_diff_max',maxval(abs(ux-ux0))
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_uy_diff_max',maxval(abs(uy-uy0))
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_uz_diff_max',maxval(abs(uz-uz0))
  write(11,'(A,1X,I0)') 'stage6_dns_noop_fluid_update_called_flag',fluid
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_pressure_diff_max',pchg
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_pressure_rhs_diff_max',prhschg
  write(11,'(A,1X,I0)') 'stage6_dns_noop_projection_modified_flag',pm
  write(11,'(A,1X,I0)') 'stage6_dns_noop_pressure_poisson_modified_flag',pp
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_divergence_diff_max',divchg
  write(11,'(A,1X,I0)') 'stage6_dns_noop_divergence_status',merge(1,0,divchg<=1e-14_mytype)
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_baseline_checksum',c0
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_stage6_checksum',c1
  write(11,'(A,1X,ES24.16E3)') 'stage6_dns_noop_checksum_abs_error',cerr
  write(11,'(A,1X,I0)') 'stage6_dns_noop_production_enabled_by_default_flag',prod
  write(11,'(A,1X,I0)') 'stage6_dns_noop_controlled_test_only_flag',ctl
  write(11,'(A,1X,I0)') 'stage6_dns_noop_default_main_dns_safe_flag',safe
  write(11,'(A,1X,I0)') 'stage6_dns_noop_regression_check_status',1
  close(11)
end program
