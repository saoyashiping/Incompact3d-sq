program fibre_stage7_config_check
  use fibre_parameters, only : mytype
  use fibre_stage7_config
  implicit none
  type(stage7_config_t) :: c
  integer :: io, valid, rej, rhs_allowed, ctest, pdns, ptwo
  integer :: s6mark, dep, modf, ppm, pm, rpc, pdc, fuc, fac, final_status
  integer :: d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12
  integer :: f_rej,f_valid,g_rej,g_valid,h_rej,h_rhs,i_rej,i_valid,i_rhs
  real(mytype) :: chg
  integer, parameter :: nx=8,ny=6,nz=5
  real(mytype) :: rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz)
  integer :: i,j,k
  logical :: ex

  open(newunit=io,file='stage7_outputs/fibre_stage7_config_check.dat',status='replace',action='write')
  inquire(file='stage6_outputs/STAGE6_CLOSED.md',exist=ex); s6mark=merge(1,0,ex); dep=s6mark
  write(io,'(A,1X,I0)') 'stage7_stage6_closed_marker_status', s6mark
  write(io,'(A,1X,I0)') 'stage7_stage6_dependency_status', dep

  call init_stage7_default_config(c)
  call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  d0=merge(1,0,c%enable_stage7); d1=merge(1,0,c%enable_nonuniform_y_ibm); d2=merge(1,0,c%enable_component_specific_layout)
  d3=merge(1,0,c%enable_real_layout_interpolation); d4=merge(1,0,c%enable_real_layout_spreading); d5=merge(1,0,c%enable_real_layout_rhs_candidate)
  d6=merge(1,0,c%enable_controlled_real_layout_test); d7=merge(1,0,c%production_dns_enabled); d8=merge(1,0,c%production_two_way_enabled)
  d9=valid; d10=rej; d11=rhs_allowed; d12=c%config_status
  write(io,'(A,1X,I0)') 'stage7_default_enable_stage7', d0
  write(io,'(A,1X,I0)') 'stage7_default_nonuniform_y_ibm_enabled', d1
  write(io,'(A,1X,I0)') 'stage7_default_component_specific_layout_enabled', d2
  write(io,'(A,1X,I0)') 'stage7_default_real_layout_interpolation_enabled', d3
  write(io,'(A,1X,I0)') 'stage7_default_real_layout_spreading_enabled', d4
  write(io,'(A,1X,I0)') 'stage7_default_real_layout_rhs_candidate_enabled', d5
  write(io,'(A,1X,I0)') 'stage7_default_controlled_test_enabled', d6
  write(io,'(A,1X,I0)') 'stage7_default_production_dns_enabled', d7
  write(io,'(A,1X,I0)') 'stage7_default_production_two_way_enabled', d8
  write(io,'(A,1X,I0)') 'stage7_default_valid_flag', d9
  write(io,'(A,1X,I0)') 'stage7_default_rejected_flag', d10
  write(io,'(A,1X,I0)') 'stage7_default_rhs_allowed_flag', d11
  write(io,'(A,1X,I0)') 'stage7_default_config_status', d12

  call init_stage7_controlled_real_layout_config(c)
  call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  c0=merge(1,0,c%enable_stage7); c1=merge(1,0,c%enable_nonuniform_y_ibm); c2=merge(1,0,c%enable_component_specific_layout)
  c3=merge(1,0,c%enable_real_layout_interpolation); c4=merge(1,0,c%enable_real_layout_spreading); c5=merge(1,0,c%enable_real_layout_rhs_candidate)
  c6=merge(1,0,c%enable_controlled_real_layout_test); c7=merge(1,0,c%production_dns_enabled); c8=merge(1,0,c%production_two_way_enabled)
  c9=valid; c10=rej; c11=rhs_allowed; c12=c%config_status
  write(io,'(A,1X,I0)') 'stage7_controlled_enable_stage7', c0
  write(io,'(A,1X,I0)') 'stage7_controlled_nonuniform_y_ibm_enabled', c1
  write(io,'(A,1X,I0)') 'stage7_controlled_component_specific_layout_enabled', c2
  write(io,'(A,1X,I0)') 'stage7_controlled_real_layout_interpolation_enabled', c3
  write(io,'(A,1X,I0)') 'stage7_controlled_real_layout_spreading_enabled', c4
  write(io,'(A,1X,I0)') 'stage7_controlled_real_layout_rhs_candidate_enabled', c5
  write(io,'(A,1X,I0)') 'stage7_controlled_test_enabled', c6
  write(io,'(A,1X,I0)') 'stage7_controlled_production_dns_enabled', c7
  write(io,'(A,1X,I0)') 'stage7_controlled_production_two_way_enabled', c8
  write(io,'(A,1X,I0)') 'stage7_controlled_valid_flag', c9
  write(io,'(A,1X,I0)') 'stage7_controlled_rejected_flag', c10
  write(io,'(A,1X,I0)') 'stage7_controlled_rhs_allowed_flag', c11
  write(io,'(A,1X,I0)') 'stage7_controlled_config_status', c12

  call init_stage7_controlled_real_layout_config(c); c%production_dns_enabled=.true.; call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  write(io,'(A,1X,I0)') 'stage7_invalid_production_dns_rejected_flag', rej
  write(io,'(A,1X,I0)') 'stage7_invalid_production_dns_valid_flag', valid
  write(io,'(A,1X,I0)') 'stage7_invalid_production_dns_rhs_allowed_flag', rhs_allowed

  call init_stage7_controlled_real_layout_config(c); c%production_two_way_enabled=.true.; call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  write(io,'(A,1X,I0)') 'stage7_invalid_production_twoway_rejected_flag', rej
  write(io,'(A,1X,I0)') 'stage7_invalid_production_twoway_valid_flag', valid
  write(io,'(A,1X,I0)') 'stage7_invalid_production_twoway_rhs_allowed_flag', rhs_allowed

  call init_stage7_controlled_real_layout_config(c); c%enable_real_layout_spreading=.false.; call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  f_rej=rej; f_valid=valid
  write(io,'(A,1X,I0)') 'stage7_invalid_interp_without_spread_rejected_flag', f_rej
  write(io,'(A,1X,I0)') 'stage7_invalid_interp_without_spread_valid_flag', f_valid

  call init_stage7_controlled_real_layout_config(c); c%enable_real_layout_interpolation=.false.; call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  g_rej=rej; g_valid=valid
  write(io,'(A,1X,I0)') 'stage7_invalid_spread_without_interp_rejected_flag', g_rej
  write(io,'(A,1X,I0)') 'stage7_invalid_spread_without_interp_valid_flag', g_valid

  call init_stage7_controlled_real_layout_config(c); c%enable_controlled_real_layout_test=.false.; call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  h_rej=rej; h_rhs=rhs_allowed
  write(io,'(A,1X,I0)') 'stage7_invalid_capability_without_controlled_test_rejected_flag', h_rej
  write(io,'(A,1X,I0)') 'stage7_invalid_capability_without_controlled_test_rhs_allowed_flag', h_rhs

  call init_stage7_controlled_real_layout_config(c); c%rho_fluid=0._mytype; call validate_stage7_config(c,valid,rej,rhs_allowed,ctest,pdns,ptwo)
  i_rej=rej; i_valid=valid; i_rhs=rhs_allowed
  write(io,'(A,1X,I0)') 'stage7_invalid_rho_rejected_flag', i_rej
  write(io,'(A,1X,I0)') 'stage7_invalid_rho_valid_flag', i_valid
  write(io,'(A,1X,I0)') 'stage7_invalid_rho_rhs_allowed_flag', i_rhs

  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k
    rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i
    rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j
  end do; end do; end do
  call init_stage7_default_config(c); call stage7_config_noop_rhs_guard(c,rhsx,rhsy,rhsz,chg,modf)
  write(io,'(A,1X,ES24.16)') 'stage7_config_noop_rhs_change_max', chg
  write(io,'(A,1X,I0)') 'stage7_config_noop_rhs_modified_flag', modf

  call stage7_config_pressure_status(ppm,pm,rpc,pdc,fuc,fac)
  write(io,'(A,1X,I0)') 'stage7_config_pressure_poisson_modified_flag', ppm
  write(io,'(A,1X,I0)') 'stage7_config_projection_modified_flag', pm
  write(io,'(A,1X,I0)') 'stage7_config_real_projection_called_flag', rpc
  write(io,'(A,1X,I0)') 'stage7_config_production_dns_called_flag', pdc
  write(io,'(A,1X,I0)') 'stage7_config_fluid_update_called_flag', fuc
  write(io,'(A,1X,I0)') 'stage7_config_fibre_advance_called_flag', fac

  final_status=merge(1,0,s6mark==1 .and. dep==1 .and. d0==0 .and. d1==0 .and. d2==0 .and. d3==0 .and. d4==0 .and. d5==0 .and. d6==0 .and. d7==0 .and. d8==0 .and. d9==1 .and. d10==0 .and. d11==0 .and. d12==1 .and. c0==1 .and. c1==1 .and. c2==1 .and. c3==1 .and. c4==1 .and. c5==1 .and. c6==1 .and. c7==0 .and. c8==0 .and. c9==1 .and. c10==0 .and. c11==1 .and. c12==1 .and. f_rej==1 .and. f_valid==0 .and. g_rej==1 .and. g_valid==0 .and. h_rej==1 .and. h_rhs==0 .and. i_rej==1 .and. i_valid==0 .and. i_rhs==0 .and. chg<=1e-14_mytype .and. modf==0 .and. ppm==0 .and. pm==0 .and. rpc==0 .and. pdc==0 .and. fuc==0 .and. fac==0)
  write(io,'(A,1X,I0)') 'stage7_config_check_status', final_status
  close(io)
end program fibre_stage7_config_check
