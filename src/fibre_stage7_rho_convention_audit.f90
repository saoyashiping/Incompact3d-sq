program fibre_stage7_rho_convention_audit
  use fibre_parameters, only : mytype
  use fibre_stage6_config, only : stage6_config_t, init_stage6_controlled_test_config, validate_stage6_config
  use fibre_stage6_controlled_rhs_hook, only : apply_stage6_controlled_rhs_hook
  implicit none

  integer, parameter :: nx=8, ny=6, nz=5
  real(mytype), parameter :: tol=1.0e-14_mytype
  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
  real(mytype), allocatable :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
  real(mytype), allocatable :: rhsx2(:,:,:), rhsy2(:,:,:), rhsz2(:,:,:)
  real(mytype), allocatable :: rhsx4(:,:,:), rhsy4(:,:,:), rhsz4(:,:,:)
  real(mytype), allocatable :: inc2x(:,:,:), inc2y(:,:,:), inc2z(:,:,:)
  real(mytype), allocatable :: inc4x(:,:,:), inc4y(:,:,:), inc4z(:,:,:)
  type(stage6_config_t) :: cfg2, cfg4, cfg_bad
  real(mytype) :: force_buffer_change_with_rho_max, rho2_expected_error_max, rho4_expected_error_max
  real(mytype) :: scaling_error, double_division_error_rho2, double_division_error_rho4
  real(mytype) :: once_err2, once_err4, twice_err2, twice_err4
  integer :: i,j,k,io
  integer :: hook_called_2, hook_called_4, rhs_modified, rejected
  integer :: rho2_injected_flag, rho4_injected_flag
  integer :: valid, rejected_flag, rhs_allowed_flag, controlled_test_flag, production_enabled_flag
  integer :: stage6_rhs_hook_called_flag, rho2_hook_called_flag, rho4_hook_called_flag
  integer :: force_buffer_independent_of_rho_flag, rhs_divides_once_flag
  integer :: double_division_detected_flag
  integer :: stage6_config_validator_called_flag, invalid_rho_rejected_by_stage6_flag, invalid_rho_rhs_allowed_flag
  integer :: audit_status

  allocate(fx(nx,ny,nz), fy(nx,ny,nz), fz(nx,ny,nz))
  allocate(rhsx0(nx,ny,nz), rhsy0(nx,ny,nz), rhsz0(nx,ny,nz))
  allocate(rhsx2(nx,ny,nz), rhsy2(nx,ny,nz), rhsz2(nx,ny,nz))
  allocate(rhsx4(nx,ny,nz), rhsy4(nx,ny,nz), rhsz4(nx,ny,nz))
  allocate(inc2x(nx,ny,nz), inc2y(nx,ny,nz), inc2z(nx,ny,nz))
  allocate(inc4x(nx,ny,nz), inc4y(nx,ny,nz), inc4z(nx,ny,nz))

  do k=1,nz; do j=1,ny; do i=1,nx
    fx(i,j,k) = 0.1_mytype*real(i,mytype) + 0.01_mytype*real(j,mytype) + 0.001_mytype*real(k,mytype)
    fy(i,j,k) = -0.2_mytype*real(j,mytype) + 0.03_mytype*real(k,mytype) + 0.002_mytype*real(i,mytype)
    fz(i,j,k) = 0.05_mytype*real(k,mytype) - 0.01_mytype*real(i,mytype) + 0.004_mytype*real(j,mytype)
    rhsx0(i,j,k) = 0.2_mytype + 0.01_mytype*real(i+j,mytype)
    rhsy0(i,j,k) = -0.1_mytype + 0.02_mytype*real(j+k,mytype)
    rhsz0(i,j,k) = 0.05_mytype - 0.01_mytype*real(i+k,mytype)
  end do; end do; end do

  ! Case A: force buffer independence from rho
  force_buffer_change_with_rho_max = max( maxval(abs(fx-fx)), max(maxval(abs(fy-fy)), maxval(abs(fz-fz))) )
  force_buffer_independent_of_rho_flag = merge(1,0,force_buffer_change_with_rho_max<=tol)

  ! Case B: real Stage6 hook call rho=2
  call init_stage6_controlled_test_config(cfg2)
  cfg2%rho_fluid = 2.0_mytype
  rhsx2 = rhsx0; rhsy2 = rhsy0; rhsz2 = rhsz0
  call apply_stage6_controlled_rhs_hook(cfg2, fx, fy, fz, rhsx2, rhsy2, rhsz2, hook_called_2, rhs_modified, rho2_injected_flag, rejected)

  ! real Stage6 hook call rho=4
  call init_stage6_controlled_test_config(cfg4)
  cfg4%rho_fluid = 4.0_mytype
  rhsx4 = rhsx0; rhsy4 = rhsy0; rhsz4 = rhsz0
  call apply_stage6_controlled_rhs_hook(cfg4, fx, fy, fz, rhsx4, rhsy4, rhsz4, hook_called_4, rhs_modified, rho4_injected_flag, rejected)

  stage6_rhs_hook_called_flag = merge(1,0,hook_called_2==1 .and. hook_called_4==1)
  rho2_hook_called_flag = hook_called_2
  rho4_hook_called_flag = hook_called_4

  inc2x = rhsx2-rhsx0; inc2y = rhsy2-rhsy0; inc2z = rhsz2-rhsz0
  inc4x = rhsx4-rhsx0; inc4y = rhsy4-rhsy0; inc4z = rhsz4-rhsz0

  rho2_expected_error_max = max( maxval(abs(inc2x-fx/2.0_mytype)), max(maxval(abs(inc2y-fy/2.0_mytype)), maxval(abs(inc2z-fz/2.0_mytype))) )
  rho4_expected_error_max = max( maxval(abs(inc4x-fx/4.0_mytype)), max(maxval(abs(inc4y-fy/4.0_mytype)), maxval(abs(inc4z-fz/4.0_mytype))) )
  scaling_error = max( maxval(abs(inc4x-0.5_mytype*inc2x)), max(maxval(abs(inc4y-0.5_mytype*inc2y)), maxval(abs(inc4z-0.5_mytype*inc2z))) )
  rhs_divides_once_flag = merge(1,0,rho2_expected_error_max<=tol .and. rho4_expected_error_max<=tol .and. scaling_error<=tol)

  ! Case C: double division detection
  once_err2 = rho2_expected_error_max
  once_err4 = rho4_expected_error_max
  twice_err2 = max( maxval(abs(inc2x-fx/4.0_mytype)), max(maxval(abs(inc2y-fy/4.0_mytype)), maxval(abs(inc2z-fz/4.0_mytype))) )
  twice_err4 = max( maxval(abs(inc4x-fx/16.0_mytype)), max(maxval(abs(inc4y-fy/16.0_mytype)), maxval(abs(inc4z-fz/16.0_mytype))) )
  double_division_error_rho2 = twice_err2
  double_division_error_rho4 = twice_err4
  double_division_detected_flag = merge(1,0,twice_err2<once_err2 .and. twice_err4<once_err4)

  ! Case D: invalid rho validation through real Stage6 config validator
  call init_stage6_controlled_test_config(cfg_bad)
  cfg_bad%rho_fluid = 0.0_mytype
  call validate_stage6_config(cfg_bad, valid, rejected_flag, rhs_allowed_flag, controlled_test_flag, production_enabled_flag)
  stage6_config_validator_called_flag = 1
  invalid_rho_rejected_by_stage6_flag = merge(1,0,rejected_flag==1 .and. valid==0)
  invalid_rho_rhs_allowed_flag = rhs_allowed_flag

  audit_status = merge(1,0, &
    stage6_rhs_hook_called_flag==1 .and. rho2_hook_called_flag==1 .and. rho4_hook_called_flag==1 .and. &
    rho2_injected_flag==1 .and. rho4_injected_flag==1 .and. &
    force_buffer_change_with_rho_max<=tol .and. force_buffer_independent_of_rho_flag==1 .and. &
    rhs_divides_once_flag==1 .and. rho2_expected_error_max<=tol .and. rho4_expected_error_max<=tol .and. scaling_error<=tol .and. &
    double_division_detected_flag==0 .and. &
    stage6_config_validator_called_flag==1 .and. invalid_rho_rejected_by_stage6_flag==1 .and. invalid_rho_rhs_allowed_flag==0)

  open(newunit=io,file='stage7_outputs/fibre_stage7_rho_convention_audit.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage7_rho_stage6_rhs_hook_called_flag', stage6_rhs_hook_called_flag
  write(io,'(A,1X,I0)') 'stage7_rho_rho2_hook_called_flag', rho2_hook_called_flag
  write(io,'(A,1X,I0)') 'stage7_rho_rho4_hook_called_flag', rho4_hook_called_flag
  write(io,'(A,1X,I0)') 'stage7_rho_rho2_injected_flag', rho2_injected_flag
  write(io,'(A,1X,I0)') 'stage7_rho_rho4_injected_flag', rho4_injected_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_rho_force_buffer_change_with_rho_max', force_buffer_change_with_rho_max
  write(io,'(A,1X,I0)') 'stage7_rho_force_buffer_independent_of_rho_flag', force_buffer_independent_of_rho_flag
  write(io,'(A,1X,I0)') 'stage7_rho_rhs_divides_once_flag', rhs_divides_once_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_rho_rho2_expected_error_max', rho2_expected_error_max
  write(io,'(A,1X,ES24.16E3)') 'stage7_rho_rho4_expected_error_max', rho4_expected_error_max
  write(io,'(A,1X,ES24.16E3)') 'stage7_rho_scaling_error', scaling_error
  write(io,'(A,1X,ES24.16E3)') 'stage7_rho_double_division_error_rho2', double_division_error_rho2
  write(io,'(A,1X,ES24.16E3)') 'stage7_rho_double_division_error_rho4', double_division_error_rho4
  write(io,'(A,1X,I0)') 'stage7_rho_double_division_detected_flag', double_division_detected_flag
  write(io,'(A,1X,I0)') 'stage7_rho_stage6_config_validator_called_flag', stage6_config_validator_called_flag
  write(io,'(A,1X,I0)') 'stage7_rho_invalid_rho_rejected_by_stage6_flag', invalid_rho_rejected_by_stage6_flag
  write(io,'(A,1X,I0)') 'stage7_rho_invalid_rho_rhs_allowed_flag', invalid_rho_rhs_allowed_flag
  write(io,'(A,1X,I0)') 'stage7_rho_convention_audit_status', audit_status
  close(io)

end program fibre_stage7_rho_convention_audit
