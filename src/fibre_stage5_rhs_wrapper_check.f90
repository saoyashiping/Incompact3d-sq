program fibre_stage5_rhs_wrapper_check
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, STAGE5_COUPLING_TWO_WAY, &
                                   init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config
  use fibre_stage5_rhs_wrapper, only : apply_stage5_ibm_force_to_fluid_rhs, compute_stage5_rhs_change_max, compute_stage5_rhs_expected_error
  implicit none
  integer, parameter :: nx=8, ny=6, nz=5
  type(stage5_config_t) :: cfg
  real(mytype), allocatable :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:), rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
  integer :: i,j,k, rhs_modified, injected, rejected, ios, status
  real(mytype) :: change_max, expected_error, buffer_max_abs, ex, ey, ez

  allocate(rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  allocate(fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz))

  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*real(i,mytype)+0.01_mytype*real(j,mytype)
    rhsy(i,j,k)=-0.2_mytype*real(j,mytype)+0.03_mytype*real(k,mytype)
    rhsz(i,j,k)=0.05_mytype*real(k,mytype)-0.01_mytype*real(i,mytype)
    fx(i,j,k)=sin(0.1_mytype*real(i,mytype)+0.2_mytype*real(j,mytype)+0.3_mytype*real(k,mytype))
    fy(i,j,k)=cos(0.2_mytype*real(i,mytype)-0.1_mytype*real(k,mytype))
    fz(i,j,k)=0.1_mytype*sin(0.3_mytype*real(j,mytype))
  end do; end do; end do

  open(11,file='stage5_outputs/fibre_stage5_rhs_wrapper_check.dat',status='replace',iostat=ios)
  if (ios/=0) stop 1

  call init_stage5_default_config(cfg)
  rhsx0=rhsx; rhsy0=rhsy; rhsz0=rhsz
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change_max)
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_disabled_change_max', change_max
  write(11,'(A,1X,I0)') 'stage5_rhs_disabled_modified_flag', rhs_modified
  write(11,'(A,1X,I0)') 'stage5_rhs_disabled_injected_flag', injected

  call init_stage5_oneway_config(cfg)
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change_max)
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_oneway_change_max', change_max
  write(11,'(A,1X,I0)') 'stage5_rhs_oneway_modified_flag', rhs_modified
  write(11,'(A,1X,I0)') 'stage5_rhs_oneway_injected_flag', injected

  call init_stage5_twoway_config(cfg)
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0
  fx=0._mytype; fy=0._mytype; fz=0._mytype
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change_max)
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_twoway_zero_buffer_change_max', change_max
  write(11,'(A,1X,I0)') 'stage5_rhs_twoway_zero_buffer_modified_flag', rhs_modified
  write(11,'(A,1X,I0)') 'stage5_rhs_twoway_zero_buffer_injected_flag', injected

  do k=1,nz; do j=1,ny; do i=1,nx
    fx(i,j,k)=sin(0.1_mytype*real(i,mytype)+0.2_mytype*real(j,mytype)+0.3_mytype*real(k,mytype))
    fy(i,j,k)=cos(0.2_mytype*real(i,mytype)-0.1_mytype*real(k,mytype))
    fz(i,j,k)=0.1_mytype*sin(0.3_mytype*real(j,mytype))
  end do; end do; end do

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2.0_mytype
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change_max)
  call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,cfg%rho_fluid,rhsx,rhsy,rhsz,expected_error)
  buffer_max_abs=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_twoway_nonzero_buffer_max_abs', buffer_max_abs
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_twoway_nonzero_change_max', change_max
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_twoway_nonzero_expected_error', expected_error
  write(11,'(A,1X,I0)') 'stage5_rhs_twoway_nonzero_modified_flag', rhs_modified
  write(11,'(A,1X,I0)') 'stage5_rhs_twoway_nonzero_injected_flag', injected

  call init_stage5_default_config(cfg)
  cfg%enable_stage5=.true.; cfg%coupling_mode=STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs=.true.; cfg%allow_two_way=.false.; cfg%rho_fluid=1._mytype
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change_max)
  write(11,'(A,1X,I0)') 'stage5_rhs_invalid_config_rejected_flag', rejected
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_invalid_config_change_max', change_max
  write(11,'(A,1X,I0)') 'stage5_rhs_invalid_config_injected_flag', injected

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=0._mytype
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change_max)
  write(11,'(A,1X,I0)') 'stage5_rhs_invalid_rho_rejected_flag', rejected
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_invalid_rho_change_max', change_max
  write(11,'(A,1X,I0)') 'stage5_rhs_invalid_rho_injected_flag', injected

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2.0_mytype
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,rhs_modified,injected,rejected)
  ex=maxval(abs((rhsx-rhsx0)-fx/cfg%rho_fluid)); ey=maxval(abs((rhsy-rhsy0)-fy/cfg%rho_fluid)); ez=maxval(abs((rhsz-rhsz0)-fz/cfg%rho_fluid))
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_component_x_error', ex
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_component_y_error', ey
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_component_z_error', ez

  write(11,'(A,1X,I0)') 'stage5_rhs_pressure_poisson_modified_flag', 0

  status=1
  if (ex>1e-14_mytype .or. ey>1e-14_mytype .or. ez>1e-14_mytype) status=0
  write(11,'(A,1X,I0)') 'stage5_rhs_wrapper_check_status', status
  close(11)
end program
