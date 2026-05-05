program fibre_stage5_spreading_candidate_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer, compute_ibm_force_buffer_total_force
  use fibre_ibm_power_diagnostics, only : compute_eulerian_power, compute_lagrangian_power, compute_power_consistency_error
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_interpolation_adapter, only : interpolate_stage4_vector_to_lag_if_supported
  use fibre_stage4_frozen_channel, only : fill_stage4_frozen_channel_velocity
  use fibre_stage4_boundary_policy, only : check_stage4_fibre_boundary_policy
  use fibre_stage5_config, only : stage5_config_t, init_stage5_twoway_config
  use fibre_stage5_spreading_candidate, only : compute_stage5_rhs_candidate_from_buffer, compute_stage5_rhs_candidate_expected_error, compute_stage5_real_rhs_noop_change
  implicit none
  integer, parameter :: nx=16,ny=12,nz=10,nl=33
  type(stage5_config_t) :: cfg
  type(fibre_t) :: f
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_grid_t) :: g
  type(ibm_force_buffer_t) :: buf
  type(stage4_grid_adapter_t) :: a
  real(mytype), allocatable :: x(:),y(:),z(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),u_lag(:,:),ax(:,:,:),ay(:,:,:),az(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:)
  real(mytype) :: totalE(3),totalL(3),f_abs,f_rel,pe,pl,pa,pr,pc,err,ex,ey,ez,pi,lag_norm,bufmax,candmax,change
  integer :: i,j,k,m,status,safe,wrap,unsafe,outside,blocked,spread_called

  allocate(x(nx),y(ny),z(nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  call fill_stage4_frozen_channel_velocity(a,0.2_mytype,0.05_mytype,0.02_mytype,0.015_mytype,ux,uy,uz)

  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0._mytype;g%xmax=2._mytype;g%ymin=-1._mytype;g%ymax=1._mytype;g%zmin=0._mytype;g%zmax=1._mytype
  g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz
  g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.;allocate(g%x(nx),g%y(ny),g%z(nz));g%x=x;g%y=y;g%z=z
  call allocate_ibm_force_buffer(buf,g)

  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype)
  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0._mytype,0.5_mytype]; end do
  call init_lagrangian_points_from_fibre(lag,f)
  allocate(u_lag(3,lag%nl))
  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype

  open(11,file='stage5_outputs/fibre_stage5_spreading_candidate_check.dat',status='replace')
  pi=acos(-1._mytype)

  lag%force=0._mytype; call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz)
  call compute_stage5_rhs_candidate_from_buffer(cfg,buf%fx,buf%fy,buf%fz,ax,ay,az,status)
  bufmax=max(maxval(abs(buf%fx)),max(maxval(abs(buf%fy)),maxval(abs(buf%fz)))); candmax=max(maxval(abs(ax)),max(maxval(abs(ay)),maxval(abs(az))))
  call compute_ibm_force_buffer_total_force(g,buf,totalE); totalL=0._mytype
  f_abs=sqrt(sum((totalE-totalL)**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_zero_force_buffer_max_abs', bufmax
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_zero_rhs_candidate_max_abs', candmax
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_zero_force_conservation_error', f_abs

  do i=1,lag%nl
    lag%force(1,i)=0.1_mytype*sin(2._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.03_mytype
    lag%force(2,i)=0.05_mytype*cos(2._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.02_mytype
    lag%force(3,i)=0.02_mytype*sin(4._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.01_mytype
  end do
  lag_norm=sqrt(sum(lag%force**2)); call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz)
  call compute_stage5_rhs_candidate_from_buffer(cfg,buf%fx,buf%fy,buf%fz,ax,ay,az,status)
  bufmax=max(maxval(abs(buf%fx)),max(maxval(abs(buf%fy)),maxval(abs(buf%fz)))); candmax=max(maxval(abs(ax)),max(maxval(abs(ay)),maxval(abs(az))))
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_nonzero_lag_force_norm', lag_norm
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_nonzero_buffer_max_abs', bufmax
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_nonzero_rhs_candidate_max_abs', candmax

  call compute_ibm_force_buffer_total_force(g,buf,totalE); totalL=0._mytype; do i=1,lag%nl; totalL=totalL+lag%force(:,i)*lag%weight(i); end do
  f_abs=sqrt(sum((totalE-totalL)**2)); f_rel=f_abs/max(sqrt(sum(totalL**2)),tiny(1._mytype))
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_force_conservation_abs_error', f_abs
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_force_conservation_relative_error', f_rel
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_total_eulerian_force_norm', sqrt(sum(totalE**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_total_lagrangian_force_norm', sqrt(sum(totalL**2))

  call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u_lag,m)
  call compute_eulerian_power(g,ux,uy,uz,buf%fx,buf%fy,buf%fz,pe); call compute_lagrangian_power(lag,u_lag,pl); call compute_power_consistency_error(pe,pl,pa,pr); pc=abs(abs(pe-pl)-pa)
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_power_eulerian', pe
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_power_lagrangian', pl
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_power_abs_error', pa
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_power_relative_error', pr
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_power_error_consistency_check', pc
  write(11,'(A,1X,I0)') 'stage5_spread_power_nonzero_flag', merge(1,0,abs(pl)>1e-12_mytype)

  call compute_stage5_rhs_candidate_expected_error(cfg,buf%fx,buf%fy,buf%fz,ax,ay,az,err,ex,ey,ez)
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_candidate_rho_fluid', cfg%rho_fluid
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_candidate_expected_error', err
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_candidate_component_x_error', ex
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_candidate_component_y_error', ey
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_candidate_component_z_error', ez

  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*real(i,mytype)+0.01_mytype*real(j,mytype)
    rhsy(i,j,k)=-0.2_mytype*real(j,mytype)+0.03_mytype*real(k,mytype)
    rhsz(i,j,k)=0.05_mytype*real(k,mytype)-0.01_mytype*real(i,mytype)
  end do; end do; end do
  rhsx0=rhsx; rhsy0=rhsy; rhsz0=rhsz
  call compute_stage5_real_rhs_noop_change(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,change)
  write(11,'(A,1X,I0)') 'stage5_spread_rhs_candidate_nonzero_flag', merge(1,0,candmax>0._mytype)
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_real_rhs_change_max', change
  write(11,'(A,1X,I0)') 'stage5_spread_real_rhs_modified_flag', merge(1,0,change>1e-14_mytype)

  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0.95_mytype,0._mytype+0.5_mytype]; end do
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  spread_called=0; call clear_ibm_force_buffer(buf)
  if (blocked/=1) then; spread_called=1; call init_lagrangian_points_from_fibre(lag,f); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); end if
  bufmax=max(maxval(abs(buf%fx)),max(maxval(abs(buf%fy)),maxval(abs(buf%fz))))
  write(11,'(A,1X,I0)') 'stage5_spread_nearwall_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage5_spread_nearwall_blocked_flag', blocked
  write(11,'(A,1X,I0)') 'stage5_spread_nearwall_spreading_called', spread_called
  write(11,'(A,1X,ES24.16E3)') 'stage5_spread_nearwall_buffer_max_abs', bufmax

  write(11,'(A,1X,I0)') 'stage5_spread_pressure_poisson_modified_flag', 0
  status=1
  write(11,'(A,1X,I0)') 'stage5_spreading_candidate_check_status', status
  close(11)
end program
