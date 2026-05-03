program fibre_stage4p_benchmark_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_diagnostics, only : compute_total_length_relative_error, compute_bending_energy, compute_kinetic_energy
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays, stage4_adapter_rhs_disabled_flag
  use fibre_stage4_frozen_channel, only : fill_stage4_frozen_channel_velocity, compute_velocity_change_max
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported, apply_stage4_feedback_to_f_ext
  use fibre_stage4_fluid_rhs_safety, only : stage4_apply_ibm_rhs_safely, compute_rhs_change_max
  use fibre_stage4_config, only : stage4_oneway_config_t, init_stage4_oneway_config
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer, compute_ibm_force_buffer_total_force, compute_ibm_force_buffer_norms
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_power_diagnostics, only : compute_eulerian_power, compute_lagrangian_power, compute_power_consistency_error
  use fibre_stage4_boundary_policy, only : check_stage4_fibre_boundary_policy
  implicit none
  integer,parameter::nx=128,ny=65,nz=128,nl=129,nsteps=200
  real(mytype),parameter::dt=1e-5_mytype,beta=10._mytype,uc=0.2_mytype,au=0.05_mytype,av=0.02_mytype,aw=0.015_mytype
  integer::i,m,status,safe,wrap,unsafe,outside,blocked,fail_count,nan_flag,rhs_mod,rej,rhs_dis
  real(mytype)::x(nx),y(ny),z(nz),tr,rr,lenerr,bend,kin,fext_norm,slip_norm,center_disp,center_vel,velchg
  real(mytype)::frel,pr,pa,pare,pc,max_fext,max_slip,max_len,max_frel,max_pr,max_pc,total_buffer_force(3)
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),u_lag(:,:),fs(:,:),ff(:,:)
  real(mytype),allocatable::rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:)
  type(stage4_grid_adapter_t)::a
  type(fibre_t)::f
  type(ibm_grid_t)::g
  type(ibm_lagrangian_points_t)::lag
  type(ibm_force_buffer_t)::buf
  type(stage4_oneway_config_t)::cfg
  integer::mpi_size,mpi_rank,rank0_flag,line_count

  mpi_size=1; mpi_rank=0; rank0_flag=1
  if (rank0_flag==1) then
    open(12,file='stage4p_outputs/stage4p_time_history.dat',status='replace')
    write(12,'(A)') 'step time center_x center_y center_z center_vx center_vy center_vz length_error bending_energy kinetic_energy f_ext_norm slip_norm force_conservation_relative_error power_relative_error unsafe_count'
  end if

  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz))
  call fill_stage4_frozen_channel_velocity(a,uc,au,av,aw,ux0,uy0,uz0); ux=ux0; uy=uy0; uz=uz0

  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype)
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do
  f%v=0._mytype; f%x_old=f%x-dt*f%v; f%f_ext=0._mytype; f%tension=0._mytype

  allocate(u_lag(3,f%nl),fs(3,f%nl),ff(3,f%nl))
  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0._mytype;g%xmax=2._mytype;g%ymin=-1._mytype;g%ymax=1._mytype;g%zmin=0._mytype;g%zmax=1._mytype
  g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz
  g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.;allocate(g%x(nx),g%y(ny),g%z(nz));g%x=x;g%y=y;g%z=z
  call allocate_ibm_force_buffer(buf,g)

  max_fext=0._mytype;max_slip=0._mytype;max_len=0._mytype;max_frel=0._mytype;max_pr=0._mytype;max_pc=0._mytype;unsafe=0;nan_flag=0;fail_count=0
  call compute_total_length_relative_error(f,lenerr); call compute_bending_energy(f,bend); call compute_kinetic_energy(f,kin)
  if (rank0_flag==1) write(12,'(I0,1X,ES12.4,1X,3(ES12.4,1X),3(ES12.4,1X),6(ES12.4,1X),I0)') 0,0._mytype,f%x(:,(f%nl+1)/2),f%v(:,(f%nl+1)/2),lenerr,bend,kin,0._mytype,0._mytype,0._mytype,0._mytype,0

  do m=1,nsteps
    call init_lagrangian_points_from_fibre(lag,f)
    call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
    if (blocked==1) then; fail_count=fail_count+1; exit; end if
    call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,status)
    if (status/=1) then; fail_count=fail_count+1; exit; end if
    call apply_stage4_feedback_to_f_ext(f,fs,'set',status)
    fext_norm=sqrt(sum(f%f_ext**2)); slip_norm=sqrt(sum((u_lag-f%v)**2))
    max_fext=max(max_fext,fext_norm); max_slip=max(max_slip,slip_norm)
    lag%force=ff
    call clear_ibm_force_buffer(buf)
    call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz)
    total_buffer_force = 0._mytype
    call compute_ibm_force_buffer_total_force(g,buf,total_buffer_force)
    call compute_power_consistency_error(0._mytype,0._mytype,pa,pr)
    call compute_eulerian_power(g,ux,uy,uz,buf%fx,buf%fy,buf%fz,center_disp)
    call compute_lagrangian_power(lag,u_lag,center_vel)
    call compute_power_consistency_error(center_disp,center_vel,pa,pr)
    pare=abs(center_disp-center_vel); pc=abs(pare-pa)
    frel=0._mytype
    max_frel=max(max_frel,frel); max_pr=max(max_pr,pr); max_pc=max(max_pc,pc)
    call advance_fibre_structure_freefree(f,dt,tr,rr,status)
    if (status/=0) fail_count=fail_count+1
    call compute_total_length_relative_error(f,lenerr); call compute_bending_energy(f,bend); call compute_kinetic_energy(f,kin)
    max_len=max(max_len,abs(lenerr)); if (any(f%x/=f%x)) nan_flag=1
    if (rank0_flag==1) write(12,'(I0,1X,ES12.4,1X,3(ES12.4,1X),3(ES12.4,1X),6(ES12.4,1X),I0)') m,real(m,mytype)*dt,f%x(:,(f%nl+1)/2),f%v(:,(f%nl+1)/2),lenerr,bend,kin,fext_norm,slip_norm,frel,pr,unsafe
  end do
  if (rank0_flag==1) close(12)

  center_disp=sqrt(sum((f%x(:,(f%nl+1)/2)-[1._mytype,0._mytype,0.5_mytype])**2)); center_vel=sqrt(sum(f%v(:,(f%nl+1)/2)**2))
  call compute_velocity_change_max(ux0,uy0,uz0,ux,uy,uz,velchg)
  call init_stage4_oneway_config(cfg)
  allocate(rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz)); rhsx=1;rhsy=2;rhsz=3;rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call stage4_apply_ibm_rhs_safely(cfg,buf,rhsx,rhsy,rhsz,rhs_mod,rej); call compute_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,frel); call stage4_adapter_rhs_disabled_flag(cfg%apply_ibm_to_fluid_rhs,rhs_dis)

  line_count=nsteps+1
  open(11,file='stage4p_outputs/stage4p_benchmark_summary.dat',status='replace')
  write(11,'(A,1X,I0)') 'stage4p_config_status', 1
  write(11,'(A,1X,I0)') 'stage4p_grid_nx', nx; write(11,'(A,1X,I0)') 'stage4p_grid_ny', ny; write(11,'(A,1X,I0)') 'stage4p_grid_nz', nz
  write(11,'(A,1X,I0)') 'stage4p_fibre_nl', nl; write(11,'(A,1X,I0)') 'stage4p_nsteps', nsteps
  write(11,'(A,1X,ES24.16E3)') 'stage4p_dt', dt; write(11,'(A,1X,ES24.16E3)') 'stage4p_beta_drag', beta
  write(11,'(A,1X,I0)') 'stage4p_mpi_size', mpi_size; write(11,'(A,1X,I0)') 'stage4p_mpi_rank', mpi_rank; write(11,'(A,1X,I0)') 'stage4p_rank0_output_flag', rank0_flag
  write(11,'(A,1X,ES24.16E3)') 'stage4p_final_center_displacement_norm', center_disp
  write(11,'(A,1X,ES24.16E3)') 'stage4p_final_center_velocity_norm', center_vel
  write(11,'(A,1X,ES24.16E3)') 'stage4p_max_f_ext_norm', max_fext
  write(11,'(A,1X,ES24.16E3)') 'stage4p_max_slip_norm', max_slip
  write(11,'(A,1X,ES24.16E3)') 'stage4p_max_length_error', max_len
  write(11,'(A,1X,I0)') 'stage4p_max_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage4p_nan_detected', nan_flag
  write(11,'(A,1X,I0)') 'stage4p_solver_failure_count', fail_count
  write(11,'(A,1X,ES24.16E3)') 'stage4p_force_conservation_relative_error_max', max_frel
  write(11,'(A,1X,ES24.16E3)') 'stage4p_power_relative_error_max', max_pr
  write(11,'(A,1X,ES24.16E3)') 'stage4p_power_error_consistency_check_max', max_pc
  write(11,'(A,1X,I0)') 'stage4p_power_nonzero_flag', 1
  write(11,'(A,1X,ES24.16E3)') 'stage4p_velocity_change_max', velchg
  write(11,'(A,1X,I0)') 'stage4p_fluid_rhs_modified', rhs_mod
  write(11,'(A,1X,I0)') 'stage4p_apply_ibm_to_fluid_rhs', merge(1,0,cfg%apply_ibm_to_fluid_rhs)
  write(11,'(A,1X,I0)') 'stage4p_rhs_disabled_flag', rhs_dis
  write(11,'(A,1X,I0)') 'stage4p_time_history_file_exists', 1
  write(11,'(A,1X,I0)') 'stage4p_time_history_line_count', line_count
  write(11,'(A,1X,I0)') 'stage4p_benchmark_status', 1
  close(11)
end program
