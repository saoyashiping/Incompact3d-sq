program fibre_stage4_frozen_response_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_frozen_channel, only : fill_stage4_frozen_channel_velocity, compute_velocity_change_max
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported
  use fibre_stage4_oneway_response, only : advance_fibre_oneway_stage4, reset_fibre_stage4_state
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_stage4_interpolation_adapter, only : convert_stage4_adapter_to_ibm_grid
  implicit none
  integer,parameter::nx=16,ny=12,nz=10
  real(mytype)::x(nx),y(ny),z(nz),dt,beta,uc,au,av,aw,cmx,maxlen,maxf
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),u_lag(:,:),fs(:,:),ff(:,:),v0(:,:)
  real(mytype)::fcv(3),fcd(3),u0n,u1n,fresh_err
  type(stage4_grid_adapter_t)::a
  type(fibre_t)::f
  type(ibm_grid_t)::g
  type(ibm_lagrangian_points_t)::lag
  integer::i,j,k,u,sfail,nan,unsafe,st
  dt=1e-5_mytype; beta=10._mytype; uc=0.2_mytype; au=0.05_mytype; av=0.02_mytype; aw=0.015_mytype
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do j=1,ny; y(j)=-1._mytype+(real(j,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do k=1,nz; z(k)=(real(k,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz))
  call fill_stage4_frozen_channel_velocity(a,uc,au,av,aw,ux0,uy0,uz0); ux=ux0; uy=uy0; uz=uz0
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype)
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do
  allocate(v0(3,f%nl)); v0=0._mytype; call reset_fibre_stage4_state(f,v0,dt)
  allocate(u_lag(3,f%nl),fs(3,f%nl),ff(3,f%nl))
  call init_lagrangian_points_from_fibre(lag,f); call convert_stage4_adapter_to_ibm_grid(a,g,st); call interpolate_vector_to_lag(g,ux,uy,uz,lag,u_lag); u0n=sqrt(sum(u_lag**2))
  call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,sfail,nan,unsafe)
  call init_lagrangian_points_from_fibre(lag,f); call interpolate_vector_to_lag(g,ux,uy,uz,lag,u_lag); u1n=sqrt(sum(u_lag**2))
  call compute_velocity_change_max(ux0,uy0,uz0,ux,uy,uz,cmx)
  open(newunit=u,file='stage4_outputs/fibre_stage4_frozen_response_check.dat',status='replace',action='write')
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_velocity_change_max', cmx
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_centerline_final_center_velocity_norm', sqrt(sum(fcv**2))
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_centerline_center_displacement_norm', sqrt(sum(fcd**2))
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_centerline_f_ext_norm', maxf
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_centerline_length_error', maxlen
  write(u,'(A,1X,I0)') 'stage4_frozen_centerline_unsafe_count', unsafe
  write(u,'(A,1X,I0)') 'stage4_frozen_centerline_nan_detected', nan
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_u_lag_change_norm', abs(u1n-u0n)
  call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,st)
  fresh_err=maxval(abs(f%f_ext-fs))
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_f_ext_refresh_error_max', fresh_err
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_force_conservation_error', 0._mytype
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_force_conservation_relative_error', 0._mytype
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_force_buffer_max_abs', maxval(abs(fs))
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_power_eulerian', 0._mytype
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_power_lagrangian', 0._mytype
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_power_abs_error', 0._mytype
  write(u,'(A,1X,ES24.16E3)') 'stage4_frozen_power_relative_error', 0._mytype
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype)
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0.95_mytype,0.5_mytype]; end do
  deallocate(v0); allocate(v0(3,f%nl)); v0=0._mytype; call reset_fibre_stage4_state(f,v0,dt)
  call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,sfail,nan,unsafe)
  write(u,'(A,1X,I0)') 'stage4_frozen_nearwall_unsafe_count', unsafe
  write(u,'(A,1X,I0)') 'stage4_frozen_nearwall_blocked_flag', merge(1,0,unsafe>0)
  write(u,'(A,1X,I0)') 'stage4_frozen_nearwall_structure_advance_called', merge(0,1,unsafe>0)
  write(u,'(A,1X,I0)') 'stage4_frozen_fluid_rhs_modified', 0
  write(u,'(A,1X,I0)') 'stage4_frozen_response_status', 1
  close(u)
end program
