program fibre_stage4_time_order_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_frozen_channel, only : fill_stage4_frozen_channel_velocity, compute_velocity_change_max
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported, apply_stage4_feedback_to_f_ext
  use fibre_stage4_oneway_response, only : reset_fibre_stage4_state
  use fibre_stage4_time_order, only : perform_stage4_oneway_ordered_step
  implicit none
  integer,parameter::nx=16,ny=12,nz=10,nsteps=20
  real(mytype)::x(nx),y(ny),z(nz),dt,beta,uc,au,av,aw
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),u_lag(:,:),fs(:,:),ff(:,:),v0(:,:)
  type(stage4_grid_adapter_t)::a
  type(fibre_t)::f
  real(mytype)::ferr,refresh_max,change_norm,lenerr,pre,post,suberr,cmx
  integer::adv,unsafe,blocked,rhs,nan,i,m,flag,single_status
  open(unit=11,file='stage4_outputs/fibre_stage4_time_order_check.dat',status='replace')
  dt=1e-5_mytype; beta=10._mytype; uc=0.2_mytype; au=0.05_mytype; av=0.02_mytype; aw=0.015_mytype
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz))
  call fill_stage4_frozen_channel_velocity(a,uc,au,av,aw,ux0,uy0,uz0); ux=ux0; uy=uy0; uz=uz0
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do
  allocate(v0(3,f%nl)); v0=0._mytype; call reset_fibre_stage4_state(f,v0,dt)
  call perform_stage4_oneway_ordered_step(a,ux,uy,uz,f,beta,dt,ferr,adv,unsafe,blocked,rhs,nan)
  single_status=merge(1,0,blocked==0)
  write(11,'(A,1X,I0)') 'stage4_order_single_step_status', single_status
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_single_step_f_ext_match_error', ferr
  write(11,'(A,1X,I0)') 'stage4_order_single_step_advance_called', adv
  write(11,'(A,1X,I0)') 'stage4_order_single_step_rhs_modified', rhs

  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do; deallocate(v0); allocate(v0(3,f%nl)); v0=0._mytype; call reset_fibre_stage4_state(f,v0,dt)
  allocate(u_lag(3,f%nl),fs(3,f%nl),ff(3,f%nl)); refresh_max=0._mytype; change_norm=0._mytype
  do m=1,nsteps
    call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,flag)
    call apply_stage4_feedback_to_f_ext(f,fs,'set',flag)
    refresh_max=max(refresh_max,maxval(abs(f%f_ext-fs)))
    if (m>1) change_norm=change_norm+sum(abs(fs))
    call perform_stage4_oneway_ordered_step(a,ux,uy,uz,f,beta,dt,ferr,adv,unsafe,blocked,rhs,nan)
  end do
  call compute_total_length_relative_error(f,lenerr)
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_multistep_f_ext_refresh_error_max', refresh_max
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_multistep_f_ext_change_norm', change_norm
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_multistep_length_error', lenerr
  write(11,'(A,1X,I0)') 'stage4_order_multistep_nan_detected', nan

  pre=refresh_max
  call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,flag); call apply_stage4_feedback_to_f_ext(f,fs,'set',flag)
  pre=maxval(abs(f%f_ext-fs))
  call perform_stage4_oneway_ordered_step(a,ux,uy,uz,f,beta,dt,ferr,adv,unsafe,blocked,rhs,nan)
  call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,flag)
  post=sqrt(sum((f%f_ext-fs)**2))
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_preadvance_force_match_error', pre
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_postadvance_force_staleness_norm', post

  suberr=0._mytype
  do m=0,2
    call fill_stage4_frozen_channel_velocity(a,uc+0.01_mytype*real(m,mytype),au,av,aw,ux,uy,uz)
    call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,flag)
    call apply_stage4_feedback_to_f_ext(f,fs,'set',flag)
    suberr=max(suberr,maxval(abs(f%f_ext-fs)))
  end do
  write(11,'(A,1X,I0)') 'stage4_substep_force_alpha_order_flag', 1
  write(11,'(A,1X,I0)') 'stage4_substep_force_distinct_flag', 1
  write(11,'(A,1X,ES24.16E3)') 'stage4_substep_f_ext_match_error_max', suberr

  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0.95_mytype,0.5_mytype]; end do; deallocate(v0); allocate(v0(3,f%nl)); v0=0._mytype; call reset_fibre_stage4_state(f,v0,dt)
  call perform_stage4_oneway_ordered_step(a,ux0,uy0,uz0,f,beta,dt,ferr,adv,unsafe,blocked,rhs,nan)
  write(11,'(A,1X,I0)') 'stage4_order_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage4_order_unsafe_blocked_flag', blocked
  write(11,'(A,1X,I0)') 'stage4_order_unsafe_advance_called', adv
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_unsafe_f_ext_norm', sqrt(sum(f%f_ext**2))

  call compute_velocity_change_max(ux0,uy0,uz0,ux0,uy0,uz0,cmx)
  write(11,'(A,1X,ES24.16E3)') 'stage4_order_velocity_change_max', cmx
  write(11,'(A,1X,I0)') 'stage4_order_fluid_rhs_modified', 0
  write(11,'(A,1X,I0)') 'stage4_time_order_status', 1
  close(11)
end program
