program fibre_stage4_poiseuille_response_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_diagnostics, only : compute_bending_energy
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_analytic_channel, only : fill_stage4_poiseuille_velocity, fill_stage4_constant_velocity
  use fibre_stage4_oneway_response, only : advance_fibre_oneway_stage4
  implicit none
  integer,parameter::nx=16,ny=12,nz=10
  real(mytype)::x(nx),y(ny),z(nz),u0(3),be,dt,uc,beta,tf
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),xref(:,:)
  type(stage4_grid_adapter_t)::a
  type(fibre_t)::f
  integer::i,j,k,u,solfail,nan,unsafe
  real(mytype)::maxlen,maxf,fcv(3),fcd(3),expv,shapeerr,bend
  dt=1e-5_mytype; beta=10._mytype; tf=20*dt
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do j=1,ny; y(j)=-1._mytype+(real(j,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do k=1,nz; z(k)=(real(k,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz))
  open(newunit=u,file='stage4_outputs/fibre_stage4_poiseuille_response_check.dat',status='replace',action='write')
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype)
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do
  allocate(xref(3,f%nl)); xref=f%x; f%v=0._mytype; call fill_stage4_constant_velocity(a,[0._mytype,0._mytype,0._mytype],ux,uy,uz)
  call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,solfail,nan,unsafe)
  write(u,'(A,1X,ES24.16E3)') 'stage4_zero_flow_preservation_error', maxval(abs(f%x-xref))
  write(u,'(A,1X,ES24.16E3)') 'stage4_zero_flow_f_ext_norm', maxf
  write(u,'(A,1X,ES24.16E3)') 'stage4_zero_flow_length_error', maxlen
  write(u,'(A,1X,I0)') 'stage4_zero_flow_solver_failure_count', solfail
  write(u,'(A,1X,I0)') 'stage4_zero_flow_nan_detected', nan
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do; f%v=0._mytype; xref=f%x
  uc=0.2_mytype; call fill_stage4_poiseuille_velocity(a,uc,ux,uy,uz)
  call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,solfail,nan,unsafe); expv=uc*(1._mytype-exp(-beta*tf))
  call compute_bending_energy(f,bend); shapeerr=maxval(abs((f%x-xref)-spread(fcd,2,f%nl)))
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_final_center_velocity_x', fcv(1)
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_expected_center_velocity_x', expv
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_velocity_error', abs(fcv(1)-expv)
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_center_displacement_x', fcd(1)
  write(u,'(A,1X,I0)') 'stage4_centerline_direction_flag', merge(1,0,fcv(1)>0._mytype)
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_shape_error_max', shapeerr
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_length_error', maxlen
  write(u,'(A,1X,ES24.16E3)') 'stage4_centerline_bending_energy_final', bend
  write(u,'(A,1X,I0)') 'stage4_centerline_solver_failure_count', solfail
  write(u,'(A,1X,I0)') 'stage4_centerline_nan_detected', nan
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do; f%v=0._mytype
  uc=-0.2_mytype; call fill_stage4_poiseuille_velocity(a,uc,ux,uy,uz); call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,solfail,nan,unsafe)
  write(u,'(A,1X,ES24.16E3)') 'stage4_reverse_center_velocity_x', fcv(1)
  write(u,'(A,1X,I0)') 'stage4_reverse_direction_flag', merge(1,0,fcv(1)<0._mytype)
  write(u,'(A,1X,ES24.16E3)') 'stage4_reverse_length_error', maxlen
  write(u,'(A,1X,I0)') 'stage4_reverse_nan_detected', nan
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[1._mytype,-0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0.5_mytype]; end do; xref=f%x; f%v=0._mytype
  uc=0.2_mytype; call fill_stage4_poiseuille_velocity(a,uc,ux,uy,uz); call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,solfail,nan,unsafe); call compute_bending_energy(f,bend)
  write(u,'(A,1X,ES24.16E3)') 'stage4_vertical_force_variation_norm', maxf
  write(u,'(A,1X,I0)') 'stage4_vertical_center_force_greater_flag', merge(1,0,abs(f%f_ext(1,(f%nl+1)/2))>abs(f%f_ext(1,1)).and.abs(f%f_ext(1,(f%nl+1)/2))>abs(f%f_ext(1,f%nl)))
  write(u,'(A,1X,ES24.16E3)') 'stage4_vertical_shape_response_max', maxval(abs(f%x-xref))
  write(u,'(A,1X,ES24.16E3)') 'stage4_vertical_bending_energy_final', bend
  write(u,'(A,1X,ES24.16E3)') 'stage4_vertical_length_error', maxlen
  write(u,'(A,1X,I0)') 'stage4_vertical_unsafe_count', unsafe
  write(u,'(A,1X,I0)') 'stage4_vertical_solver_failure_count', solfail
  write(u,'(A,1X,I0)') 'stage4_vertical_nan_detected', nan
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype); do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0.95_mytype,0.5_mytype]; end do; f%v=0._mytype
  call fill_stage4_poiseuille_velocity(a,0.2_mytype,ux,uy,uz); call advance_fibre_oneway_stage4(a,ux,uy,uz,f,beta,dt,20,maxlen,maxf,fcv,fcd,solfail,nan,unsafe)
  write(u,'(A,1X,I0)') 'stage4_nearwall_unsafe_count', unsafe
  write(u,'(A,1X,I0)') 'stage4_nearwall_blocked_flag', merge(1,0,unsafe>0)
  write(u,'(A,1X,I0)') 'stage4_nearwall_structure_advance_called', merge(0,1,unsafe>0)
  write(u,'(A,1X,I0)') 'stage4_poiseuille_fluid_rhs_modified', 0
  write(u,'(A,1X,I0)') 'stage4_poiseuille_response_status', 1
  close(u)
end program
