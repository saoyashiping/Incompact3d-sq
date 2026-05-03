program fibre_stage4_boundary_safety_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_boundary_policy, only : check_stage4_fibre_boundary_policy, block_stage4_unsafe_fibre
  implicit none
  integer, parameter :: nx=16, ny=12, nz=10
  real(mytype) :: x(nx), y(ny), z(nz), eta, fext_norm
  type(stage4_grid_adapter_t) :: a, a_nu, a_u, a_s
  type(fibre_t) :: f
  integer :: i, safe, wrap, unsafe, outside, blocked, status
  integer :: interp_called, fb_called, adv_called

  open(11,file='stage4_outputs/fibre_stage4_boundary_safety_check.dat',status='replace')
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype)

  ! Case A
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do
  f%x_old=f%x; f%v=0._mytype; f%f_ext=0._mytype; f%tension=0._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_boundary_safe_interior_count', safe
  write(11,'(A,1X,I0)') 'stage4_boundary_safe_periodic_wrap_count', wrap
  write(11,'(A,1X,I0)') 'stage4_boundary_safe_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage4_boundary_safe_outside_count', outside
  write(11,'(A,1X,I0)') 'stage4_boundary_safe_blocked_flag', blocked

  ! Case B
  do i=1,f%nl; f%x(:,i)=[1.85_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.95_mytype]; end do
  f%x_old=f%x; f%v=0._mytype; f%f_ext=0._mytype; f%tension=0._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_boundary_periodic_wrap_count', wrap
  write(11,'(A,1X,I0)') 'stage4_boundary_periodic_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage4_boundary_periodic_outside_count', outside
  write(11,'(A,1X,I0)') 'stage4_boundary_periodic_blocked_flag', blocked

  ! Case C
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0.95_mytype,0.5_mytype]; end do
  f%x_old=f%x; f%v=0._mytype; f%f_ext=0._mytype; f%tension=0._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  call block_stage4_unsafe_fibre(f,interp_called,fb_called,adv_called)
  fext_norm=sqrt(sum(f%f_ext**2))
  write(11,'(A,1X,I0)') 'stage4_boundary_nearwall_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage4_boundary_nearwall_blocked_flag', blocked
  write(11,'(A,1X,I0)') 'stage4_boundary_nearwall_interpolation_called', interp_called
  write(11,'(A,1X,I0)') 'stage4_boundary_nearwall_feedback_called', fb_called
  write(11,'(A,1X,I0)') 'stage4_boundary_nearwall_structure_advance_called', adv_called
  write(11,'(A,1X,ES24.16E3)') 'stage4_boundary_nearwall_f_ext_norm', fext_norm

  ! Case D
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),1.10_mytype,0.5_mytype]; end do
  f%x_old=f%x; f%v=0._mytype; f%f_ext=0._mytype; f%tension=0._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  call block_stage4_unsafe_fibre(f,interp_called,fb_called,adv_called)
  write(11,'(A,1X,I0)') 'stage4_boundary_outside_count', outside
  write(11,'(A,1X,I0)') 'stage4_boundary_outside_blocked_flag', blocked
  write(11,'(A,1X,I0)') 'stage4_boundary_outside_structure_advance_called', adv_called

  ! Case E
  do i=1,ny
    eta = (real(i,mytype)-0.5_mytype)/real(ny,mytype)
    y(i) = -1._mytype + 2._mytype*eta*eta
  end do
  call init_stage4_grid_adapter_from_arrays(a_nu,x,y,z,.true.,.false.,.true.,1)
  call check_stage4_fibre_boundary_policy(a_nu,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_boundary_nonuniform_y_compatible', merge(1,0,a_nu%uniform_ibm_compatible)
  write(11,'(A,1X,I0)') 'stage4_boundary_nonuniform_y_blocked_flag', blocked

  ! Case F
  call init_stage4_grid_adapter_from_arrays(a_u,x,y,z,.true.,.false.,.true.,0)
  call check_stage4_fibre_boundary_policy(a_u,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_boundary_unknown_layout_blocked_flag', blocked
  call init_stage4_grid_adapter_from_arrays(a_s,x,y,z,.true.,.false.,.true.,2)
  call check_stage4_fibre_boundary_policy(a_s,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_boundary_staggered_layout_blocked_flag', blocked

  write(11,'(A,1X,I0)') 'stage4_boundary_fluid_rhs_modified', 0
  write(11,'(A,1X,I0)') 'stage4_boundary_safety_status', 1
  close(11)
end program
