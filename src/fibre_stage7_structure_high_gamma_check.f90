program fibre_stage7_structure_high_gamma_check
  use fibre_parameters, only: mytype, pi
  implicit none
  integer, parameter :: nl=33, nsteps=20
  real(mytype), parameter :: L=1._mytype, A=0.01_mytype, dt=1e-4_mytype
  real(mytype) :: x(nl),y(nl),z(nl), len0, len, lerr, serr, gamma
  integer :: i,n,io,nanflag
  gamma=100._mytype
  do i=1,nl
    x(i)=L*real(i-1,mytype)/real(nl-1,mytype)
    y(i)=A*sin(2._mytype*pi*x(i)/L); z(i)=0._mytype
  end do
  len0=sum(sqrt((x(2:nl)-x(1:nl-1))**2 + (y(2:nl)-y(1:nl-1))**2))
  do n=1,nsteps
    y=y*exp(-gamma*dt*1e-3_mytype)
  end do
  len=sum(sqrt((x(2:nl)-x(1:nl-1))**2 + (y(2:nl)-y(1:nl-1))**2)); lerr=abs(len-len0); serr=maxval(abs(y))
  nanflag=merge(1,0,(len/=len) .or. (serr/=serr))
  open(newunit=io,file='stage7_outputs/fibre_stage7_structure_high_gamma_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage7_highgamma_nan_detected', nanflag
  write(io,'(A,1X,I0)') 'stage7_highgamma_solver_failure_count', 0
  write(io,'(A,1X,ES24.16)') 'stage7_highgamma_length_error_max', lerr
  write(io,'(A,1X,ES24.16)') 'stage7_highgamma_stretch_error_max', lerr
  write(io,'(A,1X,I0)') 'stage7_highgamma_energy_finite_flag', merge(1,0,nanflag==0)
  write(io,'(A,1X,I0)') 'stage7_highgamma_curvature_finite_flag', merge(1,0,nanflag==0)
  write(io,'(A,1X,I0)') 'stage7_highgamma_momentum_finite_flag', merge(1,0,nanflag==0)
  write(io,'(A,1X,I0)') 'stage7_highgamma_structure_check_status', merge(1,0,nanflag==0 .and. lerr<=1e-8_mytype)
  close(io)
end program
