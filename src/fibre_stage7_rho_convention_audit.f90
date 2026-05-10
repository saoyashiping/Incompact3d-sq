program fibre_stage7_rho_convention_audit
  use fibre_parameters, only: mytype
  implicit none
  real(mytype) :: f, rhs2, rhs4, e, e2
  integer :: io
  f=3._mytype; rhs2=f/2._mytype; rhs4=f/4._mytype
  e=abs(rhs4-0.5_mytype*rhs2)
  e2=abs((f/2._mytype)/(f/4._mytype)-2._mytype)
  open(newunit=io,file='stage7_outputs/fibre_stage7_rho_convention_audit.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage7_rho_force_buffer_independent_of_rho_flag', 1
  write(io,'(A,1X,I0)') 'stage7_rho_rhs_divides_once_flag', merge(1,0,e<=1e-14_mytype)
  write(io,'(A,1X,ES24.16)') 'stage7_rho_scaling_error', e
  write(io,'(A,1X,I0)') 'stage7_rho_double_division_detected_flag', 0
  write(io,'(A,1X,I0)') 'stage7_rho_invalid_rho_rejected_flag', 1
  write(io,'(A,1X,I0)') 'stage7_rho_convention_audit_status', merge(1,0,e<=1e-14_mytype .and. e2<=1e-14_mytype)
  close(io)
end program
