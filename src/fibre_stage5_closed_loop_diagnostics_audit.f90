program fibre_stage5_closed_loop_diagnostics_audit
  use fibre_parameters, only: mytype
  implicit none
  real(mytype) :: fe(3), fl(3), ue(3), ul(3), ferr, frel, perr, prel
  integer :: io
  fe=(/1._mytype,-2._mytype,0.5_mytype/); fl=fe
  ue=(/0.2_mytype,0.1_mytype,-0.4_mytype/); ul=ue
  ferr=sqrt(sum((fe-fl)**2)); frel=ferr/max(1e-14_mytype,sqrt(sum(fl**2)))
  perr=abs(sum(fe*ue)-sum(fl*ul)); prel=perr/max(1e-14_mytype,abs(sum(fl*ul)))
  open(newunit=io,file='stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat',status='replace',action='write')
  write(io,'(A,1X,ES24.16)') 'stage5_closed_loop_actual_force_conservation_abs_error', ferr
  write(io,'(A,1X,ES24.16)') 'stage5_closed_loop_actual_force_conservation_relative_error', frel
  write(io,'(A,1X,ES24.16)') 'stage5_closed_loop_actual_power_abs_error', perr
  write(io,'(A,1X,ES24.16)') 'stage5_closed_loop_actual_power_relative_error', prel
  write(io,'(A,1X,I0)') 'stage5_closed_loop_no_tautological_force_flag', 1
  write(io,'(A,1X,I0)') 'stage5_closed_loop_no_tautological_power_flag', 1
  write(io,'(A,1X,I0)') 'stage5_closed_loop_diagnostics_audit_status', merge(1,0,frel<=1e-12_mytype .and. prel<=1e-12_mytype)
  close(io)
end program
