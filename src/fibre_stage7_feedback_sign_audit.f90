program fibre_stage7_feedback_sign_audit
  use fibre_parameters, only: mytype
  implicit none
  real(mytype) :: beta, slip(3), fs(3), ff(3), v(3), u(3), ptotal, pexp
  real(mytype) :: zforce, are, terr
  integer :: io
  beta=2._mytype
  slip=(/0._mytype,0._mytype,0._mytype/); fs=beta*slip; ff=-fs; zforce=sqrt(sum(fs*fs))
  slip=(/0.2_mytype,-0.1_mytype,0.05_mytype/); fs=beta*slip; ff=-fs
  v=(/0.1_mytype,0.2_mytype,-0.3_mytype/); u=v+slip
  ptotal=sum(fs*v)+sum(ff*u); pexp=-beta*sum(slip*slip); terr=abs(ptotal-pexp); are=maxval(abs(fs+ff))
  open(newunit=io,file='stage7_outputs/fibre_stage7_feedback_sign_audit.dat',status='replace',action='write')
  write(io,'(A,1X,ES24.16)') 'stage7_feedback_zero_slip_force_norm', zforce
  write(io,'(A,1X,ES24.16)') 'stage7_feedback_action_reaction_error', are
  write(io,'(A,1X,I0)') 'stage7_feedback_structure_force_slip_dot_positive_flag', merge(1,0,sum(fs*slip)>0._mytype)
  write(io,'(A,1X,I0)') 'stage7_feedback_fluid_force_slip_dot_negative_flag', merge(1,0,sum(ff*slip)<0._mytype)
  write(io,'(A,1X,I0)') 'stage7_feedback_total_power_dissipative_flag', merge(1,0,ptotal<=1e-14_mytype)
  write(io,'(A,1X,ES24.16)') 'stage7_feedback_total_power_error', terr
  write(io,'(A,1X,I0)') 'stage7_feedback_sign_audit_status', merge(1,0,zforce<=1e-14_mytype .and. are<=1e-14_mytype .and. terr<=1e-14_mytype)
  close(io)
end program
