program fibre_stage5_rhs_audit_check
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t, init_stage5_twoway_config, validate_stage5_config
  use fibre_stage5_rhs_audit, only : stage5_rhs_audit_t, STAGE5_RHS_FORM_ACCELERATION, STAGE5_RHS_INSERT_BEFORE_PROJECTION, &
                                     STAGE5_RHS_INSERT_AFTER_PROJECTION, init_stage5_rhs_audit_policy, stage5_rhs_audit_noop, &
                                     write_stage5_rhs_audit_report
  implicit none

  type(stage5_config_t) :: cfg
  type(stage5_rhs_audit_t) :: audit
  integer :: valid, rejected, two_way_enabled, rhs_allowed, rhs_modified, report_status
  integer :: report_exists, check_status
  real(mytype) :: rhs_change_max
  real(mytype), allocatable :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
  integer :: i, j, k, ios
  logical :: exists_flag

  allocate(rhsx(8,6,5), rhsy(8,6,5), rhsz(8,6,5))
  do k=1,5
    do j=1,6
      do i=1,8
        rhsx(i,j,k)=0.1_mytype*real(i,mytype)+0.01_mytype*real(j,mytype)
        rhsy(i,j,k)=-0.2_mytype*real(j,mytype)+0.03_mytype*real(k,mytype)
        rhsz(i,j,k)=0.05_mytype*real(k,mytype)-0.01_mytype*real(i,mytype)
      end do
    end do
  end do

  call init_stage5_twoway_config(cfg)
  call init_stage5_rhs_audit_policy(audit, cfg)

  open(11,file='stage5_outputs/fibre_stage5_rhs_audit_check.dat',status='replace',iostat=ios)
  if (ios /= 0) stop 1

  write(11,'(A,1X,I0)') 'stage5_rhs_convention_status', audit%audit_status
  write(11,'(A,1X,I0)') 'stage5_rhs_uses_acceleration_form', merge(1,0,audit%rhs_form==STAGE5_RHS_FORM_ACCELERATION)
  write(11,'(A,1X,I0)') 'stage5_rhs_requires_divide_by_rho', audit%requires_divide_by_rho
  write(11,'(A,1X,I0)') 'stage5_rhs_rho_positive_flag', merge(1,0,cfg%rho_fluid>0._mytype)

  write(11,'(A,1X,I0)') 'stage5_rhs_location_identified', merge(1,0,audit%insertion_order==STAGE5_RHS_INSERT_BEFORE_PROJECTION)
  write(11,'(A,1X,I0)') 'stage5_rhs_insert_before_projection_flag', merge(1,0,audit%insertion_order==STAGE5_RHS_INSERT_BEFORE_PROJECTION)
  write(11,'(A,1X,I0)') 'stage5_rhs_insert_after_projection_flag', merge(1,0,audit%insertion_order==STAGE5_RHS_INSERT_AFTER_PROJECTION)
  write(11,'(A,1X,I0)') 'stage5_rhs_pressure_poisson_modified_flag', audit%pressure_poisson_modified

  write(11,'(A,1X,I0)') 'stage5_rhs_rk_substep_policy_status', merge(1,0,audit%current_substep_velocity_required==1 .and. audit%stale_velocity_forbidden==1 .and. audit%substep_force_recompute_required==1)
  write(11,'(A,1X,I0)') 'stage5_rhs_current_substep_velocity_required', audit%current_substep_velocity_required
  write(11,'(A,1X,I0)') 'stage5_rhs_stale_velocity_forbidden_flag', audit%stale_velocity_forbidden
  write(11,'(A,1X,I0)') 'stage5_rhs_substep_force_recompute_required', audit%substep_force_recompute_required

  call stage5_rhs_audit_noop(cfg, rhsx, rhsy, rhsz, rhs_change_max, rhs_modified)
  write(11,'(A,1X,ES24.16E3)') 'stage5_rhs_audit_noop_rhs_change_max', rhs_change_max
  write(11,'(A,1X,I0)') 'stage5_rhs_audit_rhs_modified_flag', rhs_modified

  call validate_stage5_config(cfg, valid, rejected, two_way_enabled, rhs_allowed)
  call stage5_rhs_audit_noop(cfg, rhsx, rhsy, rhsz, rhs_change_max, rhs_modified)
  write(11,'(A,1X,I0)') 'stage5_rhs_twoway_config_rhs_allowed_flag', rhs_allowed
  write(11,'(A,1X,I0)') 'stage5_rhs_twoway_config_rhs_modified_in_5_1_flag', rhs_modified

  call write_stage5_rhs_audit_report('stage5_outputs/stage5_rhs_audit_report.md', audit, cfg, report_status)
  inquire(file='stage5_outputs/stage5_rhs_audit_report.md', exist=exists_flag)
  report_exists = merge(1,0,exists_flag)
  write(11,'(A,1X,I0)') 'stage5_rhs_audit_report_exists', report_exists
  write(11,'(A,1X,I0)') 'stage5_rhs_audit_status', report_status

  check_status = 1
  if (audit%audit_status /= 1) check_status = 0
  if (rhs_change_max > 1.0e-14_mytype) check_status = 0
  if (rhs_modified /= 0) check_status = 0
  if (rhs_allowed /= 1) check_status = 0
  if (report_exists /= 1 .or. report_status /= 1) check_status = 0
  write(11,'(A,1X,I0)') 'stage5_rhs_audit_check_status', check_status

  close(11)
  deallocate(rhsx, rhsy, rhsz)
end program fibre_stage5_rhs_audit_check
