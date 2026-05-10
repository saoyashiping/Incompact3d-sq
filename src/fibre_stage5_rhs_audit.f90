module fibre_stage5_rhs_audit
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t
  implicit none

  integer, parameter :: STAGE5_RHS_FORM_ACCELERATION = 1
  integer, parameter :: STAGE5_RHS_FORM_FORCE_DENSITY = 2

  integer, parameter :: STAGE5_RHS_INSERT_BEFORE_PROJECTION = 1
  integer, parameter :: STAGE5_RHS_INSERT_AFTER_PROJECTION = 2

  type stage5_rhs_audit_t
    integer :: rhs_form
    integer :: insertion_order
    integer :: requires_divide_by_rho
    integer :: pressure_poisson_modified
    integer :: current_substep_velocity_required
    integer :: stale_velocity_forbidden
    integer :: substep_force_recompute_required
    integer :: audit_status
  end type stage5_rhs_audit_t

contains

  subroutine init_stage5_rhs_audit_policy(audit, config)
    type(stage5_rhs_audit_t), intent(out) :: audit
    type(stage5_config_t), intent(in) :: config

    audit%rhs_form = STAGE5_RHS_FORM_ACCELERATION
    audit%insertion_order = STAGE5_RHS_INSERT_BEFORE_PROJECTION
    audit%requires_divide_by_rho = 1
    audit%pressure_poisson_modified = 0
    audit%current_substep_velocity_required = 1
    audit%stale_velocity_forbidden = 1
    audit%substep_force_recompute_required = 1
    audit%audit_status = merge(1,0,config%rho_fluid > 0._mytype)
  end subroutine init_stage5_rhs_audit_policy

  subroutine stage5_rhs_audit_noop(config, rhsx, rhsy, rhsz, rhs_change_max, rhs_modified_flag)
    type(stage5_config_t), intent(in) :: config
    real(mytype), intent(inout) :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
    real(mytype), intent(out) :: rhs_change_max
    integer, intent(out) :: rhs_modified_flag
    real(mytype), allocatable :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)

    allocate(rhsx0(size(rhsx,1),size(rhsx,2),size(rhsx,3)))
    allocate(rhsy0(size(rhsy,1),size(rhsy,2),size(rhsy,3)))
    allocate(rhsz0(size(rhsz,1),size(rhsz,2),size(rhsz,3)))
    rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz

    rhs_change_max = max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
    rhs_modified_flag = 0

    if (config%rho_fluid <= 0._mytype) rhs_change_max = rhs_change_max

    deallocate(rhsx0, rhsy0, rhsz0)
  end subroutine stage5_rhs_audit_noop

  subroutine write_stage5_rhs_audit_report(filename, audit, config, status)
    character(len=*), intent(in) :: filename
    type(stage5_rhs_audit_t), intent(in) :: audit
    type(stage5_config_t), intent(in) :: config
    integer, intent(out) :: status
    integer :: ios

    status = 0
    open(21, file=filename, status='replace', iostat=ios)
    if (ios /= 0) return

    write(21,'(A)') '# Stage 5.1 RHS Audit Report'
    write(21,'(A)') 'RHS convention: acceleration form'
    write(21,'(A)') 'IBM RHS increment: f_ibm / rho_fluid'
    write(21,'(A,1X,ES24.16E3)') 'rho_fluid:', config%rho_fluid
    write(21,'(A)') 'Insertion order: during RHS assembly before pressure projection'
    write(21,'(A)') 'Pressure Poisson equation: not directly modified by IBM force'
    write(21,'(A)') 'RK substep policy: recompute IBM force using current substep velocity'
    write(21,'(A)') 'Stale velocity policy: forbidden'
    write(21,'(A)') 'Stage 5.1 mode: audit only, no RHS injection'
    write(21,'(A,1X,I0)') 'Audit status:', audit%audit_status

    close(21)
    status = 1
  end subroutine write_stage5_rhs_audit_report

end module fibre_stage5_rhs_audit
