module fibre_stage6_rhs_audit
  use fibre_parameters, only : mytype
  use fibre_stage6_config, only : stage6_config_t
  implicit none
  integer, parameter :: STAGE6_RHS_FORM_ACCELERATION = 1
  integer, parameter :: STAGE6_RHS_FORM_FORCE_DENSITY = 2
  integer, parameter :: STAGE6_INSERT_MOMENTUM_RHS = 1
  integer, parameter :: STAGE6_INSERT_PRESSURE_POISSON = 2
  integer, parameter :: STAGE6_INSERT_AFTER_PROJECTION = 3
  type stage6_rhs_audit_t
    integer :: rhs_form,insertion_target,insertion_before_projection,pressure_poisson_direct_modify
    integer :: requires_divide_by_rho,current_substep_rhs_required,current_substep_force_required
    integer :: stale_force_forbidden,force_recompute_each_substep,audit_status
  end type
contains
  subroutine init_stage6_rhs_audit_policy(config,audit)
    type(stage6_config_t), intent(in) :: config
    type(stage6_rhs_audit_t), intent(out) :: audit
    audit%rhs_form=STAGE6_RHS_FORM_ACCELERATION; audit%insertion_target=STAGE6_INSERT_MOMENTUM_RHS
    audit%insertion_before_projection=1; audit%pressure_poisson_direct_modify=0; audit%requires_divide_by_rho=1
    audit%current_substep_rhs_required=1; audit%current_substep_force_required=1; audit%stale_force_forbidden=1
    audit%force_recompute_each_substep=1; audit%audit_status=merge(1,0,config%rho_fluid>0._mytype)
  end subroutine
  subroutine stage6_rhs_audit_noop(config,rhsx,rhsy,rhsz,rhs_change_max,rhs_modified_flag)
    type(stage6_config_t), intent(in) :: config
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype), intent(out) :: rhs_change_max
    integer, intent(out) :: rhs_modified_flag
    real(mytype), allocatable :: x0(:,:,:),y0(:,:,:),z0(:,:,:)
    allocate(x0(size(rhsx,1),size(rhsx,2),size(rhsx,3)),y0(size(rhsy,1),size(rhsy,2),size(rhsy,3)),z0(size(rhsz,1),size(rhsz,2),size(rhsz,3)))
    x0=rhsx; y0=rhsy; z0=rhsz
    rhs_change_max=max(maxval(abs(rhsx-x0)),max(maxval(abs(rhsy-y0)),maxval(abs(rhsz-z0)))); rhs_modified_flag=0
    deallocate(x0,y0,z0)
  end subroutine
  subroutine write_stage6_rhs_audit_report(filename,audit,config,status)
    character(len=*), intent(in) :: filename
    type(stage6_rhs_audit_t), intent(in) :: audit
    type(stage6_config_t), intent(in) :: config
    integer, intent(out) :: status
    integer :: u,ios
    open(newunit=u,file=filename,status='replace',iostat=ios); if (ios/=0) then; status=0; return; end if
    write(u,'(A)') '# Stage 6.1 RHS Audit Report'
    write(u,'(A)') 'Stage 6 mode: audit only, no RHS injection'
    write(u,'(A)') 'RHS convention: acceleration form'
    write(u,'(A)') 'IBM RHS increment: f_ibm / rho_fluid'
    write(u,'(A,1X,ES12.4)') 'rho_fluid:',config%rho_fluid
    write(u,'(A)') 'Insertion target: momentum RHS arrays for u/v/w'
    write(u,'(A)') 'Insertion order: before pressure projection'
    write(u,'(A)') 'Pressure Poisson equation: not directly modified'
    write(u,'(A)') 'RK policy: recompute IBM force using current substep velocity/state'
    write(u,'(A)') 'Stale force policy: forbidden'
    write(u,'(A)') 'Production two-way enabled by default: false'
    close(u); status=1
  end subroutine
end module
