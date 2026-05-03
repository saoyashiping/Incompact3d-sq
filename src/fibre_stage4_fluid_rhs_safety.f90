module fibre_stage4_fluid_rhs_safety

  use fibre_parameters, only : mytype
  use fibre_stage4_config, only : stage4_oneway_config_t
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t

  implicit none
  private

  public :: stage4_should_apply_ibm_to_fluid_rhs
  public :: stage4_apply_ibm_rhs_safely
  public :: compute_rhs_change_max

contains

  subroutine stage4_should_apply_ibm_to_fluid_rhs(config, allowed, rejected_flag)
    type(stage4_oneway_config_t), intent(in) :: config
    integer, intent(out) :: allowed
    integer, intent(out) :: rejected_flag

    allowed = 0
    rejected_flag = 0

    if (config%coupling_mode == 2 .or. config%apply_ibm_to_fluid_rhs) then
      rejected_flag = 1
      allowed = 0
      return
    end if

    allowed = 0
  end subroutine stage4_should_apply_ibm_to_fluid_rhs

  subroutine stage4_apply_ibm_rhs_safely(config, buffer, rhsx, rhsy, rhsz, rhs_modified_flag, rejected_flag)
    type(stage4_oneway_config_t), intent(in) :: config
    type(ibm_force_buffer_t), intent(in) :: buffer
    real(mytype), intent(inout) :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
    integer, intent(out) :: rhs_modified_flag
    integer, intent(out) :: rejected_flag

    integer :: allowed

    call stage4_should_apply_ibm_to_fluid_rhs(config, allowed, rejected_flag)

    rhs_modified_flag = 0

    if (.not. buffer%is_allocated) then
      rejected_flag = 1
    end if

    ! Stage 4.7 safety rule: never alter main fluid RHS here.
    ! rhsx/rhsy/rhsz are intentionally left unchanged.
  end subroutine stage4_apply_ibm_rhs_safely

  subroutine compute_rhs_change_max(rhsx0, rhsy0, rhsz0, rhsx1, rhsy1, rhsz1, change_max)
    real(mytype), intent(in) :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
    real(mytype), intent(in) :: rhsx1(:,:,:), rhsy1(:,:,:), rhsz1(:,:,:)
    real(mytype), intent(out) :: change_max

    change_max = max( maxval(abs(rhsx1-rhsx0)), max(maxval(abs(rhsy1-rhsy0)), maxval(abs(rhsz1-rhsz0))) )
  end subroutine compute_rhs_change_max

end module fibre_stage4_fluid_rhs_safety
