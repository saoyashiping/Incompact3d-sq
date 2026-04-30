module fibre_stage4_feedback_adapter

  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_external_force, only : clear_fibre_external_force, set_fibre_external_force, add_fibre_external_force
  use fibre_ibm_types, only : ibm_lagrangian_points_t
  use fibre_ibm_feedback, only : compute_ibm_feedback_forces
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  use fibre_stage4_interpolation_adapter, only : stage4_adapter_can_use_uniform_collocated_ibm, interpolate_stage4_vector_to_lag_if_supported

  implicit none
  private
  public :: compute_stage4_feedback_if_supported
  public :: apply_stage4_feedback_to_f_ext
contains

  subroutine compute_stage4_feedback_if_supported(adapter, ux, uy, uz, fibre, beta_drag, u_lag, force_on_structure, force_on_fluid, status)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: beta_drag
    real(mytype), intent(out) :: u_lag(:,:)
    real(mytype), intent(out) :: force_on_structure(:,:), force_on_fluid(:,:)
    integer, intent(out) :: status

    type(ibm_lagrangian_points_t) :: lag
    integer :: can_use, i

    call stage4_adapter_can_use_uniform_collocated_ibm(adapter, can_use)
    if (can_use /= 1) then
      status = 0
      u_lag = 0._mytype
      force_on_structure = 0._mytype
      force_on_fluid = 0._mytype
      return
    end if

    lag%nl = fibre%nl
    allocate(lag%x(3,lag%nl), lag%v(3,lag%nl), lag%force(3,lag%nl), lag%weight(lag%nl))
    lag%x = fibre%x
    lag%v = fibre%v
    lag%force = 0._mytype
    lag%weight = fibre%ds
    if (lag%nl >= 1) lag%weight(1) = 0.5_mytype * fibre%ds
    if (lag%nl >= 2) lag%weight(lag%nl) = 0.5_mytype * fibre%ds

    call interpolate_stage4_vector_to_lag_if_supported(adapter, ux, uy, uz, lag, u_lag, status)
    if (status /= 1) then
      force_on_structure = 0._mytype
      force_on_fluid = 0._mytype
      return
    end if

    call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
    status = 1

    deallocate(lag%x, lag%v, lag%force, lag%weight)
  end subroutine

  subroutine apply_stage4_feedback_to_f_ext(fibre, force_on_structure, mode, status)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: force_on_structure(:,:)
    character(len=*), intent(in) :: mode
    integer, intent(out) :: status

    if (trim(mode) == 'set') then
      call set_fibre_external_force(fibre, force_on_structure)
      status = 1
    else if (trim(mode) == 'add') then
      call add_fibre_external_force(fibre, force_on_structure)
      status = 1
    else if (trim(mode) == 'clear') then
      call clear_fibre_external_force(fibre)
      status = 1
    else
      status = 0
    end if
  end subroutine

end module fibre_stage4_feedback_adapter
