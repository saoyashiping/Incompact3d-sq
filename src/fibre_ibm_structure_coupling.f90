module fibre_ibm_structure_coupling

  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_external_force, only : clear_fibre_external_force, set_fibre_external_force, add_fibre_external_force
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_ibm_feedback, only : compute_ibm_feedback_forces

  implicit none
  private

  public :: apply_ibm_structure_force_to_fibre
  public :: compute_and_apply_ibm_structure_force

contains

  subroutine apply_ibm_structure_force_to_fibre(fibre, force_on_structure, mode)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: force_on_structure(:,:)
    character(len=*), intent(in) :: mode

    if (size(force_on_structure, 1) /= 3 .or. size(force_on_structure, 2) /= fibre%nl) then
      error stop 'apply_ibm_structure_force_to_fibre: force_on_structure shape mismatch with fibre%nl'
    end if

    select case (trim(mode))
    case ('set')
      call set_fibre_external_force(fibre, force_on_structure)
    case ('add')
      call add_fibre_external_force(fibre, force_on_structure)
    case ('clear')
      call clear_fibre_external_force(fibre)
    case default
      error stop 'apply_ibm_structure_force_to_fibre: unsupported mode'
    end select
  end subroutine apply_ibm_structure_force_to_fibre

  subroutine compute_and_apply_ibm_structure_force(grid, ux, uy, uz, fibre, beta_drag, mode, &
                                                   u_lag, force_on_structure, force_on_fluid)
    type(ibm_grid_t), intent(in) :: grid
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: beta_drag
    character(len=*), intent(in) :: mode
    real(mytype), intent(out) :: u_lag(:,:), force_on_structure(:,:), force_on_fluid(:,:)

    type(ibm_lagrangian_points_t) :: lag

    integer :: l

    if (beta_drag <= 0._mytype) then
      error stop 'compute_and_apply_ibm_structure_force: beta_drag must be > 0'
    end if
    if (size(u_lag, 1) /= 3 .or. size(u_lag, 2) /= fibre%nl) then
      error stop 'compute_and_apply_ibm_structure_force: u_lag shape mismatch with fibre%nl'
    end if
    if (size(force_on_structure, 1) /= 3 .or. size(force_on_structure, 2) /= fibre%nl) then
      error stop 'compute_and_apply_ibm_structure_force: force_on_structure shape mismatch with fibre%nl'
    end if
    if (size(force_on_fluid, 1) /= 3 .or. size(force_on_fluid, 2) /= fibre%nl) then
      error stop 'compute_and_apply_ibm_structure_force: force_on_fluid shape mismatch with fibre%nl'
    end if

    lag%nl = fibre%nl
    allocate(lag%x(3, lag%nl), lag%v(3, lag%nl), lag%force(3, lag%nl), lag%weight(lag%nl))

    lag%x = fibre%x
    lag%v = fibre%v
    lag%force = 0._mytype
    lag%weight = fibre%ds
    if (lag%nl >= 1) lag%weight(1) = 0.5_mytype * fibre%ds
    if (lag%nl >= 2) lag%weight(lag%nl) = 0.5_mytype * fibre%ds

    do l = 1, lag%nl
      if (lag%weight(l) <= 0._mytype) then
        error stop 'compute_and_apply_ibm_structure_force: lag weights must be > 0'
      end if
    end do

    call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
    call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
    call apply_ibm_structure_force_to_fibre(fibre, force_on_structure, mode)

    deallocate(lag%x, lag%v, lag%force, lag%weight)
  end subroutine compute_and_apply_ibm_structure_force

end module fibre_ibm_structure_coupling
