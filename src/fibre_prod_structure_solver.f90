module fibre_prod_structure_solver
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_all_finite, &
                               fibre_prod_state_reset_forces
  use fibre_prod_bending_solver, only : fibre_prod_bending_compute_force
  use fibre_prod_tension_solver, only : fibre_prod_tension_force_diagnostic_candidate
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_structure_compute_forces
  public :: fibre_prod_structure_step_explicit
  public :: fibre_prod_structure_kinetic_energy
  public :: fibre_prod_structure_bending_energy_candidate
  public :: fibre_prod_structure_validate_state

contains

  subroutine fibre_prod_structure_compute_forces(state, gamma, ds, external_force, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: ds
    real(dp), intent(in), optional :: external_force(:, :, :)
    integer, intent(out) :: status
    real(dp), allocatable :: bend(:, :)
    real(dp), allocatable :: tension(:, :)
    integer :: i_fibre
    integer :: ierr

    status = validate_parameters(gamma=gamma, ds=ds)
    if (status /= 0) return
    if (.not. fibre_prod_structure_validate_state(state)) then
      status = 10
      return
    end if
    if (present(external_force)) then
      if (.not. force_shape_matches_state(state, external_force)) then
        status = 11
        return
      end if
      if (.not. all(ieee_is_finite(external_force))) then
        status = 12
        return
      end if
    end if

    call fibre_prod_state_reset_forces(state)
    allocate(bend(state%nnode, 3), tension(state%nnode, 3))
    do i_fibre = 1, state%nfibre
      call fibre_prod_bending_compute_force(state%x(i_fibre, :, :), gamma, ds, bend, ierr)
      if (ierr /= 0) then
        status = 20 + ierr
        exit
      end if
      call fibre_prod_tension_force_diagnostic_candidate(state%x(i_fibre, :, :), ds, 0.0_dp, tension, ierr)
      if (ierr /= 0) then
        status = 30 + ierr
        exit
      end if
      state%f_total(i_fibre, :, :) = bend + tension
      if (present(external_force)) state%f_total(i_fibre, :, :) = state%f_total(i_fibre, :, :) + external_force(i_fibre, :, :)
    end do
    deallocate(bend, tension)
  end subroutine fibre_prod_structure_compute_forces

  subroutine fibre_prod_structure_step_explicit(state, dt, rho_tilde, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: rho_tilde
    integer, intent(out) :: status

    status = 0
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      status = 1
    else if (.not. ieee_is_finite(rho_tilde) .or. rho_tilde <= 0.0_dp) then
      status = 2
    else if (.not. fibre_prod_structure_validate_state(state)) then
      status = 3
    else if (.not. all(ieee_is_finite(state%f_total))) then
      status = 4
    else
      state%a = state%f_total / rho_tilde
      state%x = state%x + dt * state%v + 0.5_dp * dt * dt * state%a
      state%v = state%v + dt * state%a
    end if
  end subroutine fibre_prod_structure_step_explicit

  pure real(dp) function fibre_prod_structure_kinetic_energy(state, rho_tilde) result(energy)
    type(fibre_prod_state_type), intent(in) :: state
    real(dp), intent(in) :: rho_tilde

    if (rho_tilde > 0.0_dp .and. fibre_prod_structure_validate_state(state)) then
      energy = 0.5_dp * rho_tilde * sum(state%v * state%v)
    else
      energy = huge(1.0_dp)
    end if
  end function fibre_prod_structure_kinetic_energy

  pure real(dp) function fibre_prod_structure_bending_energy_candidate(state, gamma, ds) result(energy)
    type(fibre_prod_state_type), intent(in) :: state
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: ds
    integer :: i_fibre
    integer :: i_node
    real(dp) :: curvature(3)

    energy = huge(1.0_dp)
    if (validate_parameters(gamma=gamma, ds=ds) /= 0) return
    if (.not. fibre_prod_structure_validate_state(state)) return

    energy = 0.0_dp
    do i_fibre = 1, state%nfibre
      do i_node = 2, state%nnode - 1
        curvature = (state%x(i_fibre, i_node - 1, :) - 2.0_dp * state%x(i_fibre, i_node, :) + &
                     state%x(i_fibre, i_node + 1, :)) / (ds * ds)
        energy = energy + 0.5_dp * gamma * sum(curvature * curvature) * ds
      end do
    end do
  end function fibre_prod_structure_bending_energy_candidate

  pure logical function fibre_prod_structure_validate_state(state) result(valid)
    type(fibre_prod_state_type), intent(in) :: state

    valid = state%nfibre >= 1 .and. state%nnode >= 5 .and. state%ndim == 3
    valid = valid .and. fibre_prod_state_all_finite(state)
  end function fibre_prod_structure_validate_state

  pure integer function validate_parameters(dt, rho_tilde, gamma, ds) result(status)
    real(dp), intent(in), optional :: dt
    real(dp), intent(in), optional :: rho_tilde
    real(dp), intent(in), optional :: gamma
    real(dp), intent(in), optional :: ds

    status = 0
    if (present(dt)) then
      if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) status = 1
    end if
    if (present(rho_tilde) .and. status == 0) then
      if (.not. ieee_is_finite(rho_tilde) .or. rho_tilde <= 0.0_dp) status = 2
    end if
    if (present(gamma) .and. status == 0) then
      if (.not. ieee_is_finite(gamma) .or. gamma < 0.0_dp) status = 3
    end if
    if (present(ds) .and. status == 0) then
      if (.not. ieee_is_finite(ds) .or. ds <= 0.0_dp) status = 4
    end if
  end function validate_parameters

  pure logical function force_shape_matches_state(state, force) result(matches)
    type(fibre_prod_state_type), intent(in) :: state
    real(dp), intent(in) :: force(:, :, :)

    matches = size(force, 1) == state%nfibre .and. size(force, 2) == state%nnode .and. size(force, 3) == 3
  end function force_shape_matches_state

end module fibre_prod_structure_solver
