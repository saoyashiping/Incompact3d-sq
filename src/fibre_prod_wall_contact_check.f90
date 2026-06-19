program fibre_prod_wall_contact_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, &
                               fibre_prod_state_destroy, fibre_prod_state_is_allocated, &
                               fibre_prod_state_all_finite
  use fibre_prod_wall_geometry, only : fibre_prod_wall_geometry_type, fibre_prod_wall_geometry_init, &
                                       fibre_prod_wall_geometry_validate, fibre_prod_wall_gap_point
  use fibre_prod_wall_contact, only : fibre_prod_wall_contact_state_type, &
                                      fibre_prod_wall_contact_evaluate, &
                                      fibre_prod_wall_contact_force_candidate
  use fibre_prod_wall_contact_diagnostics, only : fibre_prod_wall_contact_force_norm, &
                                                  fibre_prod_wall_contact_energy_candidate, &
                                                  fibre_prod_wall_contact_diagnostics_finite
  use fibre_prod_structure_solver, only : fibre_prod_structure_compute_forces, &
                                          fibre_prod_structure_step_explicit
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_state_type) :: state
  type(fibre_prod_wall_geometry_type) :: geometry
  type(fibre_prod_wall_contact_state_type) :: contact_state
  real(dp), allocatable :: contact_force(:, :, :)
  real(dp) :: lower_gap
  real(dp) :: upper_gap
  real(dp) :: min_gap
  real(dp) :: k_wall
  real(dp) :: c_wall
  integer :: status

  k_wall = 100.0_dp
  c_wall = 5.0_dp
  call fibre_prod_wall_geometry_init(geometry, 0.0_dp, 1.0_dp, 0.05_dp, 0.05_dp, status)
  if (status /= 0 .or. fibre_prod_wall_geometry_validate(geometry) /= 0) error stop 1
  call fibre_prod_state_allocate(state, 1, 5, status)
  if (status /= 0 .or. .not. fibre_prod_state_is_allocated(state)) error stop 2
  allocate(contact_force(state%nfibre, state%nnode, 3))

  call set_straight_y(state, 0.50_dp)
  call run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
  if (status /= 0 .or. contact_state%min_gap <= 0.0_dp) error stop 3
  if (contact_state%near_wall_warning .or. contact_state%penetration_detected) error stop 4
  if (fibre_prod_wall_contact_force_norm(contact_force(1, :, :)) /= 0.0_dp) error stop 5
  if (.not. fibre_prod_wall_contact_diagnostics_finite(contact_state)) error stop 6

  call set_straight_y(state, 0.08_dp)
  call run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
  if (status /= 0 .or. contact_state%min_gap <= 0.0_dp) error stop 7
  if (.not. contact_state%near_wall_warning .or. contact_state%penetration_detected) error stop 8
  if (fibre_prod_wall_contact_force_norm(contact_force(1, :, :)) /= 0.0_dp) error stop 9

  call set_straight_y(state, 0.03_dp)
  state%v(1, :, 2) = -0.1_dp
  call run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
  if (status /= 0 .or. .not. contact_state%penetration_detected) error stop 10
  if (contact_state%max_penetration <= 0.0_dp .or. contact_force(1, 3, 2) <= 0.0_dp) error stop 11
  if (.not. all(contact_force == contact_force)) error stop 12
  if (fibre_prod_wall_contact_force_norm(contact_force(1, :, :)) <= 0.0_dp) error stop 12
  call structure_step_with_contact(state, contact_force, status)
  if (status /= 0 .or. state%a(1, 3, 2) <= 0.0_dp) error stop 13

  call set_straight_y(state, 0.97_dp)
  state%v(1, :, 2) = 0.1_dp
  call run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
  if (status /= 0 .or. .not. contact_state%penetration_detected) error stop 14
  if (contact_state%max_penetration <= 0.0_dp .or. contact_force(1, 3, 2) >= 0.0_dp) error stop 15
  call structure_step_with_contact(state, contact_force, status)
  if (status /= 0 .or. state%a(1, 3, 2) >= 0.0_dp) error stop 16

  call set_straight_y(state, 0.03_dp)
  state%v(1, :, 2) = -0.5_dp
  call run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
  if (status /= 0) error stop 17
  lower_gap = contact_force(1, 3, 2)
  state%v(1, :, 2) = 0.5_dp
  call run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
  if (status /= 0 .or. contact_force(1, 3, 2) < 0.0_dp .or. contact_force(1, 3, 2) > lower_gap) error stop 18

  call fibre_prod_wall_geometry_init(geometry, 1.0_dp, 0.0_dp, 0.05_dp, 0.05_dp, status)
  if (status == 0) error stop 19
  call fibre_prod_wall_geometry_init(geometry, 0.0_dp, 1.0_dp, -0.05_dp, 0.05_dp, status)
  if (status == 0) error stop 20
  call fibre_prod_wall_geometry_init(geometry, 0.0_dp, 1.0_dp, 0.05_dp, -0.05_dp, status)
  if (status == 0) error stop 21
  call fibre_prod_wall_geometry_init(geometry, 0.0_dp, 1.0_dp, 0.05_dp, 0.05_dp, status)
  if (status /= 0) error stop 22
  call run_contact_force(state, geometry, -k_wall, c_wall, contact_force, contact_state, status)
  if (status == 0) error stop 23
  state%x(1, 3, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
  call fibre_prod_wall_contact_evaluate(geometry, state%x(1, :, :), contact_state, status)
  if (status == 0) error stop 24
  call set_straight_y(state, -10.0_dp)
  call fibre_prod_wall_contact_evaluate(geometry, state%x(1, :, :), contact_state, status)
  if (status /= 0 .or. .not. contact_state%penetration_detected) error stop 25
  if (fibre_prod_wall_contact_energy_candidate(contact_state, k_wall) < 0.0_dp) error stop 26
  call fibre_prod_wall_gap_point(geometry, state%x(1, 3, 2), lower_gap, upper_gap, min_gap, status)
  if (status /= 0 .or. .not. (min_gap == min_gap)) error stop 27

  if (.not. fibre_prod_state_all_finite(state)) error stop 28
  call fibre_prod_state_destroy(state)
  if (fibre_prod_state_is_allocated(state)) error stop 29
  deallocate(contact_force)
  print *, 'R8_FIBRE_PROD_WALL_CONTACT_CHECK PASS'

contains

  subroutine set_straight_y(state, y_value)
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(in) :: y_value
    integer :: inode

    state%x = 0.0_dp
    state%v = 0.0_dp
    state%a = 0.0_dp
    do inode = 1, state%nnode
      state%x(1, inode, 1) = 0.1_dp * real(inode - 1, dp)
      state%x(1, inode, 2) = y_value
      state%x(1, inode, 3) = 0.5_dp
    end do
  end subroutine set_straight_y

  subroutine run_contact_force(state, geometry, k_wall, c_wall, contact_force, contact_state, status)
    type(fibre_prod_state_type), intent(in) :: state
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: k_wall
    real(dp), intent(in) :: c_wall
    real(dp), intent(out) :: contact_force(:, :, :)
    type(fibre_prod_wall_contact_state_type), intent(out) :: contact_state
    integer, intent(out) :: status

    contact_force = 0.0_dp
    call fibre_prod_wall_contact_force_candidate(geometry, state%x(1, :, :), state%v(1, :, :), &
                                                 k_wall, c_wall, contact_force(1, :, :), &
                                                 contact_state, status)
  end subroutine run_contact_force

  subroutine structure_step_with_contact(state, contact_force, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(in) :: contact_force(:, :, :)
    integer, intent(out) :: status

    call fibre_prod_structure_compute_forces(state, 0.0_dp, 0.1_dp, contact_force, status)
    if (status /= 0) return
    call fibre_prod_structure_step_explicit(state, 1.0e-3_dp, 1.0_dp, status)
  end subroutine structure_step_with_contact

end program fibre_prod_wall_contact_check
