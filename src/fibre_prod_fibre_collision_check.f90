program fibre_prod_fibre_collision_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, &
                               fibre_prod_state_destroy, fibre_prod_state_is_allocated, &
                               fibre_prod_state_all_finite
  use fibre_prod_collision_geometry, only : fibre_prod_distance_segment_segment
  use fibre_prod_fibre_collision, only : fibre_prod_collision_state_type, &
                                         fibre_prod_collision_evaluate_two_fibres, &
                                         fibre_prod_collision_force_candidate_two_fibres
  use fibre_prod_collision_diagnostics, only : fibre_prod_collision_force_norm, &
                                               fibre_prod_collision_action_reaction_residual, &
                                               fibre_prod_collision_energy_candidate, &
                                               fibre_prod_collision_diagnostics_finite
  use fibre_prod_structure_solver, only : fibre_prod_structure_compute_forces, &
                                          fibre_prod_structure_step_explicit
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_state_type) :: fibre1
  type(fibre_prod_state_type) :: fibre2
  type(fibre_prod_collision_state_type) :: collision_state
  real(dp), allocatable :: force1(:, :)
  real(dp), allocatable :: force2(:, :)
  real(dp) :: radius
  real(dp) :: warning_gap
  real(dp) :: k_collision
  real(dp) :: c_collision
  real(dp) :: closing_norm
  real(dp) :: separating_norm
  real(dp) :: segment_distance
  integer :: status

  radius = 0.05_dp
  warning_gap = 0.02_dp
  k_collision = 100.0_dp
  c_collision = 5.0_dp
  call fibre_prod_state_allocate(fibre1, 1, 5, status)
  if (status /= 0 .or. .not. fibre_prod_state_is_allocated(fibre1)) error stop 1
  call fibre_prod_state_allocate(fibre2, 1, 5, status)
  if (status /= 0 .or. .not. fibre_prod_state_is_allocated(fibre2)) error stop 2
  allocate(force1(fibre1%nnode, 3), force2(fibre2%nnode, 3))

  call set_straight_fibre(fibre1, 0.0_dp, 0.0_dp)
  call set_straight_fibre(fibre2, 0.0_dp, 0.50_dp)
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status /= 0 .or. collision_state%min_gap <= 0.0_dp) error stop 3
  if (collision_state%near_contact_warning .or. collision_state%overlap_detected) error stop 4
  if (fibre_prod_collision_force_norm(force1, force2) /= 0.0_dp) error stop 5
  if (.not. fibre_prod_collision_diagnostics_finite(collision_state)) error stop 6

  call set_straight_fibre(fibre1, 0.0_dp, 0.0_dp)
  call set_straight_fibre(fibre2, 0.0_dp, 0.115_dp)
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status /= 0 .or. collision_state%min_gap <= 0.0_dp) error stop 7
  if (.not. collision_state%near_contact_warning .or. collision_state%overlap_detected) error stop 8
  if (fibre_prod_collision_force_norm(force1, force2) /= 0.0_dp) error stop 9

  call set_straight_fibre(fibre1, 0.0_dp, 0.0_dp)
  call set_straight_fibre(fibre2, 0.0_dp, 0.08_dp)
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status /= 0 .or. .not. collision_state%overlap_detected) error stop 10
  if (collision_state%max_penetration <= 0.0_dp) error stop 11
  if (fibre_prod_collision_force_norm(force1, force2) <= 0.0_dp) error stop 12
  if (fibre_prod_collision_action_reaction_residual(force1, force2) > 1.0e-12_dp) error stop 13
  if (force1(1, 2) >= 0.0_dp .or. force2(1, 2) <= 0.0_dp) error stop 14

  segment_distance = fibre_prod_distance_segment_segment(fibre1%x(1, 1, :), fibre1%x(1, 2, :), &
                                                         fibre2%x(1, 1, :), fibre2%x(1, 2, :), status)
  if (status /= 0 .or. segment_distance < 0.0_dp .or. .not. (segment_distance == segment_distance)) error stop 15
  if (collision_state%closest_pair%fibre_i /= 1 .or. collision_state%closest_pair%fibre_j /= 2) error stop 16
  if (collision_state%closest_pair%node_i < 1 .or. collision_state%closest_pair%node_j < 1) error stop 17

  fibre1%v(1, :, 2) = 0.5_dp
  fibre2%v(1, :, 2) = -0.5_dp
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status /= 0) error stop 18
  closing_norm = fibre_prod_collision_force_norm(force1, force2)
  fibre1%v(1, :, 2) = -0.5_dp
  fibre2%v(1, :, 2) = 0.5_dp
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status /= 0) error stop 19
  separating_norm = fibre_prod_collision_force_norm(force1, force2)
  if (closing_norm <= separating_norm) error stop 20
  if (force1(1, 2) > 0.0_dp .or. force2(1, 2) < 0.0_dp) error stop 21

  call structure_step_with_collision(fibre1, force1, status)
  if (status /= 0 .or. fibre1%a(1, 1, 2) >= 0.0_dp) error stop 22
  call structure_step_with_collision(fibre2, force2, status)
  if (status /= 0 .or. fibre2%a(1, 1, 2) <= 0.0_dp) error stop 23
  if (.not. fibre_prod_state_all_finite(fibre1) .or. .not. fibre_prod_state_all_finite(fibre2)) error stop 24
  if (fibre_prod_collision_energy_candidate(collision_state, k_collision) < 0.0_dp) error stop 25

  call run_collision(fibre1, fibre2, -radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status == 0) error stop 26
  call run_collision(fibre1, fibre2, radius, -warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status == 0) error stop 27
  call run_collision(fibre1, fibre2, radius, warning_gap, -k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status == 0) error stop 28
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, -c_collision, &
                     force1, force2, collision_state, status)
  if (status == 0) error stop 29
  call set_straight_fibre(fibre1, 0.0_dp, 0.0_dp)
  call set_straight_fibre(fibre2, 0.0_dp, 0.08_dp)
  fibre1%x(1, 3, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
  call fibre_prod_collision_evaluate_two_fibres(fibre1%x(1, :, :), fibre2%x(1, :, :), radius, radius, &
                                                warning_gap, collision_state, status)
  if (status == 0) error stop 30
  call set_straight_fibre(fibre1, 0.0_dp, 0.0_dp)
  fibre1%x(1, 2, :) = fibre1%x(1, 1, :)
  call run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                     force1, force2, collision_state, status)
  if (status /= 0 .or. .not. fibre_prod_collision_diagnostics_finite(collision_state)) error stop 31

  call fibre_prod_state_destroy(fibre1)
  call fibre_prod_state_destroy(fibre2)
  if (fibre_prod_state_is_allocated(fibre1) .or. fibre_prod_state_is_allocated(fibre2)) error stop 32
  deallocate(force1, force2)
  print *, 'R9_FIBRE_PROD_FIBRE_COLLISION_CHECK PASS'

contains

  subroutine set_straight_fibre(state, x_offset, y_value)
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(in) :: x_offset
    real(dp), intent(in) :: y_value
    integer :: inode

    state%x = 0.0_dp
    state%v = 0.0_dp
    state%a = 0.0_dp
    do inode = 1, state%nnode
      state%x(1, inode, 1) = x_offset + 0.1_dp * real(inode - 1, dp)
      state%x(1, inode, 2) = y_value
      state%x(1, inode, 3) = 0.0_dp
    end do
  end subroutine set_straight_fibre

  subroutine run_collision(fibre1, fibre2, radius, warning_gap, k_collision, c_collision, &
                           force1, force2, collision_state, status)
    type(fibre_prod_state_type), intent(in) :: fibre1
    type(fibre_prod_state_type), intent(in) :: fibre2
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: warning_gap
    real(dp), intent(in) :: k_collision
    real(dp), intent(in) :: c_collision
    real(dp), intent(out) :: force1(:, :)
    real(dp), intent(out) :: force2(:, :)
    type(fibre_prod_collision_state_type), intent(out) :: collision_state
    integer, intent(out) :: status

    call fibre_prod_collision_force_candidate_two_fibres(fibre1%x(1, :, :), fibre1%v(1, :, :), &
                                                         fibre2%x(1, :, :), fibre2%v(1, :, :), &
                                                         radius, radius, warning_gap, &
                                                         k_collision, c_collision, &
                                                         force1, force2, collision_state, status)
  end subroutine run_collision

  subroutine structure_step_with_collision(fibre, force, status)
    type(fibre_prod_state_type), intent(inout) :: fibre
    real(dp), intent(in) :: force(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: external_force(:, :, :)

    allocate(external_force(1, fibre%nnode, 3))
    external_force(1, :, :) = force
    call fibre_prod_structure_compute_forces(fibre, 0.0_dp, 0.1_dp, external_force, status)
    if (status == 0) call fibre_prod_structure_step_explicit(fibre, 1.0e-3_dp, 1.0_dp, status)
    deallocate(external_force)
  end subroutine structure_step_with_collision

end program fibre_prod_fibre_collision_check
