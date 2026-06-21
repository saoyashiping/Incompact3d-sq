program fibre_prod_structure_dry_step_check
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                               fibre_prod_state_get_structure_coordinates, &
                               fibre_prod_state_get_structure_velocity_or_zero
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init
  use fibre_prod_runtime_bridge, only : fibre_prod_runtime_bridge_type, &
                                        fibre_prod_runtime_bridge_init_from_rhs, &
                                        fibre_prod_runtime_bridge_apply_lambda0_noop, &
                                        fibre_prod_runtime_bridge_finalize
  use fibre_prod_state_velocity_attachment, only : fibre_prod_state_velocity_attachment_type, &
                                                   fibre_prod_state_velocity_attachment_init, &
                                                   fibre_prod_state_velocity_attachment_set_points, &
                                                   fibre_prod_state_velocity_attachment_sample, &
                                                   fibre_prod_state_velocity_attachment_attach_to_state, &
                                                   fibre_prod_state_velocity_attachment_finalize
  use fibre_prod_hydro_input_candidate, only : fibre_prod_hydro_input_candidate_type, &
                                               fibre_prod_hydro_input_candidate_init, &
                                               fibre_prod_hydro_input_candidate_compute, &
                                               fibre_prod_hydro_input_candidate_attach_to_state, &
                                               fibre_prod_hydro_input_candidate_finalize
  use fibre_prod_structure_input_handoff, only : fibre_prod_structure_input_handoff_type, &
                                                 fibre_prod_structure_input_handoff_init, &
                                                 fibre_prod_structure_input_handoff_from_candidate, &
                                                 fibre_prod_structure_input_handoff_attach_to_state, &
                                                 fibre_prod_structure_input_handoff_finalize
  use fibre_prod_structure_dry_step, only : fibre_prod_structure_dry_step_type, &
                                            fibre_prod_structure_dry_step_init, &
                                            fibre_prod_structure_dry_step_predict, &
                                            fibre_prod_structure_dry_step_check_bounded, &
                                            fibre_prod_structure_dry_step_finalize
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: n = 9
  integer, parameter :: nnode = 3
  real(dp), parameter :: beta = 2.5_dp
  real(dp), parameter :: dt = 1.0e-4_dp
  real(dp), parameter :: rho_eff = 2.0_dp
  real(dp), parameter :: tol = 1.0e-12_dp
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_state_type) :: state
  type(fibre_prod_state_velocity_attachment_type) :: attach
  type(fibre_prod_hydro_input_candidate_type) :: candidate
  type(fibre_prod_structure_input_handoff_type) :: handoff
  type(fibre_prod_structure_dry_step_type) :: dry
  type(fibre_prod_runtime_bridge_type) :: bridge
  type(fibre_prod_runtime_config_type) :: config
  type(fibre_prod_force_buffer_type) :: force_buffer
  real(dp) :: x(n), y(n), z(n)
  real(dp) :: ux(n, n, n), uy(n, n, n), uz(n, n, n)
  real(dp) :: ux0(n, n, n), uy0(n, n, n), uz0(n, n, n)
  real(dp) :: rhs_x(n, n, n), rhs_y(n, n, n), rhs_z(n, n, n)
  real(dp) :: rhs_x0(n, n, n), rhs_y0(n, n, n), rhs_z0(n, n, n)
  real(dp) :: points(nnode, 3), state_x(nnode, 3), structure_u(nnode, 3)
  real(dp) :: x0(1, nnode, 3), v0(1, nnode, 3), a0(1, nnode, 3)
  real(dp) :: sampled0(nnode, 3), hydro0(nnode, 3), input0(nnode, 3)
  real(dp) :: expected_dx(nnode, 3), expected_x(nnode, 3), expected_u(nnode, 3)
  real(dp) :: bad_input(2, 3), huge_input(nnode, 3)
  integer :: status

  call initialize_grid_coordinates(x, y, z)
  call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, n, 1, n, 1, n, .false., .false., .false., status)
  if (status /= 0) error stop 1
  call initialize_velocity_field(grid, ux, uy, uz)
  ux0 = ux; uy0 = uy; uz0 = uz
  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  rhs_x0 = rhs_x; rhs_y0 = rhs_y; rhs_z0 = rhs_z

  call setup_state_pipeline(status)
  if (status /= 0) error stop 2
  x0 = state%x; v0 = state%v; a0 = state%a
  sampled0 = state%sampled_u
  hydro0 = state%hydro_force_candidate
  input0 = state%structure_input_force
  call fibre_prod_state_get_structure_coordinates(state, state_x, status)
  if (status /= 0) error stop 3
  call fibre_prod_state_get_structure_velocity_or_zero(state, structure_u, status)
  if (status /= 0) error stop 4

  call fibre_prod_structure_dry_step_init(dry, nnode, dt, rho_eff, status)
  if (status /= 0) error stop 5
  if (.not. allocated(dry%x_trial) .or. .not. allocated(dry%u_trial) .or. .not. allocated(dry%dx_trial)) error stop 6
  if (any(shape(dry%x_trial) /= [nnode, 3]) .or. any(shape(dry%u_trial) /= [nnode, 3]) .or. &
      any(shape(dry%dx_trial) /= [nnode, 3])) error stop 7
  call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
  if (status /= 0) error stop 8
  expected_dx = dt * structure_u + 0.5_dp * dt * dt * state%structure_input_force / rho_eff
  expected_x = state_x + expected_dx
  expected_u = structure_u + dt * state%structure_input_force / rho_eff
  if (maxval(abs(dry%dx_trial - expected_dx)) > tol) error stop 9
  if (maxval(abs(dry%x_trial - expected_x)) > tol) error stop 10
  if (maxval(abs(dry%u_trial - expected_u)) > tol) error stop 11
  if (abs(dry%max_abs_dx_trial - maxval(abs(expected_dx))) > tol) error stop 12
  if (abs(dry%max_abs_u_trial - maxval(abs(expected_u))) > tol) error stop 13
  if (abs(dry%sum_abs_dx_trial - sum(abs(expected_dx))) > tol) error stop 14
  call fibre_prod_structure_dry_step_check_bounded(dry, 1.0e-2_dp, status)
  if (status /= 0) error stop 15

  huge_input = 1.0e8_dp
  call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, huge_input, status)
  if (status /= 0) error stop 16
  call fibre_prod_structure_dry_step_check_bounded(dry, 1.0e-8_dp, status)
  if (status == 0) error stop 17
  call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
  if (status /= 0) error stop 18

  if (maxval(abs(state%x - x0)) > 0.0_dp) error stop 19
  if (maxval(abs(state%v - v0)) > 0.0_dp) error stop 20
  if (maxval(abs(state%a - a0)) > 0.0_dp) error stop 21
  if (maxval(abs(state%sampled_u - sampled0)) > 0.0_dp) error stop 22
  if (maxval(abs(state%hydro_force_candidate - hydro0)) > 0.0_dp) error stop 23
  if (maxval(abs(state%structure_input_force - input0)) > 0.0_dp) error stop 24
  if (.not. same_field(ux, ux0) .or. .not. same_field(uy, uy0) .or. .not. same_field(uz, uz0)) error stop 25
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 26

  call fibre_prod_force_buffer_allocate(force_buffer, grid, status)
  if (status /= 0) error stop 27
  call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
  if (status /= 0) error stop 28
  if (any(force_buffer%fx /= 0.0_dp) .or. any(force_buffer%fy /= 0.0_dp) .or. &
      any(force_buffer%fz /= 0.0_dp)) error stop 29

  call fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 30
  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.; config%lambda_fsi = 0.0_dp; config%penalty_beta = 2.0_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 31
  call fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 32
  call run_pipeline_and_dry_step(status)
  if (status /= 0) error stop 33
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 34

  config%lambda_fsi = 1.0e-3_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 35
  call run_pipeline_and_dry_step(status)
  if (status /= 0) error stop 36
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 37
  if (any(force_buffer%fx /= 0.0_dp) .or. any(force_buffer%fy /= 0.0_dp) .or. &
      any(force_buffer%fz /= 0.0_dp)) error stop 38

  bad_input = 0.0_dp
  call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, bad_input, status)
  if (status == 0) error stop 39

  call fibre_prod_structure_dry_step_finalize(dry)
  call fibre_prod_structure_input_handoff_finalize(handoff)
  call fibre_prod_hydro_input_candidate_finalize(candidate)
  call fibre_prod_state_velocity_attachment_finalize(attach)
  call fibre_prod_runtime_bridge_finalize(bridge)
  call fibre_prod_force_buffer_destroy(force_buffer)
  call fibre_prod_state_destroy(state)
  call fibre_prod_grid_destroy(grid)
  print *, 'P0_8_STRUCTURE_DRY_STEP_CHECK PASS'
contains

  subroutine setup_state_pipeline(status)
    integer, intent(out) :: status
    structure_u = 0.0_dp
    call fibre_prod_state_allocate(state, 1, nnode, status)
    if (status /= 0) return
    points(1, :) = [0.25_dp, 0.50_dp, 0.50_dp]
    points(2, :) = [0.50_dp, 0.50_dp, 0.50_dp]
    points(3, :) = [0.75_dp, 0.50_dp, 0.50_dp]
    state%x(1, :, :) = points
    call fibre_prod_state_velocity_attachment_init(attach, nnode, status)
    if (status /= 0) return
    call fibre_prod_state_velocity_attachment_set_points(attach, state%x(1, :, :), status)
    if (status /= 0) return
    call fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    if (status /= 0) return
    call fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_init(candidate, nnode, beta, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_compute(candidate, state%sampled_u, structure_u, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_attach_to_state(candidate, state, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_init(handoff, nnode, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
  end subroutine setup_state_pipeline

  subroutine run_pipeline_and_dry_step(status)
    integer, intent(out) :: status
    call fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    if (status /= 0) return
    call fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_compute(candidate, state%sampled_u, structure_u, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_attach_to_state(candidate, state, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    if (status /= 0) return
    call fibre_prod_state_get_structure_coordinates(state, state_x, status)
    if (status /= 0) return
    call fibre_prod_state_get_structure_velocity_or_zero(state, structure_u, status)
    if (status /= 0) return
    call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
  end subroutine run_pipeline_and_dry_step

  subroutine initialize_grid_coordinates(x, y, z)
    real(dp), intent(out) :: x(:), y(:), z(:)
    integer :: i
    do i = 1, size(x)
      x(i) = real(i - 1, dp) / real(size(x) - 1, dp)
      y(i) = real(i - 1, dp) / real(size(y) - 1, dp)
      z(i) = real(i - 1, dp) / real(size(z) - 1, dp)
    end do
  end subroutine initialize_grid_coordinates

  subroutine initialize_velocity_field(grid, ux, uy, uz)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(out) :: ux(:, :, :), uy(:, :, :), uz(:, :, :)
    integer :: i, j, k
    real(dp) :: point(3), velocity(3)
    do k = 1, size(ux, 3); do j = 1, size(ux, 2); do i = 1, size(ux, 1)
      point = [grid%x(i), grid%y(j), grid%z(k)]
      velocity = analytic_velocity(point)
      ux(i, j, k) = velocity(1); uy(i, j, k) = velocity(2); uz(i, j, k) = velocity(3)
    end do; end do; end do
  end subroutine initialize_velocity_field

  subroutine initialize_rhs(rhs_x, rhs_y, rhs_z)
    real(dp), intent(out) :: rhs_x(:, :, :), rhs_y(:, :, :), rhs_z(:, :, :)
    integer :: i, j, k
    do k = 1, size(rhs_x, 3); do j = 1, size(rhs_x, 2); do i = 1, size(rhs_x, 1)
      rhs_x(i, j, k) = real(i + j + k, dp)
      rhs_y(i, j, k) = -real(i + 2*j + k, dp)
      rhs_z(i, j, k) = real(i - j + 2*k, dp)
    end do; end do; end do
  end subroutine initialize_rhs

  pure function analytic_velocity(point) result(velocity)
    real(dp), intent(in) :: point(3)
    real(dp) :: velocity(3)
    velocity(1) = 1.0_dp + point(1) + 2.0_dp * point(2) + 3.0_dp * point(3)
    velocity(2) = -2.0_dp + 2.0_dp * point(1) - point(2) + 0.5_dp * point(3)
    velocity(3) = 0.25_dp + 0.5_dp * point(1) + point(2) - point(3)
  end function analytic_velocity

  logical function same_field(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :), rhs(:, :, :)
    matches = all(lhs == rhs)
  end function same_field
end program fibre_prod_structure_dry_step_check
