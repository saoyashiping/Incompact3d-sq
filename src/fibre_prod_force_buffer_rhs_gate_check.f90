program fibre_prod_force_buffer_rhs_gate_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_reset_to_zero, fibre_prod_force_buffer_destroy
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity, fibre_prod_state_attach_hydro_force_candidate, &
                              fibre_prod_state_attach_structure_input_force, fibre_prod_state_attach_structure_u
  use fibre_prod_force_buffer_rhs_gate, only : fibre_prod_force_buffer_rhs_gate_type, &
                                               fibre_prod_force_buffer_rhs_gate_init, &
                                               fibre_prod_force_buffer_rhs_gate_apply, &
                                               fibre_prod_force_buffer_rhs_gate_check_linear_response, &
                                               fibre_prod_force_buffer_rhs_gate_finalize
  implicit none

  integer, parameter :: dp = real64
  integer :: status
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_grid_type) :: small_grid
  type(fibre_prod_force_buffer_type) :: force_buffer
  type(fibre_prod_force_buffer_type) :: zero_buffer
  type(fibre_prod_force_buffer_type) :: bad_shape_buffer
  type(fibre_prod_force_buffer_type) :: bad_value_buffer
  type(fibre_prod_force_buffer_type) :: unallocated_buffer
  type(fibre_prod_force_buffer_rhs_gate_type) :: gate
  type(fibre_prod_runtime_config_type) :: config
  type(fibre_prod_state_type) :: state
  real(dp) :: coord(3), small_coord(2)
  real(dp) :: rhs_x(3,3,3), rhs_y(3,3,3), rhs_z(3,3,3)
  real(dp) :: rhs0_x(3,3,3), rhs0_y(3,3,3), rhs0_z(3,3,3)
  real(dp) :: inc_lambda1(3,3,3), inc_lambda2(3,3,3), inc_beta1(3,3,3), inc_beta2(3,3,3)
  real(dp) :: ux(3,3,3), uy(3,3,3), uz(3,3,3), ux0(3,3,3), uy0(3,3,3), uz0(3,3,3)
  real(dp) :: x0(3,3), structure_u0(3,3), sampled0(3,3), hydro0(3,3), input0(3,3)
  integer :: i

  coord = [0.0_dp, 0.5_dp, 1.0_dp]
  small_coord = [0.0_dp, 1.0_dp]
  call fibre_prod_grid_init_from_coordinates(grid, coord, coord, coord, 1, 3, 1, 3, 1, 3, .false., .false., .false., status)
  call require(status == 0, 'grid init failed')
  call fibre_prod_grid_init_from_coordinates(small_grid, small_coord, small_coord, small_coord, 1, 2, 1, 2, 1, 2, &
                                             .false., .false., .false., status)
  call require(status == 0, 'small grid init failed')
  call fibre_prod_force_buffer_allocate(force_buffer, grid, status)
  call require(status == 0, 'force buffer allocation failed')
  call fill_nonuniform_force_buffer(force_buffer)
  call fibre_prod_force_buffer_allocate(zero_buffer, grid, status)
  call require(status == 0, 'zero buffer allocation failed')
  call fibre_prod_force_buffer_allocate(bad_shape_buffer, small_grid, status)
  call require(status == 0, 'bad shape buffer allocation failed')
  call fibre_prod_force_buffer_allocate(bad_value_buffer, grid, status)
  call require(status == 0, 'bad value buffer allocation failed')
  call fill_nonuniform_force_buffer(bad_value_buffer)
  bad_value_buffer%fx(1,1,1) = ieee_value(0.0_dp, ieee_quiet_nan)

  call setup_state(state, x0, structure_u0, sampled0, hydro0, input0, status)
  call require(status == 0, 'state setup failed')
  ux = reshape([(real(i, dp), i=1,27)], shape(ux)); uy = ux + 10.0_dp; uz = ux - 10.0_dp
  ux0 = ux; uy0 = uy; uz0 = uz
  rhs_x = ux * 0.1_dp; rhs_y = uy * 0.1_dp; rhs_z = uz * 0.1_dp
  rhs0_x = rhs_x; rhs0_y = rhs_y; rhs0_z = rhs_z

  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.; config%lambda_fsi = 0.0_dp; config%penalty_beta = 2.0_dp
  call fibre_prod_force_buffer_rhs_gate_init(gate, 3, 3, 3, config%lambda_fsi, config%penalty_beta, status)
  call require(status == 0, 'lambda0 gate init failed')
  call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
  call require(status == 0 .and. gate%lambda_zero_noop, 'lambda0 gate apply failed')
  call require(all(rhs_x == rhs0_x) .and. all(rhs_y == rhs0_y) .and. all(rhs_z == rhs0_z), 'lambda0 changed RHS')

  config%lambda_fsi = 1.0e-3_dp; config%penalty_beta = 2.0_dp
  rhs_x = rhs0_x; rhs_y = rhs0_y; rhs_z = rhs0_z
  call fibre_prod_force_buffer_rhs_gate_init(gate, 3, 3, 3, config%lambda_fsi, config%penalty_beta, status)
  call require(status == 0, 'small-lambda gate init failed')
  call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
  call require(status == 0 .and. gate%applied, 'small-lambda gate apply failed')
  inc_lambda1 = rhs_x - rhs0_x
  call require(maxval(abs(inc_lambda1)) > 0.0_dp, 'small-lambda increment is zero')
  call require(maxval(abs(inc_lambda1 - config%lambda_fsi * config%penalty_beta * force_buffer%fx)) <= 1.0e-12_dp, &
               'x increment does not match force buffer scale')
  call require(gate%max_abs_rhs_increment < 1.0e3_dp .and. gate%measured_scale_error <= 1.0e-10_dp, &
               'small-lambda diagnostics failed')
  call require(maxval(inc_lambda1) - minval(inc_lambda1) > 0.0_dp, 'increment became uniform')

  config%lambda_fsi = 2.0e-3_dp; config%penalty_beta = 2.0_dp
  rhs_x = rhs0_x; rhs_y = rhs0_y; rhs_z = rhs0_z
  call fibre_prod_force_buffer_rhs_gate_init(gate, 3, 3, 3, config%lambda_fsi, config%penalty_beta, status)
  call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
  call require(status == 0, 'lambda2 gate apply failed')
  inc_lambda2 = rhs_x - rhs0_x
  call fibre_prod_force_buffer_rhs_gate_check_linear_response(gate, inc_lambda1, inc_lambda2, 1.0e-3_dp, 2.0e-3_dp, status)
  call require(status == 0, 'lambda linear response failed')

  config%lambda_fsi = 1.0e-3_dp; config%penalty_beta = 1.0_dp
  rhs_x = rhs0_x; rhs_y = rhs0_y; rhs_z = rhs0_z
  call fibre_prod_force_buffer_rhs_gate_init(gate, 3, 3, 3, config%lambda_fsi, config%penalty_beta, status)
  call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
  call require(status == 0, 'beta1 gate apply failed')
  inc_beta1 = rhs_x - rhs0_x
  config%penalty_beta = 2.0_dp
  rhs_x = rhs0_x; rhs_y = rhs0_y; rhs_z = rhs0_z
  call fibre_prod_force_buffer_rhs_gate_init(gate, 3, 3, 3, config%lambda_fsi, config%penalty_beta, status)
  call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
  call require(status == 0, 'beta2 gate apply failed')
  inc_beta2 = rhs_x - rhs0_x
  call require(maxval(abs(inc_beta2 - 2.0_dp * inc_beta1)) <= 1.0e-12_dp, 'penalty-beta linear response failed')

  rhs_x = rhs0_x; rhs_y = rhs0_y; rhs_z = rhs0_z
  call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, zero_buffer, status)
  call require(status == 0 .and. all(rhs_x == rhs0_x) .and. all(rhs_y == rhs0_y) .and. all(rhs_z == rhs0_z), &
               'zero force buffer changed RHS')
  call assert_fail_closed(gate, config, rhs0_x, rhs0_y, rhs0_z, unallocated_buffer, 'unallocated buffer')
  call assert_fail_closed(gate, config, rhs0_x, rhs0_y, rhs0_z, bad_shape_buffer, 'shape mismatch buffer')
  call assert_fail_closed(gate, config, rhs0_x, rhs0_y, rhs0_z, bad_value_buffer, 'nonfinite buffer')

  call require(all(ux == ux0) .and. all(uy == uy0) .and. all(uz == uz0), 'velocity fields changed')
  call require(all(state%x(1,:,:) == x0), 'coordinates changed')
  call require(all(state%structure_u == structure_u0), 'structure velocity changed')
  call require(all(state%sampled_u == sampled0), 'sampled velocity changed')
  call require(all(state%hydro_force_candidate == hydro0), 'hydro candidate changed')
  call require(all(state%structure_input_force == input0), 'structure input changed')

  print *, 'P0_11_DIAGNOSTIC: force-buffer RHS gate only, lambda0 noop and small-lambda bounded response verified'
  print *, 'P0_11_FORCE_BUFFER_RHS_GATE_CHECK PASS'

  call fibre_prod_force_buffer_rhs_gate_finalize(gate)
  call fibre_prod_force_buffer_destroy(force_buffer)
  call fibre_prod_force_buffer_destroy(zero_buffer)
  call fibre_prod_force_buffer_destroy(bad_shape_buffer)
  call fibre_prod_force_buffer_destroy(bad_value_buffer)
  call fibre_prod_grid_destroy(grid)
  call fibre_prod_grid_destroy(small_grid)
  call fibre_prod_state_destroy(state)
contains
  subroutine fill_nonuniform_force_buffer(buffer)
    type(fibre_prod_force_buffer_type), intent(inout) :: buffer
    integer :: i, j, k

    do k = 1, buffer%nz_local
      do j = 1, buffer%ny_local
        do i = 1, buffer%nx_local
          buffer%fx(i,j,k) = real(i + 2*j + 3*k, dp)
          buffer%fy(i,j,k) = real(2*i - j + k, dp)
          buffer%fz(i,j,k) = -real(i + j + k, dp)
        end do
      end do
    end do
  end subroutine fill_nonuniform_force_buffer

  subroutine setup_state(state, x0, structure_u0, sampled0, hydro0, input0, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(out) :: x0(3,3), structure_u0(3,3), sampled0(3,3), hydro0(3,3), input0(3,3)
    integer, intent(out) :: status

    x0(1,:) = [0.25_dp, 0.5_dp, 0.5_dp]
    x0(2,:) = [0.50_dp, 0.5_dp, 0.5_dp]
    x0(3,:) = [0.75_dp, 0.5_dp, 0.5_dp]
    structure_u0 = 0.0_dp
    sampled0(1,:) = [1.0_dp, 0.5_dp, -0.25_dp]
    sampled0(2,:) = [1.5_dp, 0.25_dp, -0.50_dp]
    sampled0(3,:) = [2.0_dp, 0.0_dp, -0.75_dp]
    hydro0 = 2.0_dp * sampled0
    input0 = hydro0
    call fibre_prod_state_allocate(state, 1, 3, status)
    if (status /= 0) return
    state%x(1,:,:) = x0
    call fibre_prod_state_attach_structure_u(state, structure_u0, status)
    if (status /= 0) return
    call fibre_prod_state_attach_sampled_velocity(state, sampled0, status)
    if (status /= 0) return
    call fibre_prod_state_attach_hydro_force_candidate(state, hydro0, status)
    if (status /= 0) return
    call fibre_prod_state_attach_structure_input_force(state, input0, status)
  end subroutine setup_state

  subroutine assert_fail_closed(gate, config, rhs0_x, rhs0_y, rhs0_z, buffer, label)
    type(fibre_prod_force_buffer_rhs_gate_type), intent(inout) :: gate
    type(fibre_prod_runtime_config_type), intent(in) :: config
    real(dp), intent(in) :: rhs0_x(:, :, :), rhs0_y(:, :, :), rhs0_z(:, :, :)
    type(fibre_prod_force_buffer_type), intent(in) :: buffer
    character(len=*), intent(in) :: label
    real(dp) :: local_x(size(rhs0_x,1), size(rhs0_x,2), size(rhs0_x,3))
    real(dp) :: local_y(size(rhs0_y,1), size(rhs0_y,2), size(rhs0_y,3))
    real(dp) :: local_z(size(rhs0_z,1), size(rhs0_z,2), size(rhs0_z,3))
    integer :: local_status

    local_x = rhs0_x; local_y = rhs0_y; local_z = rhs0_z
    call fibre_prod_force_buffer_rhs_gate_apply(gate, config, local_x, local_y, local_z, buffer, local_status)
    call require(local_status /= 0, trim(label)//' did not fail')
    call require(all(local_x == rhs0_x) .and. all(local_y == rhs0_y) .and. all(local_z == rhs0_z), &
                 trim(label)//' changed RHS')
  end subroutine assert_fail_closed

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      print *, 'P0_11_FORCE_BUFFER_RHS_GATE_CHECK FAIL: ', trim(message)
      error stop 1
    end if
  end subroutine require
end program fibre_prod_force_buffer_rhs_gate_check
