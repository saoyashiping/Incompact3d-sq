module fibre_prod_force_buffer_rhs_gate
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_validate, &
                                        fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init, fibre_prod_main_hook_apply_force_buffer
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_is_finite, &
                                          fibre_prod_force_buffer_allocate, fibre_prod_force_buffer_destroy
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity, fibre_prod_state_get_structure_coordinates, &
                              fibre_prod_state_get_structure_velocity_or_zero, fibre_prod_state_get_structure_input_force
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
                                            fibre_prod_structure_dry_step_init, fibre_prod_structure_dry_step_predict, &
                                            fibre_prod_structure_dry_step_finalize
  use fibre_prod_structure_commit_gate, only : fibre_prod_structure_commit_gate_type, &
                                               fibre_prod_structure_commit_gate_init, &
                                               fibre_prod_structure_commit_gate_set_enabled, &
                                               fibre_prod_structure_commit_gate_evaluate, &
                                               fibre_prod_structure_commit_gate_commit_to_state, &
                                               fibre_prod_structure_commit_gate_finalize
  use fibre_prod_reaction_force_candidate, only : fibre_prod_reaction_force_candidate_type, &
                                                  fibre_prod_reaction_force_candidate_init, &
                                                  fibre_prod_reaction_force_candidate_from_structure_input, &
                                                  fibre_prod_reaction_force_candidate_finalize
  use fibre_prod_reaction_spreading_buffer, only : fibre_prod_reaction_spreading_buffer_type, &
                                                   fibre_prod_reaction_spreading_buffer_init, &
                                                   fibre_prod_reaction_spreading_buffer_apply, &
                                                   fibre_prod_reaction_spreading_buffer_finalize
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_force_buffer_rhs_gate_type
  public :: fibre_prod_force_buffer_rhs_gate_init
  public :: fibre_prod_force_buffer_rhs_gate_apply
  public :: fibre_prod_force_buffer_rhs_gate_check_linear_response
  public :: fibre_prod_force_buffer_rhs_gate_finalize
  public :: fibre_prod_force_buffer_rhs_gate_env_enabled
  public :: fibre_prod_force_buffer_rhs_gate_runtime_diagnostic

  type :: fibre_prod_force_buffer_rhs_gate_type
    logical :: initialized = .false.
    logical :: applied = .false.
    logical :: lambda_zero_noop = .false.
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    real(dp) :: lambda_fsi = 0.0_dp
    real(dp) :: penalty_beta = 0.0_dp
    real(dp) :: max_abs_force_buffer = 0.0_dp
    real(dp) :: sum_abs_force_buffer = 0.0_dp
    real(dp) :: max_abs_rhs_increment = 0.0_dp
    real(dp) :: sum_abs_rhs_increment = 0.0_dp
    real(dp) :: expected_scale = 0.0_dp
    real(dp) :: measured_scale_error = 0.0_dp
  end type fibre_prod_force_buffer_rhs_gate_type

contains

  subroutine fibre_prod_force_buffer_rhs_gate_init(gate, nx, ny, nz, lambda_fsi, penalty_beta, status)
    type(fibre_prod_force_buffer_rhs_gate_type), intent(inout) :: gate
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    real(dp), intent(in) :: lambda_fsi
    real(dp), intent(in) :: penalty_beta
    integer, intent(out) :: status

    status = 0
    call fibre_prod_force_buffer_rhs_gate_finalize(gate)
    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
      status = 1
      return
    end if
    if (.not. ieee_is_finite(lambda_fsi) .or. lambda_fsi < 0.0_dp) then
      status = 2
      return
    end if
    if (.not. ieee_is_finite(penalty_beta)) then
      status = 3
      return
    end if
    gate%initialized = .true.
    gate%nx = nx
    gate%ny = ny
    gate%nz = nz
    gate%lambda_fsi = lambda_fsi
    gate%penalty_beta = penalty_beta
    gate%expected_scale = lambda_fsi * penalty_beta
    gate%lambda_zero_noop = lambda_fsi == 0.0_dp
  end subroutine fibre_prod_force_buffer_rhs_gate_init

  subroutine fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
    type(fibre_prod_force_buffer_rhs_gate_type), intent(inout) :: gate
    type(fibre_prod_runtime_config_type), intent(in) :: config
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    type(fibre_prod_force_buffer_type), intent(in) :: force_buffer
    integer, intent(out) :: status
    real(dp), allocatable :: before_x(:, :, :)
    real(dp), allocatable :: before_y(:, :, :)
    real(dp), allocatable :: before_z(:, :, :)
    real(dp), allocatable :: inc_x(:, :, :)
    real(dp), allocatable :: inc_y(:, :, :)
    real(dp), allocatable :: inc_z(:, :, :)
    real(dp) :: scale

    status = 0
    if (.not. gate%initialized) then
      status = 4
      return
    end if
    if (fibre_prod_runtime_config_validate(config) /= 0) then
      status = 5
      return
    end if
    if (.not. same_shape(rhs_x, rhs_y) .or. .not. same_shape(rhs_x, rhs_z)) then
      status = 6
      return
    end if
    if (size(rhs_x, 1) /= gate%nx .or. size(rhs_x, 2) /= gate%ny .or. size(rhs_x, 3) /= gate%nz) then
      status = 7
      return
    end if
    if (.not. force_buffer%allocated .or. .not. allocated(force_buffer%fx) .or. &
        .not. allocated(force_buffer%fy) .or. .not. allocated(force_buffer%fz)) then
      status = 8
      return
    end if
    if (.not. same_shape(rhs_x, force_buffer%fx) .or. .not. same_shape(rhs_y, force_buffer%fy) .or. &
        .not. same_shape(rhs_z, force_buffer%fz)) then
      status = 9
      return
    end if
    if (.not. fibre_prod_force_buffer_is_finite(force_buffer)) then
      status = 10
      return
    end if

    allocate(before_x(size(rhs_x,1), size(rhs_x,2), size(rhs_x,3)), &
             before_y(size(rhs_y,1), size(rhs_y,2), size(rhs_y,3)), &
             before_z(size(rhs_z,1), size(rhs_z,2), size(rhs_z,3)), &
             inc_x(size(rhs_x,1), size(rhs_x,2), size(rhs_x,3)), &
             inc_y(size(rhs_y,1), size(rhs_y,2), size(rhs_y,3)), &
             inc_z(size(rhs_z,1), size(rhs_z,2), size(rhs_z,3)), stat=status)
    if (status /= 0) then
      status = 11
      return
    end if
    before_x = rhs_x
    before_y = rhs_y
    before_z = rhs_z
    call fibre_prod_main_hook_init(status, config)
    if (status == 0) call fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, force_buffer, status)
    if (status /= 0) then
      rhs_x = before_x
      rhs_y = before_y
      rhs_z = before_z
      return
    end if

    scale = config%lambda_fsi * config%penalty_beta
    inc_x = rhs_x - before_x
    inc_y = rhs_y - before_y
    inc_z = rhs_z - before_z
    gate%max_abs_force_buffer = max(maxval(abs(force_buffer%fx)), &
                                    max(maxval(abs(force_buffer%fy)), maxval(abs(force_buffer%fz))))
    gate%sum_abs_force_buffer = sum(abs(force_buffer%fx)) + sum(abs(force_buffer%fy)) + sum(abs(force_buffer%fz))
    gate%max_abs_rhs_increment = max(maxval(abs(inc_x)), max(maxval(abs(inc_y)), maxval(abs(inc_z))))
    gate%sum_abs_rhs_increment = sum(abs(inc_x)) + sum(abs(inc_y)) + sum(abs(inc_z))
    gate%expected_scale = scale
    gate%lambda_fsi = config%lambda_fsi
    gate%penalty_beta = config%penalty_beta
    gate%lambda_zero_noop = config%lambda_fsi == 0.0_dp
    gate%measured_scale_error = max(maxval(abs(inc_x - scale * force_buffer%fx)), &
                                    max(maxval(abs(inc_y - scale * force_buffer%fy)), &
                                        maxval(abs(inc_z - scale * force_buffer%fz))))
    if (gate%lambda_zero_noop .and. gate%max_abs_rhs_increment /= 0.0_dp) status = 12
    if (status == 0 .and. gate%measured_scale_error > 1.0e-10_dp) status = 13
    if (status /= 0) then
      rhs_x = before_x
      rhs_y = before_y
      rhs_z = before_z
      return
    end if
    gate%applied = .true.
  end subroutine fibre_prod_force_buffer_rhs_gate_apply

  subroutine fibre_prod_force_buffer_rhs_gate_check_linear_response(gate, rhs_increment_1, rhs_increment_2, &
                                                                    lambda_1, lambda_2, status)
    type(fibre_prod_force_buffer_rhs_gate_type), intent(in) :: gate
    real(dp), intent(in) :: rhs_increment_1(:, :, :)
    real(dp), intent(in) :: rhs_increment_2(:, :, :)
    real(dp), intent(in) :: lambda_1
    real(dp), intent(in) :: lambda_2
    integer, intent(out) :: status
    real(dp) :: expected_ratio
    real(dp) :: max_error

    status = 0
    if (.not. gate%initialized) then
      status = 20
      return
    end if
    if (.not. same_shape(rhs_increment_1, rhs_increment_2)) then
      status = 21
      return
    end if
    if (.not. ieee_is_finite(lambda_1) .or. .not. ieee_is_finite(lambda_2) .or. &
        lambda_1 <= 0.0_dp .or. lambda_2 <= 0.0_dp) then
      status = 22
      return
    end if
    expected_ratio = lambda_2 / lambda_1
    max_error = maxval(abs(rhs_increment_2 - expected_ratio * rhs_increment_1))
    if (max_error > 1.0e-10_dp) status = 23
  end subroutine fibre_prod_force_buffer_rhs_gate_check_linear_response

  subroutine fibre_prod_force_buffer_rhs_gate_finalize(gate)
    type(fibre_prod_force_buffer_rhs_gate_type), intent(inout) :: gate

    gate%initialized = .false.
    gate%applied = .false.
    gate%lambda_zero_noop = .false.
    gate%nx = 0
    gate%ny = 0
    gate%nz = 0
    gate%lambda_fsi = 0.0_dp
    gate%penalty_beta = 0.0_dp
    gate%max_abs_force_buffer = 0.0_dp
    gate%sum_abs_force_buffer = 0.0_dp
    gate%max_abs_rhs_increment = 0.0_dp
    gate%sum_abs_rhs_increment = 0.0_dp
    gate%expected_scale = 0.0_dp
    gate%measured_scale_error = 0.0_dp
  end subroutine fibre_prod_force_buffer_rhs_gate_finalize

  logical function fibre_prod_force_buffer_rhs_gate_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_FORCE_BUFFER_RHS_GATE_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_force_buffer_rhs_gate_env_enabled

  subroutine fibre_prod_force_buffer_rhs_gate_runtime_diagnostic(ux, uy, uz, rhs_x, rhs_y, rhs_z, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_force_buffer_type) :: force_buffer
    type(fibre_prod_force_buffer_rhs_gate_type) :: gate
    type(fibre_prod_runtime_config_type) :: config
    real(dp), allocatable :: x(:), y(:), z(:)

    status = 0
    allocate(x(size(ux, 1)), y(size(ux, 2)), z(size(ux, 3)), stat=status)
    if (status /= 0) return
    call fill_unit_coordinates(x)
    call fill_unit_coordinates(y)
    call fill_unit_coordinates(z)
    call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, size(x), 1, size(y), 1, size(z), &
                                               .false., .false., .false., status)
    if (status == 0) call build_synthetic_force_buffer(grid, force_buffer, status)
    if (status == 0) call fibre_prod_runtime_config_default(config)
    if (status == 0) then
      config%enabled = .true.
      call read_env_real('FIBRE_PROD_LAMBDA', config%lambda_fsi, 0.0_dp)
      call read_env_real('FIBRE_PROD_PENALTY_BETA', config%penalty_beta, 1.0_dp)
      call fibre_prod_force_buffer_rhs_gate_init(gate, size(rhs_x, 1), size(rhs_x, 2), size(rhs_x, 3), &
                                                config%lambda_fsi, config%penalty_beta, status)
    end if
    if (status == 0) call fibre_prod_force_buffer_rhs_gate_apply(gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)
    call fibre_prod_force_buffer_rhs_gate_finalize(gate)
    call fibre_prod_force_buffer_destroy(force_buffer)
    call fibre_prod_grid_destroy(grid)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(z)) deallocate(z)
  end subroutine fibre_prod_force_buffer_rhs_gate_runtime_diagnostic

  subroutine build_synthetic_force_buffer(grid, force_buffer, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    type(fibre_prod_force_buffer_type), intent(inout) :: force_buffer
    integer, intent(out) :: status
    type(fibre_prod_state_type) :: state
    type(fibre_prod_hydro_input_candidate_type) :: hydro
    type(fibre_prod_structure_input_handoff_type) :: handoff
    type(fibre_prod_structure_dry_step_type) :: dry
    type(fibre_prod_structure_commit_gate_type) :: commit_gate
    type(fibre_prod_reaction_force_candidate_type) :: react
    type(fibre_prod_reaction_spreading_buffer_type) :: spread
    real(dp) :: points(3, 3)
    real(dp) :: sampled_u(3, 3)
    real(dp) :: structure_u(3, 3)
    real(dp) :: state_x(3, 3)
    real(dp) :: structure_input(3, 3)

    status = 0
    points(1, :) = [0.25_dp, 0.5_dp, 0.5_dp]
    points(2, :) = [0.50_dp, 0.5_dp, 0.5_dp]
    points(3, :) = [0.75_dp, 0.5_dp, 0.5_dp]
    sampled_u(1, :) = [1.0_dp, 0.5_dp, -0.25_dp]
    sampled_u(2, :) = [1.5_dp, 0.25_dp, -0.50_dp]
    sampled_u(3, :) = [2.0_dp, 0.0_dp, -0.75_dp]
    call fibre_prod_state_allocate(state, 1, 3, status)
    if (status == 0) state%x(1, :, :) = points
    if (status == 0) call fibre_prod_state_attach_sampled_velocity(state, sampled_u, status)
    if (status == 0) call fibre_prod_state_get_structure_velocity_or_zero(state, structure_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_init(hydro, 3, 2.5_dp, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_compute(hydro, state%sampled_u, structure_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_attach_to_state(hydro, state, status)
    if (status == 0) call fibre_prod_structure_input_handoff_init(handoff, 3, status)
    if (status == 0) call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status == 0) call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    if (status == 0) call fibre_prod_state_get_structure_coordinates(state, state_x, status)
    if (status == 0) call fibre_prod_structure_dry_step_init(dry, 3, 1.0e-4_dp, 2.0_dp, status)
    if (status == 0) call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
    if (status == 0) call fibre_prod_structure_commit_gate_init(commit_gate, 3, 1.0e-2_dp, status)
    if (status == 0) call fibre_prod_structure_commit_gate_set_enabled(commit_gate, .true., status)
    if (status == 0) call fibre_prod_structure_commit_gate_evaluate(commit_gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
    if (status == 0) call fibre_prod_structure_commit_gate_commit_to_state(commit_gate, state, dry%x_trial, dry%u_trial, status)
    if (status == 0) call fibre_prod_state_get_structure_input_force(state, structure_input, status)
    if (status == 0) call fibre_prod_reaction_force_candidate_init(react, 3, status)
    if (status == 0) call fibre_prod_reaction_force_candidate_from_structure_input(react, structure_input, status)
    if (status == 0) call fibre_prod_force_buffer_allocate(force_buffer, grid, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_init(spread, grid%nx_local, grid%ny_local, grid%nz_local, 3, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_apply(spread, grid, state%x(1, :, :), react%reaction_force, &
                                                                    force_buffer, status)
    call fibre_prod_reaction_spreading_buffer_finalize(spread)
    call fibre_prod_reaction_force_candidate_finalize(react)
    call fibre_prod_structure_commit_gate_finalize(commit_gate)
    call fibre_prod_structure_dry_step_finalize(dry)
    call fibre_prod_structure_input_handoff_finalize(handoff)
    call fibre_prod_hydro_input_candidate_finalize(hydro)
    call fibre_prod_state_destroy(state)
  end subroutine build_synthetic_force_buffer

  subroutine fill_unit_coordinates(coord)
    real(dp), intent(out) :: coord(:)
    integer :: i

    do i = 1, size(coord)
      coord(i) = real(i - 1, dp) / real(max(1, size(coord) - 1), dp)
    end do
  end subroutine fill_unit_coordinates

  subroutine read_env_real(name, value, default_value)
    character(len=*), intent(in) :: name
    real(dp), intent(inout) :: value
    real(dp), intent(in) :: default_value
    character(len=128) :: raw
    integer :: length
    integer :: env_status
    integer :: read_status
    real(dp) :: parsed

    value = default_value
    call get_environment_variable(name, raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    read(raw(1:min(length, len(raw))), *, iostat=read_status) parsed
    if (read_status == 0 .and. ieee_is_finite(parsed)) value = parsed
  end subroutine read_env_real

  pure logical function same_shape(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = size(lhs, 1) == size(rhs, 1) .and. size(lhs, 2) == size(rhs, 2) .and. &
              size(lhs, 3) == size(rhs, 3)
  end function same_shape

end module fibre_prod_force_buffer_rhs_gate
