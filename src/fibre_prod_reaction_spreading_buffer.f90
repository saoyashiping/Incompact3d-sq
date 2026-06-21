module fibre_prod_reaction_spreading_buffer
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_is_initialized, fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_is_finite, &
                                          fibre_prod_force_buffer_reset_to_zero, fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_destroy
  use fibre_prod_ibm_spreading, only : fibre_prod_spread_multiple_point_forces
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity, fibre_prod_state_get_structure_coordinates, &
                              fibre_prod_state_get_structure_velocity_or_zero, &
                              fibre_prod_state_get_structure_input_force
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
                                                  fibre_prod_reaction_force_candidate_finalize, &
                                                  fibre_prod_reaction_force_candidate_env_enabled
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_reaction_spreading_buffer_type
  public :: fibre_prod_reaction_spreading_buffer_init
  public :: fibre_prod_reaction_spreading_buffer_apply
  public :: fibre_prod_reaction_spreading_buffer_check_finite_bounded
  public :: fibre_prod_reaction_spreading_buffer_finalize
  public :: fibre_prod_reaction_spreading_buffer_env_enabled
  public :: fibre_prod_reaction_spreading_buffer_env_max_allowed_force
  public :: fibre_prod_reaction_spreading_buffer_runtime_diagnostic

  type :: fibre_prod_reaction_spreading_buffer_type
    logical :: initialized = .false.
    logical :: has_spread = .false.
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: nnode = 0
    real(dp) :: max_abs_force_buffer = 0.0_dp
    real(dp) :: sum_abs_force_buffer = 0.0_dp
    real(dp) :: net_lagrangian_force(3) = 0.0_dp
    real(dp) :: net_eulerian_force(3) = 0.0_dp
    real(dp) :: conservation_error(3) = 0.0_dp
  end type fibre_prod_reaction_spreading_buffer_type

contains

  subroutine fibre_prod_reaction_spreading_buffer_init(spread, nx, ny, nz, nnode, status)
    type(fibre_prod_reaction_spreading_buffer_type), intent(inout) :: spread
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer, intent(in) :: nnode
    integer, intent(out) :: status

    status = 0
    call fibre_prod_reaction_spreading_buffer_finalize(spread)
    if (nx <= 0 .or. ny <= 0 .or. nz <= 0 .or. nnode <= 0) then
      status = 1
      return
    end if
    spread%nx = nx
    spread%ny = ny
    spread%nz = nz
    spread%nnode = nnode
    spread%initialized = .true.
  end subroutine fibre_prod_reaction_spreading_buffer_init

  subroutine fibre_prod_reaction_spreading_buffer_apply(spread, grid, fibre_points, reaction_force, force_buffer, status)
    type(fibre_prod_reaction_spreading_buffer_type), intent(inout) :: spread
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: fibre_points(:, :)
    real(dp), intent(in) :: reaction_force(:, :)
    type(fibre_prod_force_buffer_type), intent(inout) :: force_buffer
    integer, intent(out) :: status
    real(dp) :: min_weight_sum

    status = 0
    if (.not. spread%initialized) then
      status = 2
      return
    end if
    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 3
      return
    end if
    if (size(fibre_points, 1) /= spread%nnode .or. size(fibre_points, 2) /= 3) then
      status = 4
      return
    end if
    if (size(reaction_force, 1) /= spread%nnode .or. size(reaction_force, 2) /= 3) then
      status = 5
      return
    end if
    if (.not. all(ieee_is_finite(fibre_points)) .or. .not. all(ieee_is_finite(reaction_force))) then
      status = 6
      return
    end if
    if (.not. force_buffer%allocated .or. .not. allocated(force_buffer%fx) .or. &
        .not. allocated(force_buffer%fy) .or. .not. allocated(force_buffer%fz)) then
      status = 7
      return
    end if
    if (force_buffer%nx_local /= spread%nx .or. force_buffer%ny_local /= spread%ny .or. &
        force_buffer%nz_local /= spread%nz) then
      status = 8
      return
    end if
    call fibre_prod_force_buffer_reset_to_zero(force_buffer, status)
    if (status /= 0) return
    call fibre_prod_spread_multiple_point_forces(grid, force_buffer, fibre_points, reaction_force, status, min_weight_sum)
    if (status /= 0) return
    if (.not. fibre_prod_force_buffer_is_finite(force_buffer)) then
      status = 9
      return
    end if
    spread%net_lagrangian_force = sum(reaction_force, dim=1)
    call compute_eulerian_net_force(grid, force_buffer, spread%net_eulerian_force)
    spread%conservation_error = spread%net_eulerian_force - spread%net_lagrangian_force
    spread%max_abs_force_buffer = max(maxval(abs(force_buffer%fx)), &
                                      max(maxval(abs(force_buffer%fy)), maxval(abs(force_buffer%fz))))
    spread%sum_abs_force_buffer = sum(abs(force_buffer%fx)) + sum(abs(force_buffer%fy)) + sum(abs(force_buffer%fz))
    spread%has_spread = .true.
  end subroutine fibre_prod_reaction_spreading_buffer_apply

  subroutine fibre_prod_reaction_spreading_buffer_check_finite_bounded(spread, max_allowed_force, status)
    type(fibre_prod_reaction_spreading_buffer_type), intent(in) :: spread
    real(dp), intent(in) :: max_allowed_force
    integer, intent(out) :: status

    status = 0
    if (.not. spread%initialized .or. .not. spread%has_spread) then
      status = 10
      return
    end if
    if (.not. ieee_is_finite(max_allowed_force) .or. max_allowed_force <= 0.0_dp) then
      status = 11
      return
    end if
    if (.not. ieee_is_finite(spread%max_abs_force_buffer) .or. &
        .not. ieee_is_finite(spread%sum_abs_force_buffer) .or. &
        .not. all(ieee_is_finite(spread%net_lagrangian_force)) .or. &
        .not. all(ieee_is_finite(spread%net_eulerian_force)) .or. &
        .not. all(ieee_is_finite(spread%conservation_error))) then
      status = 12
      return
    end if
    if (spread%max_abs_force_buffer > max_allowed_force) status = 13
  end subroutine fibre_prod_reaction_spreading_buffer_check_finite_bounded

  subroutine fibre_prod_reaction_spreading_buffer_finalize(spread)
    type(fibre_prod_reaction_spreading_buffer_type), intent(inout) :: spread

    spread%initialized = .false.
    spread%has_spread = .false.
    spread%nx = 0
    spread%ny = 0
    spread%nz = 0
    spread%nnode = 0
    spread%max_abs_force_buffer = 0.0_dp
    spread%sum_abs_force_buffer = 0.0_dp
    spread%net_lagrangian_force = 0.0_dp
    spread%net_eulerian_force = 0.0_dp
    spread%conservation_error = 0.0_dp
  end subroutine fibre_prod_reaction_spreading_buffer_finalize

  logical function fibre_prod_reaction_spreading_buffer_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_REACTION_SPREADING_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_reaction_spreading_buffer_env_enabled

  real(dp) function fibre_prod_reaction_spreading_buffer_env_max_allowed_force() result(value)
    value = read_positive_env('FIBRE_PROD_REACTION_SPREADING_MAX_FORCE', 1.0e6_dp)
  end function fibre_prod_reaction_spreading_buffer_env_max_allowed_force

  subroutine fibre_prod_reaction_spreading_buffer_runtime_diagnostic(ux, uy, uz, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_state_type) :: state
    type(fibre_prod_hydro_input_candidate_type) :: hydro
    type(fibre_prod_structure_input_handoff_type) :: handoff
    type(fibre_prod_structure_dry_step_type) :: dry
    type(fibre_prod_structure_commit_gate_type) :: gate
    type(fibre_prod_reaction_force_candidate_type) :: react
    type(fibre_prod_reaction_spreading_buffer_type) :: spread
    type(fibre_prod_force_buffer_type) :: buffer
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: points(3, 3)
    real(dp) :: sampled_u(3, 3)
    real(dp) :: structure_u(3, 3)
    real(dp) :: state_x(3, 3)
    real(dp) :: structure_input(3, 3)

    status = 0
    if (.not. fibre_prod_reaction_force_candidate_env_enabled()) return
    allocate(x(size(ux, 1)), y(size(ux, 2)), z(size(ux, 3)), stat=status)
    if (status /= 0) return
    call fill_unit_coordinates(x)
    call fill_unit_coordinates(y)
    call fill_unit_coordinates(z)
    call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, size(x), 1, size(y), 1, size(z), &
                                               .false., .false., .false., status)
    points(1, :) = [0.25_dp, 0.5_dp, 0.5_dp]
    points(2, :) = [0.50_dp, 0.5_dp, 0.5_dp]
    points(3, :) = [0.75_dp, 0.5_dp, 0.5_dp]
    sampled_u(1, :) = [1.0_dp, 0.5_dp, -0.25_dp]
    sampled_u(2, :) = [1.5_dp, 0.25_dp, -0.50_dp]
    sampled_u(3, :) = [2.0_dp, 0.0_dp, -0.75_dp]
    if (status == 0) call fibre_prod_state_allocate(state, 1, 3, status)
    if (status == 0) state%x(1, :, :) = points
    if (status == 0) call fibre_prod_state_attach_sampled_velocity(state, sampled_u, status)
    if (status == 0) call fibre_prod_state_get_structure_velocity_or_zero(state, structure_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_init(hydro, 3, 2.0_dp, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_compute(hydro, state%sampled_u, structure_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_attach_to_state(hydro, state, status)
    if (status == 0) call fibre_prod_structure_input_handoff_init(handoff, 3, status)
    if (status == 0) call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status == 0) call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    if (status == 0) call fibre_prod_state_get_structure_coordinates(state, state_x, status)
    if (status == 0) call fibre_prod_structure_dry_step_init(dry, 3, 1.0e-4_dp, 2.0_dp, status)
    if (status == 0) call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
    if (status == 0) call fibre_prod_structure_commit_gate_init(gate, 3, 1.0e-2_dp, status)
    if (status == 0) call fibre_prod_structure_commit_gate_set_enabled(gate, .true., status)
    if (status == 0) call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
    if (status == 0) call fibre_prod_structure_commit_gate_commit_to_state(gate, state, dry%x_trial, dry%u_trial, status)
    if (status == 0) call fibre_prod_state_get_structure_input_force(state, structure_input, status)
    if (status == 0) call fibre_prod_reaction_force_candidate_init(react, 3, status)
    if (status == 0) call fibre_prod_reaction_force_candidate_from_structure_input(react, structure_input, status)
    if (status == 0) call fibre_prod_force_buffer_allocate(buffer, grid, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_init(spread, grid%nx_local, grid%ny_local, grid%nz_local, 3, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_apply(spread, grid, state%x(1, :, :), react%reaction_force, &
                                                                    buffer, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_check_finite_bounded(spread, &
                                                                                   fibre_prod_reaction_spreading_buffer_env_max_allowed_force(), status)
    call fibre_prod_force_buffer_destroy(buffer)
    call fibre_prod_reaction_spreading_buffer_finalize(spread)
    call fibre_prod_reaction_force_candidate_finalize(react)
    call fibre_prod_structure_commit_gate_finalize(gate)
    call fibre_prod_structure_dry_step_finalize(dry)
    call fibre_prod_structure_input_handoff_finalize(handoff)
    call fibre_prod_hydro_input_candidate_finalize(hydro)
    call fibre_prod_state_destroy(state)
    call fibre_prod_grid_destroy(grid)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(z)) deallocate(z)
  contains
    subroutine fill_unit_coordinates(coord)
      real(dp), intent(out) :: coord(:)
      integer :: i

      do i = 1, size(coord)
        coord(i) = real(i - 1, dp) / real(max(1, size(coord) - 1), dp)
      end do
    end subroutine fill_unit_coordinates
  end subroutine fibre_prod_reaction_spreading_buffer_runtime_diagnostic

  subroutine compute_eulerian_net_force(grid, force_buffer, net_force)
    type(fibre_prod_grid_type), intent(in) :: grid
    type(fibre_prod_force_buffer_type), intent(in) :: force_buffer
    real(dp), intent(out) :: net_force(3)

    net_force(1) = sum(force_buffer%fx * grid%cell_volume)
    net_force(2) = sum(force_buffer%fy * grid%cell_volume)
    net_force(3) = sum(force_buffer%fz * grid%cell_volume)
  end subroutine compute_eulerian_net_force

  real(dp) function read_positive_env(name, default_value) result(value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: default_value
    character(len=64) :: raw
    integer :: length
    integer :: env_status
    integer :: read_status

    value = default_value
    call get_environment_variable(name, raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    read(raw(1:min(length, len(raw))), *, iostat=read_status) value
    if (read_status /= 0 .or. .not. ieee_is_finite(value) .or. value <= 0.0_dp) value = default_value
  end function read_positive_env

end module fibre_prod_reaction_spreading_buffer
