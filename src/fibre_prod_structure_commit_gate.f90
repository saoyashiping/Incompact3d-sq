module fibre_prod_structure_commit_gate
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_get_structure_coordinates, &
                              fibre_prod_state_get_structure_velocity_or_zero, &
                              fibre_prod_state_commit_structure_trial
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
                                               fibre_prod_hydro_input_candidate_finalize, &
                                               fibre_prod_hydro_input_candidate_env_beta
  use fibre_prod_structure_input_handoff, only : fibre_prod_structure_input_handoff_type, &
                                                 fibre_prod_structure_input_handoff_init, &
                                                 fibre_prod_structure_input_handoff_from_candidate, &
                                                 fibre_prod_structure_input_handoff_attach_to_state, &
                                                 fibre_prod_structure_input_handoff_finalize
  use fibre_prod_structure_dry_step, only : fibre_prod_structure_dry_step_type, &
                                            fibre_prod_structure_dry_step_init, &
                                            fibre_prod_structure_dry_step_predict, &
                                            fibre_prod_structure_dry_step_check_bounded, &
                                            fibre_prod_structure_dry_step_finalize, &
                                            fibre_prod_structure_dry_step_env_dt, &
                                            fibre_prod_structure_dry_step_env_rho_eff
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: fibre_prod_commit_reject_none = 0
  integer, parameter, public :: fibre_prod_commit_reject_disabled = 1
  integer, parameter, public :: fibre_prod_commit_reject_not_initialized = 2
  integer, parameter, public :: fibre_prod_commit_reject_shape = 3
  integer, parameter, public :: fibre_prod_commit_reject_nonfinite = 4
  integer, parameter, public :: fibre_prod_commit_reject_bound = 5
  integer, parameter, public :: fibre_prod_commit_reject_state = 6

  public :: fibre_prod_structure_commit_gate_type
  public :: fibre_prod_structure_commit_gate_init
  public :: fibre_prod_structure_commit_gate_set_enabled
  public :: fibre_prod_structure_commit_gate_evaluate
  public :: fibre_prod_structure_commit_gate_commit_to_state
  public :: fibre_prod_structure_commit_gate_finalize
  public :: fibre_prod_structure_commit_gate_env_enabled
  public :: fibre_prod_structure_commit_gate_env_max_allowed_dx
  public :: fibre_prod_structure_commit_gate_runtime_diagnostic

  type :: fibre_prod_structure_commit_gate_type
    logical :: initialized = .false.
    logical :: gate_enabled = .false.
    logical :: accepted = .false.
    logical :: committed = .false.
    logical :: rejected = .false.
    integer :: reject_code = fibre_prod_commit_reject_not_initialized
    integer :: nnode = 0
    real(dp) :: max_allowed_dx = 0.0_dp
    real(dp) :: max_abs_dx_trial = 0.0_dp
    real(dp) :: max_abs_x_trial = 0.0_dp
    real(dp) :: max_abs_u_trial = 0.0_dp
    real(dp) :: sum_abs_dx_trial = 0.0_dp
  end type fibre_prod_structure_commit_gate_type

contains

  subroutine fibre_prod_structure_commit_gate_init(gate, nnode, max_allowed_dx, status)
    type(fibre_prod_structure_commit_gate_type), intent(inout) :: gate
    integer, intent(in) :: nnode
    real(dp), intent(in) :: max_allowed_dx
    integer, intent(out) :: status

    status = 0
    call fibre_prod_structure_commit_gate_finalize(gate)
    if (nnode <= 0) then
      status = 1
      return
    end if
    if (.not. ieee_is_finite(max_allowed_dx) .or. max_allowed_dx <= 0.0_dp) then
      status = 2
      return
    end if
    gate%initialized = .true.
    gate%gate_enabled = .false.
    gate%accepted = .false.
    gate%committed = .false.
    gate%rejected = .false.
    gate%reject_code = fibre_prod_commit_reject_none
    gate%nnode = nnode
    gate%max_allowed_dx = max_allowed_dx
  end subroutine fibre_prod_structure_commit_gate_init

  subroutine fibre_prod_structure_commit_gate_set_enabled(gate, enabled, status)
    type(fibre_prod_structure_commit_gate_type), intent(inout) :: gate
    logical, intent(in) :: enabled
    integer, intent(out) :: status

    status = 0
    if (.not. gate%initialized) then
      status = fibre_prod_commit_reject_not_initialized
      gate%rejected = .true.
      gate%accepted = .false.
      gate%reject_code = fibre_prod_commit_reject_not_initialized
      return
    end if
    gate%gate_enabled = enabled
    gate%accepted = .false.
    gate%committed = .false.
    gate%rejected = .false.
    gate%reject_code = fibre_prod_commit_reject_none
  end subroutine fibre_prod_structure_commit_gate_set_enabled

  subroutine fibre_prod_structure_commit_gate_evaluate(gate, x_trial, u_trial, dx_trial, status)
    type(fibre_prod_structure_commit_gate_type), intent(inout) :: gate
    real(dp), intent(in) :: x_trial(:, :)
    real(dp), intent(in) :: u_trial(:, :)
    real(dp), intent(in) :: dx_trial(:, :)
    integer, intent(out) :: status

    status = 0
    gate%accepted = .false.
    gate%committed = .false.
    gate%rejected = .false.
    if (.not. gate%initialized) then
      status = fibre_prod_commit_reject_not_initialized
      gate%rejected = .true.
      gate%reject_code = fibre_prod_commit_reject_not_initialized
      return
    end if
    if (size(x_trial, 1) /= gate%nnode .or. size(x_trial, 2) /= 3 .or. &
        size(u_trial, 1) /= gate%nnode .or. size(u_trial, 2) /= 3 .or. &
        size(dx_trial, 1) /= gate%nnode .or. size(dx_trial, 2) /= 3) then
      status = fibre_prod_commit_reject_shape
      gate%rejected = .true.
      gate%reject_code = fibre_prod_commit_reject_shape
      return
    end if
    if (.not. all(ieee_is_finite(x_trial)) .or. .not. all(ieee_is_finite(u_trial)) .or. &
        .not. all(ieee_is_finite(dx_trial))) then
      status = fibre_prod_commit_reject_nonfinite
      gate%rejected = .true.
      gate%reject_code = fibre_prod_commit_reject_nonfinite
      return
    end if
    gate%max_abs_dx_trial = maxval(abs(dx_trial))
    gate%max_abs_x_trial = maxval(abs(x_trial))
    gate%max_abs_u_trial = maxval(abs(u_trial))
    gate%sum_abs_dx_trial = sum(abs(dx_trial))
    if (gate%max_abs_dx_trial > gate%max_allowed_dx) then
      status = fibre_prod_commit_reject_bound
      gate%rejected = .true.
      gate%reject_code = fibre_prod_commit_reject_bound
      return
    end if
    if (.not. gate%gate_enabled) then
      gate%rejected = .true.
      gate%reject_code = fibre_prod_commit_reject_disabled
      return
    end if
    gate%accepted = .true.
    gate%rejected = .false.
    gate%reject_code = fibre_prod_commit_reject_none
  end subroutine fibre_prod_structure_commit_gate_evaluate

  subroutine fibre_prod_structure_commit_gate_commit_to_state(gate, state, x_trial, u_trial, status)
    type(fibre_prod_structure_commit_gate_type), intent(inout) :: gate
    type(fibre_prod_state_type), intent(inout) :: state
    real(dp), intent(in) :: x_trial(:, :)
    real(dp), intent(in) :: u_trial(:, :)
    integer, intent(out) :: status

    status = 0
    gate%committed = .false.
    if (.not. gate%initialized) then
      status = fibre_prod_commit_reject_not_initialized
      gate%rejected = .true.
      gate%reject_code = fibre_prod_commit_reject_not_initialized
      return
    end if
    if (.not. gate%gate_enabled .or. .not. gate%accepted .or. gate%rejected) then
      if (gate%reject_code == fibre_prod_commit_reject_none) gate%reject_code = fibre_prod_commit_reject_disabled
      return
    end if
    if (size(x_trial, 1) /= gate%nnode .or. size(x_trial, 2) /= 3 .or. &
        size(u_trial, 1) /= gate%nnode .or. size(u_trial, 2) /= 3) then
      status = fibre_prod_commit_reject_shape
      gate%rejected = .true.
      gate%accepted = .false.
      gate%reject_code = fibre_prod_commit_reject_shape
      return
    end if
    call fibre_prod_state_commit_structure_trial(state, x_trial, u_trial, status)
    if (status /= 0) then
      status = fibre_prod_commit_reject_state
      gate%rejected = .true.
      gate%accepted = .false.
      gate%reject_code = fibre_prod_commit_reject_state
      return
    end if
    gate%committed = .true.
  end subroutine fibre_prod_structure_commit_gate_commit_to_state

  subroutine fibre_prod_structure_commit_gate_finalize(gate)
    type(fibre_prod_structure_commit_gate_type), intent(inout) :: gate

    gate%initialized = .false.
    gate%gate_enabled = .false.
    gate%accepted = .false.
    gate%committed = .false.
    gate%rejected = .false.
    gate%reject_code = fibre_prod_commit_reject_not_initialized
    gate%nnode = 0
    gate%max_allowed_dx = 0.0_dp
    gate%max_abs_dx_trial = 0.0_dp
    gate%max_abs_x_trial = 0.0_dp
    gate%max_abs_u_trial = 0.0_dp
    gate%sum_abs_dx_trial = 0.0_dp
  end subroutine fibre_prod_structure_commit_gate_finalize

  logical function fibre_prod_structure_commit_gate_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_STRUCTURE_DRY_COMMIT_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_structure_commit_gate_env_enabled

  real(dp) function fibre_prod_structure_commit_gate_env_max_allowed_dx() result(value)
    value = read_positive_env('FIBRE_PROD_STRUCTURE_DRY_COMMIT_MAX_DX', 1.0e-2_dp)
  end function fibre_prod_structure_commit_gate_env_max_allowed_dx

  subroutine fibre_prod_structure_commit_gate_runtime_diagnostic(ux, uy, uz, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_state_type) :: state
    type(fibre_prod_state_velocity_attachment_type) :: attach
    type(fibre_prod_hydro_input_candidate_type) :: candidate
    type(fibre_prod_structure_input_handoff_type) :: handoff
    type(fibre_prod_structure_dry_step_type) :: dry
    type(fibre_prod_structure_commit_gate_type) :: gate
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: points(3, 3)
    real(dp) :: state_x(3, 3)
    real(dp) :: state_u(3, 3)

    status = 0
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
    if (status == 0) call fibre_prod_state_allocate(state, 1, 3, status)
    if (status == 0) state%x(1, :, :) = points
    if (status == 0) call fibre_prod_state_velocity_attachment_init(attach, 3, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_set_points(attach, points, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_init(candidate, 3, &
                                                               fibre_prod_hydro_input_candidate_env_beta(), status)
    if (status == 0) call fibre_prod_state_get_structure_velocity_or_zero(state, state_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_compute(candidate, state%sampled_u, state_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_attach_to_state(candidate, state, status)
    if (status == 0) call fibre_prod_structure_input_handoff_init(handoff, 3, status)
    if (status == 0) call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status == 0) call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    if (status == 0) call fibre_prod_state_get_structure_coordinates(state, state_x, status)
    if (status == 0) call fibre_prod_state_get_structure_velocity_or_zero(state, state_u, status)
    if (status == 0) call fibre_prod_structure_dry_step_init(dry, 3, fibre_prod_structure_dry_step_env_dt(), &
                                                            fibre_prod_structure_dry_step_env_rho_eff(), status)
    if (status == 0) call fibre_prod_structure_dry_step_predict(dry, state_x, state_u, state%structure_input_force, status)
    if (status == 0) call fibre_prod_structure_dry_step_check_bounded(dry, &
                                                                      fibre_prod_structure_commit_gate_env_max_allowed_dx(), status)
    if (status == 0) call fibre_prod_structure_commit_gate_init(gate, 3, &
                                                               fibre_prod_structure_commit_gate_env_max_allowed_dx(), status)
    if (status == 0) call fibre_prod_structure_commit_gate_set_enabled(gate, &
                                                                       fibre_prod_structure_commit_gate_env_enabled(), status)
    if (status == 0) call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
    if (status == 0) call fibre_prod_structure_commit_gate_commit_to_state(gate, state, dry%x_trial, dry%u_trial, status)
    call fibre_prod_structure_commit_gate_finalize(gate)
    call fibre_prod_structure_dry_step_finalize(dry)
    call fibre_prod_structure_input_handoff_finalize(handoff)
    call fibre_prod_hydro_input_candidate_finalize(candidate)
    call fibre_prod_state_velocity_attachment_finalize(attach)
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
  end subroutine fibre_prod_structure_commit_gate_runtime_diagnostic

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

end module fibre_prod_structure_commit_gate
