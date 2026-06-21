module fibre_prod_structure_dry_step
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_get_structure_coordinates, &
                              fibre_prod_state_get_structure_velocity_or_zero
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
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_structure_dry_step_type
  public :: fibre_prod_structure_dry_step_init
  public :: fibre_prod_structure_dry_step_predict
  public :: fibre_prod_structure_dry_step_check_bounded
  public :: fibre_prod_structure_dry_step_finalize
  public :: fibre_prod_structure_dry_step_env_enabled
  public :: fibre_prod_structure_dry_step_env_dt
  public :: fibre_prod_structure_dry_step_env_rho_eff
  public :: fibre_prod_structure_dry_step_env_max_allowed_dx
  public :: fibre_prod_structure_dry_step_runtime_diagnostic

  type :: fibre_prod_structure_dry_step_type
    logical :: initialized = .false.
    logical :: has_trial = .false.
    integer :: nnode = 0
    real(dp) :: dt = 0.0_dp
    real(dp) :: rho_eff = 1.0_dp
    real(dp), allocatable :: x_trial(:, :)
    real(dp), allocatable :: u_trial(:, :)
    real(dp), allocatable :: dx_trial(:, :)
    real(dp) :: max_abs_dx_trial = 0.0_dp
    real(dp) :: max_abs_u_trial = 0.0_dp
    real(dp) :: sum_abs_dx_trial = 0.0_dp
  end type fibre_prod_structure_dry_step_type

contains

  subroutine fibre_prod_structure_dry_step_init(dry, nnode, dt, rho_eff, status)
    type(fibre_prod_structure_dry_step_type), intent(inout) :: dry
    integer, intent(in) :: nnode
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: rho_eff
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    call fibre_prod_structure_dry_step_finalize(dry)
    if (nnode <= 0) then
      status = 1
      return
    end if
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      status = 2
      return
    end if
    if (.not. ieee_is_finite(rho_eff) .or. rho_eff <= 0.0_dp) then
      status = 3
      return
    end if
    allocate(dry%x_trial(nnode, 3), dry%u_trial(nnode, 3), dry%dx_trial(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 4
      call fibre_prod_structure_dry_step_finalize(dry)
      return
    end if
    dry%nnode = nnode
    dry%dt = dt
    dry%rho_eff = rho_eff
    dry%x_trial = 0.0_dp
    dry%u_trial = 0.0_dp
    dry%dx_trial = 0.0_dp
    dry%max_abs_dx_trial = 0.0_dp
    dry%max_abs_u_trial = 0.0_dp
    dry%sum_abs_dx_trial = 0.0_dp
    dry%has_trial = .false.
    dry%initialized = .true.
  end subroutine fibre_prod_structure_dry_step_init

  subroutine fibre_prod_structure_dry_step_predict(dry, x, structure_u, structure_input_force, status)
    type(fibre_prod_structure_dry_step_type), intent(inout) :: dry
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: structure_u(:, :)
    real(dp), intent(in) :: structure_input_force(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. dry%initialized .or. .not. allocated(dry%x_trial) .or. .not. allocated(dry%u_trial) .or. &
        .not. allocated(dry%dx_trial)) then
      status = 5
      return
    end if
    if (size(x, 1) /= dry%nnode .or. size(x, 2) /= 3) then
      status = 6
      return
    end if
    if (size(structure_u, 1) /= dry%nnode .or. size(structure_u, 2) /= 3) then
      status = 7
      return
    end if
    if (size(structure_input_force, 1) /= dry%nnode .or. size(structure_input_force, 2) /= 3) then
      status = 8
      return
    end if
    if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(structure_u)) .or. &
        .not. all(ieee_is_finite(structure_input_force))) then
      status = 9
      return
    end if

    dry%dx_trial = dry%dt * structure_u + 0.5_dp * dry%dt * dry%dt * structure_input_force / dry%rho_eff
    dry%x_trial = x + dry%dx_trial
    dry%u_trial = structure_u + dry%dt * structure_input_force / dry%rho_eff
    dry%max_abs_dx_trial = maxval(abs(dry%dx_trial))
    dry%max_abs_u_trial = maxval(abs(dry%u_trial))
    dry%sum_abs_dx_trial = sum(abs(dry%dx_trial))
    dry%has_trial = .true.
  end subroutine fibre_prod_structure_dry_step_predict

  subroutine fibre_prod_structure_dry_step_check_bounded(dry, max_allowed_dx, status)
    type(fibre_prod_structure_dry_step_type), intent(in) :: dry
    real(dp), intent(in) :: max_allowed_dx
    integer, intent(out) :: status

    status = 0
    if (.not. dry%initialized .or. .not. dry%has_trial .or. .not. allocated(dry%x_trial) .or. &
        .not. allocated(dry%u_trial) .or. .not. allocated(dry%dx_trial)) then
      status = 10
      return
    end if
    if (.not. ieee_is_finite(max_allowed_dx) .or. max_allowed_dx <= 0.0_dp) then
      status = 11
      return
    end if
    if (.not. all(ieee_is_finite(dry%x_trial)) .or. .not. all(ieee_is_finite(dry%u_trial)) .or. &
        .not. all(ieee_is_finite(dry%dx_trial))) then
      status = 12
      return
    end if
    if (dry%max_abs_dx_trial > max_allowed_dx) status = 13
  end subroutine fibre_prod_structure_dry_step_check_bounded

  subroutine fibre_prod_structure_dry_step_finalize(dry)
    type(fibre_prod_structure_dry_step_type), intent(inout) :: dry

    if (allocated(dry%x_trial)) deallocate(dry%x_trial)
    if (allocated(dry%u_trial)) deallocate(dry%u_trial)
    if (allocated(dry%dx_trial)) deallocate(dry%dx_trial)
    dry%initialized = .false.
    dry%has_trial = .false.
    dry%nnode = 0
    dry%dt = 0.0_dp
    dry%rho_eff = 1.0_dp
    dry%max_abs_dx_trial = 0.0_dp
    dry%max_abs_u_trial = 0.0_dp
    dry%sum_abs_dx_trial = 0.0_dp
  end subroutine fibre_prod_structure_dry_step_finalize

  logical function fibre_prod_structure_dry_step_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_STRUCTURE_DRY_STEP_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_structure_dry_step_env_enabled

  real(dp) function fibre_prod_structure_dry_step_env_dt() result(value)
    value = read_positive_env('FIBRE_PROD_STRUCTURE_DRY_DT', 1.0e-4_dp)
  end function fibre_prod_structure_dry_step_env_dt

  real(dp) function fibre_prod_structure_dry_step_env_rho_eff() result(value)
    value = read_positive_env('FIBRE_PROD_STRUCTURE_DRY_RHO_EFF', 1.0_dp)
  end function fibre_prod_structure_dry_step_env_rho_eff

  real(dp) function fibre_prod_structure_dry_step_env_max_allowed_dx() result(value)
    value = read_positive_env('FIBRE_PROD_STRUCTURE_DRY_MAX_DX', 1.0e-2_dp)
  end function fibre_prod_structure_dry_step_env_max_allowed_dx

  subroutine fibre_prod_structure_dry_step_runtime_diagnostic(ux, uy, uz, status)
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
                                                                      fibre_prod_structure_dry_step_env_max_allowed_dx(), status)
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
  end subroutine fibre_prod_structure_dry_step_runtime_diagnostic

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

end module fibre_prod_structure_dry_step
