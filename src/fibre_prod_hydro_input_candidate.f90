module fibre_prod_hydro_input_candidate
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity, &
                              fibre_prod_state_attach_hydro_force_candidate
  use fibre_prod_state_velocity_attachment, only : fibre_prod_state_velocity_attachment_type, &
                                                   fibre_prod_state_velocity_attachment_init, &
                                                   fibre_prod_state_velocity_attachment_set_points, &
                                                   fibre_prod_state_velocity_attachment_sample, &
                                                   fibre_prod_state_velocity_attachment_attach_to_state, &
                                                   fibre_prod_state_velocity_attachment_finalize
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_hydro_input_candidate_type
  public :: fibre_prod_hydro_input_candidate_init
  public :: fibre_prod_hydro_input_candidate_compute
  public :: fibre_prod_hydro_input_candidate_attach_to_state
  public :: fibre_prod_hydro_input_candidate_finalize
  public :: fibre_prod_hydro_input_candidate_env_enabled
  public :: fibre_prod_hydro_input_candidate_env_beta
  public :: fibre_prod_hydro_input_candidate_runtime_diagnostic

  type :: fibre_prod_hydro_input_candidate_type
    logical :: initialized = .false.
    logical :: has_candidate = .false.
    integer :: nnode = 0
    real(dp) :: beta_hydro = 1.0_dp
    real(dp), allocatable :: relative_u(:, :)
    real(dp), allocatable :: candidate_force(:, :)
    real(dp) :: max_abs_relative_u = 0.0_dp
    real(dp) :: max_abs_candidate_force = 0.0_dp
    real(dp) :: sum_abs_candidate_force = 0.0_dp
  end type fibre_prod_hydro_input_candidate_type

contains

  subroutine fibre_prod_hydro_input_candidate_init(candidate, nnode, beta_hydro, status)
    type(fibre_prod_hydro_input_candidate_type), intent(inout) :: candidate
    integer, intent(in) :: nnode
    real(dp), intent(in) :: beta_hydro
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    call fibre_prod_hydro_input_candidate_finalize(candidate)
    if (nnode <= 0) then
      status = 1
      return
    end if
    if (.not. ieee_is_finite(beta_hydro) .or. beta_hydro < 0.0_dp) then
      status = 2
      return
    end if
    allocate(candidate%relative_u(nnode, 3), candidate%candidate_force(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 3
      call fibre_prod_hydro_input_candidate_finalize(candidate)
      return
    end if
    candidate%nnode = nnode
    candidate%beta_hydro = beta_hydro
    candidate%relative_u = 0.0_dp
    candidate%candidate_force = 0.0_dp
    candidate%max_abs_relative_u = 0.0_dp
    candidate%max_abs_candidate_force = 0.0_dp
    candidate%sum_abs_candidate_force = 0.0_dp
    candidate%has_candidate = .false.
    candidate%initialized = .true.
  end subroutine fibre_prod_hydro_input_candidate_init

  subroutine fibre_prod_hydro_input_candidate_compute(candidate, sampled_u, structure_u, status)
    type(fibre_prod_hydro_input_candidate_type), intent(inout) :: candidate
    real(dp), intent(in) :: sampled_u(:, :)
    real(dp), intent(in) :: structure_u(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. candidate%initialized .or. .not. allocated(candidate%relative_u) .or. &
        .not. allocated(candidate%candidate_force)) then
      status = 4
      return
    end if
    if (size(sampled_u, 1) /= candidate%nnode .or. size(sampled_u, 2) /= 3) then
      status = 5
      return
    end if
    if (size(structure_u, 1) /= candidate%nnode .or. size(structure_u, 2) /= 3) then
      status = 6
      return
    end if
    if (.not. all(ieee_is_finite(sampled_u)) .or. .not. all(ieee_is_finite(structure_u))) then
      status = 7
      return
    end if

    candidate%relative_u = sampled_u - structure_u
    candidate%candidate_force = candidate%beta_hydro * candidate%relative_u
    candidate%max_abs_relative_u = maxval(abs(candidate%relative_u))
    candidate%max_abs_candidate_force = maxval(abs(candidate%candidate_force))
    candidate%sum_abs_candidate_force = sum(abs(candidate%candidate_force))
    candidate%has_candidate = .true.
  end subroutine fibre_prod_hydro_input_candidate_compute

  subroutine fibre_prod_hydro_input_candidate_attach_to_state(candidate, state, status)
    type(fibre_prod_hydro_input_candidate_type), intent(in) :: candidate
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(out) :: status

    status = 0
    if (.not. candidate%initialized .or. .not. candidate%has_candidate .or. &
        .not. allocated(candidate%candidate_force)) then
      status = 8
      return
    end if
    call fibre_prod_state_attach_hydro_force_candidate(state, candidate%candidate_force, status)
  end subroutine fibre_prod_hydro_input_candidate_attach_to_state

  subroutine fibre_prod_hydro_input_candidate_finalize(candidate)
    type(fibre_prod_hydro_input_candidate_type), intent(inout) :: candidate

    if (allocated(candidate%relative_u)) deallocate(candidate%relative_u)
    if (allocated(candidate%candidate_force)) deallocate(candidate%candidate_force)
    candidate%initialized = .false.
    candidate%has_candidate = .false.
    candidate%nnode = 0
    candidate%beta_hydro = 1.0_dp
    candidate%max_abs_relative_u = 0.0_dp
    candidate%max_abs_candidate_force = 0.0_dp
    candidate%sum_abs_candidate_force = 0.0_dp
  end subroutine fibre_prod_hydro_input_candidate_finalize

  logical function fibre_prod_hydro_input_candidate_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_HYDRO_INPUT_CANDIDATE_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_hydro_input_candidate_env_enabled

  real(dp) function fibre_prod_hydro_input_candidate_env_beta() result(beta)
    character(len=64) :: raw
    integer :: length
    integer :: env_status
    integer :: read_status

    beta = 1.0_dp
    call get_environment_variable('FIBRE_PROD_HYDRO_BETA', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    read(raw(1:min(length, len(raw))), *, iostat=read_status) beta
    if (read_status /= 0 .or. .not. ieee_is_finite(beta) .or. beta < 0.0_dp) beta = 1.0_dp
  end function fibre_prod_hydro_input_candidate_env_beta

  subroutine fibre_prod_hydro_input_candidate_runtime_diagnostic(ux, uy, uz, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_state_type) :: state
    type(fibre_prod_state_velocity_attachment_type) :: attach
    type(fibre_prod_hydro_input_candidate_type) :: candidate
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: structure_u(3, 3)
    real(dp) :: points(3, 3)

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
    structure_u = 0.0_dp
    if (status == 0) call fibre_prod_state_allocate(state, 1, 3, status)
    if (status == 0) state%x(1, :, :) = points
    if (status == 0) call fibre_prod_state_velocity_attachment_init(attach, 3, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_set_points(attach, points, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_init(candidate, 3, &
                                                               fibre_prod_hydro_input_candidate_env_beta(), status)
    if (status == 0) call fibre_prod_hydro_input_candidate_compute(candidate, state%sampled_u, structure_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_attach_to_state(candidate, state, status)
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
  end subroutine fibre_prod_hydro_input_candidate_runtime_diagnostic

end module fibre_prod_hydro_input_candidate
