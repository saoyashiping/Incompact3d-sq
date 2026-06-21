module fibre_prod_structure_input_handoff
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_structure_input_force
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
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_structure_input_handoff_type
  public :: fibre_prod_structure_input_handoff_init
  public :: fibre_prod_structure_input_handoff_from_candidate
  public :: fibre_prod_structure_input_handoff_attach_to_state
  public :: fibre_prod_structure_input_handoff_finalize
  public :: fibre_prod_structure_input_handoff_env_enabled
  public :: fibre_prod_structure_input_handoff_runtime_diagnostic

  type :: fibre_prod_structure_input_handoff_type
    logical :: initialized = .false.
    logical :: has_input = .false.
    integer :: nnode = 0
    real(dp), allocatable :: structure_input_force(:, :)
    real(dp) :: max_abs_input_force = 0.0_dp
    real(dp) :: sum_abs_input_force = 0.0_dp
  end type fibre_prod_structure_input_handoff_type

contains

  subroutine fibre_prod_structure_input_handoff_init(handoff, nnode, status)
    type(fibre_prod_structure_input_handoff_type), intent(inout) :: handoff
    integer, intent(in) :: nnode
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    call fibre_prod_structure_input_handoff_finalize(handoff)
    if (nnode <= 0) then
      status = 1
      return
    end if
    allocate(handoff%structure_input_force(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 2
      call fibre_prod_structure_input_handoff_finalize(handoff)
      return
    end if
    handoff%nnode = nnode
    handoff%structure_input_force = 0.0_dp
    handoff%max_abs_input_force = 0.0_dp
    handoff%sum_abs_input_force = 0.0_dp
    handoff%has_input = .false.
    handoff%initialized = .true.
  end subroutine fibre_prod_structure_input_handoff_init

  subroutine fibre_prod_structure_input_handoff_from_candidate(handoff, candidate_force, status)
    type(fibre_prod_structure_input_handoff_type), intent(inout) :: handoff
    real(dp), intent(in) :: candidate_force(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. handoff%initialized .or. .not. allocated(handoff%structure_input_force)) then
      status = 3
      return
    end if
    if (size(candidate_force, 1) /= handoff%nnode .or. size(candidate_force, 2) /= 3) then
      status = 4
      return
    end if
    if (.not. all(ieee_is_finite(candidate_force))) then
      status = 5
      return
    end if
    handoff%structure_input_force = candidate_force
    handoff%max_abs_input_force = maxval(abs(handoff%structure_input_force))
    handoff%sum_abs_input_force = sum(abs(handoff%structure_input_force))
    handoff%has_input = .true.
  end subroutine fibre_prod_structure_input_handoff_from_candidate

  subroutine fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    type(fibre_prod_structure_input_handoff_type), intent(in) :: handoff
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(out) :: status

    status = 0
    if (.not. handoff%initialized .or. .not. handoff%has_input .or. &
        .not. allocated(handoff%structure_input_force)) then
      status = 6
      return
    end if
    call fibre_prod_state_attach_structure_input_force(state, handoff%structure_input_force, status)
  end subroutine fibre_prod_structure_input_handoff_attach_to_state

  subroutine fibre_prod_structure_input_handoff_finalize(handoff)
    type(fibre_prod_structure_input_handoff_type), intent(inout) :: handoff

    if (allocated(handoff%structure_input_force)) deallocate(handoff%structure_input_force)
    handoff%initialized = .false.
    handoff%has_input = .false.
    handoff%nnode = 0
    handoff%max_abs_input_force = 0.0_dp
    handoff%sum_abs_input_force = 0.0_dp
  end subroutine fibre_prod_structure_input_handoff_finalize

  logical function fibre_prod_structure_input_handoff_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_STRUCTURE_INPUT_HANDOFF_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_structure_input_handoff_env_enabled

  subroutine fibre_prod_structure_input_handoff_runtime_diagnostic(ux, uy, uz, status)
    real(dp), intent(in) :: ux(:, :, :)
    real(dp), intent(in) :: uy(:, :, :)
    real(dp), intent(in) :: uz(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_state_type) :: state
    type(fibre_prod_state_velocity_attachment_type) :: attach
    type(fibre_prod_hydro_input_candidate_type) :: candidate
    type(fibre_prod_structure_input_handoff_type) :: handoff
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: points(3, 3)
    real(dp) :: structure_u(3, 3)

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
    if (status == 0) call fibre_prod_structure_input_handoff_init(handoff, 3, status)
    if (status == 0) call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status == 0) call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
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
  end subroutine fibre_prod_structure_input_handoff_runtime_diagnostic

end module fibre_prod_structure_input_handoff
