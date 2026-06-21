module fibre_prod_reaction_force_candidate
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_reaction_force_candidate_type
  public :: fibre_prod_reaction_force_candidate_init
  public :: fibre_prod_reaction_force_candidate_from_structure_input
  public :: fibre_prod_reaction_force_candidate_finalize
  public :: fibre_prod_reaction_force_candidate_env_enabled

  type :: fibre_prod_reaction_force_candidate_type
    logical :: initialized = .false.
    logical :: has_reaction_force = .false.
    integer :: nnode = 0
    real(dp), allocatable :: reaction_force(:, :)
    real(dp) :: max_abs_reaction_force = 0.0_dp
    real(dp) :: sum_abs_reaction_force = 0.0_dp
    real(dp) :: net_reaction_force(3) = 0.0_dp
  end type fibre_prod_reaction_force_candidate_type

contains

  subroutine fibre_prod_reaction_force_candidate_init(candidate, nnode, status)
    type(fibre_prod_reaction_force_candidate_type), intent(inout) :: candidate
    integer, intent(in) :: nnode
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    call fibre_prod_reaction_force_candidate_finalize(candidate)
    if (nnode <= 0) then
      status = 1
      return
    end if
    allocate(candidate%reaction_force(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 2
      call fibre_prod_reaction_force_candidate_finalize(candidate)
      return
    end if
    candidate%nnode = nnode
    candidate%reaction_force = 0.0_dp
    candidate%max_abs_reaction_force = 0.0_dp
    candidate%sum_abs_reaction_force = 0.0_dp
    candidate%net_reaction_force = 0.0_dp
    candidate%has_reaction_force = .false.
    candidate%initialized = .true.
  end subroutine fibre_prod_reaction_force_candidate_init

  subroutine fibre_prod_reaction_force_candidate_from_structure_input(candidate, structure_input_force, status)
    type(fibre_prod_reaction_force_candidate_type), intent(inout) :: candidate
    real(dp), intent(in) :: structure_input_force(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. candidate%initialized .or. .not. allocated(candidate%reaction_force)) then
      status = 3
      return
    end if
    if (size(structure_input_force, 1) /= candidate%nnode .or. size(structure_input_force, 2) /= 3) then
      status = 4
      return
    end if
    if (.not. all(ieee_is_finite(structure_input_force))) then
      status = 5
      return
    end if
    candidate%reaction_force = -structure_input_force
    candidate%max_abs_reaction_force = maxval(abs(candidate%reaction_force))
    candidate%sum_abs_reaction_force = sum(abs(candidate%reaction_force))
    candidate%net_reaction_force = sum(candidate%reaction_force, dim=1)
    candidate%has_reaction_force = .true.
  end subroutine fibre_prod_reaction_force_candidate_from_structure_input

  subroutine fibre_prod_reaction_force_candidate_finalize(candidate)
    type(fibre_prod_reaction_force_candidate_type), intent(inout) :: candidate

    if (allocated(candidate%reaction_force)) deallocate(candidate%reaction_force)
    candidate%initialized = .false.
    candidate%has_reaction_force = .false.
    candidate%nnode = 0
    candidate%max_abs_reaction_force = 0.0_dp
    candidate%sum_abs_reaction_force = 0.0_dp
    candidate%net_reaction_force = 0.0_dp
  end subroutine fibre_prod_reaction_force_candidate_finalize

  logical function fibre_prod_reaction_force_candidate_env_enabled() result(enabled)
    character(len=64) :: raw
    integer :: length
    integer :: env_status

    enabled = .false.
    call get_environment_variable('FIBRE_PROD_REACTION_FORCE_CANDIDATE_ENABLE', raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function fibre_prod_reaction_force_candidate_env_enabled

end module fibre_prod_reaction_force_candidate
