module fibre_prod_state
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  public :: fibre_prod_state_type
  public :: fibre_prod_state_allocate
  public :: fibre_prod_state_init_straight
  public :: fibre_prod_state_reset_forces
  public :: fibre_prod_state_destroy
  public :: fibre_prod_state_is_allocated
  public :: fibre_prod_state_all_finite
  public :: fibre_prod_state_segment_length_residual
  public :: fibre_prod_state_total_force_norm
  public :: fibre_prod_state_allocate_sampled_velocity
  public :: fibre_prod_state_attach_sampled_velocity
  public :: fibre_prod_state_allocate_hydro_force_candidate
  public :: fibre_prod_state_attach_hydro_force_candidate
  public :: fibre_prod_state_allocate_structure_input_force
  public :: fibre_prod_state_attach_structure_input_force
  public :: fibre_prod_state_get_structure_coordinates
  public :: fibre_prod_state_get_structure_velocity_or_zero

  type :: fibre_prod_state_type
    integer :: nfibre = 0
    integer :: nnode = 0
    integer :: ndim = 3
    real(real64), allocatable :: x(:, :, :)
    real(real64), allocatable :: v(:, :, :)
    real(real64), allocatable :: a(:, :, :)
    real(real64), allocatable :: tension(:, :)
    real(real64), allocatable :: f_fs(:, :, :)
    real(real64), allocatable :: f_wall(:, :, :)
    real(real64), allocatable :: f_coll(:, :, :)
    real(real64), allocatable :: f_total(:, :, :)
    real(real64), allocatable :: sampled_u(:, :)
    logical :: has_sampled_velocity = .false.
    real(real64), allocatable :: hydro_force_candidate(:, :)
    logical :: has_hydro_force_candidate = .false.
    real(real64), allocatable :: structure_input_force(:, :)
    logical :: has_structure_input_force = .false.
  end type fibre_prod_state_type

contains

  subroutine fibre_prod_state_allocate(state, nfibre, nnode, stat)
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(in) :: nfibre
    integer, intent(in) :: nnode
    integer, intent(out), optional :: stat
    integer :: ierr

    ierr = 0
    call fibre_prod_state_destroy(state)

    if (nfibre < 1 .or. nnode < 2) then
      ierr = 1
    else
      state%nfibre = nfibre
      state%nnode = nnode
      state%ndim = 3

      allocate(state%x(nfibre, nnode, 3), stat=ierr)
      if (ierr == 0) allocate(state%v(nfibre, nnode, 3), stat=ierr)
      if (ierr == 0) allocate(state%a(nfibre, nnode, 3), stat=ierr)
      if (ierr == 0) allocate(state%tension(nfibre, nnode - 1), stat=ierr)
      if (ierr == 0) allocate(state%f_fs(nfibre, nnode, 3), stat=ierr)
      if (ierr == 0) allocate(state%f_wall(nfibre, nnode, 3), stat=ierr)
      if (ierr == 0) allocate(state%f_coll(nfibre, nnode, 3), stat=ierr)
      if (ierr == 0) allocate(state%f_total(nfibre, nnode, 3), stat=ierr)
    end if

    if (ierr == 0) then
      state%x = 0.0_real64
      state%v = 0.0_real64
      state%a = 0.0_real64
      state%tension = 0.0_real64
      call fibre_prod_state_reset_forces(state)
    else
      call fibre_prod_state_destroy(state)
    end if

    if (present(stat)) stat = ierr
  end subroutine fibre_prod_state_allocate

  subroutine fibre_prod_state_init_straight(state, origin, direction, segment_length, spacing, stat)
    type(fibre_prod_state_type), intent(inout) :: state
    real(real64), intent(in) :: origin(3)
    real(real64), intent(in) :: direction(3)
    real(real64), intent(in) :: segment_length
    real(real64), intent(in), optional :: spacing
    integer, intent(out), optional :: stat
    real(real64) :: norm_dir
    real(real64) :: unit_dir(3)
    real(real64) :: node_spacing
    integer :: i_fibre
    integer :: i_node
    integer :: ierr

    ierr = 0
    if (.not. fibre_prod_state_is_allocated(state)) then
      ierr = 1
    else if (segment_length <= 0.0_real64) then
      ierr = 2
    else
      node_spacing = segment_length
      if (present(spacing)) node_spacing = spacing
      if (node_spacing <= 0.0_real64) then
        ierr = 3
      else
        norm_dir = sqrt(sum(direction * direction))
        if (norm_dir <= 0.0_real64) then
          ierr = 4
        else
          unit_dir = direction / norm_dir
          do i_fibre = 1, state%nfibre
            do i_node = 1, state%nnode
              state%x(i_fibre, i_node, :) = origin + real(i_node - 1, real64) * node_spacing * unit_dir
            end do
          end do
          state%v = 0.0_real64
          state%a = 0.0_real64
          state%tension = 0.0_real64
          call fibre_prod_state_reset_forces(state)
        end if
      end if
    end if

    if (present(stat)) stat = ierr
  end subroutine fibre_prod_state_init_straight

  subroutine fibre_prod_state_reset_forces(state)
    type(fibre_prod_state_type), intent(inout) :: state

    if (allocated(state%f_fs)) state%f_fs = 0.0_real64
    if (allocated(state%f_wall)) state%f_wall = 0.0_real64
    if (allocated(state%f_coll)) state%f_coll = 0.0_real64
    if (allocated(state%f_total)) state%f_total = 0.0_real64
  end subroutine fibre_prod_state_reset_forces

  subroutine fibre_prod_state_destroy(state)
    type(fibre_prod_state_type), intent(inout) :: state

    if (allocated(state%x)) deallocate(state%x)
    if (allocated(state%v)) deallocate(state%v)
    if (allocated(state%a)) deallocate(state%a)
    if (allocated(state%tension)) deallocate(state%tension)
    if (allocated(state%f_fs)) deallocate(state%f_fs)
    if (allocated(state%f_wall)) deallocate(state%f_wall)
    if (allocated(state%f_coll)) deallocate(state%f_coll)
    if (allocated(state%f_total)) deallocate(state%f_total)
    if (allocated(state%sampled_u)) deallocate(state%sampled_u)
    if (allocated(state%hydro_force_candidate)) deallocate(state%hydro_force_candidate)
    if (allocated(state%structure_input_force)) deallocate(state%structure_input_force)
    state%nfibre = 0
    state%nnode = 0
    state%ndim = 3
    state%has_sampled_velocity = .false.
    state%has_hydro_force_candidate = .false.
    state%has_structure_input_force = .false.
  end subroutine fibre_prod_state_destroy


  subroutine fibre_prod_state_allocate_sampled_velocity(state, nnode, status)
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(in) :: nnode
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    ierr = 0
    if (nnode <= 0) then
      status = 1
      return
    end if
    if (state%nnode /= 0 .and. state%nnode /= nnode) then
      status = 2
      return
    end if
    if (allocated(state%sampled_u)) deallocate(state%sampled_u)
    allocate(state%sampled_u(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 3
      state%has_sampled_velocity = .false.
      return
    end if
    state%sampled_u = 0.0_real64
    state%has_sampled_velocity = .false.
  end subroutine fibre_prod_state_allocate_sampled_velocity

  subroutine fibre_prod_state_attach_sampled_velocity(state, sampled_u, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(real64), intent(in) :: sampled_u(:, :)
    integer, intent(out) :: status

    status = 0
    if (size(sampled_u, 2) /= 3) then
      status = 4
      return
    end if
    if (state%nnode /= 0 .and. size(sampled_u, 1) /= state%nnode) then
      status = 5
      return
    end if
    if (.not. allocated(state%sampled_u)) then
      call fibre_prod_state_allocate_sampled_velocity(state, size(sampled_u, 1), status)
      if (status /= 0) return
    else if (size(state%sampled_u, 1) /= size(sampled_u, 1) .or. size(state%sampled_u, 2) /= 3) then
      status = 6
      return
    end if
    state%sampled_u = sampled_u
    state%has_sampled_velocity = .true.
  end subroutine fibre_prod_state_attach_sampled_velocity



  subroutine fibre_prod_state_allocate_hydro_force_candidate(state, nnode, status)
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(in) :: nnode
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    ierr = 0
    if (nnode <= 0) then
      status = 1
      return
    end if
    if (state%nnode /= 0 .and. state%nnode /= nnode) then
      status = 2
      return
    end if
    if (allocated(state%hydro_force_candidate)) deallocate(state%hydro_force_candidate)
    allocate(state%hydro_force_candidate(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 3
      state%has_hydro_force_candidate = .false.
      return
    end if
    state%hydro_force_candidate = 0.0_real64
    state%has_hydro_force_candidate = .false.
  end subroutine fibre_prod_state_allocate_hydro_force_candidate

  subroutine fibre_prod_state_attach_hydro_force_candidate(state, candidate_force, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(real64), intent(in) :: candidate_force(:, :)
    integer, intent(out) :: status

    status = 0
    if (size(candidate_force, 2) /= 3) then
      status = 4
      return
    end if
    if (state%nnode /= 0 .and. size(candidate_force, 1) /= state%nnode) then
      status = 5
      return
    end if
    if (.not. allocated(state%hydro_force_candidate)) then
      call fibre_prod_state_allocate_hydro_force_candidate(state, size(candidate_force, 1), status)
      if (status /= 0) return
    else if (size(state%hydro_force_candidate, 1) /= size(candidate_force, 1) .or. &
             size(state%hydro_force_candidate, 2) /= 3) then
      status = 6
      return
    end if
    state%hydro_force_candidate = candidate_force
    state%has_hydro_force_candidate = .true.
  end subroutine fibre_prod_state_attach_hydro_force_candidate



  subroutine fibre_prod_state_allocate_structure_input_force(state, nnode, status)
    type(fibre_prod_state_type), intent(inout) :: state
    integer, intent(in) :: nnode
    integer, intent(out) :: status
    integer :: ierr

    status = 0
    ierr = 0
    if (nnode <= 0) then
      status = 1
      return
    end if
    if (state%nnode /= 0 .and. state%nnode /= nnode) then
      status = 2
      return
    end if
    if (allocated(state%structure_input_force)) deallocate(state%structure_input_force)
    allocate(state%structure_input_force(nnode, 3), stat=ierr)
    if (ierr /= 0) then
      status = 3
      state%has_structure_input_force = .false.
      return
    end if
    state%structure_input_force = 0.0_real64
    state%has_structure_input_force = .false.
  end subroutine fibre_prod_state_allocate_structure_input_force

  subroutine fibre_prod_state_attach_structure_input_force(state, structure_input_force, status)
    type(fibre_prod_state_type), intent(inout) :: state
    real(real64), intent(in) :: structure_input_force(:, :)
    integer, intent(out) :: status

    status = 0
    if (size(structure_input_force, 2) /= 3) then
      status = 4
      return
    end if
    if (state%nnode /= 0 .and. size(structure_input_force, 1) /= state%nnode) then
      status = 5
      return
    end if
    if (.not. allocated(state%structure_input_force)) then
      call fibre_prod_state_allocate_structure_input_force(state, size(structure_input_force, 1), status)
      if (status /= 0) return
    else if (size(state%structure_input_force, 1) /= size(structure_input_force, 1) .or. &
             size(state%structure_input_force, 2) /= 3) then
      status = 6
      return
    end if
    state%structure_input_force = structure_input_force
    state%has_structure_input_force = .true.
  end subroutine fibre_prod_state_attach_structure_input_force



  subroutine fibre_prod_state_get_structure_coordinates(state, x, status)
    type(fibre_prod_state_type), intent(in) :: state
    real(real64), intent(out) :: x(:, :)
    integer, intent(out) :: status

    status = 0
    if (.not. allocated(state%x) .or. state%nfibre < 1 .or. state%nnode <= 0) then
      status = 1
      return
    end if
    if (size(x, 1) /= state%nnode .or. size(x, 2) /= 3) then
      status = 2
      return
    end if
    x = state%x(1, :, :)
  end subroutine fibre_prod_state_get_structure_coordinates

  subroutine fibre_prod_state_get_structure_velocity_or_zero(state, structure_u, status)
    type(fibre_prod_state_type), intent(in) :: state
    real(real64), intent(out) :: structure_u(:, :)
    integer, intent(out) :: status

    status = 0
    if (state%nnode <= 0) then
      status = 1
      return
    end if
    if (size(structure_u, 1) /= state%nnode .or. size(structure_u, 2) /= 3) then
      status = 2
      return
    end if
    if (allocated(state%v) .and. state%nfibre >= 1) then
      structure_u = state%v(1, :, :)
    else
      structure_u = 0.0_real64
    end if
  end subroutine fibre_prod_state_get_structure_velocity_or_zero

  pure logical function fibre_prod_state_is_allocated(state) result(is_allocated)
    type(fibre_prod_state_type), intent(in) :: state

    is_allocated = allocated(state%x) .and. allocated(state%v) .and. &
                   allocated(state%a) .and. allocated(state%tension) .and. &
                   allocated(state%f_fs) .and. allocated(state%f_wall) .and. &
                   allocated(state%f_coll) .and. allocated(state%f_total)
  end function fibre_prod_state_is_allocated

  pure logical function fibre_prod_state_all_finite(state) result(all_finite)
    type(fibre_prod_state_type), intent(in) :: state

    all_finite = fibre_prod_state_is_allocated(state)
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%x))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%v))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%a))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%tension))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%f_fs))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%f_wall))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%f_coll))
    if (all_finite) all_finite = all_finite .and. all(ieee_is_finite(state%f_total))
  end function fibre_prod_state_all_finite

  pure real(real64) function fibre_prod_state_segment_length_residual(state, expected_length) result(residual)
    type(fibre_prod_state_type), intent(in) :: state
    real(real64), intent(in) :: expected_length
    integer :: i_fibre
    integer :: i_segment
    real(real64) :: dx(3)
    real(real64) :: length

    residual = huge(1.0_real64)
    if (.not. allocated(state%x)) return
    if (state%nfibre < 1 .or. state%nnode < 2) return

    residual = 0.0_real64
    do i_fibre = 1, state%nfibre
      do i_segment = 1, state%nnode - 1
        dx = state%x(i_fibre, i_segment + 1, :) - state%x(i_fibre, i_segment, :)
        length = sqrt(sum(dx * dx))
        residual = max(residual, abs(length - expected_length))
      end do
    end do
  end function fibre_prod_state_segment_length_residual

  pure real(real64) function fibre_prod_state_total_force_norm(state) result(force_norm)
    type(fibre_prod_state_type), intent(in) :: state

    if (allocated(state%f_total)) then
      force_norm = sqrt(sum(state%f_total * state%f_total))
    else
      force_norm = huge(1.0_real64)
    end if
  end function fibre_prod_state_total_force_norm

end module fibre_prod_state
