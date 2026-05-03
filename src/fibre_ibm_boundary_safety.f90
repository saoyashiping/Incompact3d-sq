module fibre_ibm_boundary_safety

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_delta, only : delta_weight_3d

  implicit none
  private

  type, public :: ibm_point_safety_t
    integer :: nl = 0
    integer, allocatable :: status(:)
    logical, allocatable :: inside_domain(:)
    logical, allocatable :: unsafe_near_wall(:)
    logical, allocatable :: periodic_wrap_needed(:)
    real(mytype), allocatable :: min_wall_distance(:)
  end type ibm_point_safety_t

  public :: init_point_safety
  public :: destroy_point_safety
  public :: check_ibm_point_boundary_safety
  public :: summarize_ibm_boundary_safety
  public :: assert_no_unsafe_ibm_points

contains

  subroutine init_point_safety(safety, nl)
    type(ibm_point_safety_t), intent(inout) :: safety
    integer, intent(in) :: nl

    call destroy_point_safety(safety)

    safety%nl = nl
    allocate(safety%status(nl), safety%inside_domain(nl), safety%unsafe_near_wall(nl), &
             safety%periodic_wrap_needed(nl), safety%min_wall_distance(nl))

    safety%status = -1
    safety%inside_domain = .false.
    safety%unsafe_near_wall = .false.
    safety%periodic_wrap_needed = .false.
    safety%min_wall_distance = huge(1._mytype)
  end subroutine init_point_safety

  subroutine destroy_point_safety(safety)
    type(ibm_point_safety_t), intent(inout) :: safety

    if (allocated(safety%status)) deallocate(safety%status)
    if (allocated(safety%inside_domain)) deallocate(safety%inside_domain)
    if (allocated(safety%unsafe_near_wall)) deallocate(safety%unsafe_near_wall)
    if (allocated(safety%periodic_wrap_needed)) &
      deallocate(safety%periodic_wrap_needed)
    if (allocated(safety%min_wall_distance)) &
      deallocate(safety%min_wall_distance)
    safety%nl = 0
  end subroutine destroy_point_safety

  subroutine check_ibm_point_boundary_safety(grid, lag, safety)
    type(ibm_grid_t), intent(in) :: grid
    type(ibm_lagrangian_points_t), intent(in) :: lag
    type(ibm_point_safety_t), intent(inout) :: safety

    integer :: l
    real(mytype) :: x, y, z, d_ymin, d_ymax
    logical :: outside_nonperiodic, unsafe_y, wrap_needed

    if (.not. allocated(safety%status) .or. safety%nl /= lag%nl) then
      call init_point_safety(safety, lag%nl)
    else
      safety%status = -1
      safety%inside_domain = .false.
      safety%unsafe_near_wall = .false.
      safety%periodic_wrap_needed = .false.
      safety%min_wall_distance = huge(1._mytype)
    end if

    do l = 1, lag%nl
      x = lag%x(1, l)
      y = lag%x(2, l)
      z = lag%x(3, l)

      outside_nonperiodic = .false.
      unsafe_y = .false.
      wrap_needed = .false.

      if (.not. grid%periodic_x) then
        if (x < grid%xmin .or. x > grid%xmax) outside_nonperiodic = .true.
      else
        if ((x - grid%xmin) < 2._mytype * grid%dx .or. &
            (grid%xmax - x) < 2._mytype * grid%dx) then
          wrap_needed = .true.
        end if
      end if

      if (.not. grid%periodic_z) then
        if (z < grid%zmin .or. z > grid%zmax) outside_nonperiodic = .true.
      else
        if ((z - grid%zmin) < 2._mytype * grid%dz .or. &
            (grid%zmax - z) < 2._mytype * grid%dz) then
          wrap_needed = .true.
        end if
      end if

      if (.not. grid%periodic_y) then
        if (y < grid%ymin .or. y > grid%ymax) then
          outside_nonperiodic = .true.
          safety%inside_domain(l) = .false.
          safety%min_wall_distance(l) = -1._mytype
        else
          d_ymin = y - grid%ymin
          d_ymax = grid%ymax - y
          safety%inside_domain(l) = .true.
          safety%min_wall_distance(l) = min(d_ymin, d_ymax)
          if (d_ymin < 2._mytype * grid%dy .or. d_ymax < 2._mytype * grid%dy) then
            unsafe_y = .true.
          end if
        end if
      else
        safety%inside_domain(l) = .true.
        safety%min_wall_distance(l) = grid%ymax - grid%ymin
      end if

      if (outside_nonperiodic) then
        safety%status(l) = 3
      else if (unsafe_y) then
        safety%status(l) = 2
        safety%unsafe_near_wall(l) = .true.
      else if (wrap_needed) then
        safety%status(l) = 1
        safety%periodic_wrap_needed(l) = .true.
      else
        safety%status(l) = 0
      end if

      if (safety%status(l) <= 2 .and. .not. safety%inside_domain(l) .and. grid%periodic_y) then
        safety%inside_domain(l) = .true.
      end if
      if (safety%status(l) == 1) safety%periodic_wrap_needed(l) = .true.
    end do
  end subroutine check_ibm_point_boundary_safety

  subroutine summarize_ibm_boundary_safety(safety, safe_count, periodic_wrap_count, unsafe_near_wall_count, outside_count)
    type(ibm_point_safety_t), intent(in) :: safety
    integer, intent(out) :: safe_count, periodic_wrap_count, unsafe_near_wall_count, outside_count

    safe_count = count(safety%status == 0)
    periodic_wrap_count = count(safety%status == 1)
    unsafe_near_wall_count = count(safety%status == 2)
    outside_count = count(safety%status == 3)
  end subroutine summarize_ibm_boundary_safety

  subroutine assert_no_unsafe_ibm_points(safety, unsafe_count)
    type(ibm_point_safety_t), intent(in) :: safety
    integer, intent(out) :: unsafe_count

    unsafe_count = count(safety%status == 2) + count(safety%status == 3)
  end subroutine assert_no_unsafe_ibm_points

end module fibre_ibm_boundary_safety
