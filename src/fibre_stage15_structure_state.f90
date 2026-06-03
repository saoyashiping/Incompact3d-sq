module fibre_stage15_structure_state
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: n_points = 0
  real(mytype), allocatable :: x_f(:,:)
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: a_f(:,:)
  real(mytype), allocatable :: t_f(:)

  integer :: allocation_status = 0
  integer :: initialization_status = 0
  integer :: clear_status = 0
  integer :: validation_status = 0
  integer :: npts_status = 0
  integer :: x_finite_status = 0
  integer :: v_finite_status = 0
  integer :: optional_a_or_rhs_finite_status = 0
  integer :: optional_tension_finite_status = 0
  integer :: structure_advance_count = 0
  integer :: bending_solve_count = 0
  integer :: tension_solve_count = 0
  integer :: position_time_update_count = 0
  integer :: velocity_time_update_count = 0
  integer :: no_fluid_rhs_modification_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: diagnostic_write_count = 0
  integer :: diagnostic_counter_status = 0
  integer :: final_status = 0

  public :: stage15_structure_state_allocate
  public :: stage15_structure_state_initialize
  public :: stage15_structure_state_clear
  public :: stage15_structure_state_is_allocated
  public :: stage15_structure_state_validate
  public :: stage15_structure_state_write_diagnostics
  public :: stage15_structure_state_finalize
  public :: stage15_structure_state_get_count
  public :: stage15_structure_state_get_status_values

contains

  subroutine stage15_structure_state_allocate(requested_points)
    integer, intent(in) :: requested_points
    integer :: alloc_status

    call stage15_structure_state_finalize()

    if (requested_points <= 0) then
      call update_status_values()
      return
    end if

    n_points = requested_points
    allocate(x_f(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_structure_state_finalize()
      return
    end if

    allocate(v_f(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_structure_state_finalize()
      return
    end if

    allocate(a_f(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_structure_state_finalize()
      return
    end if

    allocate(t_f(n_points), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_structure_state_finalize()
      return
    end if

    x_f(:, :) = 0.0_mytype
    v_f(:, :) = 0.0_mytype
    a_f(:, :) = 0.0_mytype
    t_f(:) = 0.0_mytype
    initialization_status = 0
    clear_status = 0
    validation_status = 0
    call update_status_values()
  end subroutine stage15_structure_state_allocate

  subroutine stage15_structure_state_initialize()
    integer :: i

    if (.not. stage15_structure_state_is_allocated()) then
      initialization_status = 0
      call update_status_values()
      return
    end if

    do i = 1, n_points
      if (n_points > 1) then
        x_f(i, 1) = real(i - 1, mytype) / real(n_points - 1, mytype)
      else
        x_f(i, 1) = 0.0_mytype
      end if
      x_f(i, 2) = 0.0_mytype
      x_f(i, 3) = 0.0_mytype
    end do
    v_f(:, :) = 0.0_mytype
    a_f(:, :) = 0.0_mytype
    t_f(:) = 0.0_mytype
    initialization_status = 1
    call stage15_structure_state_validate()
  end subroutine stage15_structure_state_initialize

  subroutine stage15_structure_state_clear()
    if (.not. stage15_structure_state_is_allocated()) then
      clear_status = 0
      call update_status_values()
      return
    end if

    x_f(:, :) = 0.0_mytype
    v_f(:, :) = 0.0_mytype
    a_f(:, :) = 0.0_mytype
    t_f(:) = 0.0_mytype
    clear_status = 1
    call stage15_structure_state_validate()
  end subroutine stage15_structure_state_clear

  logical function stage15_structure_state_is_allocated()
    stage15_structure_state_is_allocated = allocated(x_f) .and. allocated(v_f) .and. &
                                           allocated(a_f) .and. allocated(t_f)
  end function stage15_structure_state_is_allocated

  subroutine stage15_structure_state_validate()
    call update_status_values()
  end subroutine stage15_structure_state_validate

  integer function stage15_structure_state_get_count()
    stage15_structure_state_get_count = n_points
  end function stage15_structure_state_get_count

  subroutine stage15_structure_state_get_status_values(allocation_out, initialization_out, clear_out, validation_out, &
                                                       npts_out, x_finite_out, v_finite_out, a_finite_out, &
                                                       tension_finite_out, structure_advance_count_out, &
                                                       bending_solve_count_out, tension_solve_count_out, &
                                                       position_time_update_count_out, velocity_time_update_count_out, &
                                                       no_fluid_rhs_modification_out, &
                                                       no_pressure_projection_modification_out, &
                                                       no_poisson_modification_out, &
                                                       no_rk3_channel_forcing_modification_out, &
                                                       diagnostic_write_count_out, diagnostic_counter_out, final_out)
    integer, intent(out) :: allocation_out
    integer, intent(out) :: initialization_out
    integer, intent(out) :: clear_out
    integer, intent(out) :: validation_out
    integer, intent(out) :: npts_out
    integer, intent(out) :: x_finite_out
    integer, intent(out) :: v_finite_out
    integer, intent(out) :: a_finite_out
    integer, intent(out) :: tension_finite_out
    integer, intent(out) :: structure_advance_count_out
    integer, intent(out) :: bending_solve_count_out
    integer, intent(out) :: tension_solve_count_out
    integer, intent(out) :: position_time_update_count_out
    integer, intent(out) :: velocity_time_update_count_out
    integer, intent(out) :: no_fluid_rhs_modification_out
    integer, intent(out) :: no_pressure_projection_modification_out
    integer, intent(out) :: no_poisson_modification_out
    integer, intent(out) :: no_rk3_channel_forcing_modification_out
    integer, intent(out) :: diagnostic_write_count_out
    integer, intent(out) :: diagnostic_counter_out
    integer, intent(out) :: final_out

    call update_status_values()
    allocation_out = allocation_status
    initialization_out = initialization_status
    clear_out = clear_status
    validation_out = validation_status
    npts_out = n_points
    x_finite_out = x_finite_status
    v_finite_out = v_finite_status
    a_finite_out = optional_a_or_rhs_finite_status
    tension_finite_out = optional_tension_finite_status
    structure_advance_count_out = structure_advance_count
    bending_solve_count_out = bending_solve_count
    tension_solve_count_out = tension_solve_count
    position_time_update_count_out = position_time_update_count
    velocity_time_update_count_out = velocity_time_update_count
    no_fluid_rhs_modification_out = no_fluid_rhs_modification_status
    no_pressure_projection_modification_out = no_pressure_projection_modification_status
    no_poisson_modification_out = no_poisson_modification_status
    no_rk3_channel_forcing_modification_out = no_rk3_channel_forcing_modification_status
    diagnostic_write_count_out = diagnostic_write_count
    diagnostic_counter_out = diagnostic_counter_status
    final_out = final_status
  end subroutine stage15_structure_state_get_status_values

  subroutine stage15_structure_state_write_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    call update_status_values()
    write(unit_id,'(A,1X,I0)') 'allocation_status', allocation_status
    write(unit_id,'(A,1X,I0)') 'initialization_status', initialization_status
    write(unit_id,'(A,1X,I0)') 'clear_status', clear_status
    write(unit_id,'(A,1X,I0)') 'validation_status', validation_status
    write(unit_id,'(A,1X,I0)') 'npts', n_points
    write(unit_id,'(A,1X,I0)') 'x_finite_status', x_finite_status
    write(unit_id,'(A,1X,I0)') 'v_finite_status', v_finite_status
    write(unit_id,'(A,1X,I0)') 'optional_a_or_rhs_finite_status', optional_a_or_rhs_finite_status
    write(unit_id,'(A,1X,I0)') 'optional_tension_finite_status', optional_tension_finite_status
    write(unit_id,'(A,1X,I0)') 'structure_advance_count', structure_advance_count
    write(unit_id,'(A,1X,I0)') 'bending_solve_count', bending_solve_count
    write(unit_id,'(A,1X,I0)') 'tension_solve_count', tension_solve_count
    write(unit_id,'(A,1X,I0)') 'position_time_update_count', position_time_update_count
    write(unit_id,'(A,1X,I0)') 'velocity_time_update_count', velocity_time_update_count
    write(unit_id,'(A,1X,I0)') 'no_fluid_rhs_modification_status', no_fluid_rhs_modification_status
    write(unit_id,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit_id,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit_id,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit_id,'(A,1X,I0)') 'diagnostic_write_count', diagnostic_write_count + 1
    write(unit_id,'(A,1X,I0)') 'diagnostic_counter_status', 1
    diagnostic_write_count = diagnostic_write_count + 1
    diagnostic_counter_status = 1
    call update_status_values()
    write(unit_id,'(A,1X,I0)') 'final_status', final_status
  end subroutine stage15_structure_state_write_diagnostics

  subroutine stage15_structure_state_finalize()
    if (allocated(x_f)) deallocate(x_f)
    if (allocated(v_f)) deallocate(v_f)
    if (allocated(a_f)) deallocate(a_f)
    if (allocated(t_f)) deallocate(t_f)
    n_points = 0
    allocation_status = 0
    initialization_status = 0
    clear_status = 0
    validation_status = 0
    npts_status = 0
    x_finite_status = 0
    v_finite_status = 0
    optional_a_or_rhs_finite_status = 0
    optional_tension_finite_status = 0
    structure_advance_count = 0
    bending_solve_count = 0
    tension_solve_count = 0
    position_time_update_count = 0
    velocity_time_update_count = 0
    no_fluid_rhs_modification_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    diagnostic_write_count = 0
    diagnostic_counter_status = 0
    final_status = 0
  end subroutine stage15_structure_state_finalize

  subroutine update_status_values()
    allocation_status = merge(1, 0, stage15_structure_state_is_allocated())
    npts_status = merge(1, 0, n_points > 0)

    if (allocation_status == 1) then
      x_finite_status = merge(1, 0, all_finite_rank2(x_f))
      v_finite_status = merge(1, 0, all_finite_rank2(v_f))
      optional_a_or_rhs_finite_status = merge(1, 0, all_finite_rank2(a_f))
      optional_tension_finite_status = merge(1, 0, all_finite_rank1(t_f))
    else
      x_finite_status = 0
      v_finite_status = 0
      optional_a_or_rhs_finite_status = 0
      optional_tension_finite_status = 0
    end if

    validation_status = merge(1, 0, allocation_status == 1 .and. npts_status == 1 .and. &
                              x_finite_status == 1 .and. v_finite_status == 1 .and. &
                              optional_a_or_rhs_finite_status == 1 .and. optional_tension_finite_status == 1)

    final_status = merge(1, 0, allocation_status == 1 .and. initialization_status == 1 .and. &
                         clear_status == 1 .and. validation_status == 1 .and. &
                         structure_advance_count == 0 .and. bending_solve_count == 0 .and. &
                         tension_solve_count == 0 .and. position_time_update_count == 0 .and. &
                         velocity_time_update_count == 0 .and. &
                         no_fluid_rhs_modification_status == 1 .and. &
                         no_pressure_projection_modification_status == 1 .and. &
                         no_poisson_modification_status == 1 .and. &
                         no_rk3_channel_forcing_modification_status == 1)
  end subroutine update_status_values

  logical function all_finite_rank2(values)
    real(mytype), intent(in) :: values(:,:)

    all_finite_rank2 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank2

  logical function all_finite_rank1(values)
    real(mytype), intent(in) :: values(:)

    all_finite_rank1 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank1

end module fibre_stage15_structure_state
