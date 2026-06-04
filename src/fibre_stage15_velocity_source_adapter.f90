module fibre_stage15_velocity_source_adapter
  use fibre_stage15_structure_state, only : stage15_structure_state_is_allocated, &
                                            stage15_structure_state_set_velocity, &
                                            stage15_structure_state_get_velocity
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  real(mytype), allocatable :: prescribed_velocity_reference(:,:)
  real(mytype), allocatable :: structure_owned_velocity_snapshot(:,:)

  integer :: n_points = 0
  integer :: structure_owned_velocity_status = 0
  integer :: prescribed_velocity_reference_status = 0
  integer :: velocity_source_adapter_status = 0
  integer :: velocity_equivalence_status = 0
  integer :: feedback_force_equivalence_status = 0
  integer :: zero_slip_status = 0
  integer :: finite_value_status = 0
  integer :: structure_advance_count = 0
  integer :: bending_solve_count = 0
  integer :: tension_solve_count = 0
  integer :: position_time_update_count = 0
  integer :: velocity_time_update_count = 0
  integer :: no_fluid_rhs_modification_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  real(mytype) :: max_velocity_source_diff = huge(1.0_mytype)
  real(mytype) :: max_feedback_force_diff = huge(1.0_mytype)
  integer :: diagnostic_write_count = 0
  integer :: final_status = 0

  public :: stage15_velocity_source_use_structure_owned
  public :: stage15_velocity_source_initialize_from_prescribed
  public :: stage15_velocity_source_compare_with_prescribed
  public :: stage15_velocity_source_record_feedback_force_diff
  public :: stage15_velocity_source_record_zero_slip
  public :: stage15_velocity_source_get_vf
  public :: stage15_velocity_source_get_status_values
  public :: stage15_velocity_source_write_diagnostics
  public :: stage15_velocity_source_finalize

contains

  logical function stage15_velocity_source_use_structure_owned()
    stage15_velocity_source_use_structure_owned = structure_owned_velocity_status == 1 .and. &
                                                  velocity_source_adapter_status == 1
  end function stage15_velocity_source_use_structure_owned

  subroutine stage15_velocity_source_initialize_from_prescribed(v_prescribed)
    real(mytype), intent(in) :: v_prescribed(:,:)
    integer :: alloc_status
    integer :: set_status
    integer :: get_status

    call stage15_velocity_source_finalize()

    if (size(v_prescribed, 1) <= 0 .or. size(v_prescribed, 2) /= 3) then
      call update_status_values()
      return
    end if

    n_points = size(v_prescribed, 1)
    allocate(prescribed_velocity_reference(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_velocity_source_finalize()
      return
    end if
    allocate(structure_owned_velocity_snapshot(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_velocity_source_finalize()
      return
    end if

    prescribed_velocity_reference(:, :) = v_prescribed(:, :)
    structure_owned_velocity_snapshot(:, :) = 0.0_mytype
    prescribed_velocity_reference_status = merge(1, 0, all_finite_rank2(prescribed_velocity_reference))

    if (stage15_structure_state_is_allocated() .and. prescribed_velocity_reference_status == 1) then
      call stage15_structure_state_set_velocity(prescribed_velocity_reference, set_status)
      call stage15_structure_state_get_velocity(structure_owned_velocity_snapshot, get_status)
      structure_owned_velocity_status = merge(1, 0, set_status == 1 .and. get_status == 1 .and. &
                                              all_finite_rank2(structure_owned_velocity_snapshot))
    else
      structure_owned_velocity_status = 0
    end if

    call stage15_velocity_source_compare_with_prescribed(0.0_mytype)
  end subroutine stage15_velocity_source_initialize_from_prescribed

  subroutine stage15_velocity_source_compare_with_prescribed(tolerance)
    real(mytype), intent(in) :: tolerance

    max_velocity_source_diff = huge(1.0_mytype)
    velocity_equivalence_status = 0
    if (allocated(prescribed_velocity_reference) .and. allocated(structure_owned_velocity_snapshot)) then
      if (size(prescribed_velocity_reference, 1) == size(structure_owned_velocity_snapshot, 1) .and. &
          size(prescribed_velocity_reference, 2) == size(structure_owned_velocity_snapshot, 2)) then
        max_velocity_source_diff = maxval(abs(structure_owned_velocity_snapshot - prescribed_velocity_reference))
        velocity_equivalence_status = merge(1, 0, max_velocity_source_diff <= tolerance)
      end if
    end if
    call update_status_values()
  end subroutine stage15_velocity_source_compare_with_prescribed

  subroutine stage15_velocity_source_record_feedback_force_diff(force_diff, tolerance)
    real(mytype), intent(in) :: force_diff
    real(mytype), intent(in) :: tolerance

    if (is_finite(force_diff)) then
      max_feedback_force_diff = force_diff
    else
      max_feedback_force_diff = huge(1.0_mytype)
    end if
    feedback_force_equivalence_status = merge(1, 0, max_feedback_force_diff <= tolerance)
    call update_status_values()
  end subroutine stage15_velocity_source_record_feedback_force_diff

  subroutine stage15_velocity_source_record_zero_slip(zero_slip_error, tolerance)
    real(mytype), intent(in) :: zero_slip_error
    real(mytype), intent(in) :: tolerance

    zero_slip_status = merge(1, 0, is_finite(zero_slip_error) .and. zero_slip_error <= tolerance)
    call update_status_values()
  end subroutine stage15_velocity_source_record_zero_slip

  subroutine stage15_velocity_source_get_vf(v_out, get_status)
    real(mytype), intent(out) :: v_out(:,:)
    integer, intent(out), optional :: get_status
    integer :: local_status

    local_status = 0
    if (allocated(structure_owned_velocity_snapshot)) then
      if (size(v_out, 1) == n_points .and. size(v_out, 2) == 3) then
        v_out(:, :) = structure_owned_velocity_snapshot(:, :)
        local_status = merge(1, 0, all_finite_rank2(v_out))
      else
        v_out(:, :) = 0.0_mytype
      end if
    else
      v_out(:, :) = 0.0_mytype
    end if
    if (present(get_status)) get_status = local_status
  end subroutine stage15_velocity_source_get_vf

  subroutine stage15_velocity_source_get_status_values(structure_owned_velocity_out, prescribed_reference_out, &
                                                       adapter_out, npts_out, velocity_diff_out, force_diff_out, &
                                                       velocity_equivalence_out, force_equivalence_out, &
                                                       zero_slip_out, finite_value_out, structure_advance_count_out, &
                                                       bending_solve_count_out, tension_solve_count_out, &
                                                       position_time_update_count_out, velocity_time_update_count_out, &
                                                       no_fluid_rhs_modification_out, &
                                                       no_pressure_projection_modification_out, &
                                                       no_poisson_modification_out, &
                                                       no_rk3_channel_forcing_modification_out, final_out)
    integer, intent(out) :: structure_owned_velocity_out
    integer, intent(out) :: prescribed_reference_out
    integer, intent(out) :: adapter_out
    integer, intent(out) :: npts_out
    real(mytype), intent(out) :: velocity_diff_out
    real(mytype), intent(out) :: force_diff_out
    integer, intent(out) :: velocity_equivalence_out
    integer, intent(out) :: force_equivalence_out
    integer, intent(out) :: zero_slip_out
    integer, intent(out) :: finite_value_out
    integer, intent(out) :: structure_advance_count_out
    integer, intent(out) :: bending_solve_count_out
    integer, intent(out) :: tension_solve_count_out
    integer, intent(out) :: position_time_update_count_out
    integer, intent(out) :: velocity_time_update_count_out
    integer, intent(out) :: no_fluid_rhs_modification_out
    integer, intent(out) :: no_pressure_projection_modification_out
    integer, intent(out) :: no_poisson_modification_out
    integer, intent(out) :: no_rk3_channel_forcing_modification_out
    integer, intent(out) :: final_out

    call update_status_values()
    structure_owned_velocity_out = structure_owned_velocity_status
    prescribed_reference_out = prescribed_velocity_reference_status
    adapter_out = velocity_source_adapter_status
    npts_out = n_points
    velocity_diff_out = max_velocity_source_diff
    force_diff_out = max_feedback_force_diff
    velocity_equivalence_out = velocity_equivalence_status
    force_equivalence_out = feedback_force_equivalence_status
    zero_slip_out = zero_slip_status
    finite_value_out = finite_value_status
    structure_advance_count_out = structure_advance_count
    bending_solve_count_out = bending_solve_count
    tension_solve_count_out = tension_solve_count
    position_time_update_count_out = position_time_update_count
    velocity_time_update_count_out = velocity_time_update_count
    no_fluid_rhs_modification_out = no_fluid_rhs_modification_status
    no_pressure_projection_modification_out = no_pressure_projection_modification_status
    no_poisson_modification_out = no_poisson_modification_status
    no_rk3_channel_forcing_modification_out = no_rk3_channel_forcing_modification_status
    final_out = final_status
  end subroutine stage15_velocity_source_get_status_values

  subroutine stage15_velocity_source_write_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    call update_status_values()
    diagnostic_write_count = diagnostic_write_count + 1
    write(unit_id,'(A,1X,I0)') 'stage15_structure_owned_velocity_status', structure_owned_velocity_status
    write(unit_id,'(A,1X,I0)') 'prescribed_velocity_reference_status', prescribed_velocity_reference_status
    write(unit_id,'(A,1X,I0)') 'velocity_source_adapter_status', velocity_source_adapter_status
    write(unit_id,'(A,1X,I0)') 'npts', n_points
    write(unit_id,'(A,1X,ES24.16)') 'max_velocity_source_diff', max_velocity_source_diff
    write(unit_id,'(A,1X,ES24.16)') 'max_feedback_force_diff', max_feedback_force_diff
    write(unit_id,'(A,1X,I0)') 'velocity_equivalence_status', velocity_equivalence_status
    write(unit_id,'(A,1X,I0)') 'feedback_force_equivalence_status', feedback_force_equivalence_status
    write(unit_id,'(A,1X,I0)') 'zero_slip_status', zero_slip_status
    write(unit_id,'(A,1X,I0)') 'finite_value_status', finite_value_status
    write(unit_id,'(A,1X,I0)') 'structure_advance_count', structure_advance_count
    write(unit_id,'(A,1X,I0)') 'bending_solve_count', bending_solve_count
    write(unit_id,'(A,1X,I0)') 'tension_solve_count', tension_solve_count
    write(unit_id,'(A,1X,I0)') 'position_time_update_count', position_time_update_count
    write(unit_id,'(A,1X,I0)') 'velocity_time_update_count', velocity_time_update_count
    write(unit_id,'(A,1X,I0)') 'no_fluid_rhs_modification_status', no_fluid_rhs_modification_status
    write(unit_id,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit_id,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit_id,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit_id,'(A,1X,I0)') 'diagnostic_write_count', diagnostic_write_count
    write(unit_id,'(A,1X,I0)') 'final_status', final_status
  end subroutine stage15_velocity_source_write_diagnostics

  subroutine stage15_velocity_source_finalize()
    if (allocated(prescribed_velocity_reference)) deallocate(prescribed_velocity_reference)
    if (allocated(structure_owned_velocity_snapshot)) deallocate(structure_owned_velocity_snapshot)
    n_points = 0
    structure_owned_velocity_status = 0
    prescribed_velocity_reference_status = 0
    velocity_source_adapter_status = 0
    velocity_equivalence_status = 0
    feedback_force_equivalence_status = 0
    zero_slip_status = 0
    finite_value_status = 0
    structure_advance_count = 0
    bending_solve_count = 0
    tension_solve_count = 0
    position_time_update_count = 0
    velocity_time_update_count = 0
    no_fluid_rhs_modification_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    max_velocity_source_diff = huge(1.0_mytype)
    max_feedback_force_diff = huge(1.0_mytype)
    diagnostic_write_count = 0
    final_status = 0
  end subroutine stage15_velocity_source_finalize

  subroutine update_status_values()
    finite_value_status = merge(1, 0, finite_values_ready())
    velocity_source_adapter_status = merge(1, 0, structure_owned_velocity_status == 1 .and. &
                                          prescribed_velocity_reference_status == 1 .and. &
                                          velocity_equivalence_status == 1 .and. finite_value_status == 1)
    final_status = merge(1, 0, velocity_source_adapter_status == 1 .and. &
                         feedback_force_equivalence_status == 1 .and. zero_slip_status == 1 .and. &
                         structure_advance_count == 0 .and. bending_solve_count == 0 .and. &
                         tension_solve_count == 0 .and. position_time_update_count == 0 .and. &
                         velocity_time_update_count == 0 .and. no_fluid_rhs_modification_status == 1 .and. &
                         no_pressure_projection_modification_status == 1 .and. &
                         no_poisson_modification_status == 1 .and. &
                         no_rk3_channel_forcing_modification_status == 1)
  end subroutine update_status_values

  logical function finite_values_ready()
    finite_values_ready = .false.
    if (.not. allocated(prescribed_velocity_reference)) return
    if (.not. allocated(structure_owned_velocity_snapshot)) return
    finite_values_ready = all_finite_rank2(prescribed_velocity_reference) .and. &
                          all_finite_rank2(structure_owned_velocity_snapshot) .and. &
                          is_finite(max_velocity_source_diff) .and. is_finite(max_feedback_force_diff)
  end function finite_values_ready

  logical function all_finite_rank2(values)
    real(mytype), intent(in) :: values(:,:)

    all_finite_rank2 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank2

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = value == value .and. abs(value) < huge(1.0_mytype)
  end function is_finite

end module fibre_stage15_velocity_source_adapter
