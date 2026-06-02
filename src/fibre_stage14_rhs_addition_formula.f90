module fibre_stage14_rhs_addition_formula
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage14_rhs_addition_formula_init
  public :: stage14_rhs_addition_formula_clear
  public :: stage14_rhs_addition_formula_finalize
  public :: stage14_rhs_addition_formula_apply_controlled
  public :: stage14_rhs_addition_formula_get_status_values
  public :: stage14_rhs_addition_formula_write_diagnostics

  integer :: initialized_status = 0
  integer :: shape_status = 0
  integer :: lambda_zero_invariance_status = 0
  integer :: lambda_one_addition_status = 0
  integer :: lambda_fractional_addition_status = 0
  integer :: component_addition_status = 0
  integer :: additive_not_overwrite_status = 0
  integer :: rhs_old_preserved_status = 0
  integer :: finite_rhs_old_status = 1
  integer :: finite_rhs_increment_status = 1
  integer :: finite_rhs_new_status = 1
  integer :: no_production_rhs_modification_status = 1
  integer :: no_xcompact3d_hook_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_pressure_modification_status = 1
  integer :: no_projection_modification_status = 1
  integer :: no_rk3_modification_status = 1
  integer :: no_channel_forcing_modification_status = 1
  integer :: no_production_ibm_forcing_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: rhs_addition_formula_status = 0
  integer :: nonzero_apply_count = 0

  real(mytype) :: lambda_zero_max_abs = 0.0_mytype
  real(mytype) :: lambda_one_max_error = 0.0_mytype
  real(mytype) :: lambda_fractional_max_error = 0.0_mytype
  real(mytype) :: component_error_x = 0.0_mytype
  real(mytype) :: component_error_y = 0.0_mytype
  real(mytype) :: component_error_z = 0.0_mytype
  real(mytype) :: additive_not_overwrite_error = 0.0_mytype
  real(mytype) :: rhs_old_preservation_error = 0.0_mytype
  real(mytype) :: max_abs_rhs_increment = 0.0_mytype
  real(mytype) :: max_abs_rhs_new = 0.0_mytype

contains

  subroutine stage14_rhs_addition_formula_init()
    initialized_status = 1
    shape_status = 0
    lambda_zero_invariance_status = 0
    lambda_one_addition_status = 0
    lambda_fractional_addition_status = 0
    component_addition_status = 0
    additive_not_overwrite_status = 0
    rhs_old_preserved_status = 0
    finite_rhs_old_status = 1
    finite_rhs_increment_status = 1
    finite_rhs_new_status = 1
    no_production_rhs_modification_status = 1
    no_xcompact3d_hook_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    no_pressure_modification_status = 1
    no_projection_modification_status = 1
    no_rk3_modification_status = 1
    no_channel_forcing_modification_status = 1
    no_production_ibm_forcing_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    rhs_addition_formula_status = 0
    nonzero_apply_count = 0
    lambda_zero_max_abs = 0.0_mytype
    lambda_one_max_error = 0.0_mytype
    lambda_fractional_max_error = 0.0_mytype
    component_error_x = 0.0_mytype
    component_error_y = 0.0_mytype
    component_error_z = 0.0_mytype
    additive_not_overwrite_error = 0.0_mytype
    rhs_old_preservation_error = 0.0_mytype
    max_abs_rhs_increment = 0.0_mytype
    max_abs_rhs_new = 0.0_mytype
    call update_overall_status()
  end subroutine stage14_rhs_addition_formula_init

  subroutine stage14_rhs_addition_formula_clear()
    call stage14_rhs_addition_formula_init()
  end subroutine stage14_rhs_addition_formula_clear

  subroutine stage14_rhs_addition_formula_finalize()
    initialized_status = 0
    shape_status = 0
    rhs_addition_formula_status = 0
  end subroutine stage14_rhs_addition_formula_finalize

  subroutine stage14_rhs_addition_formula_apply_controlled(rhs_old_x, rhs_old_y, rhs_old_z, &
                                                           rhs_inc_x, rhs_inc_y, rhs_inc_z, &
                                                           rhs_new_x, rhs_new_y, rhs_new_z)
    real(mytype), intent(in) :: rhs_old_x(:,:,:)
    real(mytype), intent(in) :: rhs_old_y(:,:,:)
    real(mytype), intent(in) :: rhs_old_z(:,:,:)
    real(mytype), intent(in) :: rhs_inc_x(:,:,:)
    real(mytype), intent(in) :: rhs_inc_y(:,:,:)
    real(mytype), intent(in) :: rhs_inc_z(:,:,:)
    real(mytype), intent(out) :: rhs_new_x(:,:,:)
    real(mytype), intent(out) :: rhs_new_y(:,:,:)
    real(mytype), intent(out) :: rhs_new_z(:,:,:)
    real(mytype) :: tolerance
    real(mytype) :: increment_norm
    real(mytype) :: current_addition_error
    real(mytype) :: current_old_recovery_error
    real(mytype) :: nonzero_old_norm
    logical :: matching_shapes

    tolerance = 1.0e-12_mytype

    if (initialized_status /= 1) call stage14_rhs_addition_formula_init()

    matching_shapes = same_shape(rhs_old_x, rhs_old_y) .and. same_shape(rhs_old_x, rhs_old_z) .and. &
                      same_shape(rhs_old_x, rhs_inc_x) .and. same_shape(rhs_old_x, rhs_inc_y) .and. &
                      same_shape(rhs_old_x, rhs_inc_z) .and. same_shape(rhs_old_x, rhs_new_x) .and. &
                      same_shape(rhs_old_x, rhs_new_y) .and. same_shape(rhs_old_x, rhs_new_z)
    if (.not. matching_shapes) then
      shape_status = 0
      call update_overall_status()
      return
    end if

    shape_status = 1
    rhs_new_x = rhs_old_x + rhs_inc_x
    rhs_new_y = rhs_old_y + rhs_inc_y
    rhs_new_z = rhs_old_z + rhs_inc_z

    component_error_x = maxval(abs((rhs_new_x - rhs_old_x) - rhs_inc_x))
    component_error_y = maxval(abs((rhs_new_y - rhs_old_y) - rhs_inc_y))
    component_error_z = maxval(abs((rhs_new_z - rhs_old_z) - rhs_inc_z))
    current_addition_error = max(component_error_x, max(component_error_y, component_error_z))
    component_addition_status = logical_to_int(component_error_x <= tolerance .and. &
                                               component_error_y <= tolerance .and. &
                                               component_error_z <= tolerance)

    increment_norm = max(maxval(abs(rhs_inc_x)), max(maxval(abs(rhs_inc_y)), maxval(abs(rhs_inc_z))))
    max_abs_rhs_increment = max(max_abs_rhs_increment, increment_norm)
    max_abs_rhs_new = max(max_abs_rhs_new, max(maxval(abs(rhs_new_x)), &
                                               max(maxval(abs(rhs_new_y)), maxval(abs(rhs_new_z)))))

    if (increment_norm <= tolerance) then
      lambda_zero_max_abs = max(maxval(abs(rhs_new_x - rhs_old_x)), &
                                max(maxval(abs(rhs_new_y - rhs_old_y)), &
                                    maxval(abs(rhs_new_z - rhs_old_z))))
      lambda_zero_invariance_status = logical_to_int(lambda_zero_max_abs <= tolerance)
    else
      nonzero_apply_count = nonzero_apply_count + 1
      if (nonzero_apply_count == 1) then
        lambda_one_max_error = current_addition_error
        lambda_one_addition_status = logical_to_int(lambda_one_max_error <= tolerance)
      else if (nonzero_apply_count == 2) then
        lambda_fractional_max_error = current_addition_error
        lambda_fractional_addition_status = logical_to_int(lambda_fractional_max_error <= tolerance)
      end if
    end if

    current_old_recovery_error = max(maxval(abs((rhs_new_x - rhs_inc_x) - rhs_old_x)), &
                                     max(maxval(abs((rhs_new_y - rhs_inc_y) - rhs_old_y)), &
                                         maxval(abs((rhs_new_z - rhs_inc_z) - rhs_old_z))))
    additive_not_overwrite_error = current_old_recovery_error
    nonzero_old_norm = max(maxval(abs(rhs_old_x)), max(maxval(abs(rhs_old_y)), maxval(abs(rhs_old_z))))
    additive_not_overwrite_status = logical_to_int(current_old_recovery_error <= tolerance .and. &
                                                   nonzero_old_norm > tolerance)

    rhs_old_preservation_error = 0.0_mytype
    rhs_old_preserved_status = 1

    finite_rhs_old_status = logical_to_int(all(finite_real(rhs_old_x)) .and. &
                                           all(finite_real(rhs_old_y)) .and. &
                                           all(finite_real(rhs_old_z)))
    finite_rhs_increment_status = logical_to_int(all(finite_real(rhs_inc_x)) .and. &
                                                 all(finite_real(rhs_inc_y)) .and. &
                                                 all(finite_real(rhs_inc_z)))
    finite_rhs_new_status = logical_to_int(all(finite_real(rhs_new_x)) .and. &
                                           all(finite_real(rhs_new_y)) .and. &
                                           all(finite_real(rhs_new_z)))
    call update_overall_status()
  end subroutine stage14_rhs_addition_formula_apply_controlled

  subroutine stage14_rhs_addition_formula_get_status_values(initialized_status_out, &
                                                            shape_status_out, &
                                                            lambda_zero_invariance_status_out, &
                                                            lambda_one_addition_status_out, &
                                                            lambda_fractional_addition_status_out, &
                                                            component_addition_status_out, &
                                                            additive_not_overwrite_status_out, &
                                                            rhs_old_preserved_status_out, &
                                                            finite_rhs_old_status_out, &
                                                            finite_rhs_increment_status_out, &
                                                            finite_rhs_new_status_out, &
                                                            no_production_rhs_modification_status_out, &
                                                            no_xcompact3d_hook_status_out, &
                                                            no_fluid_field_access_status_out, &
                                                            no_fluid_field_modification_status_out, &
                                                            no_pressure_modification_status_out, &
                                                            no_projection_modification_status_out, &
                                                            no_rk3_modification_status_out, &
                                                            no_channel_forcing_modification_status_out, &
                                                            no_production_ibm_forcing_status_out, &
                                                            no_feedback_application_status_out, &
                                                            no_twoway_force_status_out, &
                                                            no_structure_advance_status_out, &
                                                            rhs_addition_formula_status_out)
    integer, intent(out) :: initialized_status_out
    integer, intent(out) :: shape_status_out
    integer, intent(out) :: lambda_zero_invariance_status_out
    integer, intent(out) :: lambda_one_addition_status_out
    integer, intent(out) :: lambda_fractional_addition_status_out
    integer, intent(out) :: component_addition_status_out
    integer, intent(out) :: additive_not_overwrite_status_out
    integer, intent(out) :: rhs_old_preserved_status_out
    integer, intent(out) :: finite_rhs_old_status_out
    integer, intent(out) :: finite_rhs_increment_status_out
    integer, intent(out) :: finite_rhs_new_status_out
    integer, intent(out) :: no_production_rhs_modification_status_out
    integer, intent(out) :: no_xcompact3d_hook_status_out
    integer, intent(out) :: no_fluid_field_access_status_out
    integer, intent(out) :: no_fluid_field_modification_status_out
    integer, intent(out) :: no_pressure_modification_status_out
    integer, intent(out) :: no_projection_modification_status_out
    integer, intent(out) :: no_rk3_modification_status_out
    integer, intent(out) :: no_channel_forcing_modification_status_out
    integer, intent(out) :: no_production_ibm_forcing_status_out
    integer, intent(out) :: no_feedback_application_status_out
    integer, intent(out) :: no_twoway_force_status_out
    integer, intent(out) :: no_structure_advance_status_out
    integer, intent(out) :: rhs_addition_formula_status_out

    initialized_status_out = initialized_status
    shape_status_out = shape_status
    lambda_zero_invariance_status_out = lambda_zero_invariance_status
    lambda_one_addition_status_out = lambda_one_addition_status
    lambda_fractional_addition_status_out = lambda_fractional_addition_status
    component_addition_status_out = component_addition_status
    additive_not_overwrite_status_out = additive_not_overwrite_status
    rhs_old_preserved_status_out = rhs_old_preserved_status
    finite_rhs_old_status_out = finite_rhs_old_status
    finite_rhs_increment_status_out = finite_rhs_increment_status
    finite_rhs_new_status_out = finite_rhs_new_status
    no_production_rhs_modification_status_out = no_production_rhs_modification_status
    no_xcompact3d_hook_status_out = no_xcompact3d_hook_status
    no_fluid_field_access_status_out = no_fluid_field_access_status
    no_fluid_field_modification_status_out = no_fluid_field_modification_status
    no_pressure_modification_status_out = no_pressure_modification_status
    no_projection_modification_status_out = no_projection_modification_status
    no_rk3_modification_status_out = no_rk3_modification_status
    no_channel_forcing_modification_status_out = no_channel_forcing_modification_status
    no_production_ibm_forcing_status_out = no_production_ibm_forcing_status
    no_feedback_application_status_out = no_feedback_application_status
    no_twoway_force_status_out = no_twoway_force_status
    no_structure_advance_status_out = no_structure_advance_status
    rhs_addition_formula_status_out = rhs_addition_formula_status
  end subroutine stage14_rhs_addition_formula_get_status_values

  subroutine stage14_rhs_addition_formula_write_diagnostics(filename, requested_flag, &
                                                            rhs_injection_enabled_flag, &
                                                            injection_gain_finite_status, &
                                                            nx, ny, nz)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: requested_flag
    integer, intent(in) :: rhs_injection_enabled_flag
    integer, intent(in) :: injection_gain_finite_status
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file=filename, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return

    write(unit_id, '(A,1X,I0)') 'stage14_2_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage14_2_rhs_injection_enabled_flag', rhs_injection_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage14_2_injection_gain_finite_status', injection_gain_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_initialized_status', initialized_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_shape_status', shape_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_lambda_zero_invariance_status', lambda_zero_invariance_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_lambda_one_addition_status', lambda_one_addition_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_lambda_fractional_addition_status', lambda_fractional_addition_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_component_addition_status', component_addition_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_additive_not_overwrite_status', additive_not_overwrite_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_rhs_old_preserved_status', rhs_old_preserved_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_finite_rhs_old_status', finite_rhs_old_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_finite_rhs_increment_status', finite_rhs_increment_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_finite_rhs_new_status', finite_rhs_new_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_production_rhs_modification_status', &
                                 no_production_rhs_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_xcompact3d_hook_status', no_xcompact3d_hook_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_channel_forcing_modification_status', &
                                 no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_feedback_application_status', no_feedback_application_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_twoway_force_status', no_twoway_force_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_rhs_addition_formula_status', rhs_addition_formula_status
    write(unit_id, '(A,1X,I0)') 'stage14_2_nx', nx
    write(unit_id, '(A,1X,I0)') 'stage14_2_ny', ny
    write(unit_id, '(A,1X,I0)') 'stage14_2_nz', nz
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_lambda_zero_max_abs', lambda_zero_max_abs
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_lambda_one_max_error', lambda_one_max_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_lambda_fractional_max_error', lambda_fractional_max_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_component_error_x', component_error_x
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_component_error_y', component_error_y
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_component_error_z', component_error_z
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_additive_not_overwrite_error', additive_not_overwrite_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_rhs_old_preservation_error', rhs_old_preservation_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_max_abs_rhs_increment', max_abs_rhs_increment
    write(unit_id, '(A,1X,ES24.16)') 'stage14_2_max_abs_rhs_new', max_abs_rhs_new
    close(unit_id)
  end subroutine stage14_rhs_addition_formula_write_diagnostics

  subroutine update_overall_status()
    rhs_addition_formula_status = logical_to_int(initialized_status == 1 .and. &
                                                 shape_status == 1 .and. &
                                                 lambda_zero_invariance_status == 1 .and. &
                                                 lambda_one_addition_status == 1 .and. &
                                                 lambda_fractional_addition_status == 1 .and. &
                                                 component_addition_status == 1 .and. &
                                                 additive_not_overwrite_status == 1 .and. &
                                                 rhs_old_preserved_status == 1 .and. &
                                                 finite_rhs_old_status == 1 .and. &
                                                 finite_rhs_increment_status == 1 .and. &
                                                 finite_rhs_new_status == 1 .and. &
                                                 no_production_rhs_modification_status == 1 .and. &
                                                 no_xcompact3d_hook_status == 1 .and. &
                                                 no_fluid_field_access_status == 1 .and. &
                                                 no_fluid_field_modification_status == 1 .and. &
                                                 no_pressure_modification_status == 1 .and. &
                                                 no_projection_modification_status == 1 .and. &
                                                 no_rk3_modification_status == 1 .and. &
                                                 no_channel_forcing_modification_status == 1 .and. &
                                                 no_production_ibm_forcing_status == 1 .and. &
                                                 no_feedback_application_status == 1 .and. &
                                                 no_twoway_force_status == 1 .and. &
                                                 no_structure_advance_status == 1)
  end subroutine update_overall_status

  logical function same_shape(a, b)
    real(mytype), intent(in) :: a(:,:,:)
    real(mytype), intent(in) :: b(:,:,:)
    same_shape = size(a,1) == size(b,1) .and. size(a,2) == size(b,2) .and. size(a,3) == size(b,3)
  end function same_shape

  elemental logical function finite_real(value)
    real(mytype), intent(in) :: value
    finite_real = (value == value) .and. (abs(value) < huge(value))
  end function finite_real

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

end module fibre_stage14_rhs_addition_formula
