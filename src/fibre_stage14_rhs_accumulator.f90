module fibre_stage14_rhs_accumulator
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage14_rhs_accumulator_init
  public :: stage14_rhs_accumulator_clear
  public :: stage14_rhs_accumulator_finalize
  public :: stage14_rhs_accumulator_is_allocated
  public :: stage14_rhs_accumulator_get_shape
  public :: stage14_rhs_accumulator_compute_from_force_density
  public :: stage14_rhs_accumulator_get_status_values
  public :: stage14_rhs_accumulator_write_diagnostics

  integer :: nx_rhs_buf = 0
  integer :: ny_rhs_buf = 0
  integer :: nz_rhs_buf = 0
  real(mytype), allocatable :: rhs_fx_cand(:,:,:)
  real(mytype), allocatable :: rhs_fy_cand(:,:,:)
  real(mytype), allocatable :: rhs_fz_cand(:,:,:)
  real(mytype), allocatable :: rhs_increment_norm(:,:,:)
  integer, allocatable :: rhs_increment_valid_flag(:,:,:)

  integer :: allocated_status = 0
  integer :: shape_status = 0
  integer :: zero_initialization_status = 0
  integer :: clear_status = 0
  integer :: lambda_zero_status = 0
  integer :: lambda_one_scaling_status = 0
  integer :: lambda_fractional_scaling_status = 0
  integer :: component_scaling_status = 0
  integer :: finite_accumulator_status = 1
  integer :: rhs_increment_norm_finite_status = 1
  integer :: rhs_increment_valid_flag_status = 1
  integer :: no_rhs_addition_status = 1
  integer :: no_production_hook_status = 1
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
  integer :: rhs_accumulator_status = 0

  real(mytype) :: lambda_zero_max_abs = 0.0_mytype
  real(mytype) :: lambda_one_max_error = 0.0_mytype
  real(mytype) :: lambda_fractional_max_error = 0.0_mytype
  real(mytype) :: component_error_x = 0.0_mytype
  real(mytype) :: component_error_y = 0.0_mytype
  real(mytype) :: component_error_z = 0.0_mytype
  real(mytype) :: max_abs_rhs_increment = 0.0_mytype
  real(mytype) :: max_abs_rhs_increment_norm_after_clear = 0.0_mytype

contains

  subroutine stage14_rhs_accumulator_init(nx, ny, nz)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer :: alloc_status

    call stage14_rhs_accumulator_finalize()

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
      shape_status = 0
      allocated_status = 0
      call update_overall_status()
      return
    end if

    allocate(rhs_fx_cand(nx,ny,nz), rhs_fy_cand(nx,ny,nz), rhs_fz_cand(nx,ny,nz), &
             rhs_increment_norm(nx,ny,nz), rhs_increment_valid_flag(nx,ny,nz), &
             stat=alloc_status)
    if (alloc_status /= 0) then
      call stage14_rhs_accumulator_finalize()
      shape_status = 0
      allocated_status = 0
      call update_overall_status()
      return
    end if

    nx_rhs_buf = nx
    ny_rhs_buf = ny
    nz_rhs_buf = nz
    shape_status = 1
    allocated_status = 1

    rhs_fx_cand = 0.0_mytype
    rhs_fy_cand = 0.0_mytype
    rhs_fz_cand = 0.0_mytype
    rhs_increment_norm = 0.0_mytype
    rhs_increment_valid_flag = 1
    zero_initialization_status = logical_to_int(maxval(abs(rhs_fx_cand)) == 0.0_mytype .and. &
                                                maxval(abs(rhs_fy_cand)) == 0.0_mytype .and. &
                                                maxval(abs(rhs_fz_cand)) == 0.0_mytype .and. &
                                                maxval(abs(rhs_increment_norm)) == 0.0_mytype .and. &
                                                minval(rhs_increment_valid_flag) == 1 .and. &
                                                maxval(rhs_increment_valid_flag) == 1)
    clear_status = 0
    call update_finite_status()
    call update_overall_status()
  end subroutine stage14_rhs_accumulator_init

  subroutine stage14_rhs_accumulator_clear()
    if (.not. stage14_rhs_accumulator_is_allocated()) then
      clear_status = 0
      call update_overall_status()
      return
    end if

    rhs_fx_cand = 0.0_mytype
    rhs_fy_cand = 0.0_mytype
    rhs_fz_cand = 0.0_mytype
    rhs_increment_norm = 0.0_mytype
    rhs_increment_valid_flag = 1
    max_abs_rhs_increment_norm_after_clear = maxval(abs(rhs_increment_norm))
    clear_status = logical_to_int(maxval(abs(rhs_fx_cand)) == 0.0_mytype .and. &
                                  maxval(abs(rhs_fy_cand)) == 0.0_mytype .and. &
                                  maxval(abs(rhs_fz_cand)) == 0.0_mytype .and. &
                                  max_abs_rhs_increment_norm_after_clear == 0.0_mytype .and. &
                                  minval(rhs_increment_valid_flag) == 1 .and. &
                                  maxval(rhs_increment_valid_flag) == 1)
    call update_finite_status()
    call update_overall_status()
  end subroutine stage14_rhs_accumulator_clear

  subroutine stage14_rhs_accumulator_finalize()
    if (allocated(rhs_fx_cand)) deallocate(rhs_fx_cand)
    if (allocated(rhs_fy_cand)) deallocate(rhs_fy_cand)
    if (allocated(rhs_fz_cand)) deallocate(rhs_fz_cand)
    if (allocated(rhs_increment_norm)) deallocate(rhs_increment_norm)
    if (allocated(rhs_increment_valid_flag)) deallocate(rhs_increment_valid_flag)
    nx_rhs_buf = 0
    ny_rhs_buf = 0
    nz_rhs_buf = 0
    allocated_status = 0
    shape_status = 0
    zero_initialization_status = 0
    clear_status = 0
    call update_overall_status()
  end subroutine stage14_rhs_accumulator_finalize

  logical function stage14_rhs_accumulator_is_allocated()
    stage14_rhs_accumulator_is_allocated = allocated(rhs_fx_cand) .and. allocated(rhs_fy_cand) .and. &
                                           allocated(rhs_fz_cand) .and. allocated(rhs_increment_norm) .and. &
                                           allocated(rhs_increment_valid_flag)
  end function stage14_rhs_accumulator_is_allocated

  subroutine stage14_rhs_accumulator_get_shape(nx, ny, nz)
    integer, intent(out) :: nx
    integer, intent(out) :: ny
    integer, intent(out) :: nz
    nx = nx_rhs_buf
    ny = ny_rhs_buf
    nz = nz_rhs_buf
  end subroutine stage14_rhs_accumulator_get_shape

  subroutine stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, lambda_14)
    real(mytype), intent(in) :: fx_cand(:,:,:)
    real(mytype), intent(in) :: fy_cand(:,:,:)
    real(mytype), intent(in) :: fz_cand(:,:,:)
    real(mytype), intent(in) :: lambda_14
    real(mytype) :: max_error
    real(mytype) :: tolerance

    tolerance = 1.0e-12_mytype

    if (.not. stage14_rhs_accumulator_is_allocated()) then
      call update_overall_status()
      return
    end if

    if (size(fx_cand,1) /= nx_rhs_buf .or. size(fx_cand,2) /= ny_rhs_buf .or. &
        size(fx_cand,3) /= nz_rhs_buf .or. size(fy_cand,1) /= nx_rhs_buf .or. &
        size(fy_cand,2) /= ny_rhs_buf .or. size(fy_cand,3) /= nz_rhs_buf .or. &
        size(fz_cand,1) /= nx_rhs_buf .or. size(fz_cand,2) /= ny_rhs_buf .or. &
        size(fz_cand,3) /= nz_rhs_buf) then
      shape_status = 0
      call update_overall_status()
      return
    end if

    rhs_fx_cand = lambda_14 * fx_cand
    rhs_fy_cand = lambda_14 * fy_cand
    rhs_fz_cand = lambda_14 * fz_cand
    rhs_increment_norm = sqrt(rhs_fx_cand * rhs_fx_cand + &
                              rhs_fy_cand * rhs_fy_cand + &
                              rhs_fz_cand * rhs_fz_cand)
    rhs_increment_valid_flag = 1
    max_abs_rhs_increment = max(maxval(abs(rhs_fx_cand)), &
                                max(maxval(abs(rhs_fy_cand)), maxval(abs(rhs_fz_cand))))

    if (lambda_14 == 0.0_mytype) then
      lambda_zero_max_abs = max_abs_rhs_increment
      lambda_zero_status = logical_to_int(lambda_zero_max_abs <= tolerance)
    end if

    if (lambda_14 == 1.0_mytype) then
      component_error_x = maxval(abs(rhs_fx_cand - fx_cand))
      component_error_y = maxval(abs(rhs_fy_cand - fy_cand))
      component_error_z = maxval(abs(rhs_fz_cand - fz_cand))
      lambda_one_max_error = max(component_error_x, max(component_error_y, component_error_z))
      lambda_one_scaling_status = logical_to_int(lambda_one_max_error <= tolerance)
      component_scaling_status = logical_to_int(component_error_x <= tolerance .and. &
                                                component_error_y <= tolerance .and. &
                                                component_error_z <= tolerance)
    end if

    if (abs(lambda_14 - 0.1_mytype) <= tolerance) then
      max_error = max(maxval(abs(rhs_fx_cand - lambda_14 * fx_cand)), &
                      max(maxval(abs(rhs_fy_cand - lambda_14 * fy_cand)), &
                          maxval(abs(rhs_fz_cand - lambda_14 * fz_cand))))
      lambda_fractional_max_error = max_error
      lambda_fractional_scaling_status = logical_to_int(lambda_fractional_max_error <= tolerance)
    end if

    call update_finite_status()
    call update_overall_status()
  end subroutine stage14_rhs_accumulator_compute_from_force_density

  subroutine stage14_rhs_accumulator_get_status_values(allocated_status_out, &
                                                       shape_status_out, &
                                                       zero_initialization_status_out, &
                                                       clear_status_out, &
                                                       lambda_zero_status_out, &
                                                       lambda_one_scaling_status_out, &
                                                       lambda_fractional_scaling_status_out, &
                                                       component_scaling_status_out, &
                                                       finite_accumulator_status_out, &
                                                       rhs_increment_norm_finite_status_out, &
                                                       rhs_increment_valid_flag_status_out, &
                                                       no_rhs_addition_status_out, &
                                                       no_production_hook_status_out, &
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
                                                       rhs_accumulator_status_out)
    integer, intent(out) :: allocated_status_out
    integer, intent(out) :: shape_status_out
    integer, intent(out) :: zero_initialization_status_out
    integer, intent(out) :: clear_status_out
    integer, intent(out) :: lambda_zero_status_out
    integer, intent(out) :: lambda_one_scaling_status_out
    integer, intent(out) :: lambda_fractional_scaling_status_out
    integer, intent(out) :: component_scaling_status_out
    integer, intent(out) :: finite_accumulator_status_out
    integer, intent(out) :: rhs_increment_norm_finite_status_out
    integer, intent(out) :: rhs_increment_valid_flag_status_out
    integer, intent(out) :: no_rhs_addition_status_out
    integer, intent(out) :: no_production_hook_status_out
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
    integer, intent(out) :: rhs_accumulator_status_out

    allocated_status_out = allocated_status
    shape_status_out = shape_status
    zero_initialization_status_out = zero_initialization_status
    clear_status_out = clear_status
    lambda_zero_status_out = lambda_zero_status
    lambda_one_scaling_status_out = lambda_one_scaling_status
    lambda_fractional_scaling_status_out = lambda_fractional_scaling_status
    component_scaling_status_out = component_scaling_status
    finite_accumulator_status_out = finite_accumulator_status
    rhs_increment_norm_finite_status_out = rhs_increment_norm_finite_status
    rhs_increment_valid_flag_status_out = rhs_increment_valid_flag_status
    no_rhs_addition_status_out = no_rhs_addition_status
    no_production_hook_status_out = no_production_hook_status
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
    rhs_accumulator_status_out = rhs_accumulator_status
  end subroutine stage14_rhs_accumulator_get_status_values

  subroutine stage14_rhs_accumulator_write_diagnostics(filename, requested_flag, &
                                                       rhs_injection_enabled_flag, &
                                                       injection_gain_finite_status)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: requested_flag
    integer, intent(in) :: rhs_injection_enabled_flag
    integer, intent(in) :: injection_gain_finite_status
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file=filename, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return

    write(unit_id, '(A,1X,I0)') 'stage14_1_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage14_1_rhs_injection_enabled_flag', rhs_injection_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage14_1_injection_gain_finite_status', injection_gain_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_allocated_status', allocated_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_shape_status', shape_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_zero_initialization_status', zero_initialization_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_clear_status', clear_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_lambda_zero_status', lambda_zero_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_lambda_one_scaling_status', lambda_one_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_lambda_fractional_scaling_status', lambda_fractional_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_component_scaling_status', component_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_finite_accumulator_status', finite_accumulator_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_rhs_increment_norm_finite_status', rhs_increment_norm_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_rhs_increment_valid_flag_status', rhs_increment_valid_flag_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_rhs_addition_status', no_rhs_addition_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_production_hook_status', no_production_hook_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_channel_forcing_modification_status', &
                                 no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_feedback_application_status', no_feedback_application_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_twoway_force_status', no_twoway_force_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_rhs_accumulator_status', rhs_accumulator_status
    write(unit_id, '(A,1X,I0)') 'stage14_1_nx', nx_rhs_buf
    write(unit_id, '(A,1X,I0)') 'stage14_1_ny', ny_rhs_buf
    write(unit_id, '(A,1X,I0)') 'stage14_1_nz', nz_rhs_buf
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_lambda_zero_max_abs', lambda_zero_max_abs
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_lambda_one_max_error', lambda_one_max_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_lambda_fractional_max_error', lambda_fractional_max_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_component_error_x', component_error_x
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_component_error_y', component_error_y
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_component_error_z', component_error_z
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_max_abs_rhs_increment', max_abs_rhs_increment
    write(unit_id, '(A,1X,ES24.16)') 'stage14_1_max_abs_rhs_increment_norm_after_clear', &
                                      max_abs_rhs_increment_norm_after_clear
    close(unit_id)
  end subroutine stage14_rhs_accumulator_write_diagnostics

  subroutine update_finite_status()
    logical :: finite_fields
    logical :: finite_norm
    logical :: valid_flags

    if (.not. stage14_rhs_accumulator_is_allocated()) then
      finite_accumulator_status = 0
      rhs_increment_norm_finite_status = 0
      rhs_increment_valid_flag_status = 0
      return
    end if

    finite_fields = all(finite_real(rhs_fx_cand)) .and. &
                    all(finite_real(rhs_fy_cand)) .and. &
                    all(finite_real(rhs_fz_cand))
    finite_norm = all(finite_real(rhs_increment_norm))
    valid_flags = all(rhs_increment_valid_flag == 0 .or. rhs_increment_valid_flag == 1)
    finite_accumulator_status = logical_to_int(finite_fields)
    rhs_increment_norm_finite_status = logical_to_int(finite_norm)
    rhs_increment_valid_flag_status = logical_to_int(valid_flags)
  end subroutine update_finite_status

  subroutine update_overall_status()
    rhs_accumulator_status = logical_to_int(allocated_status == 1 .and. &
                                            shape_status == 1 .and. &
                                            zero_initialization_status == 1 .and. &
                                            clear_status == 1 .and. &
                                            lambda_zero_status == 1 .and. &
                                            lambda_one_scaling_status == 1 .and. &
                                            lambda_fractional_scaling_status == 1 .and. &
                                            component_scaling_status == 1 .and. &
                                            finite_accumulator_status == 1 .and. &
                                            rhs_increment_norm_finite_status == 1 .and. &
                                            rhs_increment_valid_flag_status == 1 .and. &
                                            no_rhs_addition_status == 1 .and. &
                                            no_production_hook_status == 1 .and. &
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

end module fibre_stage14_rhs_accumulator
