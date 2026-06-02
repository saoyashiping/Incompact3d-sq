module fibre_stage13_production_force_density_candidate
  use decomp_2d_constants, only : mytype
  use fibre_stage13_config, only : stage13_requested, stage13_readonly_mode, stage13_spreading_readonly_mode
  use fibre_stage13_spreading_kernel, only : stage13_spreading_kernel_init, stage13_spreading_kernel_clear, &
       stage13_spreading_kernel_spread_controlled
  use fibre_stage13_volume_normalization_audit, only : stage13_volume_normalization_audit_init, &
       stage13_volume_normalization_audit_check_volumes, stage13_volume_normalization_audit_integrate_force_density
  use fibre_stage13_conservation_sign_audit, only : stage13_conservation_sign_audit_init, &
       stage13_conservation_sign_audit_compute_expected_forces, &
       stage13_conservation_sign_audit_check_integrated_force_sign
  implicit none
  private

  real(mytype), parameter :: force_conservation_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: wrong_sign_min_separation = 1.0e-8_mytype
  integer, parameter :: nx_diag = 4
  integer, parameter :: ny_diag = 4
  integer, parameter :: nz_diag = 4

  integer :: requested_flag = 0
  integer :: readonly_mode_status = 0
  integer :: spreading_readonly_status = 0
  integer :: hook_initialized_status = 0
  integer :: hook_sample_called_status = 0
  integer :: sampled_velocity_available_status = 0
  integer :: force_density_candidate_computed_status = 0
  integer :: force_density_candidate_finite_status = 0
  integer :: force_density_norm_finite_status = 0
  integer :: integrated_force_finite_status = 0
  integer :: integrated_force_conservation_status = 0
  integer :: spreading_input_sign_status = 0
  integer :: wrong_sign_rejection_status = 0
  integer :: field_modified_status = 0
  integer :: rhs_modified_status = 0
  integer :: no_rhs_injection_status = 1
  integer :: no_production_ibm_forcing_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: production_force_density_candidate_status = 0
  integer :: sample_count = 0

  real(mytype) :: sample_sum_u = 0.0_mytype
  real(mytype) :: sample_sum_v = 0.0_mytype
  real(mytype) :: sample_sum_w = 0.0_mytype
  real(mytype) :: force_density_l2 = 0.0_mytype
  real(mytype) :: integrated_force(3) = 0.0_mytype
  real(mytype) :: expected_fibre_to_fluid(3) = 0.0_mytype
  real(mytype) :: integrated_force_error_l2 = 0.0_mytype
  real(mytype) :: wrong_sign_separation = 0.0_mytype
  real(mytype) :: max_abs_force_density = 0.0_mytype
  real(mytype) :: max_abs_force_density_norm = 0.0_mytype

  public :: stage13_production_force_density_candidate_init
  public :: stage13_production_force_density_candidate_sample
  public :: stage13_production_force_density_candidate_finalize
  public :: stage13_production_force_density_candidate_get_status_values
  public :: stage13_production_force_density_candidate_write_diagnostics

contains

  subroutine stage13_production_force_density_candidate_init()
    requested_flag = merge(1, 0, stage13_requested())
    readonly_mode_status = merge(1, 0, stage13_readonly_mode())
    spreading_readonly_status = merge(1, 0, stage13_spreading_readonly_mode())
    hook_initialized_status = 1
    hook_sample_called_status = 0
    sampled_velocity_available_status = 0
    force_density_candidate_computed_status = 0
    force_density_candidate_finite_status = 0
    force_density_norm_finite_status = 0
    integrated_force_finite_status = 0
    integrated_force_conservation_status = 0
    spreading_input_sign_status = 0
    wrong_sign_rejection_status = 0
    field_modified_status = 0
    rhs_modified_status = 0
    no_rhs_injection_status = 1
    no_production_ibm_forcing_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    sample_count = 0
    sample_sum_u = 0.0_mytype
    sample_sum_v = 0.0_mytype
    sample_sum_w = 0.0_mytype
    force_density_l2 = 0.0_mytype
    integrated_force(:) = 0.0_mytype
    expected_fibre_to_fluid(:) = 0.0_mytype
    integrated_force_error_l2 = 0.0_mytype
    wrong_sign_separation = 0.0_mytype
    max_abs_force_density = 0.0_mytype
    max_abs_force_density_norm = 0.0_mytype
    call stage13_spreading_kernel_init()
    call stage13_volume_normalization_audit_init()
    call stage13_conservation_sign_audit_init()
    call update_status()
  end subroutine stage13_production_force_density_candidate_init

  subroutine stage13_production_force_density_candidate_sample(ux, uy, uz)
    real(mytype), intent(in) :: ux(:,:,:)
    real(mytype), intent(in) :: uy(:,:,:)
    real(mytype), intent(in) :: uz(:,:,:)
    integer :: i
    integer :: j
    integer :: k
    integer :: i0
    integer :: j0
    integer :: k0
    real(mytype) :: x_eul(nx_diag)
    real(mytype) :: y_eul(ny_diag)
    real(mytype) :: z_eul(nz_diag)
    real(mytype) :: volume_eul(nx_diag,ny_diag,nz_diag)
    real(mytype) :: fx_cand(nx_diag,ny_diag,nz_diag)
    real(mytype) :: fy_cand(nx_diag,ny_diag,nz_diag)
    real(mytype) :: fz_cand(nx_diag,ny_diag,nz_diag)
    real(mytype) :: force_density_norm(nx_diag,ny_diag,nz_diag)
    real(mytype) :: x_lag(1)
    real(mytype) :: y_lag(1)
    real(mytype) :: z_lag(1)
    real(mytype) :: ds_lag(1)
    real(mytype) :: f_fluid_to_fibre(1,3)
    real(mytype) :: f_fibre_to_fluid(1,3)
    real(mytype) :: expected_fluid_to_fibre(3)
    real(mytype) :: local_integrated(3)
    real(mytype) :: fibre_to_fluid_error
    real(mytype) :: wrong_sign_error
    real(mytype) :: min_volume
    real(mytype) :: max_volume
    real(mytype) :: sample_velocity(3)
    real(mytype) :: prescribed_velocity(3)
    real(mytype) :: slip(3)
    real(mytype) :: alpha

    hook_sample_called_status = 1
    if (size(ux) <= 0 .or. size(uy) <= 0 .or. size(uz) <= 0) then
      call update_status()
      return
    end if
    if (any(shape(ux) /= shape(uy)) .or. any(shape(ux) /= shape(uz))) then
      call update_status()
      return
    end if

    ! Use a fixed low-index diagnostic control point on rank 0 rather than the
    ! local subdomain centre.  The local-centre choice changes physical location
    ! when the MPI decomposition changes and regresses Stage 14.8 np=1/2/4
    ! force-density consistency.  The +2 offset avoids boundary/end-point values
    ! while remaining valid on very small local arrays.
    i0 = min(lbound(ux, 1) + 2, ubound(ux, 1))
    j0 = min(lbound(ux, 2) + 2, ubound(ux, 2))
    k0 = min(lbound(ux, 3) + 2, ubound(ux, 3))
    sample_velocity(1) = ux(i0, j0, k0)
    sample_velocity(2) = uy(i0, j0, k0)
    sample_velocity(3) = uz(i0, j0, k0)
    prescribed_velocity(:) = 0.0_mytype
    alpha = 1.0_mytype
    slip(:) = sample_velocity(:) - prescribed_velocity(:)

    f_fluid_to_fibre(1,1) = alpha * slip(1)
    f_fluid_to_fibre(1,2) = alpha * slip(2)
    f_fluid_to_fibre(1,3) = alpha * slip(3)
    f_fibre_to_fluid(:,:) = -f_fluid_to_fibre(:,:)

    do i = 1, nx_diag
      x_eul(i) = real(i - 1, mytype)
    end do
    do j = 1, ny_diag
      y_eul(j) = real(j - 1, mytype)
    end do
    do k = 1, nz_diag
      z_eul(k) = real(k - 1, mytype)
    end do
    volume_eul(:,:,:) = 1.0_mytype
    fx_cand(:,:,:) = 0.0_mytype
    fy_cand(:,:,:) = 0.0_mytype
    fz_cand(:,:,:) = 0.0_mytype
    force_density_norm(:,:,:) = 0.0_mytype
    x_lag(1) = 1.5_mytype
    y_lag(1) = 1.5_mytype
    z_lag(1) = 1.5_mytype
    ds_lag(1) = 1.0_mytype

    call stage13_volume_normalization_audit_check_volumes(volume_eul, min_volume, max_volume)
    call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
    call stage13_spreading_kernel_spread_controlled(x_lag, y_lag, z_lag, ds_lag, f_fibre_to_fluid, &
         x_eul, y_eul, z_eul, volume_eul, fx_cand, fy_cand, fz_cand, force_density_norm)
    call stage13_volume_normalization_audit_integrate_force_density(fx_cand, fy_cand, fz_cand, volume_eul, &
         local_integrated)
    call stage13_conservation_sign_audit_compute_expected_forces(f_fluid_to_fibre, f_fibre_to_fluid, ds_lag, &
         expected_fluid_to_fibre, expected_fibre_to_fluid)
    call stage13_conservation_sign_audit_check_integrated_force_sign(local_integrated, expected_fluid_to_fibre, &
         expected_fibre_to_fluid, fibre_to_fluid_error, wrong_sign_error)

    sample_count = sample_count + 1
    sample_sum_u = sample_sum_u + sample_velocity(1)
    sample_sum_v = sample_sum_v + sample_velocity(2)
    sample_sum_w = sample_sum_w + sample_velocity(3)
    force_density_l2 = sqrt(sum(fx_cand * fx_cand + fy_cand * fy_cand + fz_cand * fz_cand))
    integrated_force(:) = local_integrated(:)
    integrated_force_error_l2 = fibre_to_fluid_error
    wrong_sign_separation = wrong_sign_error - fibre_to_fluid_error
    max_abs_force_density = max(maxval(abs(fx_cand)), maxval(abs(fy_cand)), maxval(abs(fz_cand)))
    max_abs_force_density_norm = maxval(abs(force_density_norm))

    sampled_velocity_available_status = merge(1, 0, all_finite_vector(sample_velocity))
    force_density_candidate_computed_status = 1
    force_density_candidate_finite_status = merge(1, 0, all_finite_rank3(fx_cand) .and. &
         all_finite_rank3(fy_cand) .and. all_finite_rank3(fz_cand))
    force_density_norm_finite_status = merge(1, 0, all_finite_rank3(force_density_norm))
    integrated_force_finite_status = merge(1, 0, all_finite_vector(local_integrated))
    integrated_force_conservation_status = merge(1, 0, integrated_force_error_l2 <= force_conservation_abs_tol)
    spreading_input_sign_status = merge(1, 0, vector_l2(f_fibre_to_fluid(1,:) + f_fluid_to_fibre(1,:)) <= &
         force_conservation_abs_tol)
    wrong_sign_rejection_status = merge(1, 0, integrated_force_conservation_status == 1 .and. &
         wrong_sign_separation >= wrong_sign_min_separation)
    field_modified_status = 0
    rhs_modified_status = 0
    call update_status()
  end subroutine stage13_production_force_density_candidate_sample

  subroutine stage13_production_force_density_candidate_finalize()
    call update_status()
    if (rank0_write_allowed()) call stage13_production_force_density_candidate_write_diagnostics()
  end subroutine stage13_production_force_density_candidate_finalize

  subroutine stage13_production_force_density_candidate_get_status_values(requested_out, readonly_out, &
       spreading_readonly_out, initialized_out, sample_called_out, sampled_velocity_out, candidate_computed_out, &
       candidate_finite_out, norm_finite_out, integrated_finite_out, conservation_out, spreading_sign_out, &
       wrong_sign_out, field_modified_out, rhs_modified_out, no_rhs_out, no_ibm_out, no_feedback_out, &
       no_twoway_out, no_structure_out, candidate_status_out)
    integer, intent(out) :: requested_out
    integer, intent(out) :: readonly_out
    integer, intent(out) :: spreading_readonly_out
    integer, intent(out) :: initialized_out
    integer, intent(out) :: sample_called_out
    integer, intent(out) :: sampled_velocity_out
    integer, intent(out) :: candidate_computed_out
    integer, intent(out) :: candidate_finite_out
    integer, intent(out) :: norm_finite_out
    integer, intent(out) :: integrated_finite_out
    integer, intent(out) :: conservation_out
    integer, intent(out) :: spreading_sign_out
    integer, intent(out) :: wrong_sign_out
    integer, intent(out) :: field_modified_out
    integer, intent(out) :: rhs_modified_out
    integer, intent(out) :: no_rhs_out
    integer, intent(out) :: no_ibm_out
    integer, intent(out) :: no_feedback_out
    integer, intent(out) :: no_twoway_out
    integer, intent(out) :: no_structure_out
    integer, intent(out) :: candidate_status_out

    call update_status()
    requested_out = requested_flag
    readonly_out = readonly_mode_status
    spreading_readonly_out = spreading_readonly_status
    initialized_out = hook_initialized_status
    sample_called_out = hook_sample_called_status
    sampled_velocity_out = sampled_velocity_available_status
    candidate_computed_out = force_density_candidate_computed_status
    candidate_finite_out = force_density_candidate_finite_status
    norm_finite_out = force_density_norm_finite_status
    integrated_finite_out = integrated_force_finite_status
    conservation_out = integrated_force_conservation_status
    spreading_sign_out = spreading_input_sign_status
    wrong_sign_out = wrong_sign_rejection_status
    field_modified_out = field_modified_status
    rhs_modified_out = rhs_modified_status
    no_rhs_out = no_rhs_injection_status
    no_ibm_out = no_production_ibm_forcing_status
    no_feedback_out = no_feedback_application_status
    no_twoway_out = no_twoway_force_status
    no_structure_out = no_structure_advance_status
    candidate_status_out = production_force_density_candidate_status
  end subroutine stage13_production_force_density_candidate_get_status_values

  subroutine stage13_production_force_density_candidate_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat'
    end if

    call update_status()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit,'(A,1X,I0)') 'stage13_6_requested_flag', requested_flag
    write(io_unit,'(A,1X,I0)') 'stage13_6_readonly_mode_status', readonly_mode_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_spreading_readonly_status', spreading_readonly_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_hook_initialized_status', hook_initialized_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_hook_sample_called_status', hook_sample_called_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_sampled_velocity_available_status', sampled_velocity_available_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_force_density_candidate_computed_status', &
         force_density_candidate_computed_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_force_density_candidate_finite_status', &
         force_density_candidate_finite_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_force_density_norm_finite_status', force_density_norm_finite_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_integrated_force_finite_status', integrated_force_finite_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_integrated_force_conservation_status', &
         integrated_force_conservation_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_spreading_input_sign_status', spreading_input_sign_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_wrong_sign_rejection_status', wrong_sign_rejection_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_field_modified_status', field_modified_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_rhs_modified_status', rhs_modified_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_no_rhs_injection_status', no_rhs_injection_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_no_feedback_application_status', no_feedback_application_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_no_twoway_force_status', no_twoway_force_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_no_structure_advance_status', no_structure_advance_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_production_force_density_candidate_status', &
         production_force_density_candidate_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_sample_count', sample_count
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_sample_sum_u', sample_sum_u
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_sample_sum_v', sample_sum_v
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_sample_sum_w', sample_sum_w
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_force_density_l2', force_density_l2
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_integrated_force_x', integrated_force(1)
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_integrated_force_y', integrated_force(2)
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_integrated_force_z', integrated_force(3)
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_expected_fibre_to_fluid_x', expected_fibre_to_fluid(1)
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_expected_fibre_to_fluid_y', expected_fibre_to_fluid(2)
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_expected_fibre_to_fluid_z', expected_fibre_to_fluid(3)
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_integrated_force_error_l2', integrated_force_error_l2
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_wrong_sign_separation', wrong_sign_separation
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_max_abs_force_density', max_abs_force_density
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_max_abs_force_density_norm', max_abs_force_density_norm
    close(io_unit)
  end subroutine stage13_production_force_density_candidate_write_diagnostics

  subroutine update_status()
    production_force_density_candidate_status = merge(1, 0, requested_flag == 1 .and. &
      readonly_mode_status == 1 .and. spreading_readonly_status == 1 .and. hook_initialized_status == 1 .and. &
      hook_sample_called_status == 1 .and. sampled_velocity_available_status == 1 .and. &
      force_density_candidate_computed_status == 1 .and. force_density_candidate_finite_status == 1 .and. &
      force_density_norm_finite_status == 1 .and. integrated_force_finite_status == 1 .and. &
      integrated_force_conservation_status == 1 .and. spreading_input_sign_status == 1 .and. &
      wrong_sign_rejection_status == 1 .and. field_modified_status == 0 .and. rhs_modified_status == 0 .and. &
      no_rhs_injection_status == 1 .and. no_production_ibm_forcing_status == 1 .and. &
      no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
      no_structure_advance_status == 1)
  end subroutine update_status

  real(mytype) function vector_l2(values)
    real(mytype), intent(in) :: values(3)

    vector_l2 = sqrt(sum(values * values))
  end function vector_l2

  logical function all_finite_vector(values)
    real(mytype), intent(in) :: values(:)
    integer :: i

    all_finite_vector = .true.
    do i = 1, size(values)
      if (.not. is_finite(values(i))) all_finite_vector = .false.
    end do
  end function all_finite_vector

  logical function all_finite_rank3(values)
    real(mytype), intent(in) :: values(:,:,:)
    integer :: i
    integer :: j
    integer :: k

    all_finite_rank3 = .true.
    do k = lbound(values, 3), ubound(values, 3)
      do j = lbound(values, 2), ubound(values, 2)
        do i = lbound(values, 1), ubound(values, 1)
          if (.not. is_finite(values(i,j,k))) all_finite_rank3 = .false.
        end do
      end do
    end do
  end function all_finite_rank3

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(value))
  end function is_finite

  logical function rank0_write_allowed()
    character(len=32) :: value
    integer :: status
    integer :: ios
    integer :: rank_value

    rank0_write_allowed = .true.
    call get_environment_variable('OMPI_COMM_WORLD_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('PMI_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('MPI_RANK', value=value, status=status)
    if (status == 0) then
      read(value, *, iostat=ios) rank_value
      if (ios == 0) rank0_write_allowed = (rank_value == 0)
    end if
  end function rank0_write_allowed

end module fibre_stage13_production_force_density_candidate
