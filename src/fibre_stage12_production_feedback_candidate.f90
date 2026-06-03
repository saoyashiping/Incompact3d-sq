module fibre_stage12_production_feedback_candidate
  use decomp_2d_constants, only : mytype
  use fibre_stage12_config, only : stage12_requested, stage12_readonly_mode, stage12_get_feedback_gain
  implicit none
  private

  real(mytype), parameter :: diagnostics_abs_tol = 1.0e-12_mytype

  integer :: requested_flag = 0
  integer :: readonly_mode_status = 0
  integer :: hook_initialized_status = 0
  integer :: hook_sample_called_status = 0
  integer :: sampled_velocity_available_status = 0
  integer :: force_candidate_computed_status = 0
  integer :: force_candidate_finite_status = 0
  integer :: force_norm_finite_status = 0
  integer :: power_diagnostics_finite_status = 0
  integer :: action_reaction_status = 0
  integer :: pair_power_consistency_status = 0
  integer :: field_modified_status = 0
  integer :: rhs_modified_status = 0
  integer :: no_eulerian_force_density_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: production_feedback_candidate_status = 0
  integer :: sample_count = 0

  real(mytype) :: sample_sum_u = 0.0_mytype
  real(mytype) :: sample_sum_v = 0.0_mytype
  real(mytype) :: sample_sum_w = 0.0_mytype
  real(mytype) :: force_l2 = 0.0_mytype
  real(mytype) :: p_slip = 0.0_mytype
  real(mytype) :: p_structure = 0.0_mytype
  real(mytype) :: p_fluid = 0.0_mytype
  real(mytype) :: p_pair = 0.0_mytype
  real(mytype) :: pair_power_consistency_error = 0.0_mytype
  real(mytype) :: action_reaction_error = 0.0_mytype

  public :: stage12_production_feedback_candidate_init
  public :: stage12_production_feedback_candidate_sample
  public :: stage12_production_feedback_candidate_finalize
  public :: stage12_production_feedback_candidate_get_status_values
  public :: stage12_production_feedback_candidate_write_diagnostics

contains

  subroutine stage12_production_feedback_candidate_init()
    requested_flag = merge(1, 0, stage12_requested())
    readonly_mode_status = merge(1, 0, stage12_readonly_mode())
    hook_initialized_status = 1
    hook_sample_called_status = 0
    sampled_velocity_available_status = 0
    force_candidate_computed_status = 0
    force_candidate_finite_status = 0
    force_norm_finite_status = 0
    power_diagnostics_finite_status = 0
    action_reaction_status = 0
    pair_power_consistency_status = 0
    field_modified_status = 0
    rhs_modified_status = 0
    no_eulerian_force_density_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    sample_count = 0
    sample_sum_u = 0.0_mytype
    sample_sum_v = 0.0_mytype
    sample_sum_w = 0.0_mytype
    force_l2 = 0.0_mytype
    p_slip = 0.0_mytype
    p_structure = 0.0_mytype
    p_fluid = 0.0_mytype
    p_pair = 0.0_mytype
    pair_power_consistency_error = 0.0_mytype
    action_reaction_error = 0.0_mytype
    call update_status()
  end subroutine stage12_production_feedback_candidate_init

  subroutine stage12_production_feedback_candidate_sample(ux, uy, uz)
    real(mytype), intent(in) :: ux(:,:,:)
    real(mytype), intent(in) :: uy(:,:,:)
    real(mytype), intent(in) :: uz(:,:,:)
    integer :: i0
    integer :: j0
    integer :: k0
    real(mytype) :: u_f(3)
    real(mytype) :: v_f(3)
    real(mytype) :: slip(3)
    real(mytype) :: f_fluid_to_fibre(3)
    real(mytype) :: f_fibre_to_fluid(3)
    real(mytype) :: alpha
    real(mytype) :: local_force_l2
    real(mytype) :: local_p_slip
    real(mytype) :: local_p_structure
    real(mytype) :: local_p_fluid
    real(mytype) :: local_p_pair
    real(mytype) :: local_action_reaction_error
    real(mytype) :: local_pair_power_error

    hook_sample_called_status = 1
    if (size(ux) <= 0 .or. size(uy) <= 0 .or. size(uz) <= 0) then
      call update_status()
      return
    end if
    if (any(shape(ux) /= shape(uy)) .or. any(shape(ux) /= shape(uz))) then
      call update_status()
      return
    end if

    i0 = (lbound(ux, 1) + ubound(ux, 1)) / 2
    j0 = (lbound(ux, 2) + ubound(ux, 2)) / 2
    k0 = (lbound(ux, 3) + ubound(ux, 3)) / 2

    u_f(1) = ux(i0, j0, k0)
    u_f(2) = uy(i0, j0, k0)
    u_f(3) = uz(i0, j0, k0)
    v_f(:) = 0.0_mytype
    alpha = stage12_get_feedback_gain()

    slip(:) = u_f(:) - v_f(:)
    f_fluid_to_fibre(:) = alpha * slip(:)
    f_fibre_to_fluid(:) = -f_fluid_to_fibre(:)
    local_force_l2 = sqrt(sum(f_fluid_to_fibre(:) * f_fluid_to_fibre(:)))
    local_p_slip = dot_product(f_fluid_to_fibre(:), slip(:))
    local_p_structure = dot_product(f_fluid_to_fibre(:), v_f(:))
    local_p_fluid = dot_product(f_fibre_to_fluid(:), u_f(:))
    local_p_pair = local_p_structure + local_p_fluid
    local_pair_power_error = abs(local_p_pair + local_p_slip)
    local_action_reaction_error = maxval(abs(f_fluid_to_fibre(:) + f_fibre_to_fluid(:)))

    sample_count = sample_count + 1
    sample_sum_u = sample_sum_u + u_f(1)
    sample_sum_v = sample_sum_v + u_f(2)
    sample_sum_w = sample_sum_w + u_f(3)
    force_l2 = force_l2 + local_force_l2
    p_slip = p_slip + local_p_slip
    p_structure = p_structure + local_p_structure
    p_fluid = p_fluid + local_p_fluid
    p_pair = p_structure + p_fluid
    pair_power_consistency_error = max(pair_power_consistency_error, local_pair_power_error)
    action_reaction_error = max(action_reaction_error, local_action_reaction_error)

    sampled_velocity_available_status = merge(1, 0, all_finite(u_f))
    force_candidate_computed_status = 1
    force_candidate_finite_status = merge(1, 0, all_finite(f_fluid_to_fibre) .and. all_finite(f_fibre_to_fluid))
    force_norm_finite_status = merge(1, 0, is_finite(local_force_l2))
    power_diagnostics_finite_status = merge(1, 0, is_finite(local_p_slip) .and. &
                                            is_finite(local_p_structure) .and. &
                                            is_finite(local_p_fluid) .and. is_finite(local_p_pair))
    action_reaction_status = merge(1, 0, local_action_reaction_error <= diagnostics_abs_tol)
    pair_power_consistency_status = merge(1, 0, local_pair_power_error <= diagnostics_abs_tol)
    field_modified_status = 0
    rhs_modified_status = 0
    call update_status()
  end subroutine stage12_production_feedback_candidate_sample

  subroutine stage12_production_feedback_candidate_finalize()
    call update_status()
    if (rank0_write_allowed()) call stage12_production_feedback_candidate_write_diagnostics()
  end subroutine stage12_production_feedback_candidate_finalize

  subroutine stage12_production_feedback_candidate_get_status_values(requested_out, readonly_out, initialized_out, &
                                                                     sample_called_out, sampled_velocity_out, &
                                                                     force_computed_out, force_finite_out, &
                                                                     force_norm_finite_out, power_finite_out, &
                                                                     action_reaction_out, pair_power_out, &
                                                                     field_modified_out, rhs_modified_out, &
                                                                     no_eulerian_force_density_out, &
                                                                     no_rhs_injection_out, no_ibm_spreading_out, &
                                                                     no_feedback_application_out, no_twoway_force_out, &
                                                                     no_structure_advance_out, candidate_status_out)
    integer, intent(out) :: requested_out
    integer, intent(out) :: readonly_out
    integer, intent(out) :: initialized_out
    integer, intent(out) :: sample_called_out
    integer, intent(out) :: sampled_velocity_out
    integer, intent(out) :: force_computed_out
    integer, intent(out) :: force_finite_out
    integer, intent(out) :: force_norm_finite_out
    integer, intent(out) :: power_finite_out
    integer, intent(out) :: action_reaction_out
    integer, intent(out) :: pair_power_out
    integer, intent(out) :: field_modified_out
    integer, intent(out) :: rhs_modified_out
    integer, intent(out) :: no_eulerian_force_density_out
    integer, intent(out) :: no_rhs_injection_out
    integer, intent(out) :: no_ibm_spreading_out
    integer, intent(out) :: no_feedback_application_out
    integer, intent(out) :: no_twoway_force_out
    integer, intent(out) :: no_structure_advance_out
    integer, intent(out) :: candidate_status_out

    call update_status()
    requested_out = requested_flag
    readonly_out = readonly_mode_status
    initialized_out = hook_initialized_status
    sample_called_out = hook_sample_called_status
    sampled_velocity_out = sampled_velocity_available_status
    force_computed_out = force_candidate_computed_status
    force_finite_out = force_candidate_finite_status
    force_norm_finite_out = force_norm_finite_status
    power_finite_out = power_diagnostics_finite_status
    action_reaction_out = action_reaction_status
    pair_power_out = pair_power_consistency_status
    field_modified_out = field_modified_status
    rhs_modified_out = rhs_modified_status
    no_eulerian_force_density_out = no_eulerian_force_density_status
    no_rhs_injection_out = no_rhs_injection_status
    no_ibm_spreading_out = no_ibm_spreading_status
    no_feedback_application_out = no_feedback_application_status
    no_twoway_force_out = no_twoway_force_status
    no_structure_advance_out = no_structure_advance_status
    candidate_status_out = production_feedback_candidate_status
  end subroutine stage12_production_feedback_candidate_get_status_values

  subroutine stage12_production_feedback_candidate_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat'
    end if

    call update_status()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit, '(A,1X,I0)') 'stage12_6_requested_flag', requested_flag
    write(io_unit, '(A,1X,I0)') 'stage12_6_readonly_mode_status', readonly_mode_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_hook_initialized_status', hook_initialized_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_hook_sample_called_status', hook_sample_called_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_sampled_velocity_available_status', sampled_velocity_available_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_force_candidate_computed_status', force_candidate_computed_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_force_candidate_finite_status', force_candidate_finite_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_force_norm_finite_status', force_norm_finite_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_power_diagnostics_finite_status', power_diagnostics_finite_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_action_reaction_status', action_reaction_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_pair_power_consistency_status', pair_power_consistency_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_field_modified_status', field_modified_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_rhs_modified_status', rhs_modified_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_no_eulerian_force_density_status', no_eulerian_force_density_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_no_rhs_injection_status', no_rhs_injection_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_no_ibm_spreading_status', no_ibm_spreading_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_no_feedback_application_status', no_feedback_application_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_no_twoway_force_status', no_twoway_force_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_no_structure_advance_status', no_structure_advance_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_production_feedback_candidate_status', &
                                production_feedback_candidate_status
    write(io_unit, '(A,1X,I0)') 'stage12_6_sample_count', sample_count
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_sample_sum_u', sample_sum_u
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_sample_sum_v', sample_sum_v
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_sample_sum_w', sample_sum_w
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_force_l2', force_l2
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_p_slip', p_slip
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_p_structure', p_structure
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_p_fluid', p_fluid
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_p_pair', p_pair
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_pair_power_consistency_error', &
                                        pair_power_consistency_error
    write(io_unit, '(A,1X,ES24.16E3)') 'stage12_6_action_reaction_error', action_reaction_error
    close(io_unit)
  end subroutine stage12_production_feedback_candidate_write_diagnostics

  subroutine update_status()
    production_feedback_candidate_status = merge(1, 0, requested_flag == 1 .and. readonly_mode_status == 1 .and. &
      hook_initialized_status == 1 .and. hook_sample_called_status == 1 .and. &
      sampled_velocity_available_status == 1 .and. force_candidate_computed_status == 1 .and. &
      force_candidate_finite_status == 1 .and. force_norm_finite_status == 1 .and. &
      power_diagnostics_finite_status == 1 .and. action_reaction_status == 1 .and. &
      pair_power_consistency_status == 1 .and. field_modified_status == 0 .and. rhs_modified_status == 0 .and. &
      no_eulerian_force_density_status == 1 .and. no_rhs_injection_status == 1 .and. &
      no_ibm_spreading_status == 1 .and. no_feedback_application_status == 1 .and. &
      no_twoway_force_status == 1 .and. no_structure_advance_status == 1)
  end subroutine update_status

  logical function all_finite(values)
    real(mytype), intent(in) :: values(:)
    integer :: i

    all_finite = .true.
    do i = 1, size(values)
      if (.not. is_finite(values(i))) all_finite = .false.
    end do
  end function all_finite

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

end module fibre_stage12_production_feedback_candidate
