program fibre_prod_p0_closure_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_rhs_adapter, only : fibre_prod_rhs_adapter_apply
  use fibre_prod_main_diagnostics, only : fibre_prod_rhs_signature_type
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init, fibre_prod_main_hook_apply
  use fibre_prod_velocity_bridge, only : fibre_prod_velocity_bridge_env_enabled
  use fibre_prod_state_velocity_attachment, only : fibre_prod_state_velocity_attachment_env_enabled
  use fibre_prod_hydro_input_candidate, only : fibre_prod_hydro_input_candidate_env_enabled
  use fibre_prod_structure_input_handoff, only : fibre_prod_structure_input_handoff_env_enabled
  use fibre_prod_structure_dry_step, only : fibre_prod_structure_dry_step_env_enabled
  use fibre_prod_structure_commit_gate, only : fibre_prod_structure_commit_gate_env_enabled
  use fibre_prod_reaction_force_candidate, only : fibre_prod_reaction_force_candidate_env_enabled
  use fibre_prod_reaction_spreading_buffer, only : fibre_prod_reaction_spreading_buffer_env_enabled
  use fibre_prod_force_buffer_rhs_gate, only : fibre_prod_force_buffer_rhs_gate_env_enabled
  use fibre_prod_synthetic_closed_loop, only : fibre_prod_synthetic_closed_loop_type, &
                                               fibre_prod_synthetic_closed_loop_env_enabled, &
                                               fibre_prod_synthetic_closed_loop_init, &
                                               fibre_prod_synthetic_closed_loop_run, &
                                               fibre_prod_synthetic_closed_loop_get_signature, &
                                               fibre_prod_synthetic_closed_loop_finalize
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_runtime_config_type) :: config
  type(fibre_prod_rhs_signature_type) :: before
  type(fibre_prod_rhs_signature_type) :: after
  type(fibre_prod_synthetic_closed_loop_type) :: loop_a
  type(fibre_prod_synthetic_closed_loop_type) :: loop_b
  real(dp) :: rhs_x(3,3,3), rhs_y(3,3,3), rhs_z(3,3,3)
  real(dp) :: rhs0_x(3,3,3), rhs0_y(3,3,3), rhs0_z(3,3,3)
  real(dp) :: ux(3,3,3), uy(3,3,3), uz(3,3,3)
  real(dp) :: ux0(3,3,3), uy0(3,3,3), uz0(3,3,3)
  real(dp) :: sig_a(16), sig_b(16)
  integer :: modified_cells
  integer :: status

  call require(.not. fibre_prod_velocity_bridge_env_enabled(), 'velocity sampling env gate defaulted on')
  call require(.not. fibre_prod_state_velocity_attachment_env_enabled(), 'state velocity env gate defaulted on')
  call require(.not. fibre_prod_hydro_input_candidate_env_enabled(), 'hydro candidate env gate defaulted on')
  call require(.not. fibre_prod_structure_input_handoff_env_enabled(), 'structure handoff env gate defaulted on')
  call require(.not. fibre_prod_structure_dry_step_env_enabled(), 'structure dry-step env gate defaulted on')
  call require(.not. fibre_prod_structure_commit_gate_env_enabled(), 'structure commit env gate defaulted on')
  call require(.not. fibre_prod_reaction_force_candidate_env_enabled(), 'reaction candidate env gate defaulted on')
  call require(.not. fibre_prod_reaction_spreading_buffer_env_enabled(), 'reaction spreading env gate defaulted on')
  call require(.not. fibre_prod_force_buffer_rhs_gate_env_enabled(), 'force-buffer RHS gate defaulted on')
  call require(.not. fibre_prod_synthetic_closed_loop_env_enabled(), 'synthetic closed-loop env gate defaulted on')

  call fill_fields(rhs0_x, rhs0_y, rhs0_z, ux0, uy0, uz0)
  rhs_x = rhs0_x; rhs_y = rhs0_y; rhs_z = rhs0_z
  ux = ux0; uy = uy0; uz = uz0

  call fibre_prod_runtime_config_default(config)
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
  call require(status == 0, 'disabled no-fibre RHS adapter path returned nonzero status')
  call require(modified_cells == 0, 'disabled no-fibre RHS adapter modified cells')
  call require(all(rhs_x == rhs0_x) .and. all(rhs_y == rhs0_y) .and. all(rhs_z == rhs0_z), &
               'disabled no-fibre RHS adapter changed RHS')
  call require(all(ux == ux0) .and. all(uy == uy0) .and. all(uz == uz0), 'no-fibre default changed velocity')

  call fibre_prod_main_hook_init(status, config)
  call require(status == 0, 'main hook disabled config init failed')
  call fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status)
  call require(status == 0, 'main hook disabled no-fibre apply failed')
  call require(all(rhs_x == rhs0_x) .and. all(rhs_y == rhs0_y) .and. all(rhs_z == rhs0_z), &
               'main hook disabled no-fibre apply changed RHS')

  config%enabled = .true.
  config%lambda_fsi = 0.0_dp
  call fibre_prod_rhs_adapter_apply(config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
  call require(status == 0, 'lambda0 no-contamination adapter path failed')
  call require(all(rhs_x == rhs0_x) .and. all(rhs_y == rhs0_y) .and. all(rhs_z == rhs0_z), &
               'lambda0 no-contamination adapter path changed RHS')

  call run_signature(loop_a, sig_a)
  call run_signature(loop_b, sig_b)
  call require(all(ieee_is_finite(sig_a)) .and. all(ieee_is_finite(sig_b)), 'P0 closure signatures nonfinite')
  call require(maxval(abs(sig_a - sig_b)) <= 1.0e-12_dp, 'P0 closure signature is not repeatable')

  print '(A,ES24.16)', 'P0_13_SIGNATURE sampled_norm=', sig_a(1)
  print '(A,ES24.16)', 'P0_13_SIGNATURE force_buffer_norm=', sig_a(11)
  print '(A,ES24.16)', 'P0_13_SIGNATURE lambda0_rhs_increment_norm=', sig_a(12)
  print *, 'P0_13_P0_CLOSURE_CHECK PASS'

  call fibre_prod_synthetic_closed_loop_finalize(loop_a)
  call fibre_prod_synthetic_closed_loop_finalize(loop_b)

contains

  subroutine run_signature(loop, signature)
    type(fibre_prod_synthetic_closed_loop_type), intent(inout) :: loop
    real(dp), intent(out) :: signature(16)
    integer :: local_status

    call fibre_prod_synthetic_closed_loop_init(loop, 4, 4, 4, 5, local_status)
    call require(local_status == 0, 'closure synthetic init failed')
    call fibre_prod_synthetic_closed_loop_run(loop, 0.0_dp, 2.0_dp, 1.0_dp, 1.0e-4_dp, 1.0_dp, local_status)
    call require(local_status == 0, 'closure synthetic lambda0 run failed')
    call fibre_prod_synthetic_closed_loop_get_signature(loop, signature, local_status)
    call require(local_status == 0, 'closure synthetic signature failed')
  end subroutine run_signature

  subroutine fill_fields(rhs_x, rhs_y, rhs_z, ux, uy, uz)
    real(dp), intent(out) :: rhs_x(:,:,:), rhs_y(:,:,:), rhs_z(:,:,:)
    real(dp), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer :: i, j, k

    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          rhs_x(i,j,k) = real(i + j + k, dp) * 1.0e-2_dp
          rhs_y(i,j,k) = real(i - j + k, dp) * 2.0e-2_dp
          rhs_z(i,j,k) = real(i + 2*j - k, dp) * 3.0e-2_dp
          ux(i,j,k) = real(i, dp) * 0.1_dp
          uy(i,j,k) = real(j, dp) * 0.2_dp
          uz(i,j,k) = real(k, dp) * 0.3_dp
        end do
      end do
    end do
  end subroutine fill_fields

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      print *, 'P0_13_P0_CLOSURE_CHECK FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine require

end program fibre_prod_p0_closure_check
