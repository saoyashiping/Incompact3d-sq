program fibre_prod_structure_solver_check
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, &
                               fibre_prod_state_init_straight, fibre_prod_state_destroy, &
                               fibre_prod_state_is_allocated, fibre_prod_state_all_finite
  use fibre_prod_bending_solver, only : fibre_prod_bending_compute_force, fibre_prod_bending_force_norm
  use fibre_prod_tension_solver, only : fibre_prod_tension_segment_length_residual, &
                                        fibre_prod_tension_max_stretch_error
  use fibre_prod_structure_solver, only : fibre_prod_structure_compute_forces, &
                                          fibre_prod_structure_step_explicit, &
                                          fibre_prod_structure_kinetic_energy, &
                                          fibre_prod_structure_bending_energy_candidate
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_state_type) :: state
  real(dp), allocatable :: external_force(:, :, :)
  real(dp) :: origin(3)
  real(dp) :: direction(3)
  real(dp) :: ds
  real(dp) :: dt
  real(dp) :: rho_tilde
  real(dp) :: gamma
  real(dp) :: residual
  real(dp) :: stretch_error
  real(dp) :: kinetic_energy
  real(dp) :: bending_energy
  real(dp) :: x_before
  real(dp), allocatable :: bend(:, :)
  integer :: status

  ds = 0.125_dp
  dt = 1.0e-3_dp
  rho_tilde = 2.0_dp
  gamma = 0.5_dp
  origin = [0.0_dp, 0.0_dp, 0.0_dp]
  direction = [1.0_dp, 0.0_dp, 0.0_dp]

  call fibre_prod_state_allocate(state, 1, 7, status)
  if (status /= 0 .or. .not. fibre_prod_state_is_allocated(state)) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: allocation', status
    error stop 1
  end if

  call fibre_prod_state_init_straight(state, origin, direction, ds, stat=status)
  if (status /= 0) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: straight init', status
    error stop 2
  end if

  allocate(bend(state%nnode, 3))
  call fibre_prod_bending_compute_force(state%x(1, :, :), gamma, ds, bend, status)
  if (status /= 0 .or. fibre_prod_bending_force_norm(bend) > 1.0e-12_dp) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: straight bending', status, &
             fibre_prod_bending_force_norm(bend)
    error stop 3
  end if

  state%x(1, 4, 2) = 1.0e-3_dp
  call fibre_prod_bending_compute_force(state%x(1, :, :), gamma, ds, bend, status)
  if (status /= 0 .or. fibre_prod_bending_force_norm(bend) >= huge(1.0_dp)) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: perturbed bending finite', status
    error stop 4
  end if
  call fibre_prod_state_init_straight(state, origin, direction, ds, stat=status)
  if (status /= 0) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: straight reinit', status
    error stop 5
  end if

  call fibre_prod_structure_compute_forces(state, gamma, ds, status=status)
  if (status /= 0) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: compute no-force', status
    error stop 4
  end if
  x_before = state%x(1, 4, 1)
  call fibre_prod_structure_step_explicit(state, dt, rho_tilde, status)
  if (status /= 0 .or. abs(state%x(1, 4, 1) - x_before) > 1.0e-12_dp) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: no-force invariance', status
    error stop 5
  end if

  allocate(external_force(state%nfibre, state%nnode, 3))
  external_force = 0.0_dp
  external_force(:, :, 2) = 0.25_dp
  call fibre_prod_structure_compute_forces(state, gamma, ds, external_force, status)
  if (status /= 0) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: compute external', status
    error stop 6
  end if
  call fibre_prod_structure_step_explicit(state, dt, rho_tilde, status)
  if (status /= 0 .or. state%v(1, 4, 2) <= 0.0_dp .or. state%x(1, 4, 2) <= 0.0_dp) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: external response', status, &
             state%v(1, 4, 2), state%x(1, 4, 2)
    error stop 7
  end if

  if (.not. fibre_prod_state_all_finite(state)) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: non-finite state'
    error stop 8
  end if

  residual = fibre_prod_tension_segment_length_residual(state%x(1, :, :), ds, status)
  stretch_error = fibre_prod_tension_max_stretch_error(state%x(1, :, :), ds, status)
  if (status /= 0 .or. residual < 0.0_dp .or. stretch_error < 0.0_dp) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: tension diagnostics', status, residual, stretch_error
    error stop 9
  end if

  kinetic_energy = fibre_prod_structure_kinetic_energy(state, rho_tilde)
  bending_energy = fibre_prod_structure_bending_energy_candidate(state, gamma, ds)
  if (kinetic_energy < 0.0_dp .or. bending_energy < 0.0_dp) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: energy diagnostics', kinetic_energy, bending_energy
    error stop 10
  end if

  call fibre_prod_structure_step_explicit(state, -dt, rho_tilde, status)
  if (status == 0) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: invalid dt accepted'
    error stop 11
  end if

  call fibre_prod_structure_compute_forces(state, gamma, -ds, status=status)
  if (status == 0) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: invalid ds accepted'
    error stop 12
  end if

  call fibre_prod_state_destroy(state)
  if (fibre_prod_state_is_allocated(state)) then
    print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK FAIL: destroy allocation state'
    error stop 13
  end if

  deallocate(bend)
  deallocate(external_force)
  print *, 'R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK PASS'
end program fibre_prod_structure_solver_check
