program fibre_no_force_long_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_geometry, only : fibre_center_of_mass
  use fibre_diagnostics, only : compute_bending_energy, compute_kinetic_energy
  use fibre_diagnostics, only : compute_curvature_diagnostics, compute_total_length_relative_error
  use fibre_diagnostics, only : compute_linear_momentum
  use fibre_structure_integrator, only : advance_fibre_structure_freefree

  implicit none

  type(fibre_t) :: fibre

  integer :: i, j, step, io_unit, ierr, nsteps
  integer :: rest_long_solver_failure_count, translation_long_solver_failure_count
  integer :: curved_long_solver_failure_count, curved_long_nan_detected
  integer :: curved_long_endpoint_fixed_constraint_detected

  real(mytype) :: dt, s, pi, theta0, theta, s_mid, ds_seg
  real(mytype) :: tension_residual, relative_tension_residual

  real(mytype) :: rest_long_preservation_error_max
  real(mytype) :: rest_long_final_bending_energy, rest_long_final_kinetic_energy
  real(mytype) :: rest_long_final_max_curvature, rest_long_final_length_error
  real(mytype) :: rest_long_momentum_norm, rest_long_center_drift_norm

  real(mytype) :: translation_long_center_error_norm, translation_long_velocity_error_norm
  real(mytype) :: translation_long_shape_error_max, translation_long_bending_energy
  real(mytype) :: translation_long_length_error, translation_long_momentum_error_norm

  real(mytype) :: curved_long_initial_bending_energy, curved_long_final_bending_energy
  real(mytype) :: curved_long_max_total_energy
  real(mytype) :: curved_long_initial_length_error, curved_long_final_length_error
  real(mytype) :: curved_long_max_length_error
  real(mytype) :: curved_long_final_momentum_norm, curved_long_center_drift_norm
  real(mytype) :: curved_long_final_momentum_relative
  real(mytype) :: curved_long_final_max_curvature

  real(mytype) :: center0(3), center_final(3), center_expected(3)
  real(mytype) :: u(3), momentum(3), momentum_expected(3), shift(3)
  real(mytype) :: left_before(3), right_before(3), left_after(3), right_after(3)
  real(mytype) :: max_curv_dummy, eb, ek, total_energy, length_error, momentum_scale

  pi = acos(-1.0_mytype)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  dt = 1.0e-4_mytype
  nsteps = 1000
  center0 = fibre_center_of_mass(fibre)
  rest_long_solver_failure_count = 0

  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) rest_long_solver_failure_count = rest_long_solver_failure_count + 1
  end do

  rest_long_preservation_error_max = 0._mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    rest_long_preservation_error_max = max(rest_long_preservation_error_max, &
                                           sqrt((fibre%x(1, i) - s)**2 + fibre%x(2, i)**2 + fibre%x(3, i)**2))
  end do

  call compute_bending_energy(fibre, rest_long_final_bending_energy)
  call compute_kinetic_energy(fibre, rest_long_final_kinetic_energy)
  call compute_curvature_diagnostics(fibre, rest_long_final_max_curvature, max_curv_dummy)
  call compute_total_length_relative_error(fibre, rest_long_final_length_error)
  call compute_linear_momentum(fibre, momentum)
  rest_long_momentum_norm = sqrt(sum(momentum**2))
  center_final = fibre_center_of_mass(fibre)
  rest_long_center_drift_norm = sqrt(sum((center_final - center0)**2))

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  dt = 1.0e-4_mytype
  nsteps = 1000
  u = [0.2_mytype, -0.1_mytype, 0.05_mytype]
  center0 = fibre_center_of_mass(fibre)

  do i = 1, fibre%nl
    fibre%x_old(:, i) = fibre%x(:, i) - u * dt
    fibre%v(:, i) = u
  end do

  translation_long_solver_failure_count = 0
  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) translation_long_solver_failure_count = translation_long_solver_failure_count + 1
  end do

  center_final = fibre_center_of_mass(fibre)
  center_expected = center0 + real(nsteps, mytype) * dt * u
  translation_long_center_error_norm = sqrt(sum((center_final - center_expected)**2))

  translation_long_velocity_error_norm = 0._mytype
  translation_long_shape_error_max = 0._mytype
  shift = real(nsteps, mytype) * dt * u
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    translation_long_velocity_error_norm = max(translation_long_velocity_error_norm, sqrt(sum((fibre%v(:, i) - u)**2)))
    translation_long_shape_error_max = max(translation_long_shape_error_max, sqrt((fibre%x(1, i) - (s + shift(1)))**2 + &
                                                                                   (fibre%x(2, i) - shift(2))**2 + &
                                                                                   (fibre%x(3, i) - shift(3))**2))
  end do

  call compute_bending_energy(fibre, translation_long_bending_energy)
  call compute_total_length_relative_error(fibre, translation_long_length_error)
  call compute_linear_momentum(fibre, momentum)
  momentum_expected = fibre%rho_tilde * fibre%length * u
  translation_long_momentum_error_norm = sqrt(sum((momentum - momentum_expected)**2))

  call fibre_init_straight_free_free(fibre, nl=65, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  theta0 = 0.05_mytype
  ds_seg = fibre%ds
  fibre%x(:, 1) = 0._mytype
  do j = 1, fibre%nl - 1
    s_mid = (real(j, mytype) - 0.5_mytype) * ds_seg
    theta = theta0 * (sin(pi * s_mid / fibre%length) + 0.2_mytype * sin(2._mytype * pi * s_mid / fibre%length))
    fibre%x(1, j + 1) = fibre%x(1, j) + ds_seg * cos(theta)
    fibre%x(2, j + 1) = fibre%x(2, j) + ds_seg * sin(theta)
    fibre%x(3, j + 1) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  center0 = fibre_center_of_mass(fibre)
  left_before = fibre%x(:, 1)
  right_before = fibre%x(:, fibre%nl)

  call compute_bending_energy(fibre, curved_long_initial_bending_energy)
  call compute_total_length_relative_error(fibre, curved_long_initial_length_error)
  call compute_kinetic_energy(fibre, ek)
  curved_long_max_total_energy = curved_long_initial_bending_energy + ek
  curved_long_max_length_error = curved_long_initial_length_error

  dt = 1.0e-6_mytype
  nsteps = 1000
  curved_long_solver_failure_count = 0
  curved_long_nan_detected = 0

  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) curved_long_solver_failure_count = curved_long_solver_failure_count + 1

    if (any(ieee_is_nan(fibre%x)) .or. any(ieee_is_nan(fibre%v)) .or. any(ieee_is_nan(fibre%tension))) then
      curved_long_nan_detected = 1
    end if

    call compute_bending_energy(fibre, eb)
    call compute_kinetic_energy(fibre, ek)
    total_energy = eb + ek
    curved_long_max_total_energy = max(curved_long_max_total_energy, total_energy)

    call compute_total_length_relative_error(fibre, length_error)
    curved_long_max_length_error = max(curved_long_max_length_error, length_error)
  end do

  call compute_bending_energy(fibre, curved_long_final_bending_energy)
  call compute_total_length_relative_error(fibre, curved_long_final_length_error)
  call compute_linear_momentum(fibre, momentum)
  curved_long_final_momentum_norm = sqrt(sum(momentum**2))
  momentum_scale = max(sqrt(2.0_mytype * fibre%rho_tilde * fibre%length * curved_long_max_total_energy), 1.0e-300_mytype)
  curved_long_final_momentum_relative = curved_long_final_momentum_norm / momentum_scale
  center_final = fibre_center_of_mass(fibre)
  curved_long_center_drift_norm = sqrt(sum((center_final - center0)**2))
  call compute_curvature_diagnostics(fibre, curved_long_final_max_curvature, max_curv_dummy)

  left_after = fibre%x(:, 1)
  right_after = fibre%x(:, fibre%nl)
  curved_long_endpoint_fixed_constraint_detected = 0
  if (sqrt(sum((left_after - left_before)**2)) <= 1.0e-15_mytype .and. &
      sqrt(sum((right_after - right_before)**2)) <= 1.0e-15_mytype) then
    curved_long_endpoint_fixed_constraint_detected = 1
  end if

  open(newunit=io_unit, file='stage2_outputs/fibre_no_force_long_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'rest_long_preservation_error_max = ', rest_long_preservation_error_max
  write(io_unit, '(A,ES24.16)') 'rest_long_final_bending_energy = ', rest_long_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'rest_long_final_kinetic_energy = ', rest_long_final_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'rest_long_final_max_curvature = ', rest_long_final_max_curvature
  write(io_unit, '(A,ES24.16)') 'rest_long_final_length_error = ', rest_long_final_length_error
  write(io_unit, '(A,ES24.16)') 'rest_long_momentum_norm = ', rest_long_momentum_norm
  write(io_unit, '(A,ES24.16)') 'rest_long_center_drift_norm = ', rest_long_center_drift_norm
  write(io_unit, '(A,I0)') 'rest_long_solver_failure_count = ', rest_long_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'translation_long_center_error_norm = ', translation_long_center_error_norm
  write(io_unit, '(A,ES24.16)') 'translation_long_velocity_error_norm = ', translation_long_velocity_error_norm
  write(io_unit, '(A,ES24.16)') 'translation_long_shape_error_max = ', translation_long_shape_error_max
  write(io_unit, '(A,ES24.16)') 'translation_long_bending_energy = ', translation_long_bending_energy
  write(io_unit, '(A,ES24.16)') 'translation_long_length_error = ', translation_long_length_error
  write(io_unit, '(A,ES24.16)') 'translation_long_momentum_error_norm = ', translation_long_momentum_error_norm
  write(io_unit, '(A,I0)') 'translation_long_solver_failure_count = ', translation_long_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'curved_long_initial_bending_energy = ', curved_long_initial_bending_energy
  write(io_unit, '(A,ES24.16)') 'curved_long_final_bending_energy = ', curved_long_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'curved_long_max_total_energy = ', curved_long_max_total_energy
  write(io_unit, '(A,ES24.16)') 'curved_long_initial_length_error = ', curved_long_initial_length_error
  write(io_unit, '(A,ES24.16)') 'curved_long_final_length_error = ', curved_long_final_length_error
  write(io_unit, '(A,ES24.16)') 'curved_long_max_length_error = ', curved_long_max_length_error
  write(io_unit, '(A,ES24.16)') 'curved_long_final_momentum_norm = ', curved_long_final_momentum_norm
  write(io_unit, '(A,ES24.16)') 'curved_long_final_momentum_relative = ', curved_long_final_momentum_relative
  write(io_unit, '(A,ES24.16)') 'curved_long_center_drift_norm = ', curved_long_center_drift_norm
  write(io_unit, '(A,ES24.16)') 'curved_long_final_max_curvature = ', curved_long_final_max_curvature
  write(io_unit, '(A,I0)') 'curved_long_nan_detected = ', curved_long_nan_detected
  write(io_unit, '(A,I0)') 'curved_long_solver_failure_count = ', curved_long_solver_failure_count
  write(io_unit, '(A,I0)') 'curved_long_endpoint_fixed_constraint_detected = ', curved_long_endpoint_fixed_constraint_detected

  close(io_unit)

end program fibre_no_force_long_check
