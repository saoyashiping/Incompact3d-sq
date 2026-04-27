program fibre_structure_advance_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_geometry, only : fibre_center_of_mass
  use fibre_diagnostics, only : compute_bending_energy, compute_kinetic_energy
  use fibre_diagnostics, only : compute_curvature_diagnostics, compute_total_length_relative_error
  use fibre_diagnostics, only : compute_linear_momentum
  use fibre_implicit_bending, only : build_freefree_d4_matrix
  use fibre_structure_integrator, only : build_tension_stiffness_matrix_freefree
  use fibre_structure_integrator, only : build_structure_advance_matrix_freefree
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_tension_solver, only : solve_tension_freefree
  use fibre_bending_operator, only : compute_fibre_d4_freefree

  implicit none

  type(fibre_t) :: fibre
  type(fibre_t) :: fibre_tmp

  integer :: i, step, ierr, nsteps, io_unit
  integer :: rest_solver_failure_count, translation_solver_failure_count
  integer :: stretch_solver_failure_count, curved_solver_failure_count
  integer :: curved_nan_detected, curved_endpoint_fixed_constraint_detected

  real(mytype) :: dt, s, pi, amp, eps
  real(mytype) :: tension_residual, relative_tension_residual
  real(mytype) :: rest_preservation_error_max, rest_final_bending_energy
  real(mytype) :: rest_final_kinetic_energy, rest_final_max_curvature
  real(mytype) :: rest_final_length_error, rest_tension_norm_max

  real(mytype) :: u(3), center_initial(3), center_expected(3), center_final(3)
  real(mytype) :: translation_center_error_norm, translation_velocity_error_norm
  real(mytype) :: translation_shape_error_max, translation_final_bending_energy
  real(mytype) :: translation_final_length_error, translation_momentum_error_norm
  real(mytype) :: momentum(3), momentum_expected(3), shift(3)

  real(mytype) :: stretch_initial_length_error, stretch_final_length_error
  real(mytype) :: stretch_tension_norm

  real(mytype) :: curved_initial_bending_energy, curved_final_bending_energy
  real(mytype) :: curved_initial_length_error, curved_final_length_error
  real(mytype) :: curved_final_max_curvature, curved_final_tension_norm
  real(mytype) :: left_before(3), right_before(3), left_after(3), right_after(3)

  real(mytype) :: advance_matrix_apply_maxdiff
  real(mytype), allocatable :: amat(:,:), dtmat(:,:), d4mat(:,:)
  real(mytype), allocatable :: q(:), y_mat(:), y_exp(:), dtq(:), d4q(:)
  real(mytype), allocatable :: tension(:), a_non_tension(:,:), d4x(:,:)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  dt = 1.0e-4_mytype
  nsteps = 20
  rest_solver_failure_count = 0
  rest_tension_norm_max = 0._mytype

  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) rest_solver_failure_count = rest_solver_failure_count + 1
    rest_tension_norm_max = max(rest_tension_norm_max, sqrt(sum(fibre%tension**2)))
  end do

  rest_preservation_error_max = 0._mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    rest_preservation_error_max = max(rest_preservation_error_max, &
                                       sqrt((fibre%x(1, i) - s)**2 + fibre%x(2, i)**2 + fibre%x(3, i)**2))
  end do

  call compute_bending_energy(fibre, rest_final_bending_energy)
  call compute_kinetic_energy(fibre, rest_final_kinetic_energy)
  call compute_curvature_diagnostics(fibre, rest_final_max_curvature, s)
  call compute_total_length_relative_error(fibre, rest_final_length_error)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  dt = 1.0e-4_mytype
  nsteps = 20
  u = [0.2_mytype, -0.1_mytype, 0.05_mytype]
  center_initial = fibre_center_of_mass(fibre)

  do i = 1, fibre%nl
    fibre%x_old(:, i) = fibre%x(:, i) - u * dt
    fibre%v(:, i) = u
  end do

  translation_solver_failure_count = 0
  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) translation_solver_failure_count = translation_solver_failure_count + 1
  end do

  center_final = fibre_center_of_mass(fibre)
  center_expected = center_initial + real(nsteps, mytype) * dt * u
  translation_center_error_norm = sqrt(sum((center_final - center_expected)**2))

  translation_velocity_error_norm = 0._mytype
  translation_shape_error_max = 0._mytype
  shift = real(nsteps, mytype) * dt * u
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    translation_velocity_error_norm = max(translation_velocity_error_norm, sqrt(sum((fibre%v(:, i) - u)**2)))
    translation_shape_error_max = max(translation_shape_error_max, sqrt((fibre%x(1, i) - (s + shift(1)))**2 + &
                                                                         (fibre%x(2, i) - shift(2))**2 + &
                                                                         (fibre%x(3, i) - shift(3))**2))
  end do

  call compute_bending_energy(fibre, translation_final_bending_energy)
  call compute_total_length_relative_error(fibre, translation_final_length_error)
  call compute_linear_momentum(fibre, momentum)
  momentum_expected = fibre%rho_tilde * fibre%length * u
  translation_momentum_error_norm = sqrt(sum((momentum - momentum_expected)**2))

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  eps = 1.0e-3_mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = (1._mytype + eps) * s
    fibre%x(2, i) = 0._mytype
    fibre%x(3, i) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  call compute_total_length_relative_error(fibre, stretch_initial_length_error)

  dt = 1.0e-5_mytype
  stretch_solver_failure_count = 0
  call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
  if (ierr /= 0) stretch_solver_failure_count = stretch_solver_failure_count + 1

  call compute_total_length_relative_error(fibre, stretch_final_length_error)
  stretch_tension_norm = sqrt(sum(fibre%tension**2))

  call fibre_init_straight_free_free(fibre, nl=65, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  pi = acos(-1.0_mytype)
  amp = 0.01_mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  call compute_bending_energy(fibre, curved_initial_bending_energy)
  call compute_total_length_relative_error(fibre, curved_initial_length_error)

  left_before = fibre%x(:, 1)
  right_before = fibre%x(:, fibre%nl)

  dt = 1.0e-6_mytype
  nsteps = 20
  curved_solver_failure_count = 0
  curved_nan_detected = 0

  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) curved_solver_failure_count = curved_solver_failure_count + 1
    if (any(ieee_is_nan(fibre%x)) .or. any(ieee_is_nan(fibre%v)) .or. any(ieee_is_nan(fibre%tension))) then
      curved_nan_detected = 1
    end if
  end do

  left_after = fibre%x(:, 1)
  right_after = fibre%x(:, fibre%nl)
  curved_endpoint_fixed_constraint_detected = 0
  if (sqrt(sum((left_after - left_before)**2)) <= 1.0e-15_mytype .and. &
      sqrt(sum((right_after - right_before)**2)) <= 1.0e-15_mytype) then
    curved_endpoint_fixed_constraint_detected = 1
  end if

  call compute_bending_energy(fibre, curved_final_bending_energy)
  call compute_total_length_relative_error(fibre, curved_final_length_error)
  call compute_curvature_diagnostics(fibre, curved_final_max_curvature, s)
  curved_final_tension_norm = sqrt(sum(fibre%tension**2))

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  allocate(d4x(3, fibre%nl), a_non_tension(3, fibre%nl), tension(fibre%nl - 1))
  fibre_tmp = fibre
  call compute_fibre_d4_freefree(fibre_tmp, d4x)
  a_non_tension = -(fibre%gamma / fibre%rho_tilde) * d4x
  call solve_tension_freefree(fibre, 1.0e-4_mytype, a_non_tension, tension, tension_residual, relative_tension_residual, ierr)
  if (ierr /= 0) error stop 'solve_tension_freefree failed in advance-matrix test'

  allocate(amat(fibre%nl, fibre%nl), dtmat(fibre%nl, fibre%nl), d4mat(fibre%nl, fibre%nl))
  allocate(q(fibre%nl), y_mat(fibre%nl), y_exp(fibre%nl), dtq(fibre%nl), d4q(fibre%nl))

  call build_structure_advance_matrix_freefree(fibre, 1.0e-4_mytype, tension, amat)
  call build_tension_stiffness_matrix_freefree(fibre, tension, dtmat)
  call build_freefree_d4_matrix(fibre, d4mat)

  do i = 1, fibre%nl
    q(i) = sin(pi * real(i, mytype) / real(fibre%nl + 1, mytype))
  end do

  y_mat = matmul(amat, q)
  dtq = matmul(dtmat, q)
  d4q = matmul(d4mat, q)
  y_exp = q - (1.0e-4_mytype**2) * dtq + (1.0e-4_mytype**2) * (fibre%gamma / fibre%rho_tilde) * d4q
  advance_matrix_apply_maxdiff = maxval(abs(y_mat - y_exp))

  open(newunit=io_unit, file='stage2_outputs/fibre_structure_advance_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'rest_preservation_error_max = ', rest_preservation_error_max
  write(io_unit, '(A,ES24.16)') 'rest_final_bending_energy = ', rest_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'rest_final_kinetic_energy = ', rest_final_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'rest_final_max_curvature = ', rest_final_max_curvature
  write(io_unit, '(A,ES24.16)') 'rest_final_length_error = ', rest_final_length_error
  write(io_unit, '(A,ES24.16)') 'rest_tension_norm_max = ', rest_tension_norm_max
  write(io_unit, '(A,I0)') 'rest_solver_failure_count = ', rest_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'translation_center_error_norm = ', translation_center_error_norm
  write(io_unit, '(A,ES24.16)') 'translation_velocity_error_norm = ', translation_velocity_error_norm
  write(io_unit, '(A,ES24.16)') 'translation_shape_error_max = ', translation_shape_error_max
  write(io_unit, '(A,ES24.16)') 'translation_final_bending_energy = ', translation_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'translation_final_length_error = ', translation_final_length_error
  write(io_unit, '(A,ES24.16)') 'translation_momentum_error_norm = ', translation_momentum_error_norm
  write(io_unit, '(A,I0)') 'translation_solver_failure_count = ', translation_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'stretch_initial_length_error = ', stretch_initial_length_error
  write(io_unit, '(A,ES24.16)') 'stretch_final_length_error = ', stretch_final_length_error
  write(io_unit, '(A,ES24.16)') 'stretch_tension_norm = ', stretch_tension_norm
  write(io_unit, '(A,I0)') 'stretch_solver_failure_count = ', stretch_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'curved_initial_bending_energy = ', curved_initial_bending_energy
  write(io_unit, '(A,ES24.16)') 'curved_final_bending_energy = ', curved_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'curved_initial_length_error = ', curved_initial_length_error
  write(io_unit, '(A,ES24.16)') 'curved_final_length_error = ', curved_final_length_error
  write(io_unit, '(A,ES24.16)') 'curved_final_max_curvature = ', curved_final_max_curvature
  write(io_unit, '(A,ES24.16)') 'curved_final_tension_norm = ', curved_final_tension_norm
  write(io_unit, '(A,I0)') 'curved_nan_detected = ', curved_nan_detected
  write(io_unit, '(A,I0)') 'curved_solver_failure_count = ', curved_solver_failure_count
  write(io_unit, '(A,I0)') 'curved_endpoint_fixed_constraint_detected = ', curved_endpoint_fixed_constraint_detected

  write(io_unit, '(A,ES24.16)') 'advance_matrix_apply_maxdiff = ', advance_matrix_apply_maxdiff

  close(io_unit)

end program fibre_structure_advance_check
