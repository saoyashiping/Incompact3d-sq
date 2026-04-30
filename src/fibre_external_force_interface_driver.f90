program fibre_external_force_interface_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_external_force, only : clear_fibre_external_force, set_fibre_external_force
  use fibre_geometry, only : fibre_center_of_mass
  use fibre_diagnostics, only : compute_bending_energy, compute_kinetic_energy
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_structure_integrator, only : advance_fibre_structure_freefree

  implicit none

  type(fibre_t) :: fibre

  integer :: i, step, io_unit, ierr, nsteps
  real(mytype) :: dt, s, pi
  real(mytype) :: tension_residual, relative_tension_residual
  real(mytype) :: fext(3, 33), fconst(3)
  real(mytype) :: x0(3, 33)
  real(mytype) :: center0(3), center_final(3), center_expected(3)
  real(mytype) :: v_center(3), v_expected(3), shift(3)
  real(mytype) :: w

  real(mytype) :: zero_force_preservation_error_max
  real(mytype) :: zero_force_final_bending_energy, zero_force_final_kinetic_energy
  real(mytype) :: zero_force_final_length_error
  integer :: zero_force_solver_failure_count

  real(mytype) :: uniform_force_center_accel_error_norm, uniform_force_center_velocity_error_norm
  real(mytype) :: uniform_force_length_error, uniform_force_shape_error_max
  real(mytype) :: uniform_force_external_power_final
  integer :: uniform_force_solver_failure_count

  real(mytype) :: sinusoidal_force_net_force_norm, sinusoidal_force_final_bending_energy
  real(mytype) :: sinusoidal_force_center_drift_norm, sinusoidal_force_final_length_error
  real(mytype) :: sinusoidal_force_max_displacement, sinusoidal_force_external_power_final
  integer :: sinusoidal_force_solver_failure_count

  real(mytype) :: clear_force_norm_after_clear

  pi = acos(-1.0_mytype)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  call clear_fibre_external_force(fibre)

  dt = 1.0e-4_mytype
  nsteps = 20
  zero_force_solver_failure_count = 0

  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) zero_force_solver_failure_count = zero_force_solver_failure_count + 1
  end do

  zero_force_preservation_error_max = 0._mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    zero_force_preservation_error_max = max(zero_force_preservation_error_max, &
                                            sqrt((fibre%x(1, i) - s)**2 + fibre%x(2, i)**2 + fibre%x(3, i)**2))
  end do
  call compute_bending_energy(fibre, zero_force_final_bending_energy)
  call compute_kinetic_energy(fibre, zero_force_final_kinetic_energy)
  call compute_total_length_relative_error(fibre, zero_force_final_length_error)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  call clear_fibre_external_force(fibre)
  x0 = fibre%x
  center0 = fibre_center_of_mass(fibre)

  fconst = [0._mytype, 0.1_mytype, 0._mytype]
  do i = 1, fibre%nl
    fext(:, i) = fconst
  end do
  call set_fibre_external_force(fibre, fext)

  dt = 1.0e-5_mytype
  nsteps = 20
  uniform_force_solver_failure_count = 0
  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) uniform_force_solver_failure_count = uniform_force_solver_failure_count + 1
  end do

  center_final = fibre_center_of_mass(fibre)
  center_expected = center0 + 0.5_mytype * (fconst / fibre%rho_tilde) * (real(nsteps, mytype) * dt)**2
  uniform_force_center_accel_error_norm = sqrt(sum((center_final - center_expected)**2))

  v_center = 0._mytype
  do i = 1, fibre%nl
    v_center = v_center + fibre%v(:, i)
  end do
  v_center = v_center / real(fibre%nl, mytype)
  v_expected = (fconst / fibre%rho_tilde) * (real(nsteps, mytype) * dt)
  uniform_force_center_velocity_error_norm = sqrt(sum((v_center - v_expected)**2))

  call compute_total_length_relative_error(fibre, uniform_force_length_error)
  shift = center_final - center0
  uniform_force_shape_error_max = 0._mytype
  do i = 1, fibre%nl
    uniform_force_shape_error_max = max(uniform_force_shape_error_max, sqrt(sum((fibre%x(:, i) - (x0(:, i) + shift))**2)))
  end do

  call compute_external_power(fibre, uniform_force_external_power_final)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  call clear_fibre_external_force(fibre)
  x0 = fibre%x
  center0 = fibre_center_of_mass(fibre)

  fext = 0._mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fext(2, i) = 0.01_mytype * sin(2._mytype * pi * s / fibre%length)
  end do
  call set_fibre_external_force(fibre, fext)

  sinusoidal_force_net_force_norm = 0._mytype
  do i = 1, fibre%nl
    if (i == 1 .or. i == fibre%nl) then
      w = 0.5_mytype * fibre%ds
    else
      w = fibre%ds
    end if
    sinusoidal_force_net_force_norm = sinusoidal_force_net_force_norm + w * fext(2, i)
  end do
  sinusoidal_force_net_force_norm = abs(sinusoidal_force_net_force_norm)

  dt = 1.0e-5_mytype
  nsteps = 20
  sinusoidal_force_solver_failure_count = 0
  do step = 1, nsteps
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) sinusoidal_force_solver_failure_count = sinusoidal_force_solver_failure_count + 1
  end do

  call compute_bending_energy(fibre, sinusoidal_force_final_bending_energy)
  center_final = fibre_center_of_mass(fibre)
  sinusoidal_force_center_drift_norm = sqrt(sum((center_final - center0)**2))
  call compute_total_length_relative_error(fibre, sinusoidal_force_final_length_error)

  sinusoidal_force_max_displacement = 0._mytype
  do i = 1, fibre%nl
    sinusoidal_force_max_displacement = max(sinusoidal_force_max_displacement, sqrt(sum((fibre%x(:, i) - x0(:, i))**2)))
  end do

  call compute_external_power(fibre, sinusoidal_force_external_power_final)

  call set_fibre_external_force(fibre, fext)
  call clear_fibre_external_force(fibre)
  clear_force_norm_after_clear = sqrt(sum(fibre%f_ext**2))

  open(newunit=io_unit, file='stage2_outputs/fibre_external_force_interface_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'zero_force_preservation_error_max = ', zero_force_preservation_error_max
  write(io_unit, '(A,ES24.16)') 'zero_force_final_bending_energy = ', zero_force_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'zero_force_final_kinetic_energy = ', zero_force_final_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'zero_force_final_length_error = ', zero_force_final_length_error
  write(io_unit, '(A,I0)') 'zero_force_solver_failure_count = ', zero_force_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'uniform_force_center_accel_error_norm = ', uniform_force_center_accel_error_norm
  write(io_unit, '(A,ES24.16)') 'uniform_force_center_velocity_error_norm = ', uniform_force_center_velocity_error_norm
  write(io_unit, '(A,ES24.16)') 'uniform_force_length_error = ', uniform_force_length_error
  write(io_unit, '(A,ES24.16)') 'uniform_force_shape_error_max = ', uniform_force_shape_error_max
  write(io_unit, '(A,I0)') 'uniform_force_solver_failure_count = ', uniform_force_solver_failure_count
  write(io_unit, '(A,ES24.16)') 'uniform_force_external_power_final = ', uniform_force_external_power_final

  write(io_unit, '(A,ES24.16)') 'sinusoidal_force_net_force_norm = ', sinusoidal_force_net_force_norm
  write(io_unit, '(A,ES24.16)') 'sinusoidal_force_final_bending_energy = ', sinusoidal_force_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'sinusoidal_force_center_drift_norm = ', sinusoidal_force_center_drift_norm
  write(io_unit, '(A,ES24.16)') 'sinusoidal_force_final_length_error = ', sinusoidal_force_final_length_error
  write(io_unit, '(A,ES24.16)') 'sinusoidal_force_max_displacement = ', sinusoidal_force_max_displacement
  write(io_unit, '(A,I0)') 'sinusoidal_force_solver_failure_count = ', sinusoidal_force_solver_failure_count
  write(io_unit, '(A,ES24.16)') 'sinusoidal_force_external_power_final = ', sinusoidal_force_external_power_final

  write(io_unit, '(A,ES24.16)') 'clear_force_norm_after_clear = ', clear_force_norm_after_clear

  close(io_unit)

contains

  subroutine compute_external_power(fibre, power)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: power

    integer :: j
    real(mytype) :: wj

    power = 0._mytype
    do j = 1, fibre%nl
      if (j == 1 .or. j == fibre%nl) then
        wj = 0.5_mytype * fibre%ds
      else
        wj = fibre%ds
      end if
      power = power + wj * dot_product(fibre%f_ext(:, j), fibre%v(:, j))
    end do
  end subroutine compute_external_power

end program fibre_external_force_interface_check
