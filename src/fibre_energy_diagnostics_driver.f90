program fibre_energy_diagnostics_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_diagnostics, only : compute_curvature_diagnostics
  use fibre_diagnostics, only : compute_bending_energy
  use fibre_diagnostics, only : compute_kinetic_energy
  use fibre_diagnostics, only : compute_total_structure_energy
  use fibre_diagnostics, only : compute_center_of_mass_velocity
  use fibre_diagnostics, only : compute_linear_momentum
  use fibre_diagnostics, only : compute_end_to_end_metrics
  use fibre_diagnostics, only : compute_kinetic_energy_split

  implicit none

  type(fibre_t) :: fibre

  integer :: i, io_unit
  real(mytype) :: pi, s, amp
  real(mytype) :: r_vec(3), r_norm, stretch_ratio
  real(mytype) :: vc(3), p(3), expected_p(3)
  real(mytype) :: bending_energy, kinetic_energy, total_energy
  real(mytype) :: max_curvature, rms_curvature
  real(mytype) :: total_length_relative_error
  real(mytype) :: ek_trans, ek_rel

  real(mytype) :: straight_rest_bending_energy
  real(mytype) :: straight_rest_kinetic_energy
  real(mytype) :: straight_rest_total_energy
  real(mytype) :: straight_rest_linear_momentum_norm
  real(mytype) :: straight_rest_max_curvature
  real(mytype) :: straight_rest_stretch_ratio
  real(mytype) :: straight_rest_total_length_relative_error

  real(mytype) :: translation_bending_energy
  real(mytype) :: translation_kinetic_energy
  real(mytype) :: translation_total_energy
  real(mytype) :: translation_center_velocity(3)
  real(mytype) :: translation_linear_momentum(3)
  real(mytype) :: translation_momentum_consistency_error
  real(mytype) :: translation_relative_kinetic_energy
  real(mytype) :: translation_translational_kinetic_energy
  real(mytype) :: translation_kinetic_split_error

  real(mytype) :: sine_bending_energy
  real(mytype) :: sine_kinetic_energy
  real(mytype) :: sine_total_energy
  real(mytype) :: sine_max_curvature
  real(mytype) :: sine_rms_curvature
  real(mytype) :: sine_stretch_ratio
  real(mytype) :: sine_total_length_relative_error
  real(mytype) :: sine_end_to_end_distance

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)

  call compute_bending_energy(fibre, bending_energy)
  call compute_kinetic_energy(fibre, kinetic_energy)
  call compute_total_structure_energy(fibre, total_energy)
  call compute_linear_momentum(fibre, p)
  call compute_curvature_diagnostics(fibre, max_curvature, rms_curvature)
  call compute_end_to_end_metrics(fibre, r_vec, r_norm, stretch_ratio)
  call compute_total_length_relative_error(fibre, total_length_relative_error)

  straight_rest_bending_energy = bending_energy
  straight_rest_kinetic_energy = kinetic_energy
  straight_rest_total_energy = total_energy
  straight_rest_linear_momentum_norm = sqrt(sum(p**2))
  straight_rest_max_curvature = max_curvature
  straight_rest_stretch_ratio = stretch_ratio
  straight_rest_total_length_relative_error = total_length_relative_error

  do i = 1, fibre%nl
    fibre%v(:, i) = [0.2_mytype, -0.1_mytype, 0.05_mytype]
  end do

  call compute_bending_energy(fibre, bending_energy)
  call compute_kinetic_energy(fibre, kinetic_energy)
  call compute_total_structure_energy(fibre, total_energy)
  call compute_center_of_mass_velocity(fibre, vc)
  call compute_linear_momentum(fibre, p)
  call compute_kinetic_energy_split(fibre, ek_trans, ek_rel)

  expected_p = fibre%rho_tilde * fibre%length * vc

  translation_bending_energy = bending_energy
  translation_kinetic_energy = kinetic_energy
  translation_total_energy = total_energy
  translation_center_velocity = vc
  translation_linear_momentum = p
  translation_momentum_consistency_error = sqrt(sum((p - expected_p)**2))
  translation_relative_kinetic_energy = ek_rel
  translation_translational_kinetic_energy = ek_trans
  translation_kinetic_split_error = abs(kinetic_energy - (ek_trans + ek_rel))

  fibre%v = 0._mytype
  pi = acos(-1.0_mytype)
  amp = 0.01_mytype

  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do

  call compute_bending_energy(fibre, bending_energy)
  call compute_kinetic_energy(fibre, kinetic_energy)
  call compute_total_structure_energy(fibre, total_energy)
  call compute_curvature_diagnostics(fibre, max_curvature, rms_curvature)
  call compute_end_to_end_metrics(fibre, r_vec, r_norm, stretch_ratio)
  call compute_total_length_relative_error(fibre, total_length_relative_error)

  sine_bending_energy = bending_energy
  sine_kinetic_energy = kinetic_energy
  sine_total_energy = total_energy
  sine_max_curvature = max_curvature
  sine_rms_curvature = rms_curvature
  sine_stretch_ratio = stretch_ratio
  sine_total_length_relative_error = total_length_relative_error
  sine_end_to_end_distance = r_norm

  open(newunit=io_unit, file='stage2_outputs/fibre_energy_diagnostics_check.dat', &
       status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'straight_rest_bending_energy = ', straight_rest_bending_energy
  write(io_unit, '(A,ES24.16)') 'straight_rest_kinetic_energy = ', straight_rest_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'straight_rest_total_energy = ', straight_rest_total_energy
  write(io_unit, '(A,ES24.16)') 'straight_rest_linear_momentum_norm = ', straight_rest_linear_momentum_norm
  write(io_unit, '(A,ES24.16)') 'straight_rest_max_curvature = ', straight_rest_max_curvature
  write(io_unit, '(A,ES24.16)') 'straight_rest_stretch_ratio = ', straight_rest_stretch_ratio
  write(io_unit, '(A,ES24.16)') 'straight_rest_total_length_relative_error = ', straight_rest_total_length_relative_error

  write(io_unit, '(A,ES24.16)') 'translation_bending_energy = ', translation_bending_energy
  write(io_unit, '(A,ES24.16)') 'translation_kinetic_energy = ', translation_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'translation_total_energy = ', translation_total_energy
  write(io_unit, '(A,ES24.16)') 'translation_center_velocity_x = ', translation_center_velocity(1)
  write(io_unit, '(A,ES24.16)') 'translation_center_velocity_y = ', translation_center_velocity(2)
  write(io_unit, '(A,ES24.16)') 'translation_center_velocity_z = ', translation_center_velocity(3)
  write(io_unit, '(A,ES24.16)') 'translation_linear_momentum_x = ', translation_linear_momentum(1)
  write(io_unit, '(A,ES24.16)') 'translation_linear_momentum_y = ', translation_linear_momentum(2)
  write(io_unit, '(A,ES24.16)') 'translation_linear_momentum_z = ', translation_linear_momentum(3)
  write(io_unit, '(A,ES24.16)') 'translation_momentum_consistency_error = ', translation_momentum_consistency_error
  write(io_unit, '(A,ES24.16)') 'translation_relative_kinetic_energy = ', translation_relative_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'translation_translational_kinetic_energy = ', translation_translational_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'translation_kinetic_split_error = ', translation_kinetic_split_error

  write(io_unit, '(A,ES24.16)') 'sine_bending_energy = ', sine_bending_energy
  write(io_unit, '(A,ES24.16)') 'sine_kinetic_energy = ', sine_kinetic_energy
  write(io_unit, '(A,ES24.16)') 'sine_total_energy = ', sine_total_energy
  write(io_unit, '(A,ES24.16)') 'sine_max_curvature = ', sine_max_curvature
  write(io_unit, '(A,ES24.16)') 'sine_rms_curvature = ', sine_rms_curvature
  write(io_unit, '(A,ES24.16)') 'sine_stretch_ratio = ', sine_stretch_ratio
  write(io_unit, '(A,ES24.16)') 'sine_total_length_relative_error = ', sine_total_length_relative_error
  write(io_unit, '(A,ES24.16)') 'sine_end_to_end_distance = ', sine_end_to_end_distance

  close(io_unit)

end program fibre_energy_diagnostics_check
