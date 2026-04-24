program fibre_structure_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_diagnostics, only : fibre_diagnostics_t, fibre_compute_diagnostics
  use fibre_diagnostics, only : fibre_write_diagnostics

  implicit none

  type(fibre_t) :: fibre
  type(fibre_diagnostics_t) :: diag

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=0.0_mytype, gamma=0.0_mytype)

  call fibre_compute_diagnostics(fibre, diag)
  call fibre_write_diagnostics('fibre_geometry_check.dat', fibre, diag)

end program fibre_structure_check
