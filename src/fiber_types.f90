!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_types

  use decomp_2d_constants, only : mytype

  implicit none

  logical :: fiber_active
  integer :: fiber_nl
  real(mytype) :: fiber_length
  real(mytype), dimension(3) :: fiber_center
  real(mytype), dimension(3) :: fiber_direction
  real(mytype), allocatable, dimension(:,:) :: fiber_x
  real(mytype), allocatable, dimension(:,:) :: fiber_uinterp
  real(mytype), allocatable, dimension(:,:) :: fiber_uexact
  real(mytype), allocatable, dimension(:,:) :: fiber_uerror
  real(mytype), allocatable, dimension(:) :: fiber_sumw
  real(mytype), allocatable, dimension(:) :: fiber_quad_w
  real(mytype), allocatable, dimension(:,:) :: fiber_test_force
  real(mytype) :: fiber_interp_max_error
  logical :: interp_test_active
  integer :: interp_test_case
  logical :: interp_solver_test_active
  integer :: interp_solver_output_step
  logical :: spread_test_active
  integer :: spread_test_case
  logical :: rigid_coupling_test_active
  logical :: rigid_free_test_active
  logical :: rigid_kinematics_test_active
  logical :: rigid_two_way_test_active
  logical :: fiber_flexible_active
  logical :: fiber_flex_initialized
  logical :: fiber_flex_operator_test_active
  logical :: fiber_flex_bending_test_active
  logical :: fiber_flex_constraint_test_active
  logical :: rigid_kinematics_one_way
  logical :: rigid_kinematics_standalone
  integer :: rigid_motion_case
  integer :: rigid_free_case
  integer :: rigid_kinematics_mode
  real(mytype) :: rigid_kinematics_shear_rate
  real(mytype) :: rigid_kinematics_poiseuille_umax
  real(mytype) :: rigid_kinematics_channel_height
  real(mytype) :: rigid_kinematics_lambda
  real(mytype) :: ibm_beta
  integer :: coupling_ramp_steps
  integer :: rigid_output_interval
  integer :: free_output_interval
  integer :: rigid_kinematics_output_interval
  integer :: rigid_two_way_output_interval
  real(mytype) :: rigid_two_way_force_relaxation
  real(mytype) :: rigid_two_way_velocity_relaxation
  real(mytype) :: rigid_two_way_velocity_relaxation_x
  real(mytype) :: rigid_two_way_velocity_relaxation_y
  real(mytype) :: rigid_two_way_velocity_relaxation_z
  real(mytype) :: rigid_two_way_omega_relaxation
  real(mytype) :: rigid_two_way_force_seed_relaxation
  logical :: rigid_two_way_write_turb_series
  integer :: rigid_two_way_turb_series_interval
  real(mytype) :: rigid_two_way_min_wall_gap
  logical :: rigid_two_way_parallel_streamwise_correction
  real(mytype) :: rigid_two_way_parallel_streamwise_alpha
  real(mytype) :: rigid_two_way_parallel_cosine_threshold
  logical :: rigid_two_way_parallel_ucx_implicit
  integer :: rigid_two_way_parallel_ucx_newton_iters
  real(mytype) :: rigid_two_way_parallel_ucx_newton_relaxation
  real(mytype) :: rigid_two_way_parallel_ucx_max_increment
  real(mytype) :: rigid_two_way_parallel_ucx_fd_eps
  integer :: rigid_two_way_subiterations
  real(mytype) :: rigid_two_way_subiter_slip_tol
  logical :: rigid_two_way_subiter_verbose
  logical :: rigid_two_way_startup_fit
  logical :: rigid_two_way_initialized
  real(mytype) :: fiber_mass
  real(mytype) :: fiber_inertia_perp
  real(mytype) :: fiber_ds
  real(mytype) :: fiber_length_error_max
  real(mytype) :: fiber_inext_error_max
  real(mytype) :: fiber_bc_residual_max
  integer :: fiber_flex_operator_case
  integer :: fiber_flex_bending_case
  integer :: fiber_flex_bending_nsteps
  integer :: fiber_flex_bending_output_interval
  real(mytype) :: fiber_flex_bending_dt
  real(mytype) :: fiber_bending_gamma
  integer :: fiber_flex_bending_linear_solver
  real(mytype) :: fiber_flex_bending_cg_tol
  integer :: fiber_flex_bending_cg_maxit
  real(mytype) :: fiber_flex_bending_iter_tol
  integer :: fiber_flex_bending_iter_maxit
  real(mytype) :: fiber_flex_bending_iter_tol_effective
  integer :: fiber_flex_bending_iter_maxit_effective
  integer :: fiber_flex_constraint_case
  integer :: fiber_flex_constraint_nsteps
  integer :: fiber_flex_constraint_output_interval
  real(mytype) :: fiber_flex_constraint_dt
  real(mytype) :: fiber_flex_constraint_force_amp
  integer :: fiber_flex_constraint_outer_maxit
  real(mytype) :: fiber_flex_constraint_outer_tol_x
  real(mytype) :: fiber_flex_constraint_outer_tol_g
  real(mytype) :: fiber_flex_constraint_outer_tol_rx_rel
  logical :: fiber_flex_constraint_line_search_active
  real(mytype) :: fiber_flex_constraint_line_search_beta
  integer :: fiber_flex_constraint_line_search_max_backtracks
  logical :: fiber_flex_constraint_tension_warm_start_active
  real(mytype) :: fiber_structure_rho_tilde
  real(mytype), dimension(3) :: rigid_translation_velocity
  real(mytype), dimension(3) :: fiber_xc
  real(mytype), dimension(3) :: fiber_uc
  real(mytype), dimension(3) :: fiber_p
  real(mytype), dimension(3) :: fiber_omega
  real(mytype), dimension(3) :: fiber_force_total
  real(mytype), dimension(3) :: fiber_torque_total
  real(mytype), allocatable, dimension(:) :: fiber_s_ref
  real(mytype), allocatable, dimension(:,:) :: fiber_x_ref
  real(mytype), allocatable, dimension(:,:) :: fiber_xdot
  real(mytype), allocatable, dimension(:,:) :: fiber_x_old
  real(mytype), allocatable, dimension(:,:) :: fiber_x_nm1
  real(mytype), allocatable, dimension(:,:) :: fiber_x_stage
  real(mytype), allocatable, dimension(:,:) :: fiber_xs
  real(mytype), allocatable, dimension(:,:) :: fiber_xss
  real(mytype), allocatable, dimension(:,:) :: fiber_xssss
  real(mytype), allocatable, dimension(:,:) :: fiber_slip
  real(mytype), allocatable, dimension(:,:) :: fiber_coupling_force
  ! Legacy/transitional nodal tension containers (fiber_nl).
  ! Keep for compatibility during migration; not the final staggered unknown semantics.
  real(mytype), allocatable, dimension(:) :: fiber_tension
  real(mytype), allocatable, dimension(:) :: fiber_tension_old
  ! Future main storage for staggered tension unknowns on half-grid interfaces (fiber_nl-1).
  real(mytype), allocatable, dimension(:) :: fiber_tension_half
  real(mytype), allocatable, dimension(:) :: fiber_tension_half_old
  real(mytype), allocatable, dimension(:) :: fiber_tension_half_prevstep
  real(mytype), allocatable, dimension(:,:) :: fiber_bending_force
  real(mytype), allocatable, dimension(:,:) :: fiber_tension_force
  real(mytype), allocatable, dimension(:,:) :: fiber_hydro_force
  real(mytype), allocatable, dimension(:,:) :: fiber_struct_rhs
  real(mytype), allocatable, dimension(:) :: fiber_constraint_residual
  real(mytype), allocatable, dimension(:) :: fiber_kappa
  real(mytype), allocatable, dimension(:,:,:) :: fiber_euler_force_x
  real(mytype), allocatable, dimension(:,:,:) :: fiber_euler_force_y
  real(mytype), allocatable, dimension(:,:,:) :: fiber_euler_force_z

contains

  subroutine fiber_set_defaults()

    fiber_active = .false.
    fiber_nl = 2
    fiber_length = 1._mytype
    fiber_center = 0._mytype
    fiber_direction = (/1._mytype, 0._mytype, 0._mytype/)

    interp_test_active = .false.
    interp_test_case = 1
    fiber_interp_max_error = 0._mytype
    interp_solver_test_active = .false.
    interp_solver_output_step = 1
    spread_test_active = .false.
    spread_test_case = 1
    rigid_coupling_test_active = .false.
    rigid_free_test_active = .false.
    rigid_kinematics_test_active = .false.
    rigid_two_way_test_active = .false.
    fiber_flexible_active = .false.
    fiber_flex_initialized = .false.
    fiber_flex_operator_test_active = .false.
    fiber_flex_bending_test_active = .false.
    fiber_flex_constraint_test_active = .false.
    rigid_kinematics_one_way = .true.
    rigid_kinematics_standalone = .false.
    rigid_motion_case = 1
    rigid_free_case = 1
    rigid_kinematics_mode = 1
    rigid_kinematics_shear_rate = 1._mytype
    rigid_kinematics_poiseuille_umax = 1._mytype
    rigid_kinematics_channel_height = 2._mytype
    rigid_kinematics_lambda = 1._mytype
    ibm_beta = -1.0e3_mytype
    coupling_ramp_steps = 0
    rigid_output_interval = 1
    free_output_interval = 1
    rigid_kinematics_output_interval = 1
    rigid_two_way_output_interval = 1
    rigid_two_way_force_relaxation = 0.25_mytype
    rigid_two_way_velocity_relaxation = 0.10_mytype
    rigid_two_way_velocity_relaxation_x = 0.10_mytype
    rigid_two_way_velocity_relaxation_y = 0.10_mytype
    rigid_two_way_velocity_relaxation_z = 0.10_mytype
    rigid_two_way_omega_relaxation = 0.10_mytype
    rigid_two_way_force_seed_relaxation = 0.5_mytype
    rigid_two_way_write_turb_series = .false.
    rigid_two_way_turb_series_interval = 1
    rigid_two_way_min_wall_gap = 0._mytype
    rigid_two_way_parallel_streamwise_correction = .false.
    rigid_two_way_parallel_streamwise_alpha = 0.05_mytype
    rigid_two_way_parallel_cosine_threshold = 0.98_mytype
    rigid_two_way_parallel_ucx_implicit = .false.
    rigid_two_way_parallel_ucx_newton_iters = 2
    rigid_two_way_parallel_ucx_newton_relaxation = 0.5_mytype
    rigid_two_way_parallel_ucx_max_increment = 0.05_mytype
    rigid_two_way_parallel_ucx_fd_eps = 1.0e-4_mytype
    rigid_two_way_subiterations = 2
    rigid_two_way_subiter_slip_tol = 1.0e-3_mytype
    rigid_two_way_subiter_verbose = .false.
    rigid_two_way_startup_fit = .true.
    rigid_two_way_initialized = .false.
    fiber_mass = 1._mytype
    fiber_inertia_perp = 1._mytype
    fiber_ds = 0._mytype
    fiber_length_error_max = 0._mytype
    fiber_inext_error_max = 0._mytype
    fiber_bc_residual_max = 0._mytype
    fiber_flex_operator_case = 1
    fiber_flex_bending_case = 1
    fiber_flex_bending_nsteps = 20
    fiber_flex_bending_output_interval = 1
    fiber_flex_bending_dt = 1.0e-3_mytype
    fiber_bending_gamma = 1._mytype
    fiber_flex_bending_linear_solver = 1
    fiber_flex_bending_cg_tol = 1.0e-12_mytype
    fiber_flex_bending_cg_maxit = 500
    fiber_flex_bending_iter_tol = -1._mytype
    fiber_flex_bending_iter_maxit = 0
    fiber_flex_bending_iter_tol_effective = fiber_flex_bending_cg_tol
    fiber_flex_bending_iter_maxit_effective = fiber_flex_bending_cg_maxit
    fiber_flex_constraint_case = 1
    fiber_flex_constraint_nsteps = 50
    fiber_flex_constraint_output_interval = 1
    fiber_flex_constraint_dt = 1.0e-3_mytype
    fiber_flex_constraint_force_amp = 1.0e-3_mytype
    fiber_flex_constraint_outer_maxit = 20
    fiber_flex_constraint_outer_tol_x = 1.0e-10_mytype
    fiber_flex_constraint_outer_tol_g = 1.0e-10_mytype
    fiber_flex_constraint_outer_tol_rx_rel = 1.0e-6_mytype
    fiber_flex_constraint_line_search_active = .true.
    fiber_flex_constraint_line_search_beta = 0.5_mytype
    fiber_flex_constraint_line_search_max_backtracks = 6
    fiber_flex_constraint_tension_warm_start_active = .true.
    fiber_structure_rho_tilde = 1.0_mytype
    rigid_translation_velocity = 0._mytype
    fiber_xc = 0._mytype
    fiber_uc = 0._mytype
    fiber_p = (/1._mytype, 0._mytype, 0._mytype/)
    fiber_omega = 0._mytype
    fiber_force_total = 0._mytype
    fiber_torque_total = 0._mytype

    if (allocated(fiber_x)) deallocate(fiber_x)
    if (allocated(fiber_s_ref)) deallocate(fiber_s_ref)
    if (allocated(fiber_x_ref)) deallocate(fiber_x_ref)
    if (allocated(fiber_xdot)) deallocate(fiber_xdot)
    if (allocated(fiber_x_old)) deallocate(fiber_x_old)
    if (allocated(fiber_x_nm1)) deallocate(fiber_x_nm1)
    if (allocated(fiber_x_stage)) deallocate(fiber_x_stage)
    if (allocated(fiber_xs)) deallocate(fiber_xs)
    if (allocated(fiber_xss)) deallocate(fiber_xss)
    if (allocated(fiber_xssss)) deallocate(fiber_xssss)
    if (allocated(fiber_slip)) deallocate(fiber_slip)
    if (allocated(fiber_coupling_force)) deallocate(fiber_coupling_force)
    if (allocated(fiber_tension)) deallocate(fiber_tension)
    if (allocated(fiber_tension_old)) deallocate(fiber_tension_old)
    if (allocated(fiber_tension_half)) deallocate(fiber_tension_half)
    if (allocated(fiber_tension_half_old)) deallocate(fiber_tension_half_old)
    if (allocated(fiber_tension_half_prevstep)) deallocate(fiber_tension_half_prevstep)
    if (allocated(fiber_bending_force)) deallocate(fiber_bending_force)
    if (allocated(fiber_tension_force)) deallocate(fiber_tension_force)
    if (allocated(fiber_hydro_force)) deallocate(fiber_hydro_force)
    if (allocated(fiber_struct_rhs)) deallocate(fiber_struct_rhs)
    if (allocated(fiber_constraint_residual)) deallocate(fiber_constraint_residual)
    if (allocated(fiber_kappa)) deallocate(fiber_kappa)
    if (allocated(fiber_uinterp)) deallocate(fiber_uinterp)
    if (allocated(fiber_uexact)) deallocate(fiber_uexact)
    if (allocated(fiber_uerror)) deallocate(fiber_uerror)
    if (allocated(fiber_sumw)) deallocate(fiber_sumw)
    if (allocated(fiber_quad_w)) deallocate(fiber_quad_w)
    if (allocated(fiber_test_force)) deallocate(fiber_test_force)
    if (allocated(fiber_euler_force_x)) deallocate(fiber_euler_force_x)
    if (allocated(fiber_euler_force_y)) deallocate(fiber_euler_force_y)
    if (allocated(fiber_euler_force_z)) deallocate(fiber_euler_force_z)

  end subroutine fiber_set_defaults

end module fiber_types
