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
  real(mytype) :: fiber_mass
  real(mytype) :: fiber_inertia_perp
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
  real(mytype), allocatable, dimension(:,:) :: fiber_slip
  real(mytype), allocatable, dimension(:,:) :: fiber_coupling_force
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
    fiber_mass = 1._mytype
    fiber_inertia_perp = 1._mytype
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
    if (allocated(fiber_slip)) deallocate(fiber_slip)
    if (allocated(fiber_coupling_force)) deallocate(fiber_coupling_force)
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
