!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_init

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, &
       fiber_ds, fiber_x, fiber_s_ref, fiber_xdot, fiber_x_old, fiber_x_nm1, fiber_x_stage, fiber_tension, &
       fiber_tension_old, fiber_bending_force, fiber_tension_force, fiber_hydro_force, fiber_struct_rhs, &
       fiber_constraint_residual, fiber_kappa, fiber_length_error_max, fiber_inext_error_max, fiber_bc_residual_max

  implicit none

contains

  subroutine init_flexible_fiber_state()

    ! Step 3.1 only: initialize flexible-fiber data containers.
    ! No flexible dynamics, constraints, tension solve, or IBM coupling is implemented here.

    if (.not.fiber_active) then
      if (nrank == 0) write(*,*) 'Error: fiber_flexible_active requires fiber_active = true.'
      stop
    endif

    if (.not.fiber_flexible_active) then
      if (nrank == 0) write(*,*) 'Error: init_flexible_fiber_state called while fiber_flexible_active = false.'
      stop
    endif

    if (fiber_nl < 5) then
      if (nrank == 0) write(*,*) 'Error: flexible fiber initialization requires fiber_nl >= 5.'
      stop
    endif

    if (.not.allocated(fiber_x)) then
      if (nrank == 0) write(*,*) 'Error: fiber_x must be allocated before flexible initialization.'
      stop
    endif

    if (.not.allocated(fiber_s_ref)) then
      if (nrank == 0) write(*,*) 'Error: fiber_s_ref must be allocated before flexible initialization.'
      stop
    endif

    fiber_ds = fiber_length / real(fiber_nl - 1, mytype)

    if (allocated(fiber_x_old)) deallocate(fiber_x_old)
    if (allocated(fiber_x_nm1)) deallocate(fiber_x_nm1)
    if (allocated(fiber_x_stage)) deallocate(fiber_x_stage)
    if (allocated(fiber_xdot)) deallocate(fiber_xdot)
    if (allocated(fiber_tension)) deallocate(fiber_tension)
    if (allocated(fiber_tension_old)) deallocate(fiber_tension_old)
    if (allocated(fiber_bending_force)) deallocate(fiber_bending_force)
    if (allocated(fiber_tension_force)) deallocate(fiber_tension_force)
    if (allocated(fiber_hydro_force)) deallocate(fiber_hydro_force)
    if (allocated(fiber_struct_rhs)) deallocate(fiber_struct_rhs)
    if (allocated(fiber_constraint_residual)) deallocate(fiber_constraint_residual)
    if (allocated(fiber_kappa)) deallocate(fiber_kappa)

    allocate(fiber_x_old(3, fiber_nl))
    allocate(fiber_x_nm1(3, fiber_nl))
    allocate(fiber_x_stage(3, fiber_nl))
    allocate(fiber_xdot(3, fiber_nl))
    allocate(fiber_tension(fiber_nl))
    allocate(fiber_tension_old(fiber_nl))
    allocate(fiber_bending_force(3, fiber_nl))
    allocate(fiber_tension_force(3, fiber_nl))
    allocate(fiber_hydro_force(3, fiber_nl))
    allocate(fiber_struct_rhs(3, fiber_nl))
    allocate(fiber_constraint_residual(fiber_nl))
    allocate(fiber_kappa(fiber_nl))

    fiber_x_old = fiber_x
    fiber_x_nm1 = fiber_x
    fiber_x_stage = fiber_x
    fiber_xdot = 0._mytype
    fiber_tension = 0._mytype
    fiber_tension_old = 0._mytype
    fiber_bending_force = 0._mytype
    fiber_tension_force = 0._mytype
    fiber_hydro_force = 0._mytype
    fiber_struct_rhs = 0._mytype
    fiber_constraint_residual = 0._mytype
    fiber_kappa = 0._mytype

    fiber_length_error_max = 0._mytype
    fiber_inext_error_max = 0._mytype
    fiber_bc_residual_max = 0._mytype
    fiber_flex_initialized = .true.

  end subroutine init_flexible_fiber_state

end module fiber_flex_init
