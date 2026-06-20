module fibre_prod_fsi_coupling
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type
  use fibre_prod_state, only : fibre_prod_state_type
  use fibre_prod_fsi_config, only : fibre_prod_fsi_config_type, fibre_prod_fsi_config_validate
  use fibre_prod_fluid_surrogate, only : fibre_prod_fluid_surrogate_type, &
                                         fibre_prod_fluid_surrogate_is_finite, &
                                         fibre_prod_fluid_surrogate_apply_force_density_step
  use fibre_prod_ibm_interpolation, only : fibre_prod_ibm_interpolate_velocity
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, &
                                          fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_total_force, &
                                          fibre_prod_force_buffer_destroy
  use fibre_prod_ibm_spreading, only : fibre_prod_spread_point_force
  use fibre_prod_structure_solver, only : fibre_prod_structure_compute_forces, &
                                          fibre_prod_structure_step_explicit
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_fsi_closed_loop_step

contains

  subroutine fibre_prod_fsi_closed_loop_step(grid, fluid, state, config, status, &
                                             total_fibre_force, reaction_integral)
    type(fibre_prod_grid_type), intent(in) :: grid
    type(fibre_prod_fluid_surrogate_type), intent(inout) :: fluid
    type(fibre_prod_state_type), intent(inout) :: state
    type(fibre_prod_fsi_config_type), intent(in) :: config
    integer, intent(out) :: status
    real(dp), intent(out) :: total_fibre_force(3)
    real(dp), intent(out) :: reaction_integral(3)
    type(fibre_prod_force_buffer_type) :: reaction_buffer
    real(dp), allocatable :: external_force(:, :, :)
    real(dp) :: fluid_velocity(3)
    real(dp) :: lag_force(3)
    real(dp) :: reaction_force(3)
    real(dp) :: weight_sum
    integer :: i_fibre
    integer :: i_node
    integer :: ierr

    status = fibre_prod_fsi_config_validate(config)
    total_fibre_force = 0.0_dp
    reaction_integral = 0.0_dp
    if (status /= 0) return
    if (.not. fibre_prod_fluid_surrogate_is_finite(fluid)) then
      status = 10
      return
    end if

    allocate(external_force(state%nfibre, state%nnode, 3))
    external_force = 0.0_dp
    call fibre_prod_force_buffer_allocate(reaction_buffer, grid, ierr)
    if (ierr /= 0) then
      status = 20 + ierr
    else
      do i_fibre = 1, state%nfibre
        do i_node = 1, state%nnode
          call fibre_prod_ibm_interpolate_velocity(grid, stack_velocity(fluid), state%x(i_fibre, i_node, :), &
                                                   fluid_velocity, ierr)
          if (ierr /= 0) then
            status = 30 + ierr
            exit
          end if
          lag_force = config%lambda_fsi * config%penalty_beta * (fluid_velocity - state%v(i_fibre, i_node, :))
          if (.not. config%fsi_enabled) lag_force = 0.0_dp
          external_force(i_fibre, i_node, :) = lag_force
          total_fibre_force = total_fibre_force + lag_force
          reaction_force = -lag_force
          call fibre_prod_spread_point_force(grid, reaction_buffer, state%x(i_fibre, i_node, :), &
                                             reaction_force, ierr, weight_sum)
          if (ierr /= 0) then
            status = 40 + ierr
            exit
          end if
        end do
        if (status /= 0) exit
      end do
      if (status == 0) then
        call fibre_prod_force_buffer_total_force(reaction_buffer, grid, reaction_integral, ierr)
        if (ierr /= 0) status = 50 + ierr
      end if
      if (status == 0) then
        call fibre_prod_structure_compute_forces(state, config%gamma, config%ds, external_force, ierr)
        if (ierr /= 0) status = 60 + ierr
      end if
      if (status == 0) then
        call fibre_prod_structure_step_explicit(state, config%dt, config%rho_tilde, ierr)
        if (ierr /= 0) status = 70 + ierr
      end if
      if (status == 0) then
        call fibre_prod_fluid_surrogate_apply_force_density_step(fluid, reaction_buffer, &
                                                                 config%dt, config%lambda_fsi, ierr)
        if (ierr /= 0) status = 80 + ierr
      end if
    end if
    call fibre_prod_force_buffer_destroy(reaction_buffer)
    if (allocated(external_force)) deallocate(external_force)
  end subroutine fibre_prod_fsi_closed_loop_step

  function stack_velocity(fluid) result(velocity_field)
    type(fibre_prod_fluid_surrogate_type), intent(in) :: fluid
    real(dp), allocatable :: velocity_field(:, :, :, :)

    allocate(velocity_field(fluid%nx_local, fluid%ny_local, fluid%nz_local, 3))
    velocity_field(:, :, :, 1) = fluid%u
    velocity_field(:, :, :, 2) = fluid%v
    velocity_field(:, :, :, 3) = fluid%w
  end function stack_velocity

end module fibre_prod_fsi_coupling
