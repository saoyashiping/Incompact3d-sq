module fibre_prod_fsi_diagnostics
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_all_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type
  use fibre_prod_fluid_surrogate, only : fibre_prod_fluid_surrogate_type, &
                                         fibre_prod_fluid_surrogate_kinetic_energy, &
                                         fibre_prod_fluid_surrogate_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_fsi_force_norm
  public :: fibre_prod_fsi_action_reaction_residual
  public :: fibre_prod_fsi_structure_kinetic_energy
  public :: fibre_prod_fsi_fluid_kinetic_energy
  public :: fibre_prod_fsi_all_finite
  public :: fibre_prod_fsi_lambda_zero_signature
  public :: fibre_prod_fsi_small_lambda_signature

contains

  pure real(dp) function fibre_prod_fsi_force_norm(force) result(norm_value)
    real(dp), intent(in) :: force(3)

    if (all(ieee_is_finite(force))) then
      norm_value = sqrt(sum(force * force))
    else
      norm_value = huge(1.0_dp)
    end if
  end function fibre_prod_fsi_force_norm

  pure real(dp) function fibre_prod_fsi_action_reaction_residual(fibre_force, reaction_integral) result(residual)
    real(dp), intent(in) :: fibre_force(3)
    real(dp), intent(in) :: reaction_integral(3)

    if (all(ieee_is_finite(fibre_force)) .and. all(ieee_is_finite(reaction_integral))) then
      residual = sqrt(sum((fibre_force + reaction_integral) * (fibre_force + reaction_integral)))
    else
      residual = huge(1.0_dp)
    end if
  end function fibre_prod_fsi_action_reaction_residual

  pure real(dp) function fibre_prod_fsi_structure_kinetic_energy(state, rho_tilde) result(energy)
    type(fibre_prod_state_type), intent(in) :: state
    real(dp), intent(in) :: rho_tilde

    if (rho_tilde > 0.0_dp .and. fibre_prod_state_all_finite(state)) then
      energy = 0.5_dp * rho_tilde * sum(state%v * state%v)
    else
      energy = huge(1.0_dp)
    end if
  end function fibre_prod_fsi_structure_kinetic_energy

  pure real(dp) function fibre_prod_fsi_fluid_kinetic_energy(fluid, grid) result(energy)
    type(fibre_prod_fluid_surrogate_type), intent(in) :: fluid
    type(fibre_prod_grid_type), intent(in) :: grid

    energy = fibre_prod_fluid_surrogate_kinetic_energy(fluid, grid)
  end function fibre_prod_fsi_fluid_kinetic_energy

  pure logical function fibre_prod_fsi_all_finite(state, fluid) result(all_finite)
    type(fibre_prod_state_type), intent(in) :: state
    type(fibre_prod_fluid_surrogate_type), intent(in) :: fluid

    all_finite = fibre_prod_state_all_finite(state) .and. fibre_prod_fluid_surrogate_is_finite(fluid)
  end function fibre_prod_fsi_all_finite

  pure logical function fibre_prod_fsi_lambda_zero_signature(force_norm, velocity_delta) result(is_clean)
    real(dp), intent(in) :: force_norm
    real(dp), intent(in) :: velocity_delta

    is_clean = ieee_is_finite(force_norm) .and. ieee_is_finite(velocity_delta) .and. &
               force_norm < 1.0e-14_dp .and. velocity_delta < 1.0e-12_dp
  end function fibre_prod_fsi_lambda_zero_signature

  pure logical function fibre_prod_fsi_small_lambda_signature(force_norm, response_norm) result(is_bounded)
    real(dp), intent(in) :: force_norm
    real(dp), intent(in) :: response_norm

    is_bounded = ieee_is_finite(force_norm) .and. ieee_is_finite(response_norm) .and. &
                 force_norm > 0.0_dp .and. response_norm >= 0.0_dp .and. response_norm < 1.0_dp
  end function fibre_prod_fsi_small_lambda_signature

end module fibre_prod_fsi_diagnostics
