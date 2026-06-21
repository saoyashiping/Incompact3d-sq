module fibre_prod_collision_diagnostics
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_fibre_collision, only : fibre_prod_collision_state_type
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_collision_min_gap
  public :: fibre_prod_collision_max_penetration
  public :: fibre_prod_collision_force_norm
  public :: fibre_prod_collision_action_reaction_residual
  public :: fibre_prod_collision_energy_candidate
  public :: fibre_prod_collision_diagnostics_finite
  public :: fibre_prod_collision_summary

contains

  pure real(dp) function fibre_prod_collision_min_gap(state) result(value)
    type(fibre_prod_collision_state_type), intent(in) :: state
    value = state%min_gap
  end function fibre_prod_collision_min_gap

  pure real(dp) function fibre_prod_collision_max_penetration(state) result(value)
    type(fibre_prod_collision_state_type), intent(in) :: state
    value = state%max_penetration
  end function fibre_prod_collision_max_penetration

  pure real(dp) function fibre_prod_collision_force_norm(force_i, force_j) result(value)
    real(dp), intent(in) :: force_i(:, :)
    real(dp), intent(in) :: force_j(:, :)

    if (all(ieee_is_finite(force_i)) .and. all(ieee_is_finite(force_j))) then
      value = sqrt(sum(force_i * force_i) + sum(force_j * force_j))
    else
      value = huge(1.0_dp)
    end if
  end function fibre_prod_collision_force_norm

  pure real(dp) function fibre_prod_collision_action_reaction_residual(force_i, force_j) result(value)
    real(dp), intent(in) :: force_i(:, :)
    real(dp), intent(in) :: force_j(:, :)
    real(dp) :: total(3)

    if (all(ieee_is_finite(force_i)) .and. all(ieee_is_finite(force_j))) then
      total = sum(force_i, dim=1) + sum(force_j, dim=1)
      value = sqrt(sum(total * total))
    else
      value = huge(1.0_dp)
    end if
  end function fibre_prod_collision_action_reaction_residual

  pure real(dp) function fibre_prod_collision_energy_candidate(state, k_collision) result(energy)
    type(fibre_prod_collision_state_type), intent(in) :: state
    real(dp), intent(in) :: k_collision

    if (ieee_is_finite(state%max_penetration) .and. ieee_is_finite(k_collision) .and. &
        k_collision >= 0.0_dp) then
      energy = 0.5_dp * k_collision * state%max_penetration * state%max_penetration
    else
      energy = huge(1.0_dp)
    end if
  end function fibre_prod_collision_energy_candidate

  pure logical function fibre_prod_collision_diagnostics_finite(state) result(is_finite)
    type(fibre_prod_collision_state_type), intent(in) :: state

    is_finite = state%initialized .and. ieee_is_finite(state%min_gap) .and. &
                ieee_is_finite(state%max_penetration) .and. ieee_is_finite(state%force_norm) .and. &
                ieee_is_finite(state%closest_pair%distance) .and. ieee_is_finite(state%closest_pair%gap) .and. &
                all(ieee_is_finite(state%closest_pair%normal))
  end function fibre_prod_collision_diagnostics_finite

  function fibre_prod_collision_summary(state) result(summary)
    type(fibre_prod_collision_state_type), intent(in) :: state
    character(len=200) :: summary

    write(summary, '(A,L1,A,L1,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)') &
      'near=', state%near_contact_warning, ', overlap=', state%overlap_detected, &
      ', checked=', state%npairs_checked, ', near_pairs=', state%npairs_near, &
      ', overlap_pairs=', state%npairs_overlap, ', min_gap=', state%min_gap, &
      ', max_penetration=', state%max_penetration
  end function fibre_prod_collision_summary

end module fibre_prod_collision_diagnostics
