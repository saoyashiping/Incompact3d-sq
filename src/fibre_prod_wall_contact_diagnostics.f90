module fibre_prod_wall_contact_diagnostics
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_wall_contact, only : fibre_prod_wall_contact_state_type
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_wall_contact_count_near_nodes
  public :: fibre_prod_wall_contact_count_penetrating_nodes
  public :: fibre_prod_wall_contact_force_norm
  public :: fibre_prod_wall_contact_energy_candidate
  public :: fibre_prod_wall_contact_diagnostics_finite
  public :: fibre_prod_wall_contact_summary

contains

  pure integer function fibre_prod_wall_contact_count_near_nodes(gaps, warning_gap) result(count_value)
    real(dp), intent(in) :: gaps(:)
    real(dp), intent(in) :: warning_gap

    if (.not. all(ieee_is_finite(gaps)) .or. .not. ieee_is_finite(warning_gap)) then
      count_value = -1
    else
      count_value = count(gaps >= 0.0_dp .and. gaps <= warning_gap)
    end if
  end function fibre_prod_wall_contact_count_near_nodes

  pure integer function fibre_prod_wall_contact_count_penetrating_nodes(gaps) result(count_value)
    real(dp), intent(in) :: gaps(:)

    if (.not. all(ieee_is_finite(gaps))) then
      count_value = -1
    else
      count_value = count(gaps < 0.0_dp)
    end if
  end function fibre_prod_wall_contact_count_penetrating_nodes

  pure real(dp) function fibre_prod_wall_contact_force_norm(force) result(norm_value)
    real(dp), intent(in) :: force(:, :)

    if (all(ieee_is_finite(force))) then
      norm_value = sqrt(sum(force * force))
    else
      norm_value = huge(1.0_dp)
    end if
  end function fibre_prod_wall_contact_force_norm

  pure real(dp) function fibre_prod_wall_contact_energy_candidate(state, k_wall) result(energy)
    type(fibre_prod_wall_contact_state_type), intent(in) :: state
    real(dp), intent(in) :: k_wall

    if (ieee_is_finite(state%max_penetration) .and. ieee_is_finite(k_wall) .and. k_wall >= 0.0_dp) then
      energy = 0.5_dp * k_wall * state%max_penetration * state%max_penetration
    else
      energy = huge(1.0_dp)
    end if
  end function fibre_prod_wall_contact_energy_candidate

  pure logical function fibre_prod_wall_contact_diagnostics_finite(state) result(is_finite)
    type(fibre_prod_wall_contact_state_type), intent(in) :: state

    is_finite = ieee_is_finite(state%min_gap) .and. ieee_is_finite(state%max_penetration) .and. &
                ieee_is_finite(state%normal_force_norm)
  end function fibre_prod_wall_contact_diagnostics_finite

  function fibre_prod_wall_contact_summary(state) result(summary)
    type(fibre_prod_wall_contact_state_type), intent(in) :: state
    character(len=160) :: summary

    write(summary, '(A,L1,A,L1,A,L1,A,I0,A,ES12.4,A,ES12.4)') &
      'contact=', state%contact_active, ', warning=', state%near_wall_warning, &
      ', penetration=', state%penetration_detected, ', wall=', state%nearest_wall, &
      ', min_gap=', state%min_gap, ', max_penetration=', state%max_penetration
  end function fibre_prod_wall_contact_summary

end module fibre_prod_wall_contact_diagnostics
