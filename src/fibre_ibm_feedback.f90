module fibre_ibm_feedback

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_lagrangian_points_t

  implicit none
  private

  public :: compute_ibm_feedback_forces
  public :: compute_feedback_slip_norms
  public :: compute_feedback_pair_power
  public :: compute_feedback_expected_dissipation

contains

  subroutine compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(in) :: u_lag(:,:)
    real(mytype), intent(in) :: beta_drag
    real(mytype), intent(out) :: force_on_structure(:,:), force_on_fluid(:,:)

    integer :: l
    real(mytype) :: slip(3)

    if (beta_drag <= 0._mytype) then
      error stop 'compute_ibm_feedback_forces: beta_drag must be > 0'
    end if
    if (size(u_lag, 1) /= 3 .or. size(u_lag, 2) /= lag%nl) then
      error stop 'compute_ibm_feedback_forces: u_lag shape mismatch with lag%nl'
    end if
    if (size(force_on_structure, 1) /= 3 .or. size(force_on_structure, 2) /= lag%nl) then
      error stop 'compute_ibm_feedback_forces: force_on_structure shape mismatch with lag%nl'
    end if
    if (size(force_on_fluid, 1) /= 3 .or. size(force_on_fluid, 2) /= lag%nl) then
      error stop 'compute_ibm_feedback_forces: force_on_fluid shape mismatch with lag%nl'
    end if
    if (size(lag%v, 1) /= 3 .or. size(lag%v, 2) /= lag%nl) then
      error stop 'compute_ibm_feedback_forces: lag%v shape mismatch with lag%nl'
    end if

    do l = 1, lag%nl
      slip = u_lag(:, l) - lag%v(:, l)
      force_on_structure(:, l) = beta_drag * slip
      force_on_fluid(:, l) = -force_on_structure(:, l)
    end do

    ! force_on_structure should be passed to the Stage-2 fibre external-force interface.
    ! force_on_fluid should be spread to the Eulerian IBM force density.
  end subroutine compute_ibm_feedback_forces

  subroutine compute_feedback_slip_norms(lag, u_lag, max_slip_norm, rms_slip_norm)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(in) :: u_lag(:,:)
    real(mytype), intent(out) :: max_slip_norm, rms_slip_norm

    integer :: l
    real(mytype) :: slip_sq, weighted_slip_sq_sum, weight_sum

    if (size(u_lag, 1) /= 3 .or. size(u_lag, 2) /= lag%nl) then
      error stop 'compute_feedback_slip_norms: u_lag shape mismatch with lag%nl'
    end if
    if (size(lag%v, 1) /= 3 .or. size(lag%v, 2) /= lag%nl) then
      error stop 'compute_feedback_slip_norms: lag%v shape mismatch with lag%nl'
    end if
    if (size(lag%weight) /= lag%nl) then
      error stop 'compute_feedback_slip_norms: lag%weight size mismatch with lag%nl'
    end if

    max_slip_norm = 0._mytype
    weighted_slip_sq_sum = 0._mytype
    weight_sum = 0._mytype

    do l = 1, lag%nl
      slip_sq = sum((u_lag(:, l) - lag%v(:, l))**2)
      max_slip_norm = max(max_slip_norm, sqrt(slip_sq))
      weighted_slip_sq_sum = weighted_slip_sq_sum + lag%weight(l) * slip_sq
      weight_sum = weight_sum + lag%weight(l)
    end do

    if (weight_sum <= 0._mytype) then
      error stop 'compute_feedback_slip_norms: total lagrangian weight must be > 0'
    end if
    rms_slip_norm = sqrt(weighted_slip_sq_sum / weight_sum)
  end subroutine compute_feedback_slip_norms

  subroutine compute_feedback_pair_power(lag, u_lag, force_on_structure, force_on_fluid, &
                                         power_structure, power_fluid, power_total)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(in) :: u_lag(:,:)
    real(mytype), intent(in) :: force_on_structure(:,:), force_on_fluid(:,:)
    real(mytype), intent(out) :: power_structure, power_fluid, power_total

    integer :: l

    if (size(u_lag, 1) /= 3 .or. size(u_lag, 2) /= lag%nl) then
      error stop 'compute_feedback_pair_power: u_lag shape mismatch with lag%nl'
    end if
    if (size(force_on_structure, 1) /= 3 .or. size(force_on_structure, 2) /= lag%nl) then
      error stop 'compute_feedback_pair_power: force_on_structure shape mismatch with lag%nl'
    end if
    if (size(force_on_fluid, 1) /= 3 .or. size(force_on_fluid, 2) /= lag%nl) then
      error stop 'compute_feedback_pair_power: force_on_fluid shape mismatch with lag%nl'
    end if
    if (size(lag%v, 1) /= 3 .or. size(lag%v, 2) /= lag%nl) then
      error stop 'compute_feedback_pair_power: lag%v shape mismatch with lag%nl'
    end if
    if (size(lag%weight) /= lag%nl) then
      error stop 'compute_feedback_pair_power: lag%weight size mismatch with lag%nl'
    end if

    power_structure = 0._mytype
    power_fluid = 0._mytype

    do l = 1, lag%nl
      power_structure = power_structure + dot_product(force_on_structure(:, l), lag%v(:, l)) * lag%weight(l)
      power_fluid = power_fluid + dot_product(force_on_fluid(:, l), u_lag(:, l)) * lag%weight(l)
    end do

    power_total = power_structure + power_fluid
  end subroutine compute_feedback_pair_power

  subroutine compute_feedback_expected_dissipation(lag, u_lag, beta_drag, dissipation)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(in) :: u_lag(:,:)
    real(mytype), intent(in) :: beta_drag
    real(mytype), intent(out) :: dissipation

    integer :: l

    if (beta_drag <= 0._mytype) then
      error stop 'compute_feedback_expected_dissipation: beta_drag must be > 0'
    end if
    if (size(u_lag, 1) /= 3 .or. size(u_lag, 2) /= lag%nl) then
      error stop 'compute_feedback_expected_dissipation: u_lag shape mismatch with lag%nl'
    end if
    if (size(lag%v, 1) /= 3 .or. size(lag%v, 2) /= lag%nl) then
      error stop 'compute_feedback_expected_dissipation: lag%v shape mismatch with lag%nl'
    end if
    if (size(lag%weight) /= lag%nl) then
      error stop 'compute_feedback_expected_dissipation: lag%weight size mismatch with lag%nl'
    end if

    dissipation = 0._mytype
    do l = 1, lag%nl
      dissipation = dissipation + beta_drag * sum((u_lag(:, l) - lag%v(:, l))**2) * lag%weight(l)
    end do
  end subroutine compute_feedback_expected_dissipation

end module fibre_ibm_feedback
