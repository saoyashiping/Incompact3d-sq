module fibre_ibm_power_diagnostics

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t

  implicit none
  private

  public :: compute_eulerian_power
  public :: compute_lagrangian_power
  public :: compute_power_consistency_error

contains

  subroutine compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, power)
    type(ibm_grid_t), intent(in) :: grid
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(out) :: power

    if (size(ux, 1) /= grid%nx .or. size(ux, 2) /= grid%ny .or. size(ux, 3) /= grid%nz) then
      error stop 'compute_eulerian_power: ux shape mismatch with grid'
    end if
    if (size(uy, 1) /= grid%nx .or. size(uy, 2) /= grid%ny .or. size(uy, 3) /= grid%nz) then
      error stop 'compute_eulerian_power: uy shape mismatch with grid'
    end if
    if (size(uz, 1) /= grid%nx .or. size(uz, 2) /= grid%ny .or. size(uz, 3) /= grid%nz) then
      error stop 'compute_eulerian_power: uz shape mismatch with grid'
    end if
    if (size(fx, 1) /= grid%nx .or. size(fx, 2) /= grid%ny .or. size(fx, 3) /= grid%nz) then
      error stop 'compute_eulerian_power: fx shape mismatch with grid'
    end if
    if (size(fy, 1) /= grid%nx .or. size(fy, 2) /= grid%ny .or. size(fy, 3) /= grid%nz) then
      error stop 'compute_eulerian_power: fy shape mismatch with grid'
    end if
    if (size(fz, 1) /= grid%nx .or. size(fz, 2) /= grid%ny .or. size(fz, 3) /= grid%nz) then
      error stop 'compute_eulerian_power: fz shape mismatch with grid'
    end if

    power = sum(fx * ux + fy * uy + fz * uz) * grid%cell_volume
  end subroutine compute_eulerian_power

  subroutine compute_lagrangian_power(lag, u_lag, power)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(in) :: u_lag(:,:)
    real(mytype), intent(out) :: power

    integer :: l

    if (size(u_lag, 1) /= 3 .or. size(u_lag, 2) /= lag%nl) then
      error stop 'compute_lagrangian_power: u_lag shape mismatch with lag%nl'
    end if
    if (size(lag%force, 1) /= 3 .or. size(lag%force, 2) /= lag%nl) then
      error stop 'compute_lagrangian_power: lag%force shape mismatch with lag%nl'
    end if
    if (size(lag%weight) /= lag%nl) then
      error stop 'compute_lagrangian_power: lag%weight size mismatch with lag%nl'
    end if

    power = 0._mytype
    do l = 1, lag%nl
      power = power + dot_product(lag%force(:, l), u_lag(:, l)) * lag%weight(l)
    end do
  end subroutine compute_lagrangian_power

  subroutine compute_power_consistency_error(power_e, power_l, abs_error, rel_error)
    real(mytype), intent(in) :: power_e, power_l
    real(mytype), intent(out) :: abs_error, rel_error

    real(mytype), parameter :: tiny_value = 1.0e-300_mytype

    abs_error = abs(power_e - power_l)
    rel_error = abs_error / max(abs(power_l), tiny_value)
  end subroutine compute_power_consistency_error

end module fibre_ibm_power_diagnostics
