module fibre_ibm_spreading

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_delta, only : delta_weight_3d

  implicit none
  private

  public :: clear_eulerian_force_density
  public :: spread_lag_force_to_eulerian
  public :: compute_eulerian_total_force
  public :: compute_lagrangian_total_force
  public :: compute_force_density_norms

contains

  subroutine clear_eulerian_force_density(fx, fy, fz)
    real(mytype), intent(inout) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)

    fx = 0._mytype
    fy = 0._mytype
    fz = 0._mytype
  end subroutine clear_eulerian_force_density

  subroutine spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
    type(ibm_grid_t), intent(in) :: grid
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(inout) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)

    integer :: i, j, k, l
    real(mytype) :: w, coeff

    if (size(fx, 1) /= grid%nx .or. size(fx, 2) /= grid%ny .or. size(fx, 3) /= grid%nz) then
      error stop 'spread_lag_force_to_eulerian: fx shape mismatch with grid'
    end if
    if (size(fy, 1) /= grid%nx .or. size(fy, 2) /= grid%ny .or. size(fy, 3) /= grid%nz) then
      error stop 'spread_lag_force_to_eulerian: fy shape mismatch with grid'
    end if
    if (size(fz, 1) /= grid%nx .or. size(fz, 2) /= grid%ny .or. size(fz, 3) /= grid%nz) then
      error stop 'spread_lag_force_to_eulerian: fz shape mismatch with grid'
    end if
    if (grid%cell_volume <= 0._mytype) then
      error stop 'spread_lag_force_to_eulerian: grid%cell_volume must be > 0'
    end if
    if (size(lag%force, 1) /= 3 .or. size(lag%force, 2) /= lag%nl) then
      error stop 'spread_lag_force_to_eulerian: lag%force shape mismatch with lag%nl'
    end if
    if (size(lag%weight) /= lag%nl) then
      error stop 'spread_lag_force_to_eulerian: lag%weight size mismatch with lag%nl'
    end if

    do l = 1, lag%nl
      coeff = lag%weight(l) / grid%cell_volume
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            w = delta_weight_3d(grid, lag%x(:, l), i, j, k)
            fx(i, j, k) = fx(i, j, k) + lag%force(1, l) * coeff * w
            fy(i, j, k) = fy(i, j, k) + lag%force(2, l) * coeff * w
            fz(i, j, k) = fz(i, j, k) + lag%force(3, l) * coeff * w
          end do
        end do
      end do
    end do
  end subroutine spread_lag_force_to_eulerian

  subroutine compute_eulerian_total_force(grid, fx, fy, fz, total_force)
    type(ibm_grid_t), intent(in) :: grid
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(out) :: total_force(3)

    total_force(1) = sum(fx) * grid%cell_volume
    total_force(2) = sum(fy) * grid%cell_volume
    total_force(3) = sum(fz) * grid%cell_volume
  end subroutine compute_eulerian_total_force

  subroutine compute_lagrangian_total_force(lag, total_force)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(out) :: total_force(3)

    integer :: l

    total_force = 0._mytype
    do l = 1, lag%nl
      total_force = total_force + lag%force(:, l) * lag%weight(l)
    end do
  end subroutine compute_lagrangian_total_force

  subroutine compute_force_density_norms(fx, fy, fz, max_abs_force_density, l2_force_density)
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(out) :: max_abs_force_density, l2_force_density

    max_abs_force_density = max(maxval(abs(fx)), max(maxval(abs(fy)), maxval(abs(fz))))
    l2_force_density = sqrt(sum(fx * fx + fy * fy + fz * fz))
  end subroutine compute_force_density_norms

end module fibre_ibm_spreading
