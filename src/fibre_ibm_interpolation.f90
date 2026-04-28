module fibre_ibm_interpolation

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_delta, only : delta_weight_3d

  implicit none
  private

  public :: interpolate_scalar_to_lag
  public :: interpolate_vector_to_lag
  public :: compute_interpolation_weight_sum

contains

  subroutine interpolate_scalar_to_lag(grid, scalar, lag, q_lag)
    type(ibm_grid_t), intent(in) :: grid
    real(mytype), intent(in) :: scalar(:,:,:)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(out) :: q_lag(:)

    integer :: l, i, j, k
    real(mytype) :: w

    if (size(scalar, 1) /= grid%nx .or. size(scalar, 2) /= grid%ny .or. size(scalar, 3) /= grid%nz) then
      error stop 'interpolate_scalar_to_lag: scalar shape mismatch with grid'
    end if
    if (size(q_lag) /= lag%nl) then
      error stop 'interpolate_scalar_to_lag: q_lag size mismatch with lag%nl'
    end if

    q_lag = 0._mytype

    do l = 1, lag%nl
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            w = delta_weight_3d(grid, lag%x(:, l), i, j, k)
            q_lag(l) = q_lag(l) + scalar(i, j, k) * w
          end do
        end do
      end do
    end do
  end subroutine interpolate_scalar_to_lag

  subroutine interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
    type(ibm_grid_t), intent(in) :: grid
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(out) :: u_lag(3, lag%nl)

    ! Stage 3.2 uses a synthetic collocated cell-centered grid.
    ! Xcompact3D staggered/actual velocity locations must be handled in a later step.

    call interpolate_scalar_to_lag(grid, ux, lag, u_lag(1, :))
    call interpolate_scalar_to_lag(grid, uy, lag, u_lag(2, :))
    call interpolate_scalar_to_lag(grid, uz, lag, u_lag(3, :))
  end subroutine interpolate_vector_to_lag

  subroutine compute_interpolation_weight_sum(grid, lag, weight_sum)
    type(ibm_grid_t), intent(in) :: grid
    type(ibm_lagrangian_points_t), intent(in) :: lag
    real(mytype), intent(out) :: weight_sum(lag%nl)

    integer :: l, i, j, k

    weight_sum = 0._mytype

    do l = 1, lag%nl
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            weight_sum(l) = weight_sum(l) + delta_weight_3d(grid, lag%x(:, l), i, j, k)
          end do
        end do
      end do
    end do
  end subroutine compute_interpolation_weight_sum

end module fibre_ibm_interpolation
