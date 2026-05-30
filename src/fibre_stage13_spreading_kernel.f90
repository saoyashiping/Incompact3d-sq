module fibre_stage13_spreading_kernel
  use decomp_2d_constants, only : mytype
  implicit none
  private

  real(mytype), parameter :: weight_sum_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: force_conservation_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: zero_force_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: compact_support_tol = 0.0_mytype

  integer :: initialized_status = 0
  integer :: zero_force_spreading_status = 0
  integer :: single_point_spreading_status = 0
  integer :: compact_support_status = 1
  integer :: weight_normalization_status = 1
  integer :: nonnegative_weight_status = 1
  integer :: component_spreading_status = 1
  integer :: boundary_safe_status = 1
  integer :: integrated_force_conservation_status = 1
  integer :: finite_force_density_status = 1
  integer :: force_density_norm_finite_status = 1
  integer :: clear_status = 0
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: spreading_kernel_status = 0

  real(mytype) :: zero_force_max_abs = 0.0_mytype
  integer :: single_point_support_count = 0
  real(mytype) :: max_weight_sum_error = 0.0_mytype
  real(mytype) :: min_weight = 0.0_mytype
  real(mytype) :: integrated_force_error_x = 0.0_mytype
  real(mytype) :: integrated_force_error_y = 0.0_mytype
  real(mytype) :: integrated_force_error_z = 0.0_mytype
  real(mytype) :: integrated_force_error_l2 = 0.0_mytype
  real(mytype) :: max_abs_force_density = 0.0_mytype
  real(mytype) :: max_abs_force_density_norm_after_clear = 0.0_mytype

  public :: stage13_spreading_kernel_init
  public :: stage13_spreading_kernel_clear
  public :: stage13_spreading_kernel_spread_controlled
  public :: stage13_spreading_kernel_get_status_values
  public :: stage13_spreading_kernel_get_diagnostics
  public :: stage13_spreading_kernel_write_diagnostics
  public :: stage13_spreading_kernel_finalize

contains

  subroutine stage13_spreading_kernel_init()
    initialized_status = 1
    zero_force_spreading_status = 0
    single_point_spreading_status = 0
    compact_support_status = 1
    weight_normalization_status = 1
    nonnegative_weight_status = 1
    component_spreading_status = 1
    boundary_safe_status = 1
    integrated_force_conservation_status = 1
    finite_force_density_status = 1
    force_density_norm_finite_status = 1
    clear_status = 0
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    call reset_diagnostics()
    call update_kernel_status()
  end subroutine stage13_spreading_kernel_init

  subroutine stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
    real(mytype), intent(inout) :: fx_cand(:,:,:)
    real(mytype), intent(inout) :: fy_cand(:,:,:)
    real(mytype), intent(inout) :: fz_cand(:,:,:)
    real(mytype), intent(inout) :: force_density_norm(:,:,:)

    fx_cand(:,:,:) = 0.0_mytype
    fy_cand(:,:,:) = 0.0_mytype
    fz_cand(:,:,:) = 0.0_mytype
    force_density_norm(:,:,:) = 0.0_mytype
    max_abs_force_density_norm_after_clear = maxval(abs(force_density_norm))
    clear_status = merge(1, 0, max_abs_force_density_norm_after_clear <= zero_force_abs_tol)
    call update_kernel_status()
  end subroutine stage13_spreading_kernel_clear

  subroutine stage13_spreading_kernel_spread_controlled(x_lag, y_lag, z_lag, ds_lag, f_lag, &
                                                       x_eul, y_eul, z_eul, volume_eul, &
                                                       fx_cand, fy_cand, fz_cand, force_density_norm)
    real(mytype), intent(in) :: x_lag(:)
    real(mytype), intent(in) :: y_lag(:)
    real(mytype), intent(in) :: z_lag(:)
    real(mytype), intent(in) :: ds_lag(:)
    real(mytype), intent(in) :: f_lag(:,:)
    real(mytype), intent(in) :: x_eul(:)
    real(mytype), intent(in) :: y_eul(:)
    real(mytype), intent(in) :: z_eul(:)
    real(mytype), intent(in) :: volume_eul(:,:,:)
    real(mytype), intent(inout) :: fx_cand(:,:,:)
    real(mytype), intent(inout) :: fy_cand(:,:,:)
    real(mytype), intent(inout) :: fz_cand(:,:,:)
    real(mytype), intent(inout) :: force_density_norm(:,:,:)

    integer :: p
    integer :: i0
    integer :: j0
    integer :: k0
    integer :: di
    integer :: dj
    integer :: dk
    integer :: ii
    integer :: jj
    integer :: kk
    integer :: nlag
    integer :: support_count
    integer :: max_support_count
    integer :: valid_point_count
    logical :: found_x
    logical :: found_y
    logical :: found_z
    real(mytype) :: tx
    real(mytype) :: ty
    real(mytype) :: tz
    real(mytype) :: wx(2)
    real(mytype) :: wy(2)
    real(mytype) :: wz(2)
    real(mytype) :: weight
    real(mytype) :: weight_sum
    real(mytype) :: local_min_weight
    real(mytype) :: lag_force_x
    real(mytype) :: lag_force_y
    real(mytype) :: lag_force_z
    real(mytype) :: eul_force_x
    real(mytype) :: eul_force_y
    real(mytype) :: eul_force_z
    real(mytype) :: force_scale
    real(mytype) :: max_input_force

    call reset_diagnostics()
    ! Stage 13.3 executes multiple independent controlled spreading tests
    ! through the same module instance.  Reset all per-spreading-call
    ! result statuses here so that a successful zero-force test cannot be
    ! inherited by later non-zero-force tests, and vice versa.
    zero_force_spreading_status = 0
    single_point_spreading_status = 0
    compact_support_status = 1
    weight_normalization_status = 1
    nonnegative_weight_status = 1
    component_spreading_status = 1
    boundary_safe_status = 1
    integrated_force_conservation_status = 1
    finite_force_density_status = 1
    force_density_norm_finite_status = 1
    nlag = size(x_lag)
    max_support_count = 0
    valid_point_count = 0
    local_min_weight = huge(1.0_mytype)
    lag_force_x = 0.0_mytype
    lag_force_y = 0.0_mytype
    lag_force_z = 0.0_mytype
    max_input_force = 0.0_mytype
    boundary_safe_status = 1
    compact_support_status = 1
    weight_normalization_status = 1
    nonnegative_weight_status = 1

    if (.not. input_shapes_valid(x_lag, y_lag, z_lag, ds_lag, f_lag, x_eul, y_eul, z_eul, &
                                 volume_eul, fx_cand, fy_cand, fz_cand, force_density_norm)) then
      boundary_safe_status = 0
      call update_norm_and_status(fx_cand, fy_cand, fz_cand, force_density_norm)
      call update_kernel_status()
      return
    end if

    do p = 1, nlag
      max_input_force = max(max_input_force, abs(f_lag(p,1)), abs(f_lag(p,2)), abs(f_lag(p,3)))
      call find_bracket(x_lag(p), x_eul, i0, tx, found_x)
      call find_bracket(y_lag(p), y_eul, j0, ty, found_y)
      call find_bracket(z_lag(p), z_eul, k0, tz, found_z)
      if (.not. found_x .or. .not. found_y .or. .not. found_z) then
        boundary_safe_status = 0
        cycle
      end if

      valid_point_count = valid_point_count + 1
      wx(1) = 1.0_mytype - tx
      wx(2) = tx
      wy(1) = 1.0_mytype - ty
      wy(2) = ty
      wz(1) = 1.0_mytype - tz
      wz(2) = tz
      weight_sum = 0.0_mytype
      support_count = 0

      do dk = 0, 1
        kk = k0 + dk
        do dj = 0, 1
          jj = j0 + dj
          do di = 0, 1
            ii = i0 + di
            weight = wx(di + 1) * wy(dj + 1) * wz(dk + 1)
            weight_sum = weight_sum + weight
            local_min_weight = min(local_min_weight, weight)
            if (weight < -weight_sum_abs_tol) nonnegative_weight_status = 0
            if (abs(weight) > compact_support_tol) support_count = support_count + 1
            if (volume_eul(ii,jj,kk) <= 0.0_mytype) then
              boundary_safe_status = 0
            else
              force_scale = weight * ds_lag(p) / volume_eul(ii,jj,kk)
              fx_cand(ii,jj,kk) = fx_cand(ii,jj,kk) + f_lag(p,1) * force_scale
              fy_cand(ii,jj,kk) = fy_cand(ii,jj,kk) + f_lag(p,2) * force_scale
              fz_cand(ii,jj,kk) = fz_cand(ii,jj,kk) + f_lag(p,3) * force_scale
            end if
          end do
        end do
      end do

      max_support_count = max(max_support_count, support_count)
      max_weight_sum_error = max(max_weight_sum_error, abs(weight_sum - 1.0_mytype))
      if (support_count > 8) compact_support_status = 0
      if (abs(weight_sum - 1.0_mytype) > weight_sum_abs_tol) weight_normalization_status = 0
      lag_force_x = lag_force_x + f_lag(p,1) * ds_lag(p)
      lag_force_y = lag_force_y + f_lag(p,2) * ds_lag(p)
      lag_force_z = lag_force_z + f_lag(p,3) * ds_lag(p)
    end do

    if (local_min_weight == huge(1.0_mytype)) local_min_weight = 0.0_mytype
    min_weight = local_min_weight
    single_point_support_count = max_support_count
    if (nlag == 1 .and. max_input_force > zero_force_abs_tol .and. valid_point_count == 1 .and. &
        max_support_count > 0 .and. max_support_count <= 8) then
      single_point_spreading_status = 1
    end if

    call update_norm_and_status(fx_cand, fy_cand, fz_cand, force_density_norm)
    eul_force_x = sum(fx_cand * volume_eul)
    eul_force_y = sum(fy_cand * volume_eul)
    eul_force_z = sum(fz_cand * volume_eul)
    integrated_force_error_x = eul_force_x - lag_force_x
    integrated_force_error_y = eul_force_y - lag_force_y
    integrated_force_error_z = eul_force_z - lag_force_z
    integrated_force_error_l2 = sqrt(integrated_force_error_x * integrated_force_error_x + &
                                     integrated_force_error_y * integrated_force_error_y + &
                                     integrated_force_error_z * integrated_force_error_z)
    component_spreading_status = merge(1, 0, abs(integrated_force_error_x) <= force_conservation_abs_tol .and. &
                                            abs(integrated_force_error_y) <= force_conservation_abs_tol .and. &
                                            abs(integrated_force_error_z) <= force_conservation_abs_tol)
    integrated_force_conservation_status = merge(1, 0, integrated_force_error_l2 <= force_conservation_abs_tol)
    zero_force_max_abs = max(maxval(abs(fx_cand)), maxval(abs(fy_cand)), maxval(abs(fz_cand)))
    if (max_input_force <= zero_force_abs_tol .and. zero_force_max_abs <= zero_force_abs_tol .and. &
        integrated_force_error_l2 <= force_conservation_abs_tol) then
      zero_force_spreading_status = 1
    end if
    call update_kernel_status()
  end subroutine stage13_spreading_kernel_spread_controlled

  subroutine stage13_spreading_kernel_get_status_values(initialized_out, zero_force_out, single_point_out, &
                                                       compact_support_out, weight_normalization_out, &
                                                       nonnegative_weight_out, component_spreading_out, &
                                                       boundary_safe_out, integrated_force_conservation_out, &
                                                       finite_force_density_out, norm_finite_out, clear_out, &
                                                       no_rhs_injection_out, no_ibm_spreading_out, &
                                                       no_feedback_application_out, no_twoway_force_out, &
                                                       no_structure_advance_out, no_fluid_field_access_out, &
                                                       no_fluid_field_modification_out, kernel_status_out)
    integer, intent(out) :: initialized_out
    integer, intent(out) :: zero_force_out
    integer, intent(out) :: single_point_out
    integer, intent(out) :: compact_support_out
    integer, intent(out) :: weight_normalization_out
    integer, intent(out) :: nonnegative_weight_out
    integer, intent(out) :: component_spreading_out
    integer, intent(out) :: boundary_safe_out
    integer, intent(out) :: integrated_force_conservation_out
    integer, intent(out) :: finite_force_density_out
    integer, intent(out) :: norm_finite_out
    integer, intent(out) :: clear_out
    integer, intent(out) :: no_rhs_injection_out
    integer, intent(out) :: no_ibm_spreading_out
    integer, intent(out) :: no_feedback_application_out
    integer, intent(out) :: no_twoway_force_out
    integer, intent(out) :: no_structure_advance_out
    integer, intent(out) :: no_fluid_field_access_out
    integer, intent(out) :: no_fluid_field_modification_out
    integer, intent(out) :: kernel_status_out

    call update_kernel_status()
    initialized_out = initialized_status
    zero_force_out = zero_force_spreading_status
    single_point_out = single_point_spreading_status
    compact_support_out = compact_support_status
    weight_normalization_out = weight_normalization_status
    nonnegative_weight_out = nonnegative_weight_status
    component_spreading_out = component_spreading_status
    boundary_safe_out = boundary_safe_status
    integrated_force_conservation_out = integrated_force_conservation_status
    finite_force_density_out = finite_force_density_status
    norm_finite_out = force_density_norm_finite_status
    clear_out = clear_status
    no_rhs_injection_out = no_rhs_injection_status
    no_ibm_spreading_out = no_ibm_spreading_status
    no_feedback_application_out = no_feedback_application_status
    no_twoway_force_out = no_twoway_force_status
    no_structure_advance_out = no_structure_advance_status
    no_fluid_field_access_out = no_fluid_field_access_status
    no_fluid_field_modification_out = no_fluid_field_modification_status
    kernel_status_out = spreading_kernel_status
  end subroutine stage13_spreading_kernel_get_status_values

  subroutine stage13_spreading_kernel_get_diagnostics(zero_max_out, support_count_out, weight_error_out, &
                                                     min_weight_out, error_x_out, error_y_out, error_z_out, &
                                                     error_l2_out, max_force_out, norm_after_clear_out)
    real(mytype), intent(out) :: zero_max_out
    integer, intent(out) :: support_count_out
    real(mytype), intent(out) :: weight_error_out
    real(mytype), intent(out) :: min_weight_out
    real(mytype), intent(out) :: error_x_out
    real(mytype), intent(out) :: error_y_out
    real(mytype), intent(out) :: error_z_out
    real(mytype), intent(out) :: error_l2_out
    real(mytype), intent(out) :: max_force_out
    real(mytype), intent(out) :: norm_after_clear_out

    zero_max_out = zero_force_max_abs
    support_count_out = single_point_support_count
    weight_error_out = max_weight_sum_error
    min_weight_out = min_weight
    error_x_out = integrated_force_error_x
    error_y_out = integrated_force_error_y
    error_z_out = integrated_force_error_z
    error_l2_out = integrated_force_error_l2
    max_force_out = max_abs_force_density
    norm_after_clear_out = max_abs_force_density_norm_after_clear
  end subroutine stage13_spreading_kernel_get_diagnostics

  subroutine stage13_spreading_kernel_write_diagnostics()
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file='stage13_outputs/fibre_stage13_spreading_kernel.dat', &
         status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return
    call update_kernel_status()
    write(unit_id,'(A,1X,I0)') 'stage13_spreading_kernel_initialized_status', initialized_status
    write(unit_id,'(A,1X,I0)') 'stage13_spreading_kernel_weight_normalization_status', weight_normalization_status
    write(unit_id,'(A,1X,I0)') 'stage13_spreading_kernel_conservation_status', integrated_force_conservation_status
    write(unit_id,'(A,1X,I0)') 'stage13_spreading_kernel_status', spreading_kernel_status
    close(unit_id)
  end subroutine stage13_spreading_kernel_write_diagnostics

  subroutine stage13_spreading_kernel_finalize()
    initialized_status = 0
    clear_status = 0
    call reset_diagnostics()
    call update_kernel_status()
  end subroutine stage13_spreading_kernel_finalize

  subroutine reset_diagnostics()
    zero_force_max_abs = 0.0_mytype
    single_point_support_count = 0
    max_weight_sum_error = 0.0_mytype
    min_weight = 0.0_mytype
    integrated_force_error_x = 0.0_mytype
    integrated_force_error_y = 0.0_mytype
    integrated_force_error_z = 0.0_mytype
    integrated_force_error_l2 = 0.0_mytype
    max_abs_force_density = 0.0_mytype
    max_abs_force_density_norm_after_clear = 0.0_mytype
  end subroutine reset_diagnostics

  subroutine update_kernel_status()
    if (initialized_status == 1 .and. compact_support_status == 1 .and. &
        weight_normalization_status == 1 .and. nonnegative_weight_status == 1 .and. &
        component_spreading_status == 1 .and. boundary_safe_status == 1 .and. &
        integrated_force_conservation_status == 1 .and. finite_force_density_status == 1 .and. &
        force_density_norm_finite_status == 1 .and. no_rhs_injection_status == 1 .and. &
        no_ibm_spreading_status == 1 .and. no_feedback_application_status == 1 .and. &
        no_twoway_force_status == 1 .and. no_structure_advance_status == 1 .and. &
        no_fluid_field_access_status == 1 .and. no_fluid_field_modification_status == 1) then
      spreading_kernel_status = 1
    else
      spreading_kernel_status = 0
    end if
  end subroutine update_kernel_status

  subroutine update_norm_and_status(fx_cand, fy_cand, fz_cand, force_density_norm)
    real(mytype), intent(in) :: fx_cand(:,:,:)
    real(mytype), intent(in) :: fy_cand(:,:,:)
    real(mytype), intent(in) :: fz_cand(:,:,:)
    real(mytype), intent(inout) :: force_density_norm(:,:,:)

    force_density_norm(:,:,:) = sqrt(fx_cand * fx_cand + fy_cand * fy_cand + fz_cand * fz_cand)
    finite_force_density_status = merge(1, 0, all_finite_rank3(fx_cand) .and. all_finite_rank3(fy_cand) .and. &
                                             all_finite_rank3(fz_cand))
    force_density_norm_finite_status = merge(1, 0, all_finite_rank3(force_density_norm))
    max_abs_force_density = max(maxval(abs(fx_cand)), maxval(abs(fy_cand)), maxval(abs(fz_cand)))
  end subroutine update_norm_and_status

  logical function input_shapes_valid(x_lag, y_lag, z_lag, ds_lag, f_lag, x_eul, y_eul, z_eul, &
                                      volume_eul, fx_cand, fy_cand, fz_cand, force_density_norm)
    real(mytype), intent(in) :: x_lag(:)
    real(mytype), intent(in) :: y_lag(:)
    real(mytype), intent(in) :: z_lag(:)
    real(mytype), intent(in) :: ds_lag(:)
    real(mytype), intent(in) :: f_lag(:,:)
    real(mytype), intent(in) :: x_eul(:)
    real(mytype), intent(in) :: y_eul(:)
    real(mytype), intent(in) :: z_eul(:)
    real(mytype), intent(in) :: volume_eul(:,:,:)
    real(mytype), intent(in) :: fx_cand(:,:,:)
    real(mytype), intent(in) :: fy_cand(:,:,:)
    real(mytype), intent(in) :: fz_cand(:,:,:)
    real(mytype), intent(in) :: force_density_norm(:,:,:)

    input_shapes_valid = size(x_lag) == size(y_lag) .and. size(x_lag) == size(z_lag) .and. &
                         size(x_lag) == size(ds_lag) .and. size(f_lag, 1) == size(x_lag) .and. &
                         size(f_lag, 2) >= 3 .and. size(x_eul) >= 2 .and. size(y_eul) >= 2 .and. &
                         size(z_eul) >= 2 .and. size(volume_eul, 1) == size(x_eul) .and. &
                         size(volume_eul, 2) == size(y_eul) .and. size(volume_eul, 3) == size(z_eul) .and. &
                         all(shape(fx_cand) == shape(volume_eul)) .and. all(shape(fy_cand) == shape(volume_eul)) .and. &
                         all(shape(fz_cand) == shape(volume_eul)) .and. &
                         all(shape(force_density_norm) == shape(volume_eul)) .and. all(volume_eul > 0.0_mytype)
  end function input_shapes_valid

  subroutine find_bracket(xp, x_eul, lower_index, fraction, found)
    real(mytype), intent(in) :: xp
    real(mytype), intent(in) :: x_eul(:)
    integer, intent(out) :: lower_index
    real(mytype), intent(out) :: fraction
    logical, intent(out) :: found
    integer :: i
    integer :: n

    n = size(x_eul)
    found = .false.
    lower_index = 1
    fraction = 0.0_mytype
    if (n < 2) return
    if (xp < x_eul(1) .or. xp > x_eul(n)) return
    if (xp == x_eul(n)) then
      lower_index = n - 1
      fraction = 1.0_mytype
      found = .true.
      return
    end if
    do i = 1, n - 1
      if (xp >= x_eul(i) .and. xp <= x_eul(i + 1)) then
        lower_index = i
        if (x_eul(i + 1) == x_eul(i)) return
        fraction = (xp - x_eul(i)) / (x_eul(i + 1) - x_eul(i))
        found = .true.
        return
      end if
    end do
  end subroutine find_bracket

  logical function all_finite_rank3(values)
    real(mytype), intent(in) :: values(:,:,:)

    all_finite_rank3 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank3

end module fibre_stage13_spreading_kernel
