module fibre_stage11_grid_metadata
  use decomp_2d_constants, only : mytype
  implicit none
  private

  logical :: initialized = .false.
  integer :: nx_global = 0, ny_global = 0, nz_global = 0
  integer :: istart = 0, iend = 0, jstart = 0, jend = 0, kstart = 0, kend = 0
  real(mytype) :: xmin = 0.0_mytype, xmax = 0.0_mytype
  real(mytype) :: ymin = 0.0_mytype, ymax = 0.0_mytype
  real(mytype) :: zmin = 0.0_mytype, zmax = 0.0_mytype
  real(mytype) :: dx = 0.0_mytype, dy = 0.0_mytype, dz = 0.0_mytype
  integer :: periodic_x = 0, periodic_y = 0, periodic_z = 0
  integer :: staggered_layout_policy = 0
  integer :: nonuniform_y_policy = 0

  integer :: initialized_status = 0
  integer :: global_size_status = 0
  integer :: local_bounds_status = 0
  integer :: extents_finite_status = 0
  integer :: spacing_positive_status = 0
  integer :: periodic_flags_status = 0
  integer :: staggered_layout_status = 0
  integer :: nonuniform_y_policy_status = 0
  integer :: domain_safety_status = 0
  integer :: no_fluid_field_access_status = 1
  integer :: no_velocity_sampling_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_force_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: grid_metadata_status = 0

  public :: stage11_grid_metadata_init_uniform
  public :: stage11_grid_metadata_finalize
  public :: stage11_grid_metadata_is_initialized
  public :: stage11_grid_metadata_get_status_values
  public :: stage11_grid_metadata_point_in_domain
  public :: stage11_grid_metadata_write_diagnostics
  public :: stage11_grid_metadata_get_global_size
  public :: stage11_grid_metadata_get_local_bounds
  public :: stage11_grid_metadata_get_spacing
  public :: stage11_grid_metadata_get_periodic_flags

contains

  subroutine stage11_grid_metadata_init_uniform()
    call stage11_grid_metadata_finalize()

    nx_global = 16
    ny_global = 17
    nz_global = 16

    istart = 1
    iend = nx_global
    jstart = 1
    jend = ny_global
    kstart = 1
    kend = nz_global

    xmin = 0.0_mytype
    xmax = 2.0_mytype
    ymin = 0.0_mytype
    ymax = 1.0_mytype
    zmin = 0.0_mytype
    zmax = 1.0_mytype

    dx = (xmax - xmin) / real(max(1, nx_global - 1), mytype)
    dy = (ymax - ymin) / real(max(1, ny_global - 1), mytype)
    dz = (zmax - zmin) / real(max(1, nz_global - 1), mytype)

    periodic_x = 1
    periodic_y = 0
    periodic_z = 1

    staggered_layout_policy = 1
    nonuniform_y_policy = 0

    initialized = .true.
    call update_statuses()
  end subroutine stage11_grid_metadata_init_uniform

  subroutine stage11_grid_metadata_finalize()
    initialized = .false.
    nx_global = 0
    ny_global = 0
    nz_global = 0
    istart = 0
    iend = 0
    jstart = 0
    jend = 0
    kstart = 0
    kend = 0
    xmin = 0.0_mytype
    xmax = 0.0_mytype
    ymin = 0.0_mytype
    ymax = 0.0_mytype
    zmin = 0.0_mytype
    zmax = 0.0_mytype
    dx = 0.0_mytype
    dy = 0.0_mytype
    dz = 0.0_mytype
    periodic_x = 0
    periodic_y = 0
    periodic_z = 0
    staggered_layout_policy = 0
    nonuniform_y_policy = 0

    initialized_status = 0
    global_size_status = 0
    local_bounds_status = 0
    extents_finite_status = 0
    spacing_positive_status = 0
    periodic_flags_status = 0
    staggered_layout_status = 0
    nonuniform_y_policy_status = 0
    domain_safety_status = 0
    no_fluid_field_access_status = 1
    no_velocity_sampling_status = 1
    no_fluid_field_modification_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_force_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    grid_metadata_status = 0
  end subroutine stage11_grid_metadata_finalize

  logical function stage11_grid_metadata_is_initialized()
    stage11_grid_metadata_is_initialized = initialized
  end function stage11_grid_metadata_is_initialized

  logical function stage11_grid_metadata_point_in_domain(x, y, z)
    real(mytype), intent(in) :: x, y, z
    if (.not. initialized) then
      stage11_grid_metadata_point_in_domain = .false.
      return
    end if
    stage11_grid_metadata_point_in_domain = finite_value(x) .and. finite_value(y) .and. finite_value(z) .and. &
                                            x >= xmin .and. x <= xmax .and. &
                                            y >= ymin .and. y <= ymax .and. &
                                            z >= zmin .and. z <= zmax
  end function stage11_grid_metadata_point_in_domain

  subroutine stage11_grid_metadata_write_diagnostics(unit)
    integer, intent(in) :: unit
    write(unit,'(A,1X,I0)') 'stage11_2_grid_initialized_status', initialized_status
    write(unit,'(A,1X,I0)') 'stage11_2_global_size_status', global_size_status
    write(unit,'(A,1X,I0)') 'stage11_2_local_bounds_status', local_bounds_status
    write(unit,'(A,1X,I0)') 'stage11_2_extents_finite_status', extents_finite_status
    write(unit,'(A,1X,I0)') 'stage11_2_spacing_positive_status', spacing_positive_status
    write(unit,'(A,1X,I0)') 'stage11_2_periodic_flags_status', periodic_flags_status
    write(unit,'(A,1X,I0)') 'stage11_2_staggered_layout_status', staggered_layout_status
    write(unit,'(A,1X,I0)') 'stage11_2_nonuniform_y_policy_status', nonuniform_y_policy_status
    write(unit,'(A,1X,I0)') 'stage11_2_domain_safety_status', domain_safety_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_velocity_sampling_status', no_velocity_sampling_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_rhs_injection_status', no_rhs_injection_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_ibm_spreading_status', no_ibm_spreading_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_feedback_force_status', no_feedback_force_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_twoway_force_status', no_twoway_force_status
    write(unit,'(A,1X,I0)') 'stage11_2_no_structure_advance_status', no_structure_advance_status
    write(unit,'(A,1X,I0)') 'stage11_2_grid_metadata_status', grid_metadata_status
  end subroutine

  subroutine stage11_grid_metadata_get_global_size(nx, ny, nz)
    integer, intent(out) :: nx, ny, nz
    nx = nx_global; ny = ny_global; nz = nz_global
  end subroutine
  subroutine stage11_grid_metadata_get_local_bounds(is, ie, js, je, ks, ke)
    integer, intent(out) :: is, ie, js, je, ks, ke
    is = istart; ie = iend; js = jstart; je = jend; ks = kstart; ke = kend
  end subroutine
  subroutine stage11_grid_metadata_get_spacing(dx_out, dy_out, dz_out)
    real(mytype), intent(out) :: dx_out, dy_out, dz_out
    dx_out = dx; dy_out = dy; dz_out = dz
  end subroutine
  subroutine stage11_grid_metadata_get_periodic_flags(px, py, pz)
    integer, intent(out) :: px, py, pz
    px = periodic_x; py = periodic_y; pz = periodic_z
  end subroutine

  subroutine stage11_grid_metadata_get_status_values(out_initialized_status, out_global_size_status, &
                                                     out_local_bounds_status, out_extents_finite_status, &
                                                     out_spacing_positive_status, out_periodic_flags_status, &
                                                     out_staggered_layout_status, out_nonuniform_y_policy_status, &
                                                     out_domain_safety_status, out_no_fluid_field_access_status, &
                                                     out_no_velocity_sampling_status, out_no_fluid_field_modification_status, &
                                                     out_no_rhs_injection_status, out_no_ibm_spreading_status, &
                                                     out_no_feedback_force_status, out_no_twoway_force_status, &
                                                     out_no_structure_advance_status, out_grid_metadata_status)
    integer, intent(out) :: out_initialized_status, out_global_size_status, out_local_bounds_status
    integer, intent(out) :: out_extents_finite_status, out_spacing_positive_status, out_periodic_flags_status
    integer, intent(out) :: out_staggered_layout_status, out_nonuniform_y_policy_status, out_domain_safety_status
    integer, intent(out) :: out_no_fluid_field_access_status, out_no_velocity_sampling_status
    integer, intent(out) :: out_no_fluid_field_modification_status, out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status, out_no_feedback_force_status
    integer, intent(out) :: out_no_twoway_force_status, out_no_structure_advance_status, out_grid_metadata_status

    out_initialized_status = initialized_status
    out_global_size_status = global_size_status
    out_local_bounds_status = local_bounds_status
    out_extents_finite_status = extents_finite_status
    out_spacing_positive_status = spacing_positive_status
    out_periodic_flags_status = periodic_flags_status
    out_staggered_layout_status = staggered_layout_status
    out_nonuniform_y_policy_status = nonuniform_y_policy_status
    out_domain_safety_status = domain_safety_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_velocity_sampling_status = no_velocity_sampling_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_force_status = no_feedback_force_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_grid_metadata_status = grid_metadata_status
  end subroutine

  subroutine update_statuses()
    initialized_status = 0
    if (initialized) initialized_status = 1

    global_size_status = 0
    if (nx_global > 0 .and. ny_global > 0 .and. nz_global > 0) global_size_status = 1

    local_bounds_status = 0
    if (istart >= 1 .and. istart <= iend .and. iend <= nx_global .and. &
        jstart >= 1 .and. jstart <= jend .and. jend <= ny_global .and. &
        kstart >= 1 .and. kstart <= kend .and. kend <= nz_global) local_bounds_status = 1

    extents_finite_status = 0
    if (finite_value(xmin) .and. finite_value(xmax) .and. xmax > xmin .and. &
        finite_value(ymin) .and. finite_value(ymax) .and. ymax > ymin .and. &
        finite_value(zmin) .and. finite_value(zmax) .and. zmax > zmin) extents_finite_status = 1

    spacing_positive_status = 0
    if (finite_value(dx) .and. finite_value(dy) .and. finite_value(dz) .and. dx > 0.0_mytype .and. dy > 0.0_mytype .and. dz > 0.0_mytype) spacing_positive_status = 1

    periodic_flags_status = 0
    if ((periodic_x == 0 .or. periodic_x == 1) .and. (periodic_y == 0 .or. periodic_y == 1) .and. (periodic_z == 0 .or. periodic_z == 1)) periodic_flags_status = 1

    staggered_layout_status = 0
    if (staggered_layout_policy >= 0 .and. staggered_layout_policy <= 2) staggered_layout_status = 1

    nonuniform_y_policy_status = 0
    if (nonuniform_y_policy >= 0 .and. nonuniform_y_policy <= 2) nonuniform_y_policy_status = 1

    domain_safety_status = 0
    if (stage11_grid_metadata_point_in_domain(0.5_mytype*(xmin+xmax),0.5_mytype*(ymin+ymax),0.5_mytype*(zmin+zmax))) domain_safety_status = 1

    grid_metadata_status = 1
    if (initialized_status /= 1) grid_metadata_status = 0
    if (global_size_status /= 1) grid_metadata_status = 0
    if (local_bounds_status /= 1) grid_metadata_status = 0
    if (extents_finite_status /= 1) grid_metadata_status = 0
    if (spacing_positive_status /= 1) grid_metadata_status = 0
    if (periodic_flags_status /= 1) grid_metadata_status = 0
    if (staggered_layout_status /= 1) grid_metadata_status = 0
    if (nonuniform_y_policy_status /= 1) grid_metadata_status = 0
    if (domain_safety_status /= 1) grid_metadata_status = 0
    if (no_fluid_field_access_status /= 1) grid_metadata_status = 0
    if (no_velocity_sampling_status /= 1) grid_metadata_status = 0
    if (no_fluid_field_modification_status /= 1) grid_metadata_status = 0
    if (no_rhs_injection_status /= 1) grid_metadata_status = 0
    if (no_ibm_spreading_status /= 1) grid_metadata_status = 0
    if (no_feedback_force_status /= 1) grid_metadata_status = 0
    if (no_twoway_force_status /= 1) grid_metadata_status = 0
    if (no_structure_advance_status /= 1) grid_metadata_status = 0
  end subroutine

  logical function finite_value(value)
    real(mytype), intent(in) :: value
    finite_value = (value == value) .and. (abs(value) < huge(value))
  end function

end module fibre_stage11_grid_metadata
