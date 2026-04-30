module fibre_diagnostics

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_geometry, only : fibre_segment_length, fibre_second_derivative_internal
  use fibre_geometry, only : fibre_end_to_end_distance, fibre_center_of_mass

  implicit none
  private

  type, public :: fibre_diagnostics_t
    real(mytype) :: max_segment_length_error = 0._mytype
    real(mytype) :: min_segment_length = 0._mytype
    real(mytype) :: max_segment_length = 0._mytype
    real(mytype) :: end_to_end_distance = 0._mytype
    real(mytype) :: center_of_mass(3) = 0._mytype
    real(mytype) :: max_curvature = 0._mytype
  end type fibre_diagnostics_t

  public :: fibre_compute_diagnostics
  public :: fibre_write_diagnostics

  public :: compute_total_discrete_length
  public :: compute_rms_segment_length_error
  public :: compute_total_length_relative_error
  public :: compute_curvature_diagnostics
  public :: compute_bending_energy
  public :: compute_kinetic_energy
  public :: compute_total_structure_energy
  public :: compute_center_of_mass_velocity
  public :: compute_linear_momentum
  public :: compute_end_to_end_metrics
  public :: compute_kinetic_energy_split
  public :: compute_centroid_angular_momentum

contains

  pure function trapezoid_weight(i, nl, ds) result(w)
    integer, intent(in) :: i, nl
    real(mytype), intent(in) :: ds
    real(mytype) :: w

    if (i == 1 .or. i == nl) then
      w = 0.5_mytype * ds
    else
      w = ds
    end if
  end function trapezoid_weight

  pure function cross_product(a, b) result(c)
    real(mytype), intent(in) :: a(3), b(3)
    real(mytype) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  subroutine compute_total_discrete_length(fibre, total_length)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: total_length

    real(mytype), allocatable :: seg_length(:)

    allocate(seg_length(fibre%nl - 1))
    call fibre_segment_length(fibre, seg_length)
    total_length = sum(seg_length)
    deallocate(seg_length)
  end subroutine compute_total_discrete_length

  subroutine compute_rms_segment_length_error(fibre, rms_error)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: rms_error

    real(mytype), allocatable :: seg_length(:)

    allocate(seg_length(fibre%nl - 1))
    call fibre_segment_length(fibre, seg_length)
    rms_error = sqrt(sum((seg_length / fibre%ds - 1._mytype)**2) / real(fibre%nl - 1, mytype))
    deallocate(seg_length)
  end subroutine compute_rms_segment_length_error

  subroutine compute_total_length_relative_error(fibre, rel_error)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: rel_error

    real(mytype) :: total_length

    call compute_total_discrete_length(fibre, total_length)
    rel_error = abs(total_length / fibre%length - 1._mytype)
  end subroutine compute_total_length_relative_error

  subroutine compute_curvature_diagnostics(fibre, max_curvature, rms_curvature)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: max_curvature
    real(mytype), intent(out) :: rms_curvature

    real(mytype), allocatable :: x_ss(:,:)
    real(mytype) :: kappa, weighted_sum
    integer :: i

    allocate(x_ss(3, fibre%nl))
    call fibre_second_derivative_internal(fibre, x_ss)

    max_curvature = 0._mytype
    weighted_sum = 0._mytype

    do i = 1, fibre%nl
      kappa = sqrt(sum(x_ss(:, i)**2))
      max_curvature = max(max_curvature, kappa)
      weighted_sum = weighted_sum + trapezoid_weight(i, fibre%nl, fibre%ds) * kappa**2
    end do

    rms_curvature = sqrt(weighted_sum / fibre%length)

    deallocate(x_ss)
  end subroutine compute_curvature_diagnostics

  subroutine compute_bending_energy(fibre, bending_energy)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: bending_energy

    real(mytype), allocatable :: x_ss(:,:)
    integer :: i

    allocate(x_ss(3, fibre%nl))
    call fibre_second_derivative_internal(fibre, x_ss)

    bending_energy = 0._mytype
    do i = 1, fibre%nl
      bending_energy = bending_energy + trapezoid_weight(i, fibre%nl, fibre%ds) * sum(x_ss(:, i)**2)
    end do
    bending_energy = 0.5_mytype * fibre%gamma * bending_energy

    deallocate(x_ss)
  end subroutine compute_bending_energy

  subroutine compute_kinetic_energy(fibre, kinetic_energy)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: kinetic_energy

    integer :: i

    kinetic_energy = 0._mytype
    do i = 1, fibre%nl
      kinetic_energy = kinetic_energy + trapezoid_weight(i, fibre%nl, fibre%ds) * sum(fibre%v(:, i)**2)
    end do
    kinetic_energy = 0.5_mytype * fibre%rho_tilde * kinetic_energy
  end subroutine compute_kinetic_energy

  subroutine compute_total_structure_energy(fibre, total_energy)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: total_energy

    real(mytype) :: bending_energy, kinetic_energy

    call compute_bending_energy(fibre, bending_energy)
    call compute_kinetic_energy(fibre, kinetic_energy)
    total_energy = bending_energy + kinetic_energy
  end subroutine compute_total_structure_energy

  subroutine compute_center_of_mass_velocity(fibre, vc)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: vc(3)

    integer :: i

    vc = 0._mytype
    do i = 1, fibre%nl
      vc = vc + trapezoid_weight(i, fibre%nl, fibre%ds) * fibre%v(:, i)
    end do
    vc = vc / fibre%length
  end subroutine compute_center_of_mass_velocity

  subroutine compute_linear_momentum(fibre, p)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: p(3)

    integer :: i

    p = 0._mytype
    do i = 1, fibre%nl
      p = p + trapezoid_weight(i, fibre%nl, fibre%ds) * fibre%v(:, i)
    end do
    p = fibre%rho_tilde * p
  end subroutine compute_linear_momentum

  subroutine compute_end_to_end_metrics(fibre, r_vec, r_norm, stretch_ratio)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: r_vec(3)
    real(mytype), intent(out) :: r_norm
    real(mytype), intent(out) :: stretch_ratio

    r_vec = fibre%x(:, fibre%nl) - fibre%x(:, 1)
    r_norm = sqrt(sum(r_vec**2))
    stretch_ratio = r_norm / fibre%length
  end subroutine compute_end_to_end_metrics

  subroutine compute_kinetic_energy_split(fibre, ek_trans, ek_rel)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: ek_trans
    real(mytype), intent(out) :: ek_rel

    real(mytype) :: vc(3), dv(3)
    integer :: i

    call compute_center_of_mass_velocity(fibre, vc)

    ek_trans = 0.5_mytype * fibre%rho_tilde * fibre%length * sum(vc**2)

    ek_rel = 0._mytype
    do i = 1, fibre%nl
      dv = fibre%v(:, i) - vc
      ek_rel = ek_rel + trapezoid_weight(i, fibre%nl, fibre%ds) * sum(dv**2)
    end do
    ek_rel = 0.5_mytype * fibre%rho_tilde * ek_rel
  end subroutine compute_kinetic_energy_split

  subroutine compute_centroid_angular_momentum(fibre, lc)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: lc(3)

    real(mytype) :: x_c(3)
    integer :: i

    x_c = fibre_center_of_mass(fibre)
    lc = 0._mytype

    do i = 1, fibre%nl
      lc = lc + trapezoid_weight(i, fibre%nl, fibre%ds) * cross_product(fibre%x(:, i) - x_c, fibre%v(:, i))
    end do
    lc = fibre%rho_tilde * lc
  end subroutine compute_centroid_angular_momentum

  subroutine fibre_compute_diagnostics(fibre, diag)
    type(fibre_t), intent(in) :: fibre
    type(fibre_diagnostics_t), intent(out) :: diag

    real(mytype), allocatable :: seg_length(:)
    real(mytype) :: rms_curvature

    allocate(seg_length(fibre%nl - 1))
    call fibre_segment_length(fibre, seg_length)

    diag%max_segment_length_error = maxval(abs(seg_length / fibre%ds - 1._mytype))
    diag%min_segment_length = minval(seg_length)
    diag%max_segment_length = maxval(seg_length)
    diag%end_to_end_distance = fibre_end_to_end_distance(fibre)
    diag%center_of_mass = fibre_center_of_mass(fibre)
    call compute_curvature_diagnostics(fibre, diag%max_curvature, rms_curvature)

    deallocate(seg_length)
  end subroutine fibre_compute_diagnostics

  subroutine fibre_write_diagnostics(filename, fibre, diag)
    character(len=*), intent(in) :: filename
    type(fibre_t), intent(in) :: fibre
    type(fibre_diagnostics_t), intent(in) :: diag

    integer :: io_unit

    open(newunit=io_unit, file=filename, status='replace', action='write', form='formatted')
    write(io_unit, '(A)') '# free-free fibre geometry diagnostics'
    write(io_unit, '(A,I0)') 'nl = ', fibre%nl
    write(io_unit, '(A,ES24.16)') 'length = ', fibre%length
    write(io_unit, '(A,ES24.16)') 'ds = ', fibre%ds
    write(io_unit, '(A,ES24.16)') 'max_segment_length_error = ', diag%max_segment_length_error
    write(io_unit, '(A,ES24.16)') 'min_segment_length = ', diag%min_segment_length
    write(io_unit, '(A,ES24.16)') 'max_segment_length = ', diag%max_segment_length
    write(io_unit, '(A,ES24.16)') 'end_to_end_distance = ', diag%end_to_end_distance
    write(io_unit, '(A,3(ES24.16,1X))') 'center_of_mass = ', diag%center_of_mass
    write(io_unit, '(A,ES24.16)') 'max_curvature = ', diag%max_curvature
    close(io_unit)
  end subroutine fibre_write_diagnostics

end module fibre_diagnostics
