program fibre_bending_operator_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_bending_operator, only : compute_fibre_d4_freefree
  use fibre_bending_operator, only : compute_fibre_bending_force_freefree

  implicit none

  type(fibre_t) :: fibre

  real(mytype), allocatable :: d4x(:,:), fb(:,:)
  real(mytype), allocatable :: d4y_exact(:), d4y_num(:)

  integer :: i, io_unit
  integer :: skip
  real(mytype) :: pi, s, amp, k_wave, coeff
  real(mytype) :: straight_max_abs_d4, straight_max_abs_fb
  real(mytype) :: cubic_inner_max_abs_d4
  real(mytype) :: sine_inner_max_abs_error_d4, sine_inner_relative_l2_error_d4
  real(mytype) :: d4_left_endpoint_change_norm, d4_right_endpoint_change_norm
  real(mytype) :: x_left_before(3), x_right_before(3), x_left_after(3), x_right_after(3)
  real(mytype) :: sine_l2_error_nl33, sine_l2_error_nl65, sine_l2_error_nl129
  real(mytype) :: sine_error_ratio_33_to_65, sine_error_ratio_65_to_129

  pi = acos(-1.0_mytype)
  skip = 4

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  allocate(d4x(3, fibre%nl), fb(3, fibre%nl))

  x_left_before = fibre%x(:, 1)
  x_right_before = fibre%x(:, fibre%nl)

  call compute_fibre_d4_freefree(fibre, d4x)
  call compute_fibre_bending_force_freefree(fibre, fb)

  x_left_after = fibre%x(:, 1)
  x_right_after = fibre%x(:, fibre%nl)

  d4_left_endpoint_change_norm = sqrt(sum((x_left_after - x_left_before)**2))
  d4_right_endpoint_change_norm = sqrt(sum((x_right_after - x_right_before)**2))

  straight_max_abs_d4 = maxval(abs(d4x))
  straight_max_abs_fb = maxval(abs(fb))

  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = 0.15_mytype * s**3 - 0.2_mytype * s**2 + 0.05_mytype * s + 0.01_mytype
    fibre%x(3, i) = 0._mytype
  end do

  call compute_fibre_d4_freefree(fibre, d4x)
  cubic_inner_max_abs_d4 = maxval(abs(d4x(:, 4:fibre%nl - 3)))

  amp = 0.01_mytype
  k_wave = 2._mytype
  coeff = (k_wave * pi / fibre%length)**4

  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(k_wave * pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do

  call compute_fibre_d4_freefree(fibre, d4x)

  allocate(d4y_exact(fibre%nl), d4y_num(fibre%nl))
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    d4y_exact(i) = amp * coeff * sin(k_wave * pi * s / fibre%length)
    d4y_num(i) = d4x(2, i)
  end do

  sine_inner_max_abs_error_d4 = maxval(abs(d4y_num(1 + skip:fibre%nl - skip) - d4y_exact(1 + skip:fibre%nl - skip)))
  sine_inner_relative_l2_error_d4 = sqrt(sum((d4y_num(1 + skip:fibre%nl - skip) - d4y_exact(1 + skip:fibre%nl - skip))**2)) / &
                                    max(sqrt(sum(d4y_exact(1 + skip:fibre%nl - skip)**2)), tiny(1.0_mytype))

  deallocate(d4y_exact, d4y_num)
  deallocate(d4x, fb)

  call compute_sine_l2_error_for_nl(33, sine_l2_error_nl33)
  call compute_sine_l2_error_for_nl(65, sine_l2_error_nl65)
  call compute_sine_l2_error_for_nl(129, sine_l2_error_nl129)

  sine_error_ratio_33_to_65 = sine_l2_error_nl33 / max(sine_l2_error_nl65, tiny(1.0_mytype))
  sine_error_ratio_65_to_129 = sine_l2_error_nl65 / max(sine_l2_error_nl129, tiny(1.0_mytype))

  open(newunit=io_unit, file='stage2_outputs/fibre_bending_operator_check.dat', status='replace', action='write', form='formatted')
  write(io_unit, '(A,ES24.16)') 'straight_max_abs_d4 = ', straight_max_abs_d4
  write(io_unit, '(A,ES24.16)') 'straight_max_abs_fb = ', straight_max_abs_fb
  write(io_unit, '(A,ES24.16)') 'cubic_inner_max_abs_d4 = ', cubic_inner_max_abs_d4
  write(io_unit, '(A,ES24.16)') 'sine_inner_max_abs_error_d4 = ', sine_inner_max_abs_error_d4
  write(io_unit, '(A,ES24.16)') 'sine_inner_relative_l2_error_d4 = ', sine_inner_relative_l2_error_d4
  write(io_unit, '(A,ES24.16)') 'd4_left_endpoint_change_norm = ', d4_left_endpoint_change_norm
  write(io_unit, '(A,ES24.16)') 'd4_right_endpoint_change_norm = ', d4_right_endpoint_change_norm
  write(io_unit, '(A,ES24.16)') 'sine_l2_error_nl33 = ', sine_l2_error_nl33
  write(io_unit, '(A,ES24.16)') 'sine_l2_error_nl65 = ', sine_l2_error_nl65
  write(io_unit, '(A,ES24.16)') 'sine_l2_error_nl129 = ', sine_l2_error_nl129
  write(io_unit, '(A,ES24.16)') 'sine_error_ratio_33_to_65 = ', sine_error_ratio_33_to_65
  write(io_unit, '(A,ES24.16)') 'sine_error_ratio_65_to_129 = ', sine_error_ratio_65_to_129
  close(io_unit)

contains

  subroutine compute_sine_l2_error_for_nl(nl, l2_error)
    integer, intent(in) :: nl
    real(mytype), intent(out) :: l2_error

    type(fibre_t) :: fibre_local
    real(mytype), allocatable :: d4x_local(:,:), d4y_exact_local(:)
    real(mytype) :: s_local, amp_local, k_wave_local, coeff_local
    integer :: j

    call fibre_init_straight_free_free(fibre_local, nl=nl, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
    allocate(d4x_local(3, fibre_local%nl), d4y_exact_local(fibre_local%nl))

    amp_local = 0.01_mytype
    k_wave_local = 2._mytype
    coeff_local = (k_wave_local * pi / fibre_local%length)**4

    do j = 1, fibre_local%nl
      s_local = real(j - 1, mytype) * fibre_local%ds
      fibre_local%x(1, j) = s_local
      fibre_local%x(2, j) = amp_local * sin(k_wave_local * pi * s_local / fibre_local%length)
      fibre_local%x(3, j) = 0._mytype
      d4y_exact_local(j) = amp_local * coeff_local * sin(k_wave_local * pi * s_local / fibre_local%length)
    end do

    call compute_fibre_d4_freefree(fibre_local, d4x_local)

    l2_error = sqrt(sum((d4x_local(2, 1 + skip:fibre_local%nl - skip) - d4y_exact_local(1 + skip:fibre_local%nl - skip))**2)) / &
               max(sqrt(sum(d4y_exact_local(1 + skip:fibre_local%nl - skip)**2)), tiny(1.0_mytype))

    deallocate(d4x_local, d4y_exact_local)
  end subroutine compute_sine_l2_error_for_nl

end program fibre_bending_operator_check
