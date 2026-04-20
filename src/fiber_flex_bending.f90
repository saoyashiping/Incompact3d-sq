!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_bending

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, fiber_ds, &
       fiber_s_ref, fiber_x, fiber_flex_bending_test_active, fiber_flex_bending_case, fiber_flex_bending_nsteps, &
       fiber_flex_bending_output_interval, fiber_flex_bending_dt, fiber_bending_gamma
  use fiber_flex_ops, only : apply_free_end_d2_scalar, apply_free_end_bending_operator_scalar
  use fiber_io, only : write_fiber_flex_bending_initial_points, write_fiber_flex_bending_final_points, &
       write_fiber_flex_bending_series, write_fiber_flex_bending_summary

  implicit none

contains

  subroutine assemble_free_end_bending_matrix(bmat, gamma_b)
    real(mytype), intent(out) :: bmat(fiber_nl, fiber_nl)
    real(mytype), intent(in) :: gamma_b
    real(mytype) :: evec(fiber_nl), bcol(fiber_nl)
    integer :: j

    bmat = 0._mytype
    do j = 1, fiber_nl
      evec = 0._mytype
      evec(j) = 1._mytype
      call apply_free_end_bending_operator_scalar(evec, bcol, gamma_b)
      bmat(:,j) = bcol
    enddo
  end subroutine assemble_free_end_bending_matrix

  subroutine build_bending_backward_euler_matrix(dt_b, gamma_b, amat)
    ! Step 3.3 semantics:
    !   x^{n+1} - x^n = -dt_b * B(x^{n+1})
    ! => (I + dt_b * B) x^{n+1} = x^n
    ! Here B = -gamma*d4 (constant gamma).
    real(mytype), intent(in) :: dt_b, gamma_b
    real(mytype), intent(out) :: amat(fiber_nl, fiber_nl)
    real(mytype) :: bmat(fiber_nl, fiber_nl)
    integer :: i

    call assemble_free_end_bending_matrix(bmat, gamma_b)
    amat = dt_b * bmat
    do i = 1, fiber_nl
      amat(i,i) = amat(i,i) + 1._mytype
    enddo
  end subroutine build_bending_backward_euler_matrix

  subroutine solve_linear_system_dense(ain, b, x)
    real(mytype), intent(in) :: ain(fiber_nl, fiber_nl), b(fiber_nl)
    real(mytype), intent(out) :: x(fiber_nl)
    real(mytype) :: a(fiber_nl, fiber_nl), rhs(fiber_nl), rowtmp(fiber_nl), piv, fac, tmp
    integer :: i, j, k, p

    a = ain
    rhs = b
    x = 0._mytype

    do k = 1, fiber_nl
      p = k
      piv = abs(a(k,k))
      do i = k+1, fiber_nl
        if (abs(a(i,k)) > piv) then
          piv = abs(a(i,k))
          p = i
        endif
      enddo
      if (piv <= 1.0e-14_mytype) then
        if (nrank == 0) write(*,*) 'Error: singular matrix in solve_linear_system_dense.'
        stop
      endif
      if (p /= k) then
        rowtmp = a(k,:); a(k,:) = a(p,:); a(p,:) = rowtmp
        tmp = rhs(k); rhs(k) = rhs(p); rhs(p) = tmp
      endif
      do i = k+1, fiber_nl
        fac = a(i,k) / a(k,k)
        a(i,k:fiber_nl) = a(i,k:fiber_nl) - fac * a(k,k:fiber_nl)
        rhs(i) = rhs(i) - fac * rhs(k)
      enddo
    enddo

    do i = fiber_nl, 1, -1
      if (i < fiber_nl) then
        x(i) = (rhs(i) - sum(a(i,i+1:fiber_nl) * x(i+1:fiber_nl))) / a(i,i)
      else
        x(i) = rhs(i) / a(i,i)
      endif
    enddo
  end subroutine solve_linear_system_dense

  subroutine advance_bending_backward_euler(x_old, x_new, dt_b, gamma_b)
    real(mytype), intent(in) :: x_old(3, fiber_nl), dt_b, gamma_b
    real(mytype), intent(out) :: x_new(3, fiber_nl)
    real(mytype) :: amat(fiber_nl, fiber_nl)
    integer :: c

    call build_bending_backward_euler_matrix(dt_b, gamma_b, amat)
    do c = 1, 3
      call solve_linear_system_dense(amat, x_old(c,:), x_new(c,:))
    enddo
  end subroutine advance_bending_backward_euler

  subroutine compute_bending_energy(x_in, gamma_b, ebend)
    real(mytype), intent(in) :: x_in(3, fiber_nl), gamma_b
    real(mytype), intent(out) :: ebend
    real(mytype) :: d2(3, fiber_nl), norm2(fiber_nl)
    integer :: i

    do i = 1, 3
      call apply_free_end_d2_scalar(x_in(i,:), d2(i,:))
    enddo
    norm2 = d2(1,:)**2 + d2(2,:)**2 + d2(3,:)**2
    ebend = 0.5_mytype * gamma_b * fiber_ds * (0.5_mytype * norm2(1) + sum(norm2(2:fiber_nl-1)) + 0.5_mytype * norm2(fiber_nl))
  end subroutine compute_bending_energy

  subroutine initialize_bending_test_shape(case_id, x0)
    integer, intent(in) :: case_id
    real(mytype), intent(out) :: x0(3, fiber_nl)
    real(mytype) :: smin, q, g, amp
    integer :: l

    select case (case_id)
    case (1)
      x0 = fiber_x
    case (2)
      smin = -0.5_mytype * fiber_length
      amp = 0.1_mytype * fiber_length
      do l = 1, fiber_nl
        q = (fiber_s_ref(l) - smin) / fiber_length
        g = q**4 * (1._mytype - q)**4
        x0(1,l) = fiber_s_ref(l)
        x0(2,l) = amp * g
        x0(3,l) = 0._mytype
      enddo
    case default
      if (nrank == 0) write(*,*) 'Error: fiber_flex_bending_case must be 1 or 2.'
      stop
    end select
  end subroutine initialize_bending_test_shape

  subroutine run_flexible_bending_test()
    real(mytype), allocatable :: xold(:,:), xnew(:,:), xinit(:,:), energy(:), maxupd(:)
    real(mytype) :: e0, ecur, tnow, max_update, straight_err, final_disp
    logical :: energy_monotone
    integer :: step, outint

    if (.not.fiber_active .or. .not.fiber_flexible_active .or. .not.fiber_flex_initialized .or. &
         .not.fiber_flex_bending_test_active) then
      if (nrank == 0) write(*,*) 'Error: invalid state for run_flexible_bending_test.'
      stop
    endif

    allocate(xold(3,fiber_nl), xnew(3,fiber_nl), xinit(3,fiber_nl), &
         energy(0:fiber_flex_bending_nsteps), maxupd(0:fiber_flex_bending_nsteps))

    call initialize_bending_test_shape(fiber_flex_bending_case, xold)
    xinit = xold
    call write_fiber_flex_bending_initial_points(xold)
    call compute_bending_energy(xold, fiber_bending_gamma, e0)
    energy(0) = e0
    maxupd(0) = 0._mytype

    outint = max(1, fiber_flex_bending_output_interval)
    call write_fiber_flex_bending_series(0, 0._mytype, energy(0), maxupd(0), .true.)
    do step = 1, fiber_flex_bending_nsteps
      call advance_bending_backward_euler(xold, xnew, fiber_flex_bending_dt, fiber_bending_gamma)
      max_update = maxval(abs(xnew - xold))
      call compute_bending_energy(xnew, fiber_bending_gamma, ecur)
      energy(step) = ecur
      maxupd(step) = max_update
      tnow = real(step, mytype) * fiber_flex_bending_dt
      if (mod(step, outint) == 0 .or. step == fiber_flex_bending_nsteps) then
        call write_fiber_flex_bending_series(step, tnow, ecur, max_update, .false.)
      endif
      xold = xnew
    enddo

    call write_fiber_flex_bending_final_points(xold)

    energy_monotone = .true.
    do step = 1, fiber_flex_bending_nsteps
      if (energy(step) > energy(step-1) + 1.0e-12_mytype) energy_monotone = .false.
    enddo

    straight_err = 0._mytype
    final_disp = 0._mytype
    if (fiber_flex_bending_case == 1) then
      straight_err = maxval(abs(xold - xinit))
    else if (fiber_flex_bending_case == 2) then
      final_disp = sqrt(sum((xold - xinit)**2) / real(3*fiber_nl, mytype))
    endif

    call write_fiber_flex_bending_summary(fiber_flex_bending_case, fiber_flex_bending_dt, fiber_bending_gamma, &
         fiber_flex_bending_nsteps, energy(0), energy(fiber_flex_bending_nsteps), energy_monotone, &
         maxupd(fiber_flex_bending_nsteps), straight_err, final_disp)

    deallocate(xold, xnew, xinit, energy, maxupd)
  end subroutine run_flexible_bending_test

end module fiber_flex_bending
