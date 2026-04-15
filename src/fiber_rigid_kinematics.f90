!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_rigid_kinematics

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use param, only : dt, ifirst, ilast, xlx, zlz, gdt
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fiber_types, only : fiber_active, fiber_nl, rigid_kinematics_test_active, rigid_kinematics_one_way, &
       rigid_kinematics_mode, rigid_kinematics_output_interval, fiber_x, fiber_xc, fiber_p, fiber_omega, fiber_xdot, &
       fiber_uinterp, fiber_s_ref, fiber_length
  use fiber_interp, only : run_fiber_interp_solver_readonly

  implicit none

  logical, save :: kinematics_initialized = .false.
  real(mytype), save :: fiber_length_ref = 0._mytype
  real(mytype), save :: fiber_p0(3) = 0._mytype
  integer, save :: mode_case_id = 1
  integer, save :: last_output_itime = -1

contains

  pure real(mytype) function wrap_periodic(value, period_length)
    real(mytype), intent(in) :: value, period_length
    wrap_periodic = value
    if (period_length > 0._mytype) then
      wrap_periodic = modulo(value, period_length)
      if (wrap_periodic < 0._mytype) wrap_periodic = wrap_periodic + period_length
    endif
  end function wrap_periodic

  pure function cross_product(a, b) result(c)
    real(mytype), intent(in) :: a(3), b(3)
    real(mytype) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  integer pure function infer_case_id(mode_value, p0) result(case_id)
    integer, intent(in) :: mode_value
    real(mytype), intent(in) :: p0(3)
    real(mytype) :: angle_deg

    case_id = mode_value
    if (mode_value == 2) then
      angle_deg = acos(max(-1._mytype, min(1._mytype, p0(1)))) * 180._mytype / acos(-1._mytype)
      if (abs(angle_deg - 15._mytype) < 5._mytype) case_id = 15
      if (abs(angle_deg - 45._mytype) < 5._mytype) case_id = 45
      if (abs(angle_deg - 75._mytype) < 5._mytype) case_id = 75
    endif
  end function infer_case_id

  subroutine fill_case_tag(case_id, tag)
    integer, intent(in) :: case_id
    character(len=*), intent(out) :: tag

    select case (case_id)
    case (1)
      tag = 'case_shear'
    case (15)
      tag = 'case_poiseuille_angle15'
    case (45)
      tag = 'case_poiseuille_angle45'
    case (75)
      tag = 'case_poiseuille_angle75'
    case default
      if (nrank == 0) then
        write(*,*) 'Error: rigid_kinematics_mode must map to shear(1) or poiseuille angles.'
        write(*,*) 'Allowed values: 15/45/75 or mode=2 with those initial angles.'
      endif
      stop
    end select
  end subroutine fill_case_tag

  subroutine rigid_kinematics_init()
    real(mytype) :: pnorm
    character(len=64) :: tag_check

    if (.not.fiber_active) return
    if (.not.rigid_kinematics_test_active) return
    if (.not.rigid_kinematics_one_way) then
      if (nrank == 0) write(*,*) 'Error: rigid_kinematics_test_active requires rigid_kinematics_one_way = true.'
      stop
    endif

    pnorm = sqrt(sum(fiber_p**2))
    if (pnorm <= 0._mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_p must be non-zero for rigid kinematics test.'
      stop
    endif
    fiber_p = fiber_p / pnorm
    fiber_p0 = fiber_p
    mode_case_id = infer_case_id(rigid_kinematics_mode, fiber_p0)
    call fill_case_tag(mode_case_id, tag_check)

    call rigid_kinematics_reconstruct_points()
    if (.not.allocated(fiber_xdot)) allocate(fiber_xdot(3, fiber_nl))
    fiber_xdot = 0._mytype
    fiber_omega = 0._mytype

    fiber_length_ref = sqrt(sum((fiber_x(:,fiber_nl) - fiber_x(:,1))**2))
    kinematics_initialized = .true.
    last_output_itime = -1
  end subroutine rigid_kinematics_init

  subroutine rigid_kinematics_reconstruct_points()
    integer :: l

    do l = 1, fiber_nl
      fiber_x(1,l) = wrap_periodic(fiber_xc(1) + fiber_s_ref(l) * fiber_p(1), xlx)
      fiber_x(2,l) = fiber_xc(2) + fiber_s_ref(l) * fiber_p(2)
      fiber_x(3,l) = wrap_periodic(fiber_xc(3) + fiber_s_ref(l) * fiber_p(3), zlz)
    enddo
  end subroutine rigid_kinematics_reconstruct_points

  subroutine rigid_kinematics_step(uxe, uye, uze, time, itime, isubstep)
    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    real(mytype), intent(in) :: time
    integer, intent(in) :: itime, isubstep

    real(mytype) :: dt_stage, ucm(3), grad_along(3), pnorm
    real(mytype) :: length_error, p_norm_error

    if (.not.fiber_active) return
    if (.not.rigid_kinematics_test_active) return
    if (.not.kinematics_initialized) call rigid_kinematics_init()

    call rigid_kinematics_reconstruct_points()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    dt_stage = gdt(isubstep)
    ucm = sum(fiber_uinterp, dim=2) / real(fiber_nl, mytype)
    fiber_xc = fiber_xc + dt_stage * ucm
    fiber_xc(1) = wrap_periodic(fiber_xc(1), xlx)
    fiber_xc(3) = wrap_periodic(fiber_xc(3), zlz)

    grad_along = (fiber_uinterp(:,fiber_nl) - fiber_uinterp(:,1)) / max(fiber_length, 1.0e-12_mytype)
    fiber_omega = cross_product(fiber_p, grad_along)
    fiber_p = fiber_p + dt_stage * cross_product(fiber_omega, fiber_p)
    pnorm = sqrt(sum(fiber_p**2))
    if (pnorm > 0._mytype) fiber_p = fiber_p / pnorm

    call rigid_kinematics_reconstruct_points()
    call rigid_kinematics_diagnostics(length_error, p_norm_error)
    call rigid_kinematics_update_point_velocity(ucm)
    call write_kinematics_outputs(itime, time, length_error, p_norm_error)
    if (isubstep == size(gdt) .and. itime == ilast) then
      call write_dt_sensitivity_summary(time, length_error, p_norm_error)
    endif
  end subroutine rigid_kinematics_step

  subroutine rigid_kinematics_update_point_velocity(ucm)
    real(mytype), intent(in) :: ucm(3)
    integer :: l
    real(mytype) :: r(3)

    if (.not.allocated(fiber_xdot)) allocate(fiber_xdot(3, fiber_nl))
    do l = 1, fiber_nl
      r = fiber_x(:,l) - fiber_xc
      fiber_xdot(:,l) = ucm + cross_product(fiber_omega, r)
    enddo
  end subroutine rigid_kinematics_update_point_velocity

  subroutine rigid_kinematics_diagnostics(length_error, p_norm_error)
    real(mytype), intent(out) :: length_error, p_norm_error
    real(mytype) :: length_now

    length_now = sqrt(sum((fiber_x(:,fiber_nl) - fiber_x(:,1))**2))
    length_error = abs(length_now - fiber_length_ref)
    p_norm_error = abs(sqrt(sum(fiber_p**2)) - 1._mytype)
  end subroutine rigid_kinematics_diagnostics

  subroutine write_kinematics_outputs(itime, time, length_error, p_norm_error)
    integer, intent(in) :: itime
    real(mytype), intent(in) :: time, length_error, p_norm_error
    integer :: ifile_traj, ifile_orient
    logical :: exists_traj, exists_orient
    character(len=256) :: file_traj, file_orient
    character(len=64) :: tag

    if (nrank /= 0) return
    if (last_output_itime == itime) return
    if (.not.(itime == ifirst .or. mod(itime, max(1, rigid_kinematics_output_interval)) == 0)) return
    if (.not.all(ieee_is_finite((/time, fiber_xc(1), fiber_xc(2), fiber_xc(3), fiber_p(1), fiber_p(2), fiber_p(3), &
         fiber_omega(1), fiber_omega(2), fiber_omega(3), length_error, p_norm_error/)))) return

    call fill_case_tag(mode_case_id, tag)

    file_traj = 'xcm_trajectory_' // trim(tag) // '.dat'
    file_orient = 'orientation_series_' // trim(tag) // '.dat'

    inquire(file=trim(file_traj), exist=exists_traj)
    if (exists_traj) then
      open(newunit=ifile_traj, file=trim(file_traj), status='old', action='write', position='append')
    else
      open(newunit=ifile_traj, file=trim(file_traj), status='replace', action='write')
      write(ifile_traj,'(A)') 'itime time xcm_x xcm_y xcm_z p_x p_y p_z omega_x omega_y omega_z length_error p_norm_error'
    endif
    write(ifile_traj,'(I10,1X,12(ES24.16,1X))') itime, time, fiber_xc(1), fiber_xc(2), fiber_xc(3), &
         fiber_p(1), fiber_p(2), fiber_p(3), fiber_omega(1), fiber_omega(2), fiber_omega(3), length_error, p_norm_error
    close(ifile_traj)

    inquire(file=trim(file_orient), exist=exists_orient)
    if (exists_orient) then
      open(newunit=ifile_orient, file=trim(file_orient), status='old', action='write', position='append')
    else
      open(newunit=ifile_orient, file=trim(file_orient), status='replace', action='write')
      write(ifile_orient,'(A)') 'itime time p_x p_y p_z omega_x omega_y omega_z length_error p_norm_error'
    endif
    write(ifile_orient,'(I10,1X,10(ES24.16,1X))') itime, time, fiber_p(1), fiber_p(2), fiber_p(3), &
         fiber_omega(1), fiber_omega(2), fiber_omega(3), length_error, p_norm_error
    close(ifile_orient)

    last_output_itime = itime
  end subroutine write_kinematics_outputs

  subroutine write_dt_sensitivity_summary(time, length_error, p_norm_error)
    real(mytype), intent(in) :: time, length_error, p_norm_error
    integer :: ifile_dt
    logical :: exists_dt
    character(len=64) :: tag
    character(len=256) :: file_dt

    if (nrank /= 0) return

    call fill_case_tag(mode_case_id, tag)
    file_dt = 'dt_sensitivity_summary.dat'

    inquire(file=trim(file_dt), exist=exists_dt)
    if (exists_dt) then
      open(newunit=ifile_dt, file=trim(file_dt), status='old', action='write', position='append')
    else
      open(newunit=ifile_dt, file=trim(file_dt), status='replace', action='write')
      write(ifile_dt,'(A)') &
           'case_tag dt time_final xcm_x xcm_y xcm_z p_x p_y p_z omega_x omega_y omega_z length_error p_norm_error'
    endif
    write(ifile_dt,'(A,1X,13(ES24.16,1X))') trim(tag), dt, time, fiber_xc(1), fiber_xc(2), fiber_xc(3), &
         fiber_p(1), fiber_p(2), fiber_p(3), fiber_omega(1), fiber_omega(2), fiber_omega(3), length_error, p_norm_error
    close(ifile_dt)
  end subroutine write_dt_sensitivity_summary

end module fiber_rigid_kinematics
