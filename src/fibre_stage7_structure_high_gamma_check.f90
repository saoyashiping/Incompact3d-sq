program fibre_stage7_structure_high_gamma_check
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_tension_solver, only : solve_tension_freefree
  use fibre_bending_operator, only : compute_fibre_bending_force_freefree
  use fibre_boundary_freefree, only : fibre_boundary_residual_t, fibre_compute_freefree_boundary_residual
  use fibre_diagnostics, only : compute_curvature_diagnostics, compute_total_structure_energy, compute_linear_momentum
  implicit none

  integer, parameter :: nl=33, nsteps=20
  real(mytype), parameter :: L=1._mytype, dt=1.0e-5_mytype, gamma_high=10._mytype, theta_amp=0.02_mytype
  real(mytype), parameter :: tol_len=1.0e-8_mytype, tol_stretch=1.0e-8_mytype

  type(fibre_t) :: fibre
  type(fibre_boundary_residual_t) :: boundary_residual
  real(mytype), allocatable :: fb(:,:), a_non_tension(:,:), tension(:), ds_ref(:), seg_length(:)
  real(mytype) :: tr, rr, length_error, stretch_error, current_length, ds0, s_mid, theta, tiny_len
  real(mytype) :: max_length_error, max_stretch_error, final_energy, max_curvature, rms_curvature
  real(mytype) :: momentum(3), max_velocity_norm, initial_length, final_length, initial_segment_error_max
  integer :: i, m, step, ierr, io
  integer :: nan_detected, solver_failure_count
  integer :: real_structure_solver_called_flag, real_tension_solver_called_flag
  integer :: real_bending_path_called_flag, freefree_boundary_path_called_flag
  integer :: energy_finite_flag, curvature_finite_flag, momentum_finite_flag, status

  call fibre_init_straight_free_free(fibre, nl, L, 1.0_mytype, gamma_high)
  ds0 = L / real(nl-1, mytype)

  fibre%x(:,1) = 0._mytype
  do m=1,nl-1
    s_mid = (real(m,mytype)-0.5_mytype) * ds0
    theta = theta_amp * sin(2._mytype*acos(-1._mytype)*s_mid/L)
    fibre%x(1,m+1) = fibre%x(1,m) + ds0*cos(theta)
    fibre%x(2,m+1) = fibre%x(2,m) + ds0*sin(theta)
    fibre%x(3,m+1) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype
  fibre%f_ext = 0._mytype

  allocate(fb(3,nl), a_non_tension(3,nl), tension(nl-1), ds_ref(nl-1), seg_length(nl-1))

  do i=1,nl-1
    ds_ref(i)=sqrt(sum((fibre%x(:,i+1)-fibre%x(:,i))**2))
  end do
  initial_length = sum(ds_ref)
  initial_segment_error_max = maxval(abs(ds_ref/ds0 - 1._mytype))
  tiny_len = tiny(1._mytype)

  real_structure_solver_called_flag = 0
  real_tension_solver_called_flag = 0
  real_bending_path_called_flag = 0
  freefree_boundary_path_called_flag = 0
  solver_failure_count = 0
  nan_detected = 0
  max_length_error = 0._mytype
  max_stretch_error = 0._mytype
  max_velocity_norm = 0._mytype

  do step=1,nsteps
    call compute_fibre_bending_force_freefree(fibre, fb)
    real_bending_path_called_flag = 1
    a_non_tension = (fb + fibre%f_ext) / fibre%rho_tilde

    call solve_tension_freefree(fibre, dt, a_non_tension, tension, tr, rr, ierr)
    real_tension_solver_called_flag = 1
    if (ierr /= 0) solver_failure_count = solver_failure_count + 1

    call fibre_compute_freefree_boundary_residual(fibre, boundary_residual)
    freefree_boundary_path_called_flag = 1

    call advance_fibre_structure_freefree(fibre, dt, tr, rr, ierr)
    real_structure_solver_called_flag = 1
    if (ierr /= 0) solver_failure_count = solver_failure_count + 1

    do i=1,nl-1
      seg_length(i)=sqrt(sum((fibre%x(:,i+1)-fibre%x(:,i))**2))
    end do
    stretch_error = maxval(abs(seg_length/max(ds_ref,tiny_len) - 1._mytype))
    current_length = sum(seg_length)
    length_error = abs(current_length - initial_length) / max(initial_length, tiny_len)
    max_stretch_error = max(max_stretch_error, stretch_error)
    max_length_error = max(max_length_error, length_error)
    max_velocity_norm = max(max_velocity_norm, sqrt(maxval(sum(fibre%v**2,dim=1))))

    if (.not.all(ieee_is_finite(fibre%x)) .or. .not.all(ieee_is_finite(fibre%v)) .or. .not.all(ieee_is_finite(fibre%tension))) then
      nan_detected = 1
    end if
  end do

  final_length = current_length
  call compute_total_structure_energy(fibre, final_energy)
  call compute_curvature_diagnostics(fibre, max_curvature, rms_curvature)
  call compute_linear_momentum(fibre, momentum)

  energy_finite_flag = merge(1,0,ieee_is_finite(final_energy))
  curvature_finite_flag = merge(1,0,ieee_is_finite(max_curvature) .and. ieee_is_finite(rms_curvature))
  momentum_finite_flag = merge(1,0,all(ieee_is_finite(momentum)))

  status = merge(1,0, real_structure_solver_called_flag==1 .and. real_tension_solver_called_flag==1 .and. &
      real_bending_path_called_flag==1 .and. freefree_boundary_path_called_flag==1 .and. &
      nan_detected==0 .and. solver_failure_count==0 .and. max_length_error<=tol_len .and. max_stretch_error<=tol_stretch .and. &
      energy_finite_flag==1 .and. curvature_finite_flag==1 .and. momentum_finite_flag==1 .and. initial_segment_error_max<=1.0e-14_mytype)

  open(newunit=io,file='stage7_outputs/fibre_stage7_structure_high_gamma_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage7_highgamma_real_structure_solver_called_flag', real_structure_solver_called_flag
  write(io,'(A,1X,I0)') 'stage7_highgamma_real_tension_solver_called_flag', real_tension_solver_called_flag
  write(io,'(A,1X,I0)') 'stage7_highgamma_real_bending_path_called_flag', real_bending_path_called_flag
  write(io,'(A,1X,I0)') 'stage7_highgamma_freefree_boundary_path_called_flag', freefree_boundary_path_called_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_gamma_value', gamma_high
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_dt', dt
  write(io,'(A,1X,I0)') 'stage7_highgamma_nsteps', nsteps
  write(io,'(A,1X,I0)') 'stage7_highgamma_nl', nl
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_theta_amp', theta_amp
  write(io,'(A,1X,I0)') 'stage7_highgamma_nan_detected', nan_detected
  write(io,'(A,1X,I0)') 'stage7_highgamma_solver_failure_count', solver_failure_count
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_length_error_max', max_length_error
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_stretch_error_max', max_stretch_error
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_initial_segment_error_max', initial_segment_error_max
  write(io,'(A,1X,I0)') 'stage7_highgamma_energy_finite_flag', energy_finite_flag
  write(io,'(A,1X,I0)') 'stage7_highgamma_curvature_finite_flag', curvature_finite_flag
  write(io,'(A,1X,I0)') 'stage7_highgamma_momentum_finite_flag', momentum_finite_flag
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_initial_length', initial_length
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_final_length', final_length
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_final_energy', final_energy
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_max_curvature', max_curvature
  write(io,'(A,1X,ES24.16E3)') 'stage7_highgamma_max_velocity_norm', max_velocity_norm
  write(io,'(A,1X,I0)') 'stage7_highgamma_structure_check_status', status
  close(io)

end program fibre_stage7_structure_high_gamma_check
