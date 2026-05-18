program fibre_stage5_closed_loop_diagnostics_audit
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_spreading, only : clear_eulerian_force_density, spread_lag_force_to_eulerian, compute_eulerian_total_force, compute_lagrangian_total_force
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_ibm_power_diagnostics, only : compute_eulerian_power, compute_lagrangian_power, compute_power_consistency_error
  implicit none

  integer, parameter :: nx=8, ny=6, nz=5, nlag=4
  real(mytype), parameter :: tiny_val=1.0e-30_mytype

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag
  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
  real(mytype), allocatable :: u_lag(:,:)
  real(mytype) :: force_abs_error, force_rel_error, power_abs_error, power_rel_error
  real(mytype) :: power_e, power_l
  real(mytype) :: fe_total(3), fl_total(3)
  integer :: i, j, k, io
  integer :: eulerian_force_computed_flag, lagrangian_force_computed_flag
  integer :: eulerian_power_computed_flag, lagrangian_power_computed_flag
  integer :: used_spreading_buffer_flag, used_lagrangian_force_flag, used_interpolated_velocity_flag
  integer :: real_spreading_module_called_flag, real_interpolation_module_called_flag
  integer :: no_tautological_force_flag, no_tautological_power_flag, audit_status

  call init_synthetic_case(grid, lag, ux, uy, uz, fx, fy, fz, u_lag)

  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  real_spreading_module_called_flag = 1

  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  real_interpolation_module_called_flag = 1

  call compute_eulerian_total_force(grid, fx, fy, fz, fe_total)
  call compute_lagrangian_total_force(lag, fl_total)
  call compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, power_e)
  call compute_lagrangian_power(lag, u_lag, power_l)
  call compute_power_consistency_error(power_e, power_l, power_abs_error, power_rel_error)

  force_abs_error = sqrt(sum((fe_total-fl_total)**2))
  force_rel_error = force_abs_error / max(sqrt(sum(fl_total**2)), tiny_val)

  eulerian_force_computed_flag = merge(1,0,sum(abs(fe_total)) > 0._mytype)
  lagrangian_force_computed_flag = merge(1,0,sum(abs(fl_total)) > 0._mytype)
  eulerian_power_computed_flag = merge(1,0,abs(power_e) > 0._mytype)
  lagrangian_power_computed_flag = merge(1,0,abs(power_l) > 0._mytype)
  used_spreading_buffer_flag = merge(1,0,sum(abs(fx)+abs(fy)+abs(fz)) > 0._mytype)
  used_lagrangian_force_flag = merge(1,0,sum(abs(lag%force)) > 0._mytype)
  used_interpolated_velocity_flag = merge(1,0,sum(abs(u_lag)) > 0._mytype)

  no_tautological_force_flag = merge(1,0, eulerian_force_computed_flag==1 .and. lagrangian_force_computed_flag==1 .and. &
       used_spreading_buffer_flag==1 .and. used_lagrangian_force_flag==1 .and. &
       real_spreading_module_called_flag==1 .and. real_interpolation_module_called_flag==1)
  no_tautological_power_flag = merge(1,0, eulerian_power_computed_flag==1 .and. lagrangian_power_computed_flag==1 .and. &
       used_interpolated_velocity_flag==1 .and. used_lagrangian_force_flag==1 .and. &
       real_spreading_module_called_flag==1 .and. real_interpolation_module_called_flag==1)

  audit_status = merge(1,0, no_tautological_force_flag==1 .and. no_tautological_power_flag==1 .and. &
       force_abs_error<=1e-12_mytype .and. force_rel_error<=1e-12_mytype .and. &
       power_abs_error<=1e-12_mytype .and. power_rel_error<=1e-12_mytype)

  open(newunit=io,file='stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage5_closed_loop_eulerian_force_computed_flag', eulerian_force_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_lagrangian_force_computed_flag', lagrangian_force_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_eulerian_power_computed_flag', eulerian_power_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_lagrangian_power_computed_flag', lagrangian_power_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_used_spreading_buffer_flag', used_spreading_buffer_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_used_lagrangian_force_flag', used_lagrangian_force_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_used_interpolated_velocity_flag', used_interpolated_velocity_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_real_spreading_module_called_flag', real_spreading_module_called_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_real_interpolation_module_called_flag', real_interpolation_module_called_flag
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_force_conservation_abs_error', force_abs_error
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_force_conservation_relative_error', force_rel_error
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_power_abs_error', power_abs_error
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_power_relative_error', power_rel_error
  write(io,'(A,1X,I0)') 'stage5_closed_loop_no_tautological_force_flag', no_tautological_force_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_no_tautological_power_flag', no_tautological_power_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_diagnostics_audit_status', audit_status
  close(io)

contains
  subroutine init_synthetic_case(grid, lag, ux, uy, uz, fx, fy, fz, u_lag)
    type(ibm_grid_t), intent(out) :: grid
    type(ibm_lagrangian_points_t), intent(out) :: lag
    real(mytype), allocatable, intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(mytype), allocatable, intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), allocatable, intent(out) :: u_lag(:,:)
    integer :: ii, jj, kk

    grid%nx = nx; grid%ny = ny; grid%nz = nz
    grid%xmin = 0._mytype; grid%xmax = 1._mytype
    grid%ymin = 0._mytype; grid%ymax = 1._mytype
    grid%zmin = 0._mytype; grid%zmax = 1._mytype
    grid%dx = (grid%xmax-grid%xmin)/real(nx,mytype)
    grid%dy = (grid%ymax-grid%ymin)/real(ny,mytype)
    grid%dz = (grid%zmax-grid%zmin)/real(nz,mytype)
    grid%cell_volume = grid%dx*grid%dy*grid%dz
    grid%periodic_x = .false.; grid%periodic_y = .false.; grid%periodic_z = .false.
    allocate(grid%x(nx), grid%y(ny), grid%z(nz))
    do ii=1,nx; grid%x(ii)=grid%xmin + (real(ii,mytype)-0.5_mytype)*grid%dx; end do
    do jj=1,ny; grid%y(jj)=grid%ymin + (real(jj,mytype)-0.5_mytype)*grid%dy; end do
    do kk=1,nz; grid%z(kk)=grid%zmin + (real(kk,mytype)-0.5_mytype)*grid%dz; end do

    allocate(ux(nx,ny,nz), uy(nx,ny,nz), uz(nx,ny,nz), fx(nx,ny,nz), fy(nx,ny,nz), fz(nx,ny,nz), u_lag(3,nlag))
    do kk=1,nz
      do jj=1,ny
        do ii=1,nx
          ux(ii,jj,kk) = 0.1_mytype + 0.01_mytype*real(ii,mytype) + 0.02_mytype*real(jj,mytype)
          uy(ii,jj,kk) = -0.2_mytype + 0.03_mytype*real(jj,mytype) + 0.01_mytype*real(kk,mytype)
          uz(ii,jj,kk) = 0.05_mytype + 0.01_mytype*real(ii,mytype) - 0.02_mytype*real(kk,mytype)
        end do
      end do
    end do

    lag%nl = nlag
    allocate(lag%x(3,nlag), lag%v(3,nlag), lag%force(3,nlag), lag%weight(nlag))
    lag%x(:,1) = (/0.35_mytype, 0.30_mytype, 0.32_mytype/)
    lag%x(:,2) = (/0.45_mytype, 0.50_mytype, 0.41_mytype/)
    lag%x(:,3) = (/0.62_mytype, 0.46_mytype, 0.57_mytype/)
    lag%x(:,4) = (/0.71_mytype, 0.64_mytype, 0.69_mytype/)
    lag%force(:,1) = (/ 0.8_mytype, -0.3_mytype,  0.5_mytype/)
    lag%force(:,2) = (/-0.6_mytype,  0.4_mytype, -0.2_mytype/)
    lag%force(:,3) = (/ 0.3_mytype,  0.7_mytype, -0.1_mytype/)
    lag%force(:,4) = (/-0.2_mytype, -0.5_mytype,  0.9_mytype/)
    lag%weight = (/0.12_mytype, 0.10_mytype, 0.09_mytype, 0.11_mytype/)
    lag%v = 0._mytype
  end subroutine init_synthetic_case

end program fibre_stage5_closed_loop_diagnostics_audit
