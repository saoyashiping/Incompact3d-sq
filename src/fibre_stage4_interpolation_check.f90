program fibre_stage4_interpolation_check

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_lagrangian_points_t, ibm_grid_t
  use fibre_ibm_grid, only : destroy_ibm_grid
  use fibre_ibm_interpolation, only : compute_interpolation_weight_sum
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays, destroy_stage4_grid_adapter
  use fibre_stage4_interpolation_adapter, only : interpolate_stage4_vector_to_lag_if_supported, convert_stage4_adapter_to_ibm_grid

  implicit none

  integer, parameter :: nx=16, ny=12, nz=10
  real(mytype), parameter :: pi = 3.1415926535897932384626433832795_mytype
  real(mytype) :: x(nx), y(ny), y_nonuniform(ny), z(nz)
  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:), wsum(:), qtrue(:)
  real(mytype), allocatable :: ux_lag(:,:), ux_lag_lin(:,:), ux_lag_poi(:,:), ux_lag_per(:,:)
  type(stage4_grid_adapter_t) :: a_uniform, a_nonuniform, a_unknown, a_staggered
  type(ibm_lagrangian_points_t) :: lag3, lag2
  type(ibm_grid_t) :: grid
  integer :: i,j,k,ios, status, unit_id
  real(mytype) :: eta, max_err, temp_err, y_center, h, uc, lx, lz, q
  integer :: const_status, linear_status, poi_status, periodic_status
  integer :: nonuniform_status, unknown_status, staggered_status
  integer :: monotonicity_flag

  do i=1,nx
    x(i) = (real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype))
  end do
  do j=1,ny
    y(j) = -1._mytype + (real(j,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype))
    eta = (real(j,mytype)-0.5_mytype)/real(ny,mytype)
    y_nonuniform(j) = -1._mytype + 2._mytype*eta*eta
  end do
  do k=1,nz
    z(k) = (real(k,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype))
  end do

  call init_stage4_grid_adapter_from_arrays(a_uniform, x, y, z, .true., .false., .true., 1)
  call init_stage4_grid_adapter_from_arrays(a_nonuniform, x, y_nonuniform, z, .true., .false., .true., 1)
  call init_stage4_grid_adapter_from_arrays(a_unknown, x, y, z, .true., .false., .true., 0)
  call init_stage4_grid_adapter_from_arrays(a_staggered, x, y, z, .true., .false., .true., 2)

  lag3%nl = 3
  allocate(lag3%x(3,3), lag3%v(3,3), lag3%force(3,3), lag3%weight(3))
  lag3%x(:,1) = [0.73_mytype, 0.17_mytype, 0.41_mytype]
  lag3%x(:,2) = [1.22_mytype,-0.25_mytype, 0.62_mytype]
  lag3%x(:,3) = [1.00_mytype, 0.00_mytype, 0.50_mytype]
  lag3%v = 0._mytype; lag3%force=0._mytype; lag3%weight=1._mytype

  lag2%nl = 2
  allocate(lag2%x(3,2), lag2%v(3,2), lag2%force(3,2), lag2%weight(2))
  lag2%x(:,1) = [0.05_mytype, 0.00_mytype, 0.95_mytype]
  lag2%x(:,2) = [1.95_mytype, 0.00_mytype, 0.05_mytype]
  lag2%v = 0._mytype; lag2%force=0._mytype; lag2%weight=1._mytype

  allocate(ux(nx,ny,nz), uy(nx,ny,nz), uz(nx,ny,nz))
  allocate(ux_lag(3,lag3%nl), ux_lag_lin(3,lag3%nl), ux_lag_poi(3,lag3%nl), ux_lag_per(3,lag2%nl))

  ux = 0.2_mytype; uy=-0.1_mytype; uz=0.05_mytype
  call interpolate_stage4_vector_to_lag_if_supported(a_uniform, ux, uy, uz, lag3, ux_lag, const_status)
  max_err = maxval(abs(ux_lag(1,:)-0.2_mytype))
  max_err = max(max_err, maxval(abs(ux_lag(2,:)+0.1_mytype)))
  max_err = max(max_err, maxval(abs(ux_lag(3,:)-0.05_mytype)))

  call convert_stage4_adapter_to_ibm_grid(a_uniform, grid, status)
  allocate(wsum(lag3%nl))
  call compute_interpolation_weight_sum(grid, lag3, wsum)
  call destroy_ibm_grid(grid)

  ! linear
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=1._mytype+0.2_mytype*x(i)-0.1_mytype*y(j)+0.05_mytype*z(k)
    uy(i,j,k)=-0.3_mytype+0.4_mytype*y(j)
    uz(i,j,k)=0.7_mytype-0.2_mytype*z(k)+0.1_mytype*x(i)
  end do; end do; end do
  call interpolate_stage4_vector_to_lag_if_supported(a_uniform, ux, uy, uz, lag3, ux_lag_lin, linear_status)
  max_err = max(max_err, 0._mytype)
  temp_err = 0._mytype
  do i=1,lag3%nl
    temp_err = max(temp_err, abs(ux_lag_lin(1,i) - (1._mytype + 0.2_mytype*lag3%x(1,i)-0.1_mytype*lag3%x(2,i)+0.05_mytype*lag3%x(3,i))))
    temp_err = max(temp_err, abs(ux_lag_lin(2,i) - (-0.3_mytype + 0.4_mytype*lag3%x(2,i))))
    temp_err = max(temp_err, abs(ux_lag_lin(3,i) - (0.7_mytype - 0.2_mytype*lag3%x(3,i)+0.1_mytype*lag3%x(1,i))))
  end do

  ! poiseuille
  uc = 1._mytype; y_center=0._mytype; h=1._mytype
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=uc*(1._mytype-((y(j)-y_center)/h)**2)
    uy(i,j,k)=0._mytype; uz(i,j,k)=0._mytype
  end do; end do; end do
  lag3%x(:,1)=[1._mytype,-0.5_mytype,0.5_mytype]
  lag3%x(:,2)=[1._mytype, 0.0_mytype,0.5_mytype]
  lag3%x(:,3)=[1._mytype, 0.5_mytype,0.5_mytype]
  call interpolate_stage4_vector_to_lag_if_supported(a_uniform, ux, uy, uz, lag3, ux_lag_poi, poi_status)

  ! periodic
  lx=2._mytype; lz=1._mytype
  do k=1,nz; do j=1,ny; do i=1,nx
    q = sin(2._mytype*pi*x(i)/lx) + 0.25_mytype*cos(2._mytype*pi*z(k)/lz) + 0.1_mytype*y(j)
    ux(i,j,k)=q; uy(i,j,k)=0._mytype; uz(i,j,k)=0._mytype
  end do; end do; end do
  call interpolate_stage4_vector_to_lag_if_supported(a_uniform, ux, uy, uz, lag2, ux_lag_per, periodic_status)

  ! blocked tests
  call interpolate_stage4_vector_to_lag_if_supported(a_nonuniform, ux, uy, uz, lag2, ux_lag_per, nonuniform_status)
  call interpolate_stage4_vector_to_lag_if_supported(a_unknown, ux, uy, uz, lag2, ux_lag_per, unknown_status)
  call interpolate_stage4_vector_to_lag_if_supported(a_staggered, ux, uy, uz, lag2, ux_lag_per, staggered_status)

  monotonicity_flag = merge(1,0, ux_lag_poi(1,2) > ux_lag_poi(1,1) .and. ux_lag_poi(1,2) > ux_lag_poi(1,3))

  open(newunit=unit_id,file='stage4_outputs/fibre_stage4_interpolation_check.dat',status='replace',action='write',iostat=ios)
  write(unit_id,'(A,1X,I0)') 'stage4_const_interp_status', const_status
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_const_interp_max_error', maxval([abs(ux_lag(1,1)-0.2_mytype),abs(ux_lag(1,2)-0.2_mytype),abs(ux_lag(1,3)-0.2_mytype), &
                                                                            abs(ux_lag(2,1)+0.1_mytype),abs(ux_lag(2,2)+0.1_mytype),abs(ux_lag(2,3)+0.1_mytype), &
                                                                            abs(ux_lag(3,1)-0.05_mytype),abs(ux_lag(3,2)-0.05_mytype),abs(ux_lag(3,3)-0.05_mytype)])
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_const_weight_sum_max_error', maxval(abs(wsum-1._mytype))
  write(unit_id,'(A,1X,I0)') 'stage4_linear_interp_status', linear_status
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_linear_interp_max_error', temp_err
  write(unit_id,'(A,1X,I0)') 'stage4_poiseuille_interp_status', poi_status
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_poiseuille_interp_max_error', maxval(abs([ux_lag_poi(1,1)-0.75_mytype, ux_lag_poi(1,2)-1._mytype, ux_lag_poi(1,3)-0.75_mytype]))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_poiseuille_center_value_error', abs(ux_lag_poi(1,2)-1._mytype)
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_poiseuille_symmetry_error', abs(ux_lag_poi(1,1)-ux_lag_poi(1,3))
  write(unit_id,'(A,1X,I0)') 'stage4_poiseuille_monotonicity_flag', monotonicity_flag
  write(unit_id,'(A,1X,I0)') 'stage4_periodic_interp_status', periodic_status
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_periodic_interp_error_max', max(abs(ux_lag_per(1,1)-(sin(2._mytype*pi*0.05_mytype/lx)+0.25_mytype*cos(2._mytype*pi*0.95_mytype/lz))), &
                                                                             abs(ux_lag_per(1,2)-(sin(2._mytype*pi*1.95_mytype/lx)+0.25_mytype*cos(2._mytype*pi*0.05_mytype/lz))))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_periodic_weight_sum_max_error', maxval(abs(wsum-1._mytype))
  write(unit_id,'(A,1X,I0)') 'stage4_nonuniform_y_uniform_ibm_compatible', merge(1,0,a_nonuniform%uniform_ibm_compatible)
  write(unit_id,'(A,1X,I0)') 'stage4_nonuniform_y_blocked_flag', merge(1,0,nonuniform_status==0)
  write(unit_id,'(A,1X,I0)') 'stage4_unknown_layout_blocked_flag', merge(1,0,unknown_status==0)
  write(unit_id,'(A,1X,I0)') 'stage4_staggered_layout_blocked_flag', merge(1,0,staggered_status==0)
  write(unit_id,'(A,1X,I0)') 'stage4_interp_fluid_rhs_modified', 0
  write(unit_id,'(A,1X,I0)') 'stage4_interpolation_status', 1
  close(unit_id)

end program fibre_stage4_interpolation_check
