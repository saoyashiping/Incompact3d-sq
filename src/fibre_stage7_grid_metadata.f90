module fibre_stage7_grid_metadata
  use fibre_parameters, only : mytype, pi
  use fibre_stage7_config
  implicit none
  type stage7_channel_grid_t
    integer :: nx, ny, nz, periodic_x, periodic_y, periodic_z
    integer :: x_uniform_flag, y_uniform_flag, y_nonuniform_flag, z_uniform_flag, y_monotonic_flag, grid_valid_flag
    real(mytype) :: xmin, xmax, ymin, ymax, zmin, zmax, dx, dz
    real(mytype), allocatable :: y_center(:), y_face(:), dy_cell(:), volume_y(:)
  end type
contains
  subroutine init_stage7_nonuniform_channel_grid(grid, nx, ny, nz)
    type(stage7_channel_grid_t), intent(out) :: grid
    integer, intent(in) :: nx, ny, nz
    integer :: j, valid, rej
    real(mytype) :: eta, a
    a=1.5_mytype
    grid%nx=nx; grid%ny=ny; grid%nz=nz; grid%xmin=0; grid%xmax=2*pi; grid%zmin=0; grid%zmax=pi; grid%ymin=-1; grid%ymax=1
    grid%dx=(grid%xmax-grid%xmin)/real(nx,mytype); grid%dz=(grid%zmax-grid%zmin)/real(nz,mytype)
    grid%periodic_x=1; grid%periodic_y=0; grid%periodic_z=1
    allocate(grid%y_center(ny),grid%y_face(ny+1),grid%dy_cell(ny),grid%volume_y(ny))
    do j=1,ny+1
      eta=-1._mytype + 2._mytype*real(j-1,mytype)/real(ny,mytype)
      grid%y_face(j)=tanh(a*eta)/tanh(a)
    end do
    do j=1,ny
      grid%dy_cell(j)=grid%y_face(j+1)-grid%y_face(j); grid%y_center(j)=0.5_mytype*(grid%y_face(j)+grid%y_face(j+1)); grid%volume_y(j)=grid%dx*grid%dy_cell(j)*grid%dz
    end do
    grid%x_uniform_flag=1; grid%z_uniform_flag=1; grid%y_uniform_flag=0; grid%y_nonuniform_flag=1
    call validate_stage7_channel_grid(grid,valid,rej)
  end subroutine
  subroutine init_stage7_uniform_channel_grid(grid, nx, ny, nz)
    type(stage7_channel_grid_t), intent(out) :: grid
    integer, intent(in) :: nx, ny, nz
    integer :: j, valid, rej
    grid%nx=nx; grid%ny=ny; grid%nz=nz; grid%xmin=0; grid%xmax=2*pi; grid%zmin=0; grid%zmax=pi; grid%ymin=-1; grid%ymax=1
    grid%dx=(grid%xmax-grid%xmin)/real(nx,mytype); grid%dz=(grid%zmax-grid%zmin)/real(nz,mytype)
    grid%periodic_x=1; grid%periodic_y=0; grid%periodic_z=1
    allocate(grid%y_center(ny),grid%y_face(ny+1),grid%dy_cell(ny),grid%volume_y(ny))
    do j=1,ny+1
      grid%y_face(j)=grid%ymin + (grid%ymax-grid%ymin)*real(j-1,mytype)/real(ny,mytype)
    end do
    do j=1,ny
      grid%dy_cell(j)=grid%y_face(j+1)-grid%y_face(j); grid%y_center(j)=0.5_mytype*(grid%y_face(j)+grid%y_face(j+1)); grid%volume_y(j)=grid%dx*grid%dy_cell(j)*grid%dz
    end do
    grid%x_uniform_flag=1; grid%z_uniform_flag=1; grid%y_uniform_flag=1; grid%y_nonuniform_flag=0
    call validate_stage7_channel_grid(grid,valid,rej)
  end subroutine
  subroutine validate_stage7_channel_grid(grid, valid, rejected_flag)
    type(stage7_channel_grid_t), intent(inout) :: grid
    integer, intent(out) :: valid, rejected_flag
    integer :: j
    valid=1
    if (grid%nx<=1 .or. grid%ny<=1 .or. grid%nz<=1) valid=0
    if (grid%dx<=0 .or. grid%dz<=0) valid=0
    do j=1,grid%ny
      if (grid%dy_cell(j)<=0 .or. grid%volume_y(j)<=0) valid=0
    end do
    do j=1,grid%ny
      if (grid%y_face(j+1)<=grid%y_face(j)) valid=0
    end do
    if (abs(grid%y_face(1)-grid%ymin)>1e-12_mytype) valid=0
    if (abs(grid%y_face(grid%ny+1)-grid%ymax)>1e-12_mytype) valid=0
    if (grid%periodic_x/=1 .or. grid%periodic_z/=1 .or. grid%periodic_y/=0) valid=0
    grid%y_monotonic_flag=merge(1,0,valid==1)
    rejected_flag=merge(0,1,valid==1)
    grid%grid_valid_flag=valid
  end subroutine
  subroutine compute_stage7_total_volume(grid,total_volume)
    type(stage7_channel_grid_t), intent(in) :: grid
    real(mytype), intent(out) :: total_volume
    total_volume=sum(grid%volume_y)*real(grid%nx*grid%nz,mytype)
  end subroutine
  subroutine compute_stage7_volume_error(grid, abs_error)
    type(stage7_channel_grid_t), intent(in) :: grid
    real(mytype), intent(out) :: abs_error
    real(mytype) :: tv, ev
    call compute_stage7_total_volume(grid,tv)
    ev=(grid%xmax-grid%xmin)*(grid%ymax-grid%ymin)*(grid%zmax-grid%zmin)
    abs_error=abs(tv-ev)
  end subroutine
  subroutine stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,rhs_change_max,rhs_modified_flag)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype), intent(out) :: rhs_change_max
    integer, intent(out) :: rhs_modified_flag
    real(mytype), allocatable :: x0(:,:,:),y0(:,:,:),z0(:,:,:)
    allocate(x0(size(rhsx,1),size(rhsx,2),size(rhsx,3)),y0(size(rhsy,1),size(rhsy,2),size(rhsy,3)),z0(size(rhsz,1),size(rhsz,2),size(rhsz,3)))
    x0=rhsx; y0=rhsy; z0=rhsz
    rhs_change_max=max(maxval(abs(rhsx-x0)),max(maxval(abs(rhsy-y0)),maxval(abs(rhsz-z0))))
    rhs_modified_flag=0; deallocate(x0,y0,z0)
  end subroutine
  subroutine stage7_grid_pressure_status(pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag, production_dns_called_flag, fluid_update_called_flag, fibre_advance_called_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag, production_dns_called_flag, fluid_update_called_flag, fibre_advance_called_flag
    pressure_poisson_modified_flag=0; projection_modified_flag=0; real_projection_called_flag=0; production_dns_called_flag=0; fluid_update_called_flag=0; fibre_advance_called_flag=0
  end subroutine
end module fibre_stage7_grid_metadata
