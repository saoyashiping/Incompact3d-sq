module fibre_stage7_channel_grid_adapter
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  implicit none
contains
  subroutine init_stage7_channel_grid_from_arrays(grid,nx,ny,nz,xmin,xmax,zmin,zmax,y_face_in,y_center_in,periodic_x,periodic_z,valid,rejected)
    type(stage7_channel_grid_t),intent(out)::grid
    integer,intent(in)::nx,ny,nz,periodic_x,periodic_z
    real(mytype),intent(in)::xmin,xmax,zmin,zmax,y_face_in(ny+1),y_center_in(ny)
    integer,intent(out)::valid,rejected
    integer::j
    call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz)
    grid%xmin=xmin; grid%xmax=xmax; grid%zmin=zmin; grid%zmax=zmax
    grid%dx=(xmax-xmin)/real(nx,mytype); grid%dz=(zmax-zmin)/real(nz,mytype)
    grid%y_face=y_face_in; grid%y_center=y_center_in
    do j=1,ny
      grid%dy_cell(j)=grid%y_face(j+1)-grid%y_face(j)
      grid%volume_y(j)=grid%dx*grid%dy_cell(j)*grid%dz
    end do
    grid%periodic_x=periodic_x; grid%periodic_z=periodic_z; grid%periodic_y=0
    call validate_stage7_channel_grid(grid,valid,rejected)
  end subroutine
  subroutine init_stage7_channel_grid_from_faces(grid,nx,ny,nz,xmin,xmax,zmin,zmax,y_face_in,periodic_x,periodic_z,valid,rejected)
    type(stage7_channel_grid_t),intent(out)::grid
    integer,intent(in)::nx,ny,nz,periodic_x,periodic_z
    real(mytype),intent(in)::xmin,xmax,zmin,zmax,y_face_in(ny+1)
    integer,intent(out)::valid,rejected
    real(mytype)::yc(ny); integer::j
    do j=1,ny; yc(j)=0.5_mytype*(y_face_in(j)+y_face_in(j+1)); end do
    call init_stage7_channel_grid_from_arrays(grid,nx,ny,nz,xmin,xmax,zmin,zmax,y_face_in,yc,periodic_x,periodic_z,valid,rejected)
  end subroutine
  subroutine compare_stage7_grid_metadata(g1,g2,err_max)
    type(stage7_channel_grid_t),intent(in)::g1,g2
    real(mytype),intent(out)::err_max
    err_max=max(abs(g1%dx-g2%dx),abs(g1%dz-g2%dz))
    err_max=max(err_max,maxval(abs(g1%y_face-g2%y_face)))
    err_max=max(err_max,maxval(abs(g1%y_center-g2%y_center)))
    err_max=max(err_max,maxval(abs(g1%dy_cell-g2%dy_cell)))
    err_max=max(err_max,maxval(abs(g1%volume_y-g2%volume_y)))
  end subroutine
end module
