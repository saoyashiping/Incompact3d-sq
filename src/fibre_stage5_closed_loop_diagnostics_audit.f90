program fibre_stage5_closed_loop_diagnostics_audit
  use fibre_parameters, only: mytype
  implicit none

  integer, parameter :: nx=8, ny=6, nz=5, nlag=4
  real(mytype), parameter :: dx=0.125_mytype, dy=0.2_mytype, dz=0.25_mytype
  real(mytype), parameter :: tiny_val=1.0e-30_mytype

  real(mytype) :: dv
  real(mytype), dimension(nx,ny,nz) :: ux, uy, uz
  real(mytype), dimension(3,nx,ny,nz) :: f_eul
  real(mytype), dimension(3,nlag) :: xlag, flag, ulag
  real(mytype), dimension(nlag) :: ds
  real(mytype), dimension(3) :: fe_total, fl_total
  real(mytype) :: pe_total, pl_total
  real(mytype) :: force_abs_error, force_rel_error, power_abs_error, power_rel_error
  integer :: i, j, k, l, io
  integer :: eulerian_force_computed_flag, lagrangian_force_computed_flag
  integer :: eulerian_power_computed_flag, lagrangian_power_computed_flag
  integer :: used_spreading_buffer_flag, used_lagrangian_force_flag, used_interpolated_velocity_flag
  integer :: no_tautological_force_flag, no_tautological_power_flag, audit_status

  dv = dx*dy*dz
  f_eul = 0._mytype

  do k=1,nz
    do j=1,ny
      do i=1,nx
        ux(i,j,k) = 0.1_mytype + 0.01_mytype*real(i,mytype) + 0.02_mytype*real(j,mytype)
        uy(i,j,k) = -0.2_mytype + 0.03_mytype*real(j,mytype) + 0.01_mytype*real(k,mytype)
        uz(i,j,k) = 0.05_mytype + 0.01_mytype*real(i,mytype) - 0.02_mytype*real(k,mytype)
      end do
    end do
  end do

  xlag(:,1) = (/2.3_mytype, 2.6_mytype, 2.2_mytype/)
  xlag(:,2) = (/3.7_mytype, 3.2_mytype, 2.8_mytype/)
  xlag(:,3) = (/5.1_mytype, 4.4_mytype, 3.6_mytype/)
  xlag(:,4) = (/6.2_mytype, 2.7_mytype, 4.1_mytype/)

  flag(:,1) = (/ 0.8_mytype, -0.3_mytype,  0.5_mytype/)
  flag(:,2) = (/-0.6_mytype,  0.4_mytype, -0.2_mytype/)
  flag(:,3) = (/ 0.3_mytype,  0.7_mytype, -0.1_mytype/)
  flag(:,4) = (/-0.2_mytype, -0.5_mytype,  0.9_mytype/)
  ds = (/0.12_mytype, 0.10_mytype, 0.09_mytype, 0.11_mytype/)

  do l=1,nlag
    call spread_force_cic(nx,ny,nz,xlag(:,l),flag(:,l),ds(l),dv,f_eul)
    call interpolate_velocity_cic(nx,ny,nz,xlag(:,l),ux,uy,uz,ulag(:,l))
  end do

  fe_total = 0._mytype
  pe_total = 0._mytype
  do k=1,nz
    do j=1,ny
      do i=1,nx
        fe_total(:) = fe_total(:) + f_eul(:,i,j,k)*dv
        pe_total = pe_total + (f_eul(1,i,j,k)*ux(i,j,k) + f_eul(2,i,j,k)*uy(i,j,k) + f_eul(3,i,j,k)*uz(i,j,k))*dv
      end do
    end do
  end do

  fl_total = 0._mytype
  pl_total = 0._mytype
  do l=1,nlag
    fl_total(:) = fl_total(:) + flag(:,l)*ds(l)
    pl_total = pl_total + dot_product(flag(:,l),ulag(:,l))*ds(l)
  end do

  force_abs_error = sqrt(sum((fe_total-fl_total)**2))
  force_rel_error = force_abs_error/max(sqrt(sum(fl_total**2)),tiny_val)
  power_abs_error = abs(pe_total-pl_total)
  power_rel_error = power_abs_error/max(abs(pl_total),tiny_val)

  eulerian_force_computed_flag = merge(1,0,sum(abs(fe_total)) > 0._mytype)
  lagrangian_force_computed_flag = merge(1,0,sum(abs(fl_total)) > 0._mytype)
  eulerian_power_computed_flag = 1
  lagrangian_power_computed_flag = 1
  used_spreading_buffer_flag = merge(1,0,sum(abs(f_eul)) > 0._mytype)
  used_lagrangian_force_flag = merge(1,0,sum(abs(flag)) > 0._mytype)
  used_interpolated_velocity_flag = merge(1,0,sum(abs(ulag)) > 0._mytype)

  no_tautological_force_flag = merge(1,0, used_spreading_buffer_flag==1 .and. used_lagrangian_force_flag==1 .and. &
       eulerian_force_computed_flag==1 .and. lagrangian_force_computed_flag==1)
  no_tautological_power_flag = merge(1,0, used_interpolated_velocity_flag==1 .and. used_lagrangian_force_flag==1 .and. &
       eulerian_power_computed_flag==1 .and. lagrangian_power_computed_flag==1)

  audit_status = merge(1,0, force_abs_error<=1e-12_mytype .and. force_rel_error<=1e-12_mytype .and. &
       power_abs_error<=1e-12_mytype .and. power_rel_error<=1e-12_mytype .and. &
       no_tautological_force_flag==1 .and. no_tautological_power_flag==1)

  open(newunit=io,file='stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage5_closed_loop_eulerian_force_computed_flag', eulerian_force_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_lagrangian_force_computed_flag', lagrangian_force_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_eulerian_power_computed_flag', eulerian_power_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_lagrangian_power_computed_flag', lagrangian_power_computed_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_used_spreading_buffer_flag', used_spreading_buffer_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_used_lagrangian_force_flag', used_lagrangian_force_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_used_interpolated_velocity_flag', used_interpolated_velocity_flag
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_force_conservation_abs_error', force_abs_error
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_force_conservation_relative_error', force_rel_error
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_power_abs_error', power_abs_error
  write(io,'(A,1X,ES24.16E3)') 'stage5_closed_loop_actual_power_relative_error', power_rel_error
  write(io,'(A,1X,I0)') 'stage5_closed_loop_no_tautological_force_flag', no_tautological_force_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_no_tautological_power_flag', no_tautological_power_flag
  write(io,'(A,1X,I0)') 'stage5_closed_loop_diagnostics_audit_status', audit_status
  close(io)

contains
  subroutine spread_force_cic(nx,ny,nz,xl,fl,ds,dv,fbuf)
    integer, intent(in) :: nx,ny,nz
    real(mytype), intent(in) :: xl(3), fl(3), ds, dv
    real(mytype), intent(inout) :: fbuf(3,nx,ny,nz)
    integer :: i0,j0,k0,ii,jj,kk,i,j,k
    real(mytype) :: rx,ry,rz,wx,wy,wz,w
    i0 = floor(xl(1)); j0 = floor(xl(2)); k0 = floor(xl(3))
    rx = xl(1)-real(i0,mytype); ry = xl(2)-real(j0,mytype); rz = xl(3)-real(k0,mytype)
    do kk=0,1
      wz = merge(1._mytype-rz,rz,kk==0)
      k = k0+kk
      do jj=0,1
        wy = merge(1._mytype-ry,ry,jj==0)
        j = j0+jj
        do ii=0,1
          wx = merge(1._mytype-rx,rx,ii==0)
          i = i0+ii
          w = wx*wy*wz
          fbuf(:,i,j,k) = fbuf(:,i,j,k) + fl(:)*ds*w/dv
        end do
      end do
    end do
  end subroutine

  subroutine interpolate_velocity_cic(nx,ny,nz,xl,ux,uy,uz,ul)
    integer, intent(in) :: nx,ny,nz
    real(mytype), intent(in) :: xl(3),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz)
    real(mytype), intent(out) :: ul(3)
    integer :: i0,j0,k0,ii,jj,kk,i,j,k
    real(mytype) :: rx,ry,rz,wx,wy,wz,w
    ul = 0._mytype
    i0 = floor(xl(1)); j0 = floor(xl(2)); k0 = floor(xl(3))
    rx = xl(1)-real(i0,mytype); ry = xl(2)-real(j0,mytype); rz = xl(3)-real(k0,mytype)
    do kk=0,1
      wz = merge(1._mytype-rz,rz,kk==0)
      k = k0+kk
      do jj=0,1
        wy = merge(1._mytype-ry,ry,jj==0)
        j = j0+jj
        do ii=0,1
          wx = merge(1._mytype-rx,rx,ii==0)
          i = i0+ii
          w = wx*wy*wz
          ul(1) = ul(1) + w*ux(i,j,k)
          ul(2) = ul(2) + w*uy(i,j,k)
          ul(3) = ul(3) + w*uz(i,j,k)
        end do
      end do
    end do
  end subroutine

end program fibre_stage5_closed_loop_diagnostics_audit
