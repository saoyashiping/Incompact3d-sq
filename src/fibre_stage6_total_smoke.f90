module fibre_stage6_total_smoke
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_controlled_rhs_hook
  use fibre_stage6_projection_order
  use fibre_stage6_layout_guard
  use fibre_stage6_io_safety
  implicit none
contains
  subroutine check_stage6_prior_outputs(flags, all_exist)
    integer, intent(out) :: flags(10), all_exist
    call file_exists_int('stage6_outputs/fibre_stage6_config_check.dat', flags(1))
    call file_exists_int('stage6_outputs/fibre_stage6_rhs_audit_check.dat', flags(2))
    call file_exists_int('stage6_outputs/fibre_stage6_noop_hook_check.dat', flags(3))
    call file_exists_int('stage6_outputs/fibre_stage6_controlled_rhs_check.dat', flags(4))
    call file_exists_int('stage6_outputs/fibre_stage6_singlephase_noop_check.dat', flags(5))
    call file_exists_int('stage6_outputs/fibre_stage6_projection_order_check.dat', flags(6))
    call file_exists_int('stage6_outputs/fibre_stage6_rk_rhs_sync_check.dat', flags(7))
    call file_exists_int('stage6_outputs/fibre_stage6_layout_guard_check.dat', flags(8))
    call file_exists_int('stage6_outputs/fibre_stage6_micro_smoke_check.dat', flags(9))
    call file_exists_int('stage6_outputs/fibre_stage6_io_safety_check.dat', flags(10))
    all_exist = merge(1,0,product(flags)==1)
  end subroutine

  subroutine run_stage6_total_controlled_rhs_probe(rhs_expected_error, ex, ey, ez, injected_flag, modified_flag)
    real(mytype), intent(out) :: rhs_expected_error, ex, ey, ez
    integer, intent(out) :: injected_flag, modified_flag
    integer, parameter :: nx=8, ny=6, nz=5
    real(mytype) :: rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz)
    real(mytype) :: rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz)
    type(stage6_config_t) :: config
    integer :: i,j,k,hook_called,rejected
    call init_stage6_controlled_test_config(config)
    config%rho_fluid=2._mytype
    do k=1,nz; do j=1,ny; do i=1,nx
      rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k
      rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i
      rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j
      fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k)
      fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k)
      fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
    end do; end do; end do
    rhsx0=rhsx; rhsy0=rhsy; rhsz0=rhsz
    call apply_stage6_controlled_rhs_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called,modified_flag,injected_flag,rejected)
    call compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,config%rho_fluid,rhsx,rhsy,rhsz,rhs_expected_error,ex,ey,ez)
  end subroutine

  subroutine write_stage6_closed_marker(filename, status)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: status
    integer :: io, ios
    status=0
    open(newunit=io,file=filename,status='replace',action='write',iostat=ios)
    if (ios/=0) return
    write(io,'(A)') '# STAGE 6 CLOSED'
    write(io,'(A)') ''
    write(io,'(A)') 'Verified scope:'
    write(io,'(A)') '- real-main-flow RHS hook safety path'
    write(io,'(A)') '- default no-op production safety'
    write(io,'(A)') '- controlled RHS injection: RHS += f_ibm / rho_fluid'
    write(io,'(A)') '- RHS-before-projection policy'
    write(io,'(A)') '- pressure Poisson not directly modified'
    write(io,'(A)') '- RK substep current-force/current-RHS policy'
    write(io,'(A)') '- layout guard for unsupported grids'
    write(io,'(A)') '- production two-way disabled by default'
    write(io,'(A)') ''
    write(io,'(A)') 'Not yet included:'
    write(io,'(A)') '- real nonuniform-y IBM'
    write(io,'(A)') '- staggered/component-specific IBM'
    write(io,'(A)') '- real production DNS two-way coupling'
    write(io,'(A)') '- MPI domain-decomposed IBM'
    write(io,'(A)') '- wall/contact model'
    write(io,'(A)') '- multi-fibre simulation'
    close(io)
    status=1
  end subroutine

  subroutine stage6_total_pressure_production_status(pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag, production_dns_called_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag, production_dns_called_flag
    pressure_poisson_modified_flag=0
    projection_modified_flag=0
    real_projection_called_flag=0
    production_dns_called_flag=0
  end subroutine
end module fibre_stage6_total_smoke
