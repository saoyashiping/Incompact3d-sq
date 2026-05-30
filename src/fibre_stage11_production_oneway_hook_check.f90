program fibre_stage11_production_oneway_hook_check
  use decomp_2d_constants, only: mytype
  use fibre_stage11_config, only: stage11_config_load, stage11_requested, stage11_readonly_mode
  use fibre_stage11_production_oneway_hook, only: stage11_production_oneway_init, stage11_production_oneway_sample, stage11_production_oneway_finalize, stage11_production_oneway_get_status_values
  implicit none
  integer, parameter :: nx=8, ny=9, nz=8
  real(mytype) :: ux(nx,ny,nz), uy(nx,ny,nz), uz(nx,ny,nz), su0,sv0,sw0,su1,sv1,sw1
  integer :: i,j,k, rf,ro,a,b,c,d,e,f,g,h,imod,jrhs,k1,l,m,n,o,p,io,ios
  integer :: field_unchanged_status
  logical :: pass
  call execute_command_line('mkdir -p stage11_outputs', exitstat=ios)
  call stage11_config_load(); rf=merge(1,0,stage11_requested()); ro=merge(1,0,stage11_readonly_mode())
  do k=1,nz; do j=1,ny; do i=1,nx
    ux(i,j,k)=sin(real(i+j+k,mytype)); uy(i,j,k)=cos(real(i-j+k,mytype)); uz(i,j,k)=0.1_mytype*real(i*j*k,mytype)
  enddo; enddo; enddo
  su0=sum(ux); sv0=sum(uy); sw0=sum(uz)
  call stage11_production_oneway_init(); call stage11_production_oneway_sample(ux,uy,uz); call stage11_production_oneway_finalize()
  su1=sum(ux); sv1=sum(uy); sw1=sum(uz)
  field_unchanged_status=merge(1,0,abs(su1-su0)<1.0e-14_mytype .and. abs(sv1-sv0)<1.0e-14_mytype .and. abs(sw1-sw0)<1.0e-14_mytype)
  call stage11_production_oneway_get_status_values(a,b,c,d,e,f,g,h,imod,jrhs,k1,l,m,n,o,p)
  open(newunit=io,file='stage11_outputs/fibre_stage11_5_production_oneway_hook_check.dat',status='replace',action='write',iostat=ios)
  write(io,'(A,1X,I0)') 'stage11_5_check_requested_flag', rf
  write(io,'(A,1X,I0)') 'stage11_5_check_readonly_mode_status', ro
  write(io,'(A,1X,I0)') 'stage11_5_check_hook_initialized_status', c
  write(io,'(A,1X,I0)') 'stage11_5_check_hook_sample_called_status', d
  write(io,'(A,1X,I0)') 'stage11_5_check_sample_performed_status', e
  write(io,'(A,1X,I0)') 'stage11_5_check_sample_count_status', f
  write(io,'(A,1X,I0)') 'stage11_5_check_sampled_velocity_finite_status', g
  write(io,'(A,1X,I0)') 'stage11_5_check_field_unchanged_status', field_unchanged_status
  write(io,'(A,1X,I0)') 'stage11_5_check_rhs_modified_status', jrhs
  write(io,'(A,1X,I0)') 'stage11_5_check_no_rhs_injection_status', k1
  write(io,'(A,1X,I0)') 'stage11_5_check_no_ibm_spreading_status', l
  write(io,'(A,1X,I0)') 'stage11_5_check_no_feedback_force_status', m
  write(io,'(A,1X,I0)') 'stage11_5_check_no_twoway_force_status', n
  write(io,'(A,1X,I0)') 'stage11_5_check_no_structure_advance_status', o
  write(io,'(A,1X,I0)') 'stage11_5_check_production_oneway_hook_status', p
  close(io)
  pass = (rf==1 .and. ro==1 .and. c==1 .and. d==1 .and. e==1 .and. f==1 .and. g==1 .and. field_unchanged_status==1 .and. jrhs==0 .and. k1==1 .and. l==1 .and. m==1 .and. n==1 .and. o==1 .and. p==1)
  if (pass) then
    print *, 'STAGE 11.5 PRODUCTION ONEWAY HOOK CHECK VERDICT: PASS'
  else
    print *, 'STAGE 11.5 PRODUCTION ONEWAY HOOK CHECK VERDICT: FAIL'
    stop 1
  endif
end program
