program fibre_stage7_grid_metadata_check
  use fibre_parameters, only : mytype
  use fibre_stage7_grid_metadata
  implicit none
  type(stage7_channel_grid_t) :: g, gu, gbad
  integer :: io, flags6, flags70, v, r, modf, ppm, pm, rpc, pdc, fuc, fac, final_status
  integer :: s6dat_exists, s6smoke, s6prior, found_smoke, found_prior
  real(mytype) :: dymin, dymax, vmin, vmax, tv, ev, verr, uverr, chg, tol_geom, tol_noop
  integer :: nx,ny,nz, inv1,inv2,inv3,inv4
  integer :: uniform_y_detected_flag_saved, uniform_y_nonuniform_flag_saved
  real(mytype) :: uniform_y_volume_error_saved
  real(mytype), allocatable :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
  logical :: ex
  nx=16; ny=17; nz=12
  tol_geom=1.0e-12_mytype; tol_noop=1.0e-14_mytype
  open(newunit=io,file='stage7_outputs/fibre_stage7_grid_metadata_check.dat',status='replace',action='write')
  inquire(file='stage6_outputs/STAGE6_CLOSED.md',exist=ex); flags6=merge(1,0,ex)
  inquire(file='stage6_outputs/fibre_stage6_total_smoke_check.dat',exist=ex); s6dat_exists=merge(1,0,ex)
  call read_int_key_from_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',s6smoke,found_smoke)
  call read_int_key_from_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',s6prior,found_prior)
  if (.not.(flags6==1 .and. s6dat_exists==1 .and. found_smoke==1 .and. found_prior==1 .and. s6smoke==1 .and. s6prior==1)) flags6=0
  inquire(file='stage7_outputs/fibre_stage7_config_check.dat',exist=ex); flags70=merge(1,0,ex)
  write(io,'(A,1X,I0)') 'stage7_grid_stage6_closed_marker_status', flags6
  write(io,'(A,1X,I0)') 'stage7_grid_stage6_total_smoke_output_exists', s6dat_exists
  write(io,'(A,1X,I0)') 'stage7_grid_stage6_total_smoke_status', merge(s6smoke,0,found_smoke==1)
  write(io,'(A,1X,I0)') 'stage7_grid_stage6_all_prior_outputs_status', merge(s6prior,0,found_prior==1)
  write(io,'(A,1X,I0)') 'stage7_grid_stage7_0_dependency_status', flags70

  call init_stage7_nonuniform_channel_grid(g,nx,ny,nz)
  write(io,'(A,1X,I0)') 'stage7_grid_x_uniform_flag', g%x_uniform_flag
  write(io,'(A,1X,I0)') 'stage7_grid_z_uniform_flag', g%z_uniform_flag
  write(io,'(A,1X,I0)') 'stage7_grid_periodic_x_flag', g%periodic_x
  write(io,'(A,1X,I0)') 'stage7_grid_periodic_z_flag', g%periodic_z

  dymin=minval(g%dy_cell); dymax=maxval(g%dy_cell)
  write(io,'(A,1X,I0)') 'stage7_grid_y_nonuniform_flag', g%y_nonuniform_flag
  write(io,'(A,1X,I0)') 'stage7_grid_y_uniform_flag', g%y_uniform_flag
  write(io,'(A,1X,I0)') 'stage7_grid_y_monotonic_flag', g%y_monotonic_flag
  write(io,'(A,1X,ES24.16)') 'stage7_grid_dy_min', dymin
  write(io,'(A,1X,ES24.16)') 'stage7_grid_dy_max', dymax
  write(io,'(A,1X,I0)') 'stage7_grid_dy_min_positive_flag', merge(1,0,dymin>0._mytype)

  vmin=minval(g%volume_y); vmax=maxval(g%volume_y); call compute_stage7_total_volume(g,tv); ev=(g%xmax-g%xmin)*(g%ymax-g%ymin)*(g%zmax-g%zmin); call compute_stage7_volume_error(g,verr)
  write(io,'(A,1X,ES24.16)') 'stage7_grid_volume_min', vmin
  write(io,'(A,1X,ES24.16)') 'stage7_grid_volume_max', vmax
  write(io,'(A,1X,I0)') 'stage7_grid_volume_positive_flag', merge(1,0,vmin>0._mytype)
  write(io,'(A,1X,ES24.16)') 'stage7_grid_total_volume', tv
  write(io,'(A,1X,ES24.16)') 'stage7_grid_expected_total_volume', ev
  write(io,'(A,1X,ES24.16)') 'stage7_grid_total_volume_error', verr

  write(io,'(A,1X,ES24.16)') 'stage7_grid_dy_face_consistency_error_max', g%dy_face_consistency_error_max
  write(io,'(A,1X,I0)') 'stage7_grid_dy_face_consistency_status', g%dy_face_consistency_status
  write(io,'(A,1X,ES24.16)') 'stage7_grid_volume_formula_error_max', g%volume_formula_error_max
  write(io,'(A,1X,I0)') 'stage7_grid_volume_formula_status', g%volume_formula_status

  write(io,'(A,1X,I0)') 'stage7_grid_ymin_correct_flag', merge(1,0,abs(g%y_face(1)-g%ymin)<=1e-12_mytype)
  write(io,'(A,1X,I0)') 'stage7_grid_ymax_correct_flag', merge(1,0,abs(g%y_face(g%ny+1)-g%ymax)<=1e-12_mytype)
  write(io,'(A,1X,I0)') 'stage7_grid_wall_bounds_status', merge(1,0,abs(g%y_face(1)-g%ymin)<=1e-12_mytype .and. abs(g%y_face(g%ny+1)-g%ymax)<=1e-12_mytype)

  call init_stage7_uniform_channel_grid(gu,nx,ny,nz); call compute_stage7_volume_error(gu,uverr)
  uniform_y_detected_flag_saved=gu%y_uniform_flag
  uniform_y_nonuniform_flag_saved=gu%y_nonuniform_flag
  uniform_y_volume_error_saved=uverr
  write(io,'(A,1X,I0)') 'stage7_grid_uniform_y_detected_flag', uniform_y_detected_flag_saved
  write(io,'(A,1X,I0)') 'stage7_grid_uniform_y_nonuniform_flag', uniform_y_nonuniform_flag_saved
  write(io,'(A,1X,ES24.16)') 'stage7_grid_uniform_y_volume_error', uniform_y_volume_error_saved

  call init_stage7_nonuniform_channel_grid(gu,1,ny,nz); call validate_stage7_channel_grid(gu,v,r); inv1=r
  call init_stage7_nonuniform_channel_grid(gu,nx,ny,nz); gu%dy_cell(3)=gu%dy_cell(3)*1.25_mytype; call validate_stage7_channel_grid(gu,v,r); inv2=r
  call init_stage7_nonuniform_channel_grid(gu,nx,ny,nz); gu%periodic_y=1; call validate_stage7_channel_grid(gu,v,r); inv3=r
  call init_stage7_nonuniform_channel_grid(gbad,nx,ny,nz); gbad%volume_y(3)=gbad%volume_y(3)*1.25_mytype; call validate_stage7_channel_grid(gbad,v,r); inv4=r
  write(io,'(A,1X,I0)') 'stage7_grid_invalid_size_rejected_flag', inv1
  write(io,'(A,1X,I0)') 'stage7_grid_invalid_dy_rejected_flag', inv2
  write(io,'(A,1X,I0)') 'stage7_grid_invalid_boundary_rejected_flag', inv3
  write(io,'(A,1X,I0)') 'stage7_grid_invalid_volume_formula_rejected_flag', inv4

  allocate(rhsx(8,6,5),rhsy(8,6,5),rhsz(8,6,5))
  rhsx=1; rhsy=2; rhsz=3
  call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,chg,modf)
  write(io,'(A,1X,ES24.16)') 'stage7_grid_noop_rhs_change_max', chg
  write(io,'(A,1X,I0)') 'stage7_grid_noop_rhs_modified_flag', modf

  call stage7_grid_pressure_status(ppm,pm,rpc,pdc,fuc,fac)
  write(io,'(A,1X,I0)') 'stage7_grid_pressure_poisson_modified_flag', ppm
  write(io,'(A,1X,I0)') 'stage7_grid_projection_modified_flag', pm
  write(io,'(A,1X,I0)') 'stage7_grid_real_projection_called_flag', rpc
  write(io,'(A,1X,I0)') 'stage7_grid_production_dns_called_flag', pdc
  write(io,'(A,1X,I0)') 'stage7_grid_fluid_update_called_flag', fuc
  write(io,'(A,1X,I0)') 'stage7_grid_fibre_advance_called_flag', fac

  final_status=1
  if (flags6/=1) final_status=0
  if (s6dat_exists/=1) final_status=0
  if (s6smoke/=1 .or. found_smoke/=1) final_status=0
  if (s6prior/=1 .or. found_prior/=1) final_status=0
  if (flags70/=1) final_status=0
  if (g%x_uniform_flag/=1 .or. g%z_uniform_flag/=1) final_status=0
  if (g%periodic_x/=1 .or. g%periodic_z/=1) final_status=0
  if (g%y_nonuniform_flag/=1 .or. g%y_uniform_flag/=0 .or. g%y_monotonic_flag/=1) final_status=0
  if (dymin<=0._mytype .or. dymax<=dymin) final_status=0
  if (vmin<=0._mytype .or. vmax<vmin) final_status=0
  if (verr>tol_geom) final_status=0
  if (uniform_y_detected_flag_saved/=1 .or. uniform_y_nonuniform_flag_saved/=0) final_status=0
  if (uniform_y_volume_error_saved>tol_geom) final_status=0
  if (g%dy_face_consistency_status/=1 .or. g%dy_face_consistency_error_max>tol_geom) final_status=0
  if (g%volume_formula_status/=1 .or. g%volume_formula_error_max>tol_geom) final_status=0
  if (inv1/=1 .or. inv2/=1 .or. inv3/=1 .or. inv4/=1) final_status=0
  if (chg>tol_noop .or. modf/=0) final_status=0
  if (ppm/=0 .or. pm/=0 .or. rpc/=0 .or. pdc/=0 .or. fuc/=0 .or. fac/=0) final_status=0
  write(io,'(A,1X,I0)') 'stage7_grid_metadata_check_status', final_status
  close(io)
end program fibre_stage7_grid_metadata_check

subroutine read_int_key_from_dat(filename,key,value,found_flag)
  implicit none
  character(len=*), intent(in) :: filename, key
  integer, intent(out) :: value, found_flag
  integer :: io, ios
  character(len=512) :: line
  found_flag=0
  value=0
  open(newunit=io,file=filename,status='old',action='read',iostat=ios)
  if (ios/=0) return
  do
    read(io,'(A)',iostat=ios) line
    if (ios/=0) exit
    if (index(adjustl(line),trim(key)//' ')==1) then
      read(line(len_trim(key)+1:),*,iostat=ios) value
      if (ios==0) found_flag=1
      exit
    end if
  end do
  close(io)
end subroutine read_int_key_from_dat
