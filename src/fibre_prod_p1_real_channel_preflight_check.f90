program fibre_prod_p1_real_channel_preflight_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_p1_real_channel_preflight
  implicit none
  real(real64) :: ux(4,4,4), uy(4,4,4), uz(4,4,4)
  integer :: status
  ux=1.0_real64; uy=0.0_real64; uz=0.0_real64
  call fibre_prod_p1_real_channel_preflight_init_values(1,33,0.0_real64,0.10_real64,status); if (status/=0) error stop 1
  call fibre_prod_p1_real_channel_preflight_record_sampling(ux,uy,uz,status); if (status/=0) error stop 2
  call fibre_prod_p1_real_channel_preflight_write_diagnostics(status); if (status/=0) error stop 3
  ux(2,2,2)=ieee_value(ux(2,2,2), ieee_quiet_nan)
  call fibre_prod_p1_real_channel_preflight_record_sampling(ux,uy,uz,status); if (status==0) error stop 4
  call fibre_prod_p1_real_channel_preflight_init_values(1,33,1.0_real64,0.10_real64,status); if (status==0) error stop 5
  call fibre_prod_p1_real_channel_preflight_init_values(2,33,0.0_real64,0.10_real64,status); if (status==0) error stop 6
  call fibre_prod_p1_real_channel_preflight_init_values(1,32,0.0_real64,0.10_real64,status); if (status==0) error stop 7
  call fibre_prod_p1_real_channel_preflight_init_values(1,33,0.0_real64,0.60_real64,status); if (status==0) error stop 8
  write(*,'(A)') 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS'
end program
