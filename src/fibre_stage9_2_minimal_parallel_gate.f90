program fibre_stage9_2_minimal_parallel_gate
  use mpi
  use decomp_2d
  use decomp_2d_mpi, only: nrank, nproc
  use decomp_2d_constants, only: mytype
  implicit none

  integer :: ierr, fail_local, fail_global
  integer :: nxg, nyg, nzg, p_row_l, p_col_l
  logical :: periodic_bc(3)
  integer :: i
  integer(kind=8) :: local_cells, total_cells

  call MPI_Init(ierr)
  fail_local = 0

  nxg = 8
  nyg = 8
  nzg = 8
  p_row_l = 0
  p_col_l = 0
  periodic_bc = (/ .false., .false., .false. /)

  call decomp_2d_init(nxg, nyg, nzg, p_row_l, p_col_l, periodic_bc)

  write(*,'(A,I0,A,I0,A,3(I0,1X),A,3(I0,1X),A,3(I0,1X))') &
    '[Stage9.2] rank=', nrank, ' of ', nproc, ' xsize=', xsize(1), xsize(2), xsize(3), &
    ' xstart=', xstart(1), xstart(2), xstart(3), ' xend=', xend(1), xend(2), xend(3)

  write(*,'(A,I0,A,3(I0,1X),A,3(I0,1X),A,3(I0,1X))') &
    '[Stage9.2] rank=', nrank, ' ysize=', ysize(1), ysize(2), ysize(3), &
    ' ystart=', ystart(1), ystart(2), ystart(3), ' yend=', yend(1), yend(2), yend(3)

  write(*,'(A,I0,A,3(I0,1X),A,3(I0,1X),A,3(I0,1X))') &
    '[Stage9.2] rank=', nrank, ' zsize=', zsize(1), zsize(2), zsize(3), &
    ' zstart=', zstart(1), zstart(2), zstart(3), ' zend=', zend(1), zend(2), zend(3)

  do i = 1, 3
    if (xsize(i) < 0 .or. ysize(i) < 0 .or. zsize(i) < 0) fail_local = 1

    if (xsize(i) > 0) then
      if (xstart(i) > xend(i)) fail_local = 1
      if (xsize(i) /= xend(i) - xstart(i) + 1) fail_local = 1
    end if
    if (ysize(i) > 0) then
      if (ystart(i) > yend(i)) fail_local = 1
      if (ysize(i) /= yend(i) - ystart(i) + 1) fail_local = 1
    end if
    if (zsize(i) > 0) then
      if (zstart(i) > zend(i)) fail_local = 1
      if (zsize(i) /= zend(i) - zstart(i) + 1) fail_local = 1
    end if

    if (xstS(i) < 1 .or. xstV(i) < 1) fail_local = 1
    if (xenS(i) < xstS(i)-1) fail_local = 1
    if (xenV(i) < xstV(i)-1) fail_local = 1
  end do

  local_cells = int(max(0,xsize(1)),8) * int(max(0,xsize(2)),8) * int(max(0,xsize(3)),8)
  call MPI_Allreduce(local_cells, total_cells, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (total_cells <= 0_8) fail_local = 1

  call MPI_Allreduce(fail_local, fail_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

  if (nrank == 0) then
    if (fail_global == 0) then
      write(*,'(A)') 'STAGE 9.2 PROGRAM VERDICT: PASS'
    else
      write(*,'(A)') 'STAGE 9.2 PROGRAM VERDICT: FAIL'
    end if
  end if

  call decomp_2d_finalize()
  call MPI_Finalize(ierr)

  if (fail_global /= 0) stop 1
end program fibre_stage9_2_minimal_parallel_gate
