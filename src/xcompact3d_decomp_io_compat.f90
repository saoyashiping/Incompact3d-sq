module xcompact3d_decomp_io_compat
  use decomp_2d_constants, only : mytype
  use decomp_2d, only : DECOMP_INFO

  implicit none
  private

  integer, parameter, public :: decomp_2d_read_mode = 0
  integer, parameter, public :: decomp_2d_write_mode = 1
  integer, parameter, public :: decomp_2d_append_mode = 2

  public :: decomp_2d_init_io, decomp_2d_register_variable
  public :: decomp_2d_open_io, decomp_2d_close_io
  public :: decomp_2d_start_io, decomp_2d_end_io
  public :: decomp_2d_write_one, decomp_2d_read_one
  public :: decomp_2d_write_plane
  public :: fine_to_coarseS, fine_to_coarseV
  public :: gen_iodir_name

  interface decomp_2d_write_one
    module procedure x3d_write_one_r3_simple
    module procedure x3d_write_one_r3_legacy
  end interface

  interface decomp_2d_read_one
    module procedure x3d_read_one_r3_simple
    module procedure x3d_read_one_r3_legacy
  end interface

  interface decomp_2d_write_plane
    module procedure x3d_write_plane_r3_simple
    module procedure x3d_write_plane_r3_legacy
  end interface

  interface fine_to_coarseS
    module procedure fine_to_coarseS_r3
  end interface

  interface fine_to_coarseV
    module procedure fine_to_coarseV_r3
  end interface

contains

  subroutine decomp_2d_init_io(io_name)
    character(len=*), intent(in) :: io_name
  end subroutine

  subroutine decomp_2d_register_variable(io_name, varname, ipencil, iscalar, output2D, dtype, opt_decomp, opt_nplanes)
    character(len=*), intent(in) :: io_name, varname
    integer, intent(in) :: ipencil, iscalar, output2D
    integer, intent(in) :: dtype
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    integer, intent(in), optional :: opt_nplanes
  end subroutine

  subroutine decomp_2d_open_io(io_name, dirname, mode)
    character(len=*), intent(in) :: io_name, dirname
    integer, intent(in) :: mode
  end subroutine

  subroutine decomp_2d_close_io(io_name, dirname)
    character(len=*), intent(in) :: io_name, dirname
  end subroutine

  subroutine decomp_2d_start_io(io_name, dirname)
    character(len=*), intent(in) :: io_name, dirname
  end subroutine

  subroutine decomp_2d_end_io(io_name, dirname)
    character(len=*), intent(in) :: io_name, dirname
  end subroutine

  function gen_iodir_name(prefix, io_name) result(name)
    character(len=*), intent(in) :: prefix, io_name
    character(len=256) :: name
    name = trim(prefix) // '/' // trim(io_name)
  end function

  subroutine fine_to_coarseS_r3(ipencil, fine, coarse)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: fine(:,:,:)
    real(mytype), intent(out) :: coarse(:,:,:)
    integer :: i, j, k, ifine, jfine, kfine
    integer :: nf1, nf2, nf3, nc1, nc2, nc3

    nf1 = size(fine,1); nf2 = size(fine,2); nf3 = size(fine,3)
    nc1 = size(coarse,1); nc2 = size(coarse,2); nc3 = size(coarse,3)

    do k = 1, nc3
      kfine = min(nf3, max(1, 1 + (k-1) * max(1, nf3 / max(1,nc3))))
      do j = 1, nc2
        jfine = min(nf2, max(1, 1 + (j-1) * max(1, nf2 / max(1,nc2))))
        do i = 1, nc1
          ifine = min(nf1, max(1, 1 + (i-1) * max(1, nf1 / max(1,nc1))))
          coarse(i,j,k) = fine(ifine,jfine,kfine)
        end do
      end do
    end do
  end subroutine

  subroutine fine_to_coarseV_r3(ipencil, fine, coarse)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: fine(:,:,:)
    real(mytype), intent(out) :: coarse(:,:,:)
    call fine_to_coarseS_r3(ipencil, fine, coarse)
  end subroutine

  subroutine x3d_write_one_r3_simple(ipencil, array, filename, mode, opt_decomp, reduce_prec, opt_deferred_writes)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    logical, intent(in), optional :: opt_deferred_writes
    ! Stage 9.1 build-only compatibility path; real restart/output semantics must be audited in Stage 9.10.
  end subroutine

  subroutine x3d_write_one_r3_legacy(ipencil, array, dirname, filename, mode, io_name, opt_decomp, reduce_prec, opt_deferred_writes)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: dirname, filename
    integer, intent(in) :: mode
    character(len=*), intent(in), optional :: io_name
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    logical, intent(in), optional :: opt_deferred_writes
    character(len=:), allocatable :: fullpath

    fullpath = trim(dirname)//'/'//trim(filename)
    call x3d_write_one_r3_simple(ipencil, array, fullpath, mode, opt_decomp, reduce_prec, opt_deferred_writes)
  end subroutine

  subroutine x3d_read_one_r3_simple(ipencil, array, filename, opt_decomp, reduce_prec)
    integer, intent(in) :: ipencil
    real(mytype), intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: filename
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    ! Stage 9.1 build-only compatibility path; real restart/output semantics must be audited in Stage 9.10.
    array = 0._mytype
  end subroutine

  subroutine x3d_read_one_r3_legacy(ipencil, array, dirname, filename, io_name, opt_decomp, reduce_prec)
    integer, intent(in) :: ipencil
    real(mytype), intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: dirname, filename
    character(len=*), intent(in), optional :: io_name
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    character(len=:), allocatable :: fullpath

    fullpath = trim(dirname)//'/'//trim(filename)
    call x3d_read_one_r3_simple(ipencil, array, fullpath, opt_decomp, reduce_prec)
  end subroutine

  subroutine x3d_write_plane_r3_simple(ipencil, array, output2D, plane_index, filename)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    integer, intent(in) :: output2D
    integer, intent(in) :: plane_index
    character(len=*), intent(in) :: filename
    ! Stage 9.1 build-only compatibility path; real restart/output semantics must be audited in Stage 9.10.
  end subroutine

  subroutine x3d_write_plane_r3_legacy(ipencil, array, output2D, plane_index, dirname, filename, io_name)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    integer, intent(in) :: output2D
    integer, intent(in) :: plane_index
    character(len=*), intent(in) :: dirname, filename
    character(len=*), intent(in), optional :: io_name
    character(len=:), allocatable :: fullpath

    fullpath = trim(dirname)//'/'//trim(filename)
    call x3d_write_plane_r3_simple(ipencil, array, output2D, plane_index, fullpath)
  end subroutine

end module xcompact3d_decomp_io_compat
