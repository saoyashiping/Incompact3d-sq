!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_delta

  use decomp_2d_constants, only : mytype

  implicit none

contains

  pure real(mytype) function delta_kernel_1d(r, h)
    ! 4-point regularized kernel (Peskin type), compact support |r/h| <= 2
    real(mytype), intent(in) :: r, h
    real(mytype) :: q, aq

    q = r / h
    aq = abs(q)

    if (aq <= 1._mytype) then
      delta_kernel_1d = (3._mytype - 2._mytype * aq + sqrt(1._mytype + 4._mytype * aq - 4._mytype * aq * aq)) / (8._mytype * h)
    else if (aq <= 2._mytype) then
      delta_kernel_1d = (5._mytype - 2._mytype * aq - sqrt(-7._mytype + 12._mytype * aq - 4._mytype * aq * aq)) / (8._mytype * h)
    else
      delta_kernel_1d = 0._mytype
    endif

  end function delta_kernel_1d

  pure real(mytype) function delta_kernel_3d(dx, dy, dz, hx, hy, hz)
    real(mytype), intent(in) :: dx, dy, dz, hx, hy, hz

    delta_kernel_3d = delta_kernel_1d(dx, hx) * delta_kernel_1d(dy, hy) * delta_kernel_1d(dz, hz)

  end function delta_kernel_3d

end module fiber_delta
