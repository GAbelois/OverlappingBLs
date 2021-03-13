!-----------------------------------------------------------------------
! 
!----------------------------------------------------------------------- 
subroutine fftn_mpi(in, out)
!-----------------------------------------------------------------------
! FFT in x and y directions using MPI and 1.5D decomposition
!----------------------------------------------------------------------- 
use types, only: rprec
use param
use fft
implicit none
! ---
real(kind=rprec),    dimension(ldx, nynpy), intent(in)  :: in
complex(kind=rprec), dimension(nxhnpy, nyt), intent(out) :: out
! ---
complex(kind=rprec), dimension(lhx, nynpy)  :: in_hat_x
complex(kind=rprec), dimension(nxhnpy, nyt) :: in_hat_y

integer :: i, j
! ---
do j = 1, nynpy
    call dfftw_execute_dft_r2c(plan_xf, in(1:nxt,j), in_hat_x(:,j))       ! 1D real to complex FFT in x direction
end do

in_hat_x(lhx, :) = (0._rprec, 0._rprec)

call transpose_y_to_x(in_hat_x, in_hat_y)               ! transpose to x pencil
! ---
do i = 1, nxhnpy
    call dfftw_execute_dft(plan_ycf, in_hat_y(i,:), out(i,:))         ! 1D complex to complex FFT in y direction
end do

out(:, nyt/2+1) = (0._rprec, 0._rprec)
! ---
end subroutine fftn_mpi
!-----------------------------------------------------------------------
! 
!----------------------------------------------------------------------- 
subroutine ifftn_mpi(in, out)
!-----------------------------------------------------------------------
! Inverse FFT in x and y directions using MPI and 1.5D decomposition
! Need to do ifft in reversed order of fft.
!----------------------------------------------------------------------- 
use types, only: rprec
use param
use fft
implicit none
! ---
complex(kind=rprec), dimension(nxhnpy, nyt), intent(in) :: in
real(kind=rprec),    dimension(ldx, nynpy), intent(out) :: out
! ---
complex(kind=rprec), dimension(lhx, nynpy)  :: in_hat_x
complex(kind=rprec), dimension(nxhnpy, nyt) :: in_hat_y

integer :: i, j
! ---
do i = 1, nxhnpy
    call dfftw_execute_dft(plan_ycb, in(i,:), in_hat_y(i,:))
end do

call transpose_x_to_y(in_hat_y, in_hat_x)

do j = 1, nynpy
    call dfftw_execute_dft_c2r(plan_xb, in_hat_x(:,j), out(1:nxt,j))       ! 1D real to complex FFT in x direction
end do
! ---
end subroutine ifftn_mpi
!-----------------------------------------------------------------------
! Transpose data between x and y directions
!----------------------------------------------------------------------- 
subroutine transpose_y_to_x(in, out)
!-----------------------------------------------------------------------
!  transpose to x pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: lhx, nynpy, nxhnpy, nyt, npy, comm_ver, mpi_cprec, ierr
implicit none
! ---
complex(rprec), dimension(lhx, nynpy),  intent(in)  :: in
complex(rprec), dimension(nxhnpy, nyt), intent(out) :: out
! ---
complex(rprec), dimension(nxhnpy, nynpy, 0:npy-1) :: work1
complex(rprec), dimension(nxhnpy, nynpy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nxhnpy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(st:en, :)
    st = en
    en = st + nxhnpy - 1
end do

num = nxhnpy * nynpy
call mpi_alltoall(work1, num, mpi_cprec, work2, num, mpi_cprec, comm_ver, ierr)

st = 1
en = nynpy
do ipy = 0, npy-1
    out(:, st:en) = work2(:,:,ipy)
    st = en + 1
    en = st + nynpy - 1
end do
! ---
end subroutine transpose_y_to_x
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_x_to_y(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: lhx, nynpy, nxhnpy, nyt, npy, comm_ver, mpi_cprec, ierr
implicit none
! ---
complex(rprec), dimension(nxhnpy, nyt), intent(in)  :: in
complex(rprec), dimension(lhx, nynpy),  intent(out) :: out
! ---
complex(rprec), dimension(nxhnpy, nynpy, 0:npy-1) :: work1
complex(rprec), dimension(nxhnpy, nynpy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nynpy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(:,st:en)
    st = en + 1
    en = st + nynpy - 1
end do

num = nxhnpy * nynpy
call mpi_alltoall(work1, num, mpi_cprec, work2, num, mpi_cprec, comm_ver, ierr)

st = 1
en = nxhnpy
do ipy = 0, npy-1
    out(st:en, :) = work2(:,:,ipy)
    st = en
    en = st + nxhnpy - 1
end do
! ---
end subroutine transpose_x_to_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_y_to_x2(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt2, nynpy, nx2npy, nyt, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxt2, nynpy),  intent(in)  :: in
real(rprec), dimension(nx2npy, nyt), intent(out) :: out
! ---
real(rprec), dimension(nx2npy, nynpy, 0:npy-1) :: work1
real(rprec), dimension(nx2npy, nynpy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nx2npy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(st:en, :)
    st = en + 1
    en = st + nx2npy - 1
end do

num = nx2npy * nynpy
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nynpy
do ipy = 0, npy-1
    out(:, st:en) = work2(:,:,ipy)
    st = en + 1
    en = st + nynpy - 1
end do
! ---
end subroutine transpose_y_to_x2
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_x2_to_y(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt2, nyt2, nx2npy, ny2npy, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nx2npy, nyt2),  intent(in)  :: in
real(rprec), dimension(nxt2, ny2npy), intent(out) :: out
! ---
real(rprec), dimension(nx2npy, ny2npy, 0:npy-1) :: work1
real(rprec), dimension(nx2npy, ny2npy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = ny2npy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(:, st:en)
    st = en + 1
    en = st + ny2npy - 1
end do

num = nx2npy * ny2npy
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nx2npy
do ipy = 0, npy-1
    out(st:en, :) = work2(:,:,ipy)
    st = en + 1
    en = st + nx2npy - 1
end do
! ---
end subroutine transpose_x2_to_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_y2_to_x(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, ny2npy, nxnpy, nyt2, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxt, ny2npy),  intent(in)  :: in
real(rprec), dimension(nxnpy, nyt2), intent(out) :: out
! ---
real(rprec), dimension(nxnpy, ny2npy, 0:npy-1) :: work1
real(rprec), dimension(nxnpy, ny2npy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nxnpy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(st:en, :)
    st = en + 1
    en = st + nxnpy - 1
end do

num = nxnpy * ny2npy
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = ny2npy
do ipy = 0, npy-1
    out(:, st:en) = work2(:,:,ipy)
    st = en + 1
    en = st + ny2npy - 1
end do
! ---
end subroutine transpose_y2_to_x
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_phy_x_to_y(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, nynpy, nxnpy, nyt, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxnpy, nyt), intent(in)  :: in
real(rprec), dimension(nxt, nynpy), intent(out) :: out
! ---
real(rprec), dimension(nxnpy, nynpy, 0:npy-1) :: work1
real(rprec), dimension(nxnpy, nynpy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nynpy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(:,st:en)
    st = en + 1
    en = st + nynpy - 1
end do

num = nxnpy * nynpy
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nxnpy
do ipy = 0, npy-1
    out(st:en, :) = work2(:,:,ipy)
    st = en + 1
    en = st + nxnpy - 1
end do
! ---
end subroutine transpose_phy_x_to_y
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine transpose_phy2_x_to_y(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, nynpy, nxnpy, nyt, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxnpy, nyt+1), intent(in)  :: in
real(rprec), dimension(nxt, nynpy+1), intent(out) :: out
! ---
real(rprec), dimension(nxnpy, nynpy+1, 0:npy-1) :: work1
real(rprec), dimension(nxnpy, nynpy+1, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nynpy+1
do ipy = 0, npy-1
    work1(:,:,ipy) = in(:,st:en)
    st = en
    en = st + nynpy
end do

num = nxnpy * (nynpy+1)
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nxnpy
do ipy = 0, npy-1
    out(st:en, :) = work2(:,:,ipy)
    st = en + 1
    en = st + nxnpy - 1
end do
! ---
end subroutine transpose_phy2_x_to_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_phy_y_to_x(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, nynpy, nxnpy, nyt, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxt, nynpy), intent(in)  :: in
real(rprec), dimension(nxnpy, nyt), intent(out) :: out
! ---
real(rprec), dimension(nxnpy, nynpy, 0:npy-1) :: work1
real(rprec), dimension(nxnpy, nynpy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nxnpy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(st:en,:)
    st = en + 1
    en = st + nxnpy - 1
end do

num = nxnpy * nynpy
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nynpy
do ipy = 0, npy-1
    out(:, st:en) = work2(:,:,ipy)
    st = en + 1
    en = st + nynpy - 1
end do
! ---
end subroutine transpose_phy_y_to_x

