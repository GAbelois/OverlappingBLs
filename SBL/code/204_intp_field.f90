!-----------------------------------------------------------------------
!   Interpolate first in z, then in x and y
!-----------------------------------------------------------------------
subroutine create_fft_plan_intp()
!-----------------------------------------------------------------------
!   create fft plans for spectral interpolation   
!-----------------------------------------------------------------------
use types, only:rprec
use param, only:nxt, nyt
use intp
implicit none
include 'fftw3.f'
! ---
real(rprec), dimension(nxt/2) :: arr_nxt
real(rprec), dimension(nyt/2) :: arr_nyt
real(rprec), dimension(nxt)   :: arr_nxt2
real(rprec), dimension(nyt)   :: arr_nyt2
complex(rprec), dimension(nxt/2/2+1) :: arr_lhx
complex(rprec), dimension(nyt/2/2+1) :: arr_lhy
complex(rprec), dimension(nxt/2+1)   :: arr_lhx2
complex(rprec), dimension(nyt/2+1)   :: arr_lhy2 
! ---
call dfftw_plan_dft_r2c_1d(plan_xf_intp, nxt/2, arr_nxt, arr_lhx, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(plan_yf_intp, nyt/2, arr_nyt, arr_lhy, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(plan_x2b_intp,  nxt, arr_lhx2, arr_nxt2, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(plan_y2b_intp,  nyt, arr_lhy2, arr_nyt2, FFTW_PATIENT, FFTW_UNALIGNED)
! ---
end subroutine create_fft_plan_intp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_coarse_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use intp
use io, only:fcumulative_time
implicit none
! ---
character(80) :: fname
logical :: exst
integer(MPI_OFFSET_KIND) :: disp0, disp3d 
! ---
allocate(u_tmp(nxt/2, nynpy/2, 0:(nz-1)/2+1), v_tmp(nxt/2, nynpy/2, 0:(nz-1)/2+1), &
         w_tmp(nxt/2, nynpy/2, 1:(nz-1)/2+1), theta_tmp(nxt/2, nynpy/2, 0:(nz-1)/2+1))
! ---
gsize2(1) = nxt/2
gsize2(2) = nyt/2
gsize2(3) = (nzt-1)/2
lsize2(1) = nxt/2
lsize2(2) = nynpy/2
lsize2(3) = (nz-1)/2
start2(1) = 0
start2(2) = nynpy/2 * coordy
start2(3) = (nz-1)/2 * coordz
! ---
inquire (file=fcumulative_time, exist=exst)
if (exst) then
    open(1, file=fcumulative_time)
    read(1, *) nums
    close(1)
else
    write(*, *) "Restart file NOT found"
    nums = 0
end if

disp3d = np  * sizeof(u_tmp(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2))
write (fname, '(a,i8.8,a)') '../input/vel_tt.out'

disp0 = 0
call read_file_mpi_3d_intp(u_tmp(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_intp(v_tmp(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_intp(w_tmp(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2), fname, disp0)

write (fname, '(a,i8.8,a)') '../input/temp_tt.out'
disp0  = 0
disp3d = np  * sizeof(theta_tmp(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2))
call read_file_mpi_3d_intp(theta_tmp(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2), fname, disp0)
! ---
end subroutine read_coarse_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine intp_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use sgsmodule
use io, only:fcumulative_time
use intp
implicit none
! ---
real(rprec), dimension(nxt/2, nynpy/2, nz-1) :: u_tmp2, v_tmp2, w_tmp2, theta_tmp2          
integer :: jz
! ---
call z_intp_uv(u_tmp,     u_tmp2)
call z_intp_uv(v_tmp,     v_tmp2)
call z_intp_w( w_tmp,     w_tmp2)
call z_intp_uv(theta_tmp, theta_tmp2)
! ---                  
do jz = 1, nz-1
    call zero_padding_intp(u_tmp2(:,:,jz),    u(1:nxt,1:nynpy,jz))
    call zero_padding_intp(v_tmp2(:,:,jz),    v(1:nxt,1:nynpy,jz))
    call zero_padding_intp(w_tmp2(:,:,jz),    w(1:nxt,1:nynpy,jz))
    call zero_padding_intp(theta_tmp2(:,:,jz), theta(1:nxt, 1:nynpy,jz))
end do
! ---
end subroutine intp_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine z_intp_uv(var, var_fine)
!-----------------------------------------------------------------------
!                                              
!-----------------------------------------------------------------------
use types, only: rprec
use param
use fft 
implicit none
! --- 
real(kind=rprec), dimension(nxt/2, nynpy/2, 0:(nz-1)/2+1), intent(in)  :: var
real(kind=rprec), dimension(nxt/2, nynpy/2, nz-1),     intent(out) :: var_fine
integer :: jz
! ---
call mpi_sendrecv(var(1, 1, (nz-1)/2),  nxt/2*nynpy/2, MPI_RPREC, up,   1,  &
                  var(1, 1, 0),         nxt/2*nynpy/2, MPI_RPREC, down, 1,  &
                  comm, status, ierr)
call mpi_sendrecv(var(1, 1, 1),          nxt/2*nynpy/2, MPI_RPREC, down, 2,  &
                  var(1, 1, (nz-1)/2+1), nxt/2*nynpy/2, MPI_RPREC, up,   2,  &
                  comm, status, ierr)                  

if ( coordz == 0 ) then
    var_fine(:, :, 1) = var(:, :, 1)
    ! do jz = 1, (nz-1)/2-1
    !     var_fine(:, :, 2*jz)   = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz+1)
    !     var_fine(:, :, 2*jz+1) = 0.25_rprec*var(:, :, jz) + 0.75_rprec*var(:, :, jz+1)
    ! end do
    ! var_fine(:, :, nz-1) = 0.75_rprec*var(:, :, (nz-1)/2) + 0.25_rprec*var(:, :, (nz-1)/2+1)
    var_fine(:, :, 2) = 0.75_rprec*var(:, :, 1) + 0.25_rprec*var(:, :, 2)
    do jz = 2, (nz-1)/2
        var_fine(:, :, 2*jz)   = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz+1)
        var_fine(:, :, 2*jz-1) = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz-1) 
    end do
else if ( coordz == npz-1 ) then
    do jz = 1, (nz-1)/2-1
        var_fine(:, :, 2*jz)   = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz+1)
        var_fine(:, :, 2*jz-1) = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz-1)
    end do
    var_fine(:, :, nz-2) = 0.75_rprec*var(:, :, (nz-1)/2) + 0.25_rprec*var(:, :, (nz-1)/2-1) 
    var_fine(:, :, nz-1) = var(:, :, (nz-1)/2)
else 
    do jz = 1, (nz-1)/2
        var_fine(:, :, 2*jz)   = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz+1)
        var_fine(:, :, 2*jz-1) = 0.75_rprec*var(:, :, jz) + 0.25_rprec*var(:, :, jz-1)
    end do
end if
! ---
end subroutine z_intp_uv
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine z_intp_w(var, var_fine)
!-----------------------------------------------------------------------
!                                               
!-----------------------------------------------------------------------
use types, only: rprec
use param
use fft 
implicit none
! --- 
real(kind=rprec), dimension(nxt/2, nynpy/2, (nz-1)/2+1), intent(in)  :: var
real(kind=rprec), dimension(nxt/2, nynpy/2, nz-1),       intent(out) :: var_fine
integer :: jz
! ---
call mpi_sendrecv(var(1, 1, 1),          nxt/2*nynpy/2, MPI_RPREC, down, 3,  &
                  var(1, 1, (nz-1)/2+1), nxt/2*nynpy/2, MPI_RPREC, up,   3,  &
                  comm, status, ierr)

if ( coordz == npz-1 ) then
    var_fine(:, :, nz-1) = var(:, :, (nz-1)/2)
    var_fine(:, :, nz-2) = var(:, :, (nz-1)/2)
    do jz = 1, (nz-1)/2-1
        var_fine(:, :, 2*jz-1) = var(:, :, jz)
        var_fine(:, :, 2*jz)   = (var(:, :, jz)+var(:, :, jz+1))/2._rprec
    end do
else
    do jz = 1, (nz-1)/2
        var_fine(:, :, 2*jz-1) = var(:, :, jz)
        var_fine(:, :, 2*jz)   = (var(:, :, jz)+var(:, :, jz+1))/2._rprec
    end do
end if
! ---
end subroutine z_intp_w
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine zero_padding_intp(var, var_big)
!-----------------------------------------------------------------------
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs                                              
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: nxt, nynpy, nyt, npy
use intp 
implicit none
! --- note we're calling with 2D arrays
real(kind=rprec), dimension(nxt/2, nynpy/2), intent(in)  :: var
real(kind=rprec), dimension(nxt, nynpy),     intent(out) :: var_big
! ---
complex(kind=rprec), dimension(nxt/2/2+1, nynpy/2) :: var_hat_x
complex(kind=rprec), dimension(nxt/npy, nyt/2/2+1) :: var_hat_y
complex(kind=rprec), dimension(nxt/2+1, nynpy/2)   :: var_hat_x_big
complex(kind=rprec), dimension(nxt/npy, nyt/2+1)   :: var_hat_y_big
real(kind=rprec), dimension(nxt, nynpy/2)   :: var_x_big
real(kind=rprec), dimension(nxt/npy, nyt/2) :: var_y
real(kind=rprec), dimension(nxt/npy, nyt) :: var_y_big

real(rprec), dimension(nxt/2, nynpy/2) :: tmp

integer :: jx, jy
! --- make sure the big array is zeroed!
var_hat_x_big = 0._rprec
var_hat_y_big = 0._rprec

tmp(:, :) = var(:, :) / (nxt/2)

do jy = 1, nynpy/2
    call dfftw_execute_dft_r2c(plan_xf_intp, tmp(:,jy), var_hat_x(:,jy))       ! 1D real to complex FFT in x direction
end do

do jx = 1, nxt/2/2
    var_hat_x_big(jx, :) = var_hat_x(jx, :)
end do

do jy = 1, nynpy/2
    call dfftw_execute_dft_c2r(plan_x2b_intp, var_hat_x_big(:,jy), var_x_big(:,jy))
end do

call transpose_y_to_x2_intp(var_x_big, var_y)

var_y = var_y / (nyt/2)
do jx = 1, nxt/npy
    call dfftw_execute_dft_r2c(plan_yf_intp, var_y(jx,:), var_hat_y(jx,:))
end do

do jy = 1, nyt/2/2
    var_hat_y_big(:, jy) = var_hat_y(:, jy)
end do

do jx = 1, nxt/npy
    call dfftw_execute_dft_c2r(plan_y2b_intp, var_hat_y_big(jx,:), var_y_big(jx,:)) 
end do

call transpose_x2_to_y_intp(var_y_big, var_big)

var_big(:, :) = var_big(:, :)
! ---
end subroutine zero_padding_intp
!-----------------------------------------------------------------------
!                                              
!-----------------------------------------------------------------------
subroutine transpose_y_to_x2_intp(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, nynpy, nyt, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxt, nynpy/2),   intent(in)  :: in
real(rprec), dimension(nxt/npy, nyt/2), intent(out) :: out
! ---
real(rprec), dimension(nxt/npy, nynpy/2, 0:npy-1) :: work1
real(rprec), dimension(nxt/npy, nynpy/2, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nxt/npy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(st:en, :)
    st = en + 1
    en = st + nxt/npy - 1
end do

num = nxt/npy * nynpy/2
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nynpy/2
do ipy = 0, npy-1
    out(:, st:en) = work2(:,:,ipy)
    st = en + 1
    en = st + nynpy/2 - 1
end do
! ---
end subroutine transpose_y_to_x2_intp
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine transpose_x2_to_y_intp(in, out)
!-----------------------------------------------------------------------
!   transpose x pencil to y pencil
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, nyt, nynpy, npy, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nxt/npy, nyt),  intent(in)  :: in
real(rprec), dimension(nxt, nynpy), intent(out) :: out
! ---
real(rprec), dimension(nxt/npy, nynpy, 0:npy-1) :: work1
real(rprec), dimension(nxt/npy, nynpy, 0:npy-1) :: work2
integer :: ipy, st, en, num
! ---
st = 1
en = nynpy
do ipy = 0, npy-1
    work1(:,:,ipy) = in(:, st:en)
    st = en + 1
    en = st + nynpy - 1
end do

num = nxt/npy * nynpy
call mpi_alltoall(work1, num, mpi_rprec, work2, num, mpi_rprec, comm_ver, ierr)

st = 1
en = nxt/npy
do ipy = 0, npy-1
    out(st:en, :) = work2(:,:,ipy)
    st = en + 1
    en = st + nxt/npy - 1
end do
! ---
end subroutine transpose_x2_to_y_intp
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine read_file_mpi_3d_intp(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: nxt, nyt, nynpy, nzt, nz, np, mpi_rprec, comm, ierr
use mpi
use intp, only:gsize2, lsize2, start2
implicit none
! ---
real(rprec), dimension(1:nxt/2, 1:nynpy/2, 1:(nz-1)/2), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
logical :: exst
integer :: bufsize
! integer(MPI_OFFSET_KIND) :: disp
! ---
! disp = coordy * sizeof(var)
bufsize = size(var)
! ---
inquire (file=fn, exist=exst)
if (exst) then
    call mpi_type_create_subarray(3, gsize2, lsize2, start2, MPI_ORDER_FORTRAN,  &
                                  mpi_rprec, filetype, ierr)
    call mpi_type_commit(filetype, ierr)

    call mpi_file_open(comm, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
    call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
    call mpi_file_close(fid, ierr)
else
    write(*,*) "FILE DOES NOT EXIST"
    stop
end if
! ---
end subroutine read_file_mpi_3d_intp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_file_mpi_2d_intp(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: nxt, nyt, nynpy, nzt, nz, np, mpi_rprec, coordz, comm_ver, comm, ierr
use mpi
use intp, only:gsize2, lsize2, start2
implicit none
! ---
real(rprec), dimension(1:nxt/2, 1:nynpy/2), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
logical :: exst
integer :: bufsize
! ---
bufsize = size(var)
! ---
inquire (file=fn, exist=exst)
if (exst) then
    if ( coordz == 0 ) then
        call mpi_type_create_subarray(2, gsize2(1:2), lsize2(1:2), start2(1:2), MPI_ORDER_FORTRAN,  &
                                      mpi_rprec, filetype, ierr)
        call mpi_type_commit(filetype, ierr)

        call mpi_file_open(comm_ver, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
        call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
        call mpi_file_read_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
        call mpi_file_close(fid, ierr)
    end if
else
    write(*,*) "FILE DOES NOT EXIST"
    stop
end if
! ---
end subroutine read_file_mpi_2d_intp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
! subroutine collect_data_intp1(a_local, a_global)
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param, only: nxt, nynpy, nyt, nz, nzt, np,    &
!                  idy, idz, mpi_rprec, comm, ierr
! implicit none
! ! ---
! real(kind=rprec), dimension(nxt/2, nynpy/2, (nz-1)/2), intent(in) :: a_local
! real(kind=rprec), dimension(nxt/2, nyt/2, (nzt-1)/2),  intent(out) :: a_global
! ! ---
! real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
! integer :: ip, ipy, ipz, jx, jy, jz
! ! ---
! allocate(temp(nxt/2, nynpy/2, (nz-1)/2, 0:np-1))
! call mpi_gather(a_local, nxt/2*nynpy/2*(nz-1)/2, mpi_rprec,   &
!                 temp,    nxt/2*nynpy/2*(nz-1)/2, mpi_rprec,   &
!                 0, comm, ierr)

! do ip = 0, np - 1
!     ipy = idy(ip)
!     ipz = idz(ip)
!     do jz = 1, (nz-1)/2
!     do jy = 1, nynpy/2
!     do jx = 1, nxt/2
!         a_global(jx,ipy*nynpy/2+jy,ipz*(nz-1)/2+jz) = temp(jx, jy, jz, ip) 
!     end do
!     end do
!     end do
! end do

! deallocate(temp)
! ! ---
! end subroutine collect_data_intp1
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! subroutine collect_data_intp2(a_local, a_global)
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param, only: nxt, nynpy, nyt, nz, nzt, np,    &
!                  idy, idz, mpi_rprec, comm, ierr
! implicit none
! ! ---
! real(kind=rprec), dimension(nxt, nynpy, (nz-1)/2), intent(in)  :: a_local
! real(kind=rprec), dimension(nxt, nyt, (nzt-1)/2),  intent(out) :: a_global
! ! ---
! real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
! integer :: ip, ipy, ipz, jx, jy, jz
! ! ---
! allocate(temp(nxt, nynpy, (nz-1)/2, 0:np-1))
! call mpi_gather(a_local, nxt*nynpy*(nz-1)/2, mpi_rprec,   &
!                 temp,    nxt*nynpy*(nz-1)/2, mpi_rprec,   &
!                 0, comm, ierr)

! do ip = 0, np - 1
!     ipy = idy(ip)
!     ipz = idz(ip)
!     do jz = 1, (nz-1)/2
!     do jy = 1, nynpy
!     do jx = 1, nxt
!         a_global(jx,ipy*nynpy+jy,ipz*(nz-1)/2+jz) = temp(jx, jy, jz, ip) 
!     end do
!     end do
!     end do
! end do

! deallocate(temp)
! ! ---
! end subroutine collect_data_intp2
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! subroutine collect_data_intp3(a_local, a_global)
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param, only: nxt, nynpy, nyt, nz, nzt, np,    &
!                  idy, idz, mpi_rprec, comm, ierr
! implicit none
! ! ---
! real(kind=rprec), dimension(nxt, nynpy, nz-1), intent(in)  :: a_local
! real(kind=rprec), dimension(nxt, nyt, nzt-1),  intent(out) :: a_global
! ! ---
! real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
! integer :: ip, ipy, ipz, jx, jy, jz
! ! ---
! allocate(temp(nxt, nynpy, nz-1, 0:np-1))
! call mpi_gather(a_local, nxt*nynpy*(nz-1), mpi_rprec,   &
!                 temp,    nxt*nynpy*(nz-1), mpi_rprec,   &
!                 0, comm, ierr)

! do ip = 0, np - 1
!     ipy = idy(ip)
!     ipz = idz(ip)
!     do jz = 1, nz-1
!     do jy = 1, nynpy
!     do jx = 1, nxt
!         a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
!     end do
!     end do
!     end do
! end do

! deallocate(temp)
! ! ---
! end subroutine collect_data_intp3
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! subroutine check_field_intp1(var)
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
! use types, only: rprec
! use param
! use sim_param, only: u, v, w, p, theta, PCon
! use stokes_drift, only: ust, vst
! implicit none
! ! ---
! real(rprec), dimension(nxt/2, nynpy/2, (nz-1)/2), intent(in) :: var
! real(rprec), dimension(nxt/2, nyt/2, (nzt-1)/2) :: var_tot

! character(80) :: fname
! integer :: jx, jy, jz, ipcon
! real(rprec) :: xx, yy, zz
! ! ---
! call collect_data_intp1(var, var_tot)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/z/var_1D_check_.csv'
!     open(9,file=fname)
!     write(9, '(3a)') 'var'
!     do jz = 1, (nzt-1)/2
!        !zz = (jz-1) * z_i / (nzt-1)        
!         write(9, '(e15.7)') sum(var_tot(:, :, jz))/float(nxt/2*nyt/2)                     
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_field_intp1
! !--------------------------------------------------------------------!
! !                                               
! !--------------------------------------------------------------------!
! subroutine check_field_intp2(var)
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
! use types, only: rprec
! use param
! use sim_param, only: u, v, w, p, theta, PCon
! use stokes_drift, only: ust, vst
! implicit none
! ! ---
! real(rprec), dimension(nxt, nynpy, (nz-1)/2), intent(in) :: var
! real(rprec), dimension(nxt, nyt, (nzt-1)/2) :: var_tot

! character(80) :: fname
! integer :: jz
! ! ---
! call collect_data_intp2(var, var_tot)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/z/var2_1D_check_.csv'
!     open(9,file=fname)
!     write(9, '(11a)') 'var' 
!     do jz = 1, (nzt-1)/2
!         write(9, '(e15.7)') sum(var_tot(:, :, jz))/float(nxt*nyt)
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_field_intp2
! !--------------------------------------------------------------------!
! !                                               
! !--------------------------------------------------------------------!
! subroutine check_field_intp3(var)
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
! use types, only: rprec
! use param
! use sim_param, only: u, v, w, p, theta, PCon
! use stokes_drift, only: ust, vst
! implicit none
! ! ---
! real(rprec), dimension(nxt, nynpy, nz-1), intent(in) :: var
! real(rprec), dimension(nxt, nyt, nzt-1) :: var_tot

! character(80) :: fname
! integer :: jz
! ! ---
! call collect_data_intp3(var, var_tot)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/z/var3_1D_check_.csv'
!     open(9,file=fname)
!     write(9, '(11a)') 'var' 
!     do jz = 1, nzt-1
!         write(9, '(e15.7)') sum(var_tot(:, :, jz))/float(nxt*nyt)
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_field_intp3
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
! subroutine calc_hor_avg_intp(in, out)    
! !-----------------------------------------------------------------------
! !   calculate the horizontal-average value of the 3d field 
! !-----------------------------------------------------------------------     
! use types, only: rprec
! use param, only: nxt, nyt, nynpy, comm_ver, mpi_rprec, ierr
! use mpi
! implicit none
! ! ---
! real(rprec), dimension(nxt/2, nynpy/2), intent(in)  :: in
! real(rprec), intent(out) :: out
! ! ---
! real(rprec) :: tmp
! integer :: jz
! ! ---
! tmp = sum(in(1:nxt/2, 1:nynpy/2)) / (nxt/2 * nyt/2)
! call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ! ---    
! end subroutine calc_hor_avg_intp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
! subroutine collect_data_z_intp(a_local, a_global)
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param, only: nz, nzt, npz, idz_col, mpi_rprec, comm_col, ierr
! implicit none
! ! ---
! real(kind=rprec), dimension((nz-1)/2),  intent(in)  :: a_local
! real(kind=rprec), dimension((nzt-1)/2), intent(out) :: a_global
! ! ---
! real(kind=rprec), dimension((nz-1)/2, 0:npz-1) :: temp
! integer :: ip, ipz, jz
! ! ---
! call mpi_allgather(a_local, (nz-1)/2, mpi_rprec, temp, (nz-1)/2, mpi_rprec, comm_col, ierr)

! do ipz = 0, npz - 1
!     ip = idz_col(ipz)
!     do jz = 1, (nz-1)/2
!         a_global(ip*(nz-1)/2+jz) = temp(jz, ipz)
!     end do
! end do
! ! ---
! end subroutine collect_data_z_intp
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!