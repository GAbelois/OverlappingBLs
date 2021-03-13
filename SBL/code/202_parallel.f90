!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_3d_sp(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nxt, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(nxt, nynpy, 0:nz), intent(in)  :: a_local
real, dimension(nxt, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(nxt, nynpy, 0:nz, 0:np-1))
call mpi_gather(a_local, nxt*nynpy*(nz+1), mpi_rprec,   &
                temp,    nxt*nynpy*(nz+1), mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz-1
    do jy = 1, nynpy
    do jx = 1, nxt
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = real(temp(jx, jy, jz, ip)) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data_3d_sp
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
! subroutine collect_data_xy_sp(a_local, a_global)
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param, only: nxt, nynpy, nyt, npy, idy, mpi_rprec, comm_ver, ierr
! implicit none
! ! ---
! real(kind=rprec), dimension(nxt, nynpy), intent(in)  :: a_local
! real, dimension(nxt, nyt),  intent(out) :: a_global
! ! ---
! real(kind=rprec), dimension(:, :, :), allocatable :: temp
! integer :: ip, ipy, ipz, jx, jy, jz
! ! ---
! allocate(temp(nxt, nynpy, 0:npy-1))
! call mpi_gather(a_local, nxt*nynpy, mpi_rprec,   &
!                 temp,    nxt*nynpy, mpi_rprec,   &
!                 0, comm_ver, ierr)

! do ip = 0, npy - 1
!     ipy = idy(ip)
!     do jy = 1, nynpy
!     do jx = 1, nxt
!         a_global(jx,ipy*nynpy+jy) = real(temp(jx, jy, ip)) 
!     end do
!     end do
! end do

! deallocate(temp)
! ! ---
! end subroutine collect_data_xy_sp
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data2_3d_sp(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, nz), intent(in)  :: a_local
real, dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, nz, 0:np-1))
call mpi_gather(a_local, ldx*nynpy*nz, mpi_rprec,   &
                temp,    ldx*nynpy*nz, mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz-1
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data2_3d_sp
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(in)  :: a_local
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, 0:nz, 0:np-1))
call mpi_gather(a_local, ldx*nynpy*(nz+1), mpi_rprec,   &
                temp,    ldx*nynpy*(nz+1), mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz-1
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data2_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, nz), intent(in)  :: a_local
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, nz, 0:np-1))
call mpi_gather(a_local, ldx*nynpy*nz, mpi_rprec,   &
                temp,    ldx*nynpy*nz, mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz-1
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data2_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data3_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: nxt, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(nxt, nynpy, nz), intent(in)  :: a_local
real(kind=rprec), dimension(nxt, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(nxt, nynpy, nz, 0:np-1))
call mpi_gather(a_local, nxt*nynpy*nz, mpi_rprec,   &
                temp,    nxt*nynpy*nz, mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz-1
    do jy = 1, nynpy
    do jx = 1, nxt
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data3_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data4_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                    idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(in)  :: a_local
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, 0:nynpy+1, 0:nz, 0:np-1))
call mpi_gather(a_local, ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                temp,    ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz-1
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data4_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
! subroutine collect_data4_3d_sp(a_local, a_global)
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
!                  idy, idz, mpi_rprec, comm, ierr
! implicit none
! ! ---
! real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(in)  :: a_local
! real, dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ! ---
! real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
! integer :: ip, ipy, ipz, jx, jy, jz
! ! ---
! allocate(temp(ldx, 0:nynpy+1, 0:nz, 0:np-1))
! call mpi_gather(a_local, ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
!                 temp,    ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
!                 0, comm, ierr)

! do ip = 0, np - 1
!     ipy = idy(ip)
!     ipz = idz(ip)
!     do jz = 1, nz-1
!     do jy = 1, nynpy
!     do jx = 1, ldx
!         a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
!     end do
!     end do
!     end do
! end do

! deallocate(temp)
! ! ---
! end subroutine collect_data4_3d_sp
! !--------------------------------------------------------------------!
! !                                              
! !--------------------------------------------------------------------!
subroutine collect_data_z(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: nz, nzt, npz, idz_col, mpi_rprec, comm_col, ierr
implicit none
! ---
real(kind=rprec), dimension(nz-1),  intent(in)  :: a_local
real(kind=rprec), dimension(nzt-1), intent(out) :: a_global
! ---
real(kind=rprec), dimension(nz-1, 0:npz-1) :: temp
integer :: ip, ipz, jz
! ---
call mpi_allgather(a_local, nz-1, mpi_rprec, temp, nz-1, mpi_rprec, comm_col, ierr)

do ipz = 0, npz - 1
    ip = idz_col(ipz)
    do jz = 1, nz - 1
        a_global(ip*(nz-1)+jz) = temp(jz, ipz)
    end do
end do
! ---
end subroutine collect_data_z
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_z_sp(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: nz, nzt, npz, idz_col, mpi_rprec, comm_col, ierr
implicit none
! ---
real(kind=rprec), dimension(nz-1),  intent(in)  :: a_local
real, dimension(nzt-1), intent(out) :: a_global
! ---
real(kind=rprec), dimension(nz-1, 0:npz-1) :: temp
integer :: ip, ipz, jz
! ---
call mpi_allgather(a_local, nz-1, mpi_rprec, temp, nz-1, mpi_rprec, comm_col, ierr)

do ipz = 0, npz - 1
    ip = idz_col(ipz)
    do jz = 1, nz - 1
        a_global(ip*(nz-1)+jz) = real(temp(jz, ipz))
    end do
end do
! ---
end subroutine collect_data_z_sp
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine scatter_data(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, 0:nz, 0:np-1))

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*nynpy*(nz+1), mpi_rprec,   &
                 a_local, ldx*nynpy*(nz+1), mpi_rprec,   &
                 0, comm, ierr)

deallocate(temp)
! ---
end subroutine scatter_data
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine scatter_data2(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                    idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, nynpy, nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, nz, 0:np-1))

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*nynpy*nz, mpi_rprec,   &
                 a_local, ldx*nynpy*nz, mpi_rprec,   &
                 0, comm, ierr)

deallocate(temp)
! ---
end subroutine scatter_data2
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine scatter_data3(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real, dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, 0:nz, 0:np-1))

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*nynpy*(nz+1), mpi_rprec,   &
                 a_local, ldx*nynpy*(nz+1), mpi_rprec,   &
                 0, comm, ierr)

deallocate(temp)
! ---
end subroutine scatter_data3
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine scatter_data4(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param
implicit none
! ---
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
real(kind=rprec), dimension(:,:), allocatable :: temp_xz0, temp_xz1, temp_xz2, temp_xz3
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, 0:nynpy+1, 0:nz, 0:np-1))
allocate(temp_xz0(ldx, 0:nz), temp_xz1(ldx, 0:nz), temp_xz2(ldx, 0:nz), temp_xz3(ldx, 0:nz))
temp = 0._rprec
temp_xz0 = 0._rprec
temp_xz1 = 0._rprec
temp_xz2 = 0._rprec
temp_xz3 = 0._rprec

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                 a_local, ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                 0, comm, ierr)

temp_xz0(:, :) = a_local(:, 1, :)
call mpi_sendrecv(temp_xz0(1, 1), ldx*(nz+1), MPI_RPREC, south, 1,       &
                  temp_xz1(1, 1), ldx*(nz+1), MPI_RPREC, north, 1,       &
                  comm, status, ierr)
a_local(:,nynpy+1,:) = temp_xz1(:, :)


temp_xz2(:, :) = a_local(:, nynpy, :)
call mpi_sendrecv(temp_xz2(1, 1), ldx*(nz+1), MPI_RPREC, south, 2,       &
                  temp_xz3(1, 1), ldx*(nz+1), MPI_RPREC, north, 2,       &
                  comm, status, ierr)
a_local(:,0,:) = temp_xz3(:, :)                  
                  
deallocate(temp, temp_xz0, temp_xz1, temp_xz2, temp_xz3)
! ---
end subroutine scatter_data4
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_xh(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxhnpy, lhx, nyt, npy, coordz, idy_ver,  &
                 comm_ver, mpi_cprec, ierr
implicit none
! ---
complex(rprec), dimension(nxhnpy, nyt), intent(in) :: in
complex(rprec), dimension(lhx, nyt),   intent(out) :: out
! ---
complex(rprec), dimension(nxhnpy, nyt, 0:npy-1) :: tmp
integer :: ip, ipy, jy, num
! ---
num = nxhnpy * nyt
call mpi_allgather(in, num, mpi_cprec, tmp, num, mpi_cprec, comm_ver, ierr)

do ipy = 0, npy-1
    ip = idy_ver(ipy)
    do jy = 1, nxhnpy
        out(ip*(nxhnpy-1)+jy, :) = tmp(jy, :, ipy)
    end do
end do
! ---
end subroutine collect_data_xh
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine collect_data_y(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, npy, coordz, idy_ver, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nynpy), intent(in) :: in
real(rprec), dimension(nyt),  intent(out) :: out
! ---
real(rprec), dimension(nynpy, 0:npy-1) :: tmp
integer :: ip, ipy, jy
! ---
call mpi_allgather(in, nynpy, mpi_rprec, tmp, nynpy, mpi_rprec, comm_ver, ierr)

do ipy = 0, npy-1
    ip = idy_ver(ipy)
    do jy = 1, nynpy
        out(ip*nynpy+jy) = tmp(jy, ipy)
    end do
end do
! ---
end subroutine collect_data_y 
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine collect_data_y_sp(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, npy, coordz, idy_ver, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nynpy), intent(in) :: in
real, dimension(nyt),  intent(out) :: out
! ---
real(rprec), dimension(nynpy, 0:npy-1) :: tmp
integer :: ip, ipy, jy
! ---
call mpi_allgather(in, nynpy, mpi_rprec, tmp, nynpy, mpi_rprec, comm_ver, ierr)

do ipy = 0, npy-1
    ip = idy_ver(ipy)
    do jy = 1, nynpy
        out(ip*nynpy+jy) = real(tmp(jy, ipy))
    end do
end do
! ---
end subroutine collect_data_y_sp 
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine scatter_y(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, npy, idy_ver, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nyt),   intent(in)  :: in
real(rprec), dimension(nynpy), intent(out) :: out
! ---
real(rprec), dimension(nynpy, 0:npy-1) :: tmp
integer :: ipy, j, ind
! ---
do ipy = 0, npy-1
    do j = 1, nynpy
        !ind = ipy*nynpy + j
        !tmp(j,rank_ver(ipy)) = in(ind)
        ind = idy_ver(ipy)*nynpy + j
        tmp(j,ipy) = in(ind)    
    end do
end do

call mpi_scatter(tmp, nynpy, mpi_rprec,     &
                 out, nynpy, mpi_rprec, 0, comm_ver, ierr)
! ---
end subroutine scatter_y
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine scatter_z(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nz, nzt, npz, idz_col, comm_col, mpi_rprec, ierr
implicit none
! ---
real, dimension(nzt-1),   intent(in)  :: in
real(rprec), dimension(nz-1), intent(out) :: out
! ---
real(rprec), dimension(nz-1, 0:npz-1) :: tmp
integer :: ipz, k, ind
! ---
do ipz = 0, npz-1
    do k = 1, nz-1
        ind = idz_col(ipz)*(nz-1) + k
        tmp(k,ipz) = in(ind)    
    end do
end do

call mpi_scatter(tmp, nz-1, mpi_rprec,     &
                 out, nz-1, mpi_rprec, 0, comm_col, ierr)
! ---
end subroutine scatter_z
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine data_exchange()
!-----------------------------------------------------------------------
!   exchange data between neighboring zones
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
integer :: ipcon, jx, jy, jz, tag
real :: temp_w
! --- synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz
call mpi_sendrecv(u(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   1,  &
                  u(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 1,  &
                  comm, status, ierr)
call mpi_sendrecv(v(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   2,  &
                  v(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 2,  &
                  comm, status, ierr)
call mpi_sendrecv(w(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   3,  &
                  w(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 3,  &
                  comm, status, ierr)
call mpi_sendrecv(u(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 4,  &
                  u(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   4,  &
                  comm, status, ierr)
call mpi_sendrecv(v(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 5,  &
                  v(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   5,  &
                  comm, status, ierr)
call mpi_sendrecv(w(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 6,  &
                  w(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   6,  &
                  comm, status, ierr)
if (theta_flag) then
    call mpi_sendrecv(theta(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   7,  &
                      theta(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 7,  &
                      comm, status, ierr)
    call mpi_sendrecv(theta(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 8,  &
                      theta(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   8,  &
                      comm, status, ierr)
end if

if (PCon_FLAG) then
    tag = 9
    do ipcon = 1, npcon
        call mpi_sendrecv(PCon(1, 0, nz - 1, ipcon), ldx*(nynpy+2), MPI_RPREC, up,   tag,    &
                          PCon(1, 0, 0, ipcon),      ldx*(nynpy+2), MPI_RPREC, down, tag,    &
                          comm, status, ierr)
        tag = tag + 1
    end do
end if

if ( coordz == 0 ) then
    !--set 0-level velocities to BOGUS
    u(:, :, 0) = BOGUS
    v(:, :, 0) = BOGUS
    w(:, :, 0) = BOGUS
    if (theta_flag) theta(:, :, 0) = BOGUS
else if ( coordz == npz - 1 ) then
    u(:, :, nz) = BOGUS
    v(:, :, nz) = BOGUS
    !w(:, :, nz) = BOGUS
    if (theta_flag) theta(:, :, nz) = BOGUS
end if
! ---
end subroutine data_exchange
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine data_exchange2(a_local)
!-----------------------------------------------------------------------
!   exchange data for PCon
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:,:), allocatable :: temp_xz0, temp_xz1, temp_xz2, temp_xz3

integer :: ipcon, jx, jy, jz
! ---
allocate(temp_xz0(ldx, 0:nz), temp_xz1(ldx, 0:nz), temp_xz2(ldx, 0:nz), temp_xz3(ldx, 0:nz))
temp_xz0 = 0._rprec
temp_xz1 = 0._rprec
temp_xz2 = 0._rprec
temp_xz3 = 0._rprec

! ---
temp_xz0(:, :) = a_local(:, 1, :)
call mpi_sendrecv(temp_xz0(1, 0), ldx*(nz+1), MPI_RPREC, south, 1,       &
                  temp_xz1(1, 0), ldx*(nz+1), MPI_RPREC, north, 1,       &
                  comm, status, ierr)
a_local(:,nynpy+1,:) = temp_xz1(:, :)


temp_xz2(:, :) = a_local(:, nynpy, :)
call mpi_sendrecv(temp_xz2(1, 0), ldx*(nz+1), MPI_RPREC, north, 2,       &
                  temp_xz3(1, 0), ldx*(nz+1), MPI_RPREC, south, 2,       &
                  comm, status, ierr)
a_local(:,0,:) = temp_xz3(:, :) 

! --- synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz
call mpi_sendrecv(a_local(1, 0, nz - 1), ldx*(nynpy+2), MPI_RPREC, up,   3,    &
                  a_local(1, 0, 0),      ldx*(nynpy+2), MPI_RPREC, down, 3,    &
                  comm, status, ierr)
call mpi_sendrecv(a_local(1, 0, 1),      ldx*(nynpy+2), MPI_RPREC, down, 4,    &
                  a_local(1, 0, nz),     ldx*(nynpy+2), MPI_RPREC, up,   4,    &
                  comm, status, ierr)
! ---
deallocate(temp_xz0, temp_xz1, temp_xz2, temp_xz3)
! ---
end subroutine data_exchange2
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine data_exchange_y(a_local, y_per)
!-----------------------------------------------------------------------
!   exchange data for PCon
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
real(kind=rprec), dimension(ldx, 0:nynpy+1, nz), intent(out) :: a_local
logical, intent(in) :: y_per 
! ---
real(kind=rprec), dimension(ldx,nz) :: temp0, temp1, temp0m, temp1m
integer :: sou, nor
! ---
temp0  = 0._rprec
temp1  = 0._rprec
temp0m = 0._rprec
temp1m = 0._rprec

! ---
temp1(:, :)  = a_local(:, 1, :)
temp1m(:, :) = a_local(:, nynpy, :)

if ( idy(rank) > 0 ) then
    sou = rank - npz
else 
    sou = mpi_proc_null
end if

if ( idy(rank) < npy-1 ) then
    nor = rank + npz
else
    nor = mpi_proc_null
end if

call mpi_sendrecv(temp1(:,:),  ldx*nz, MPI_RPREC, sou,  1,  &
                  temp0m(:,:), ldx*nz, MPI_RPREC, nor,  1,  &
                  comm, status, ierr)
call mpi_sendrecv(temp1m(:,:), ldx*nz, MPI_RPREC, nor,  2,  &
                  temp0(:,:),  ldx*nz, MPI_RPREC, sou,  2,  &
                  comm, status, ierr)
! --- periodic condition
if ( y_per ) then
    if ( idy(rank) == 0 ) then
        sou = rank + (npy-1)*npz
    else
        sou = mpi_proc_null
    end if

    if ( idy(rank) == npy-1 ) then
        nor = rank - (npy-1)*npz
    else
        nor = mpi_proc_null
    end if

    call mpi_sendrecv(temp1(:,:),  ldx*nz, MPI_RPREC, sou,  1,  &
                      temp0m(:,:), ldx*nz, MPI_RPREC, nor,  1,  &
                      comm, status, ierr)
    call mpi_sendrecv(temp1m(:,:), ldx*nz, MPI_RPREC, nor,  2,  &
                      temp0(:,:),  ldx*nz, MPI_RPREC, sou,  2,  &
                      comm, status, ierr)
end if

if ( idy(rank) > 0 .or. y_per ) then
    a_local(:,0,:) = temp0(:, :)
end if

if ( idy(rank) < npy-1 .or. y_per ) then
    a_local(:,nynpy+1,:) = temp0m(:, :)
end if

! ---
end subroutine data_exchange_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine data_exchange_y2(a_local, y_per)
!-----------------------------------------------------------------------
!   exchange data for PCon
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
real(kind=rprec), dimension(nxt+2, 1:nynpy+2, nz+2), intent(out) :: a_local
logical, intent(in) :: y_per 
! ---
real(kind=rprec), dimension(nxt+2, nz+2) :: temp0, temp1, temp0m, temp1m
integer :: sou, nor
! ---
temp0  = 0._rprec
temp1  = 0._rprec
temp0m = 0._rprec
temp1m = 0._rprec

! ---
temp1(:, :)  = a_local(:, 2, :)
temp1m(:, :) = a_local(:, nynpy+1, :)

if ( idy(rank) > 0 ) then
    sou = rank - npz
else 
    sou = mpi_proc_null
end if

if ( idy(rank) < npy-1 ) then
    nor = rank + npz
else
    nor = mpi_proc_null
end if

call mpi_sendrecv(temp1(:,:),  (nxt+2)*(nz+2), MPI_RPREC, sou,  1,  &
                  temp0m(:,:), (nxt+2)*(nz+2), MPI_RPREC, nor,  1,  &
                  comm, status, ierr)
call mpi_sendrecv(temp1m(:,:), (nxt+2)*(nz+2), MPI_RPREC, nor,  2,  &
                  temp0(:,:),  (nxt+2)*(nz+2), MPI_RPREC, sou,  2,  &
                  comm, status, ierr)
! --- periodic condition
if ( y_per ) then
    if ( idy(rank) == 0 ) then
        sou = rank + (npy-1)*npz
    else
        sou = mpi_proc_null
    end if

    if ( idy(rank) == npy-1 ) then
        nor = rank - (npy-1)*npz
    else
        nor = mpi_proc_null
    end if

    call mpi_sendrecv(temp1(:,:),  (nxt+2)*(nz+2), MPI_RPREC, sou,  1,  &
                      temp0m(:,:), (nxt+2)*(nz+2), MPI_RPREC, nor,  1,  &
                      comm, status, ierr)
    call mpi_sendrecv(temp1m(:,:), (nxt+2)*(nz+2), MPI_RPREC, nor,  2,  &
                      temp0(:,:),  (nxt+2)*(nz+2), MPI_RPREC, sou,  2,  &
                      comm, status, ierr)
end if

if ( idy(rank) > 0 .or. y_per ) then
    a_local(:,1,:) = temp0(:, :)
end if

if ( idy(rank) < npy-1 .or. y_per ) then
    a_local(:,nynpy+2,:) = temp0m(:, :)
end if

! ---
end subroutine data_exchange_y2
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine avg_y(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nynpy), intent(in)  :: in
real(rprec), intent(out) :: out
! ---
real(rprec) :: tmp
integer :: ipy, j, ind
! ---
tmp = 0._rprec
do j = 1, nynpy
    tmp = tmp + in(j)
end do
call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)

out = out / nyt
! ---
end subroutine avg_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine avg_y3d(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxt, nynpy, nyt, nz, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nxt, nynpy, nz-1), intent(in)  :: in
real(rprec), dimension(nxt, nz-1), intent(out) :: out
! ---
real(rprec), dimension(nxt, nz-1) :: tmp
integer :: jx, jz
! ---
tmp = 0._rprec
do jz = 1, nz-1
do jx = 1, nxt
    tmp(jx,jz) = sum(in(jx,1:nynpy,jz)) / nyt
end do
end do
call mpi_allreduce(tmp, out, nxt*(nz-1), mpi_rprec, mpi_sum, comm_ver, ierr)
! ---
end subroutine avg_y3d
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine calc_hor_avg(in, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: ldx, nynpy, nz, nxt, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(ldx, nynpy, 0:nz), intent(in)  :: in
real(rprec), dimension(0:nz-1), intent(out) :: out
! ---
real(rprec), dimension(0:nz-1) :: tmp
integer :: jz
! ---
do jz = 0, nz-1
    tmp(jz) = sum(in(1:nxt, 1:nynpy, jz)) / (nxt * nyt)
end do
call mpi_allreduce(tmp, out, nz, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_hor_avg
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------     
subroutine calc_hor_avg2(in, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: ldx, nynpy, nz, nxt, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(ldx, nynpy, nz), intent(in)  :: in
real(rprec), dimension(nz-1), intent(out) :: out
! ---
real(rprec), dimension(nz-1) :: tmp
integer :: jz
! ---
do jz = 1, nz-1
    tmp(jz) = sum (in(1:nxt, 1:nynpy, jz)) / (nxt * nyt)
end do
call mpi_allreduce(tmp, out, nz-1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_hor_avg2
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine calc_hor_avg3(in, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: ldx, nynpy, nz, nxt, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nxt, nynpy), intent(in)  :: in
real(rprec), intent(out) :: out
! ---
real(rprec) :: tmp
integer :: jz
! ---
tmp = sum(in(1:nxt, 1:nynpy)) / (nxt * nyt)
call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_hor_avg3
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine calc_rms(in1, in2, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: nxt, nynpy, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nxt, nynpy), intent(in)  :: in1, in2
real(rprec), intent(out) :: out
! ---
real(rprec) :: tmp
! ---
tmp = sum (in1*in2) / (nxt * nyt)
call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_rms
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
subroutine write_file_mpi_3d(var, fn, disp0) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!   12/30/2018 -- Bicheng Chen (chabby@ucla.edu) -- First created
!   Modified by Chao Yan for 2d parallelization
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nynpy, nz, np, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nxt, 1:nynpy, 1:nz-1), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
integer :: bufsize
! integer(MPI_OFFSET_KIND) :: disp
! ---
bufsize = size(var)
! disp = np * sizeof(var)
! ---
call mpi_type_create_subarray(3, gsize, lsize, start, MPI_ORDER_FORTRAN,  &
                              mpi_rprec, filetype, ierr)
call mpi_type_commit(filetype, ierr)

call mpi_file_open(comm, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE,   &
                   MPI_INFO_NULL, fid, ierr)          
call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
call mpi_file_write_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
call mpi_file_close(fid, ierr)
! ---
end subroutine write_file_mpi_3d
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine write_file_mpi_3d_sp(var, fn, disp0) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nynpy, nz, np, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real, dimension(1:nxt, 1:nynpy, 1:nz-1), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
integer :: bufsize
! integer(MPI_OFFSET_KIND) :: disp
! ---
bufsize = size(var)
! disp = np * sizeof(var)
! ---
call mpi_type_create_subarray(3, gsize, lsize, start, MPI_ORDER_FORTRAN,  &
                              mpi_real, filetype, ierr)
call mpi_type_commit(filetype, ierr)

call mpi_file_open(comm, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fid, ierr)          
call mpi_file_set_view(fid, disp0, mpi_real, filetype, 'native', MPI_INFO_NULL, ierr)
call mpi_file_write_all(fid, var, bufsize, mpi_real, MPI_STATUS_IGNORE, ierr)
call mpi_file_close(fid, ierr)
! ---
end subroutine write_file_mpi_3d_sp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine write_file_mpi_3dFraction_sp(var, fn, disp0, buffersize) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!-----------------------------------------------------------------------
use types, only: rprec
use param !, only: nxt, nynpy, nz, np, mpi_rprec, comm, ierr, gsize_opt, lsize_opt, start_opt
use mpi
implicit none
! ---
real, dimension(1:nxout, 1:nynpy, 1:nz-1), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
integer, intent(in) :: buffersize
! ---
integer :: fid, filetype
! ---
call mpi_type_create_subarray(3, gsize_opt, lsize_opt, start_opt, MPI_ORDER_FORTRAN,  &
                              mpi_real, filetype, ierr)
call mpi_type_commit(filetype, ierr)
call mpi_file_open(comm, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fid, ierr)          
call mpi_file_set_view(fid, disp0, mpi_real, filetype, 'native', MPI_INFO_NULL, ierr)
call mpi_file_write_all(fid, var, buffersize, mpi_real, MPI_STATUS_IGNORE, ierr)
call mpi_file_close(fid, ierr)
! ---
end subroutine write_file_mpi_3dFraction_sp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine write_file_mpi_3dFraction(var, fn, disp0, buffersize) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!-----------------------------------------------------------------------
use types, only: rprec
use param !, only: nxt, nynpy, nz, np, mpi_rprec, comm, ierr, gsize_opt, lsize_opt, start_opt
use mpi
implicit none
! ---
real(rprec), dimension(1:nxout, 1:nynpy, 1:nz-1), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
integer, intent(in) :: buffersize
! ---
integer :: fid, filetype
! ---
call mpi_type_create_subarray(3, gsize_opt, lsize_opt, start_opt, MPI_ORDER_FORTRAN,  &
                              mpi_rprec, filetype, ierr)
call mpi_type_commit(filetype, ierr)
call mpi_file_open(comm, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fid, ierr)          
call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
call mpi_file_write_all(fid, var, buffersize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
call mpi_file_close(fid, ierr)
! ---
end subroutine write_file_mpi_3dFraction
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_file_mpi_3d(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nynpy, nz, np, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nxt, 1:nynpy, 1:nz-1), intent(in) :: var
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
    call mpi_type_create_subarray(3, gsize, lsize, start, MPI_ORDER_FORTRAN,  &
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
end subroutine read_file_mpi_3d
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine write_file_mpi_xy(var, fn, disp0) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!   12/30/2018 -- Bicheng Chen (chabby@ucla.edu) -- First created
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nynpy, nz, np, comm_ver, coordz, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nxt, 1:nynpy), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
integer :: bufsize
! ---
bufsize = size(var)
! ---
if ( coordz == 0 ) then
    call mpi_type_create_subarray(2, gsize(1:2), lsize(1:2), start(1:2), MPI_ORDER_FORTRAN,  &
                                  mpi_rprec, filetype, ierr)
    call mpi_type_commit(filetype, ierr)

    call mpi_file_open(comm_ver, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE,   &
                       MPI_INFO_NULL, fid, ierr)
    call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)   
    call mpi_file_write_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
    call mpi_file_close(fid, ierr)
end if
! ---
end subroutine write_file_mpi_xy
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine write_file_mpi_yz(var, fn, disp0) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!   12/30/2018 -- Bicheng Chen (chabby@ucla.edu) -- First created
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nynpy, nz, np, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nynpy, 1:nz-1), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
integer :: bufsize
! ---
bufsize = size(var)
! ---
call mpi_type_create_subarray(2, gsize(2:3), lsize(2:3), start(2:3), MPI_ORDER_FORTRAN,  &
                                mpi_rprec, filetype, ierr)
call mpi_type_commit(filetype, ierr)

call mpi_file_open(comm, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE,   &
                    MPI_INFO_NULL, fid, ierr)
call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)   
call mpi_file_write_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
call mpi_file_close(fid, ierr)
! ---
end subroutine write_file_mpi_yz
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine write_file_mpi_fringe(var, fn, disp0) 
!-----------------------------------------------------------------------
!   Write one 3d variable to binary file
!   12/30/2018 -- Bicheng Chen (chabby@ucla.edu) -- First created
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:
use param, only: jx_relax, nxt, nynpy, nz, np, mpi_rprec, comm, ierr, &
                 gsize_fringe, lsize_fringe, start_fringe
use mpi
implicit none
! ---
real(rprec), dimension(1:jx_relax, 1:nynpy, 1:nz-1), intent(in) :: var
character(*), intent(in) :: fn
integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ---
integer :: fid, filetype
integer :: bufsize
! ---
bufsize = size(var)
! ---
call mpi_type_create_subarray(3, gsize_fringe, lsize_fringe, start_fringe, MPI_ORDER_FORTRAN,  &
                              mpi_rprec, filetype, ierr)
call mpi_type_commit(filetype, ierr)

call mpi_file_open(comm, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE,   &
                   MPI_INFO_NULL, fid, ierr)
call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)   
call mpi_file_write_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
call mpi_file_close(fid, ierr)
! ---
end subroutine write_file_mpi_fringe
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
! subroutine write_file_mpi_xy_sp(var, fn, disp0) 
! !-----------------------------------------------------------------------
! !   Write one 3d variable to binary file
! !   12/30/2018 -- Bicheng Chen (chabby@ucla.edu) -- First created
! !-----------------------------------------------------------------------
! use types, only: rprec
! use io, only:gsize, lsize, start
! use param, only: nxt, nynpy, nz, np, comm_ver, coordz, mpi_rprec, comm, ierr
! use mpi
! implicit none
! ! ---
! real, dimension(1:nxt, 1:nynpy), intent(in) :: var
! character(*), intent(in) :: fn
! integer(MPI_OFFSET_KIND), intent(in) :: disp0
! ! ---
! integer :: fid, filetype
! integer :: bufsize
! ! ---
! bufsize = size(var)
! ! ---
! if ( coordz == 0 ) then
!     call mpi_type_create_subarray(2, gsize(1:2), lsize(1:2), start(1:2), MPI_ORDER_FORTRAN,  &
!                                   mpi_real, filetype, ierr)
!     call mpi_type_commit(filetype, ierr)

!     call mpi_file_open(comm_ver, fn, MPI_MODE_WRONLY+MPI_MODE_CREATE,   &
!                        MPI_INFO_NULL, fid, ierr)
!     call mpi_file_set_view(fid, disp0, mpi_real, filetype, 'native', MPI_INFO_NULL, ierr)   
!     call mpi_file_write_all(fid, var, bufsize, mpi_real, MPI_STATUS_IGNORE, ierr)
!     call mpi_file_close(fid, ierr)
! end if
! ! ---
! end subroutine write_file_mpi_xy_sp
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
subroutine read_file_mpi_3d_inflow(var, fn, disp0)
!-----------------------------------------------------------------------
!   inflow version
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nyt, nynpy, nz, nzt, np, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nxt/2, 1:nynpy, 1:nz-1), intent(in) :: var
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
    call mpi_type_create_subarray(3, (/nxt/2, nyt, nzt-1/), (/nxt/2, nynpy, nz-1/), start,   &
                                  MPI_ORDER_FORTRAN, mpi_rprec, filetype, ierr)
    call mpi_type_commit(filetype, ierr)

    call mpi_file_open(comm, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
    call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
    call mpi_file_close(fid, ierr)
else
    write(*,*) "INFLOW FILE DOES NOT EXIST"
    stop
end if
! ---
end subroutine read_file_mpi_3d_inflow
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_file_mpi_xy(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nynpy, nz, np, comm_ver, coordz, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nxt, 1:nynpy), intent(in) :: var
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
        call mpi_type_create_subarray(2, gsize(1:2), lsize(1:2), start(1:2), MPI_ORDER_FORTRAN,  &
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
end subroutine read_file_mpi_xy
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_file_mpi_xy_inflow(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nxt, nynpy, nyt, nz, np, comm_ver, coordz, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nxt/2, 1:nynpy), intent(in) :: var
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
        call mpi_type_create_subarray(2, (/nxt/2, nyt/), (/nxt/2, nynpy/), start(1:2),   &
                                      MPI_ORDER_FORTRAN, mpi_rprec, filetype, ierr)
        call mpi_type_commit(filetype, ierr)

        call mpi_file_open(comm_ver, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
        call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
        call mpi_file_read_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
        call mpi_file_close(fid, ierr)
    end if
else
    write(*,*) "INFLOW FILE DOES NOT EXIST"
    stop
end if
! ---
end subroutine read_file_mpi_xy_inflow
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_file_mpi_yz(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use io, only:gsize, lsize, start
use param, only: nynpy, nz, np, mpi_rprec, comm, ierr
use mpi
implicit none
! ---
real(rprec), dimension(1:nynpy, 1:nz-1), intent(in) :: var
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
    call mpi_type_create_subarray(2, gsize(2:3), lsize(2:3), start(2:3), MPI_ORDER_FORTRAN,  &
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
end subroutine read_file_mpi_yz
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_file_mpi_fringe(var, fn, disp0)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: jx_relax, nynpy, nz, np, mpi_rprec, comm, ierr, &
                 gsize_fringe, lsize_fringe, start_fringe
use mpi
implicit none
! ---
real(rprec), dimension(1:jx_relax, 1:nynpy, 1:nz-1), intent(in) :: var
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
    call mpi_type_create_subarray(3, gsize_fringe, lsize_fringe, start_fringe, MPI_ORDER_FORTRAN,  &
                                    mpi_rprec, filetype, ierr)
    call mpi_type_commit(filetype, ierr)

    call mpi_file_open(comm, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
    call mpi_file_set_view(fid, disp0, mpi_rprec, filetype, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fid, var, bufsize, mpi_rprec, MPI_STATUS_IGNORE, ierr)
    call mpi_file_close(fid, ierr)
else
    write(*,*) "FRINGE FILE DOES NOT EXIST"
    stop
end if
! ---
end subroutine read_file_mpi_fringe