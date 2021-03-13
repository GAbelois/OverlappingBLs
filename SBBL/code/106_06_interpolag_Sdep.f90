!-----------------------------------------------------------------------
! this is the w-node version
!-----------------------------------------------------------------------
subroutine interpolag_Sdep()
!-----------------------------------------------------------------------
! This subroutine computes the values of F_LM and F_MM
! at positions (x-u*dt) for use in the subroutine lagrangian
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sgsmodule
use stokes_drift, only: ust, vst
implicit none
! ---  
integer :: jx, jy, jz, jjx, jjy, jjz, jz_min, jz_max

!--saves here: force heap storage
real(rprec) :: frac_x, frac_y, frac_z
real(rprec) :: comp_x, comp_y, comp_z

integer :: addx, addy, addz
integer :: jxaddx, jyaddy, jzaddz

if (VERBOSE) write (*, *) 'started interpolag_Sdep'

! --- creates dummy arrays FF_LM and FF_MM to use in the subroutine

do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    FF_LM(jx + 1, jy + 1, jz + 1) = F_LM(jx, jy, jz)
    FF_MM(jx + 1, jy + 1, jz + 1) = F_MM(jx, jy, jz)
    FF_QN(jx + 1, jy + 1, jz + 1) = F_QN(jx, jy, jz)
    FF_NN(jx + 1, jy + 1, jz + 1) = F_NN(jx, jy, jz)
end do
end do
end do

!This is a bit like witch craft but the following lines do take care
!of all the edges including the corners
FF_LM(1, :, :) = FF_LM(nxt + 1, :, :)
FF_LM(nxt + 2, :, :) = FF_LM(2, :, :)

! FF_LM(:, 1, :) = FF_LM(:, nynpy + 1, :)
! FF_LM(:, nynpy + 2, :) = FF_LM(:, 2, :)
call data_exchange_y2(FF_LM, .true.)

!--send F_LM @ jz=nz-1 to F_LM @ jz=0'
!  i.e. FF_LM @ jz=nz to FF_LM @ jz=1'
call mpi_sendrecv(FF_LM(1, 1, nz), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   1, &
                  FF_LM(1, 1, 1),  (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 1, &
                  comm, status, ierr)
!--F_LM @ jz=nz and F_LM @ jz=1' should already be in sync (test?)
!  i.e. FF_LM @ jz=nz+1 and FF_LM @ jz=2'
  
if ( coordz == 0 ) then
    FF_LM(:, :, 1) = FF_LM(:, :, 2)
end if

!--send F_LM @ jz=2 to F_LM @ jz=nz+1'
!  i.e. FF_LM @ jz=3 to FF_LM @ jz=nz+2'
! call mpi_sendrecv(FF_LM(1, 1, 2),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 2, &
!                   FF_LM(1, 1, nz + 1), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   2, &
!                   comm, status, ierr)
call mpi_sendrecv(FF_LM(1, 1, 3),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 2, &
                  FF_LM(1, 1, nz + 2), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   2, &
                  comm, status, ierr)
  
if ( coordz == npz - 1 ) then
    ! FF_LM(:, :, nz + 1) = FF_LM(:, :, nz)
    FF_LM(:, :, nz + 2) = FF_LM(:, :, nz + 1)
end if

FF_MM(1, :, :) = FF_MM(nxt + 1, :, :)
FF_MM(nxt + 2, :, :) = FF_MM(2, :, :)

! FF_MM(:, 1, :) = FF_MM(:, nynpy + 1, :)
! FF_MM(:, nynpy + 2, :) = FF_MM(:, 2, :)
call data_exchange_y2(FF_MM, .true.)

!--send F_MM @ jz=nz-1 to F_MM @ jz=0'
!  i.e. FF_MM @ jz=nz to FF_MM @ jz=1'
call mpi_sendrecv(FF_MM(1, 1, nz), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   3, &
                  FF_MM(1, 1, 1),  (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 3, &
                  comm, status, ierr)
!--F_MM @ jz=nz and F_MM @ jz=1' should already be in sync (test?)
!  i.e. FF_MM @ jz=nz+1 and FF_MM @ jz=2'
  
if ( coordz == 0 ) then
    FF_MM(:, :, 1) = FF_MM(:, :, 2)
end if

!--send F_MM @ jz=2 to F_MM @ jz=nz+1'
!  i.e. FF_MM @ jz=3 to FF_MM @ jz=nz+2'
! call mpi_sendrecv(FF_MM(1, 1, 2),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 4, &
!                   FF_MM(1, 1, nz + 1), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   4, &
!                   comm, status, ierr)
call mpi_sendrecv(FF_MM(1, 1, 3),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 4, &
                  FF_MM(1, 1, nz + 2), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   4, &
                  comm, status, ierr)

if ( coordz == npz - 1 ) then
    ! FF_MM(:, :, nz + 1) = FF_MM(:, :, nz)
    FF_MM(:, :, nz + 2) = FF_MM(:, :, nz + 1)
end if

FF_QN(1, :, :) = FF_QN(nxt + 1, :, :)
FF_QN(nxt + 2, :, :) = FF_QN(2, :, :)

! FF_QN(:, 1, :) = FF_QN(:, nynpy + 1, :)
! FF_QN(:, nynpy + 2, :) = FF_QN(:, 2, :)
call data_exchange_y2(FF_QN, .true.)

!--send F_QN @ jz=nz-1 to F_QN @ jz=0'
!  i.e. FF_QN @ jz=nz to FF_QN @ jz=1'
call mpi_sendrecv(FF_QN(1, 1, nz), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   5, &
                  FF_QN(1, 1, 1),  (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 5, &
                  comm, status, ierr)
!--F_QN @ jz=nz and F_QN @ jz=1' should already be in sync (test?)
!  i.e. FF_QN @ jz=nz+1 and FF_QN @ jz=2'

if ( coordz == 0 ) then
    FF_QN(:, :, 1) = FF_QN(:, :, 2)
end if

!--send F_QN @ jz=2 to F_QN @ jz=nz+1'
!  i.e. FF_QN @ jz=3 to FF_QN @ jz=nz+2'
! call mpi_sendrecv(FF_QN(1, 1, 2),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 6, &
!                   FF_QN(1, 1, nz + 1), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   6, &
!                   comm, status, ierr)
call mpi_sendrecv(FF_QN(1, 1, 3),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 6, &
                  FF_QN(1, 1, nz + 2), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   6, &
                  comm, status, ierr)

if ( coordz == npz - 1 ) then
    ! FF_QN(:, :, nz + 1) = FF_QN(:, :, nz)
    FF_QN(:, :, nz + 2) = FF_QN(:, :, nz + 1)
end if

FF_NN(1, :, :) = FF_NN(nxt + 1, :, :)
FF_NN(nxt + 2, :, :) = FF_NN(2, :, :)

! FF_NN(:, 1, :) = FF_NN(:, nynpy + 1, :)
! FF_NN(:, nynpy + 2, :) = FF_NN(:, 2, :)
call data_exchange_y2(FF_NN, .true.)

!--send F_NN @ jz=nz-1 to F_NN @ jz=0'
!  i.e. FF_NN @ jz=nz to FF_NN @ jz=1'
call mpi_sendrecv(FF_NN(1, 1, nz), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   7, &
                  FF_NN(1, 1, 1),  (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 7, &
                  comm, status, ierr)
!--F_NN @ jz=nz and F_NN @ jz=1' should already be in sync (test?)
!  i.e. FF_NN @ jz=nz+1 and FF_NN @ jz=2'

if ( coordz == 0 ) then
    FF_NN(:, :, 1) = FF_NN(:, :, 2)
end if

!--send F_NN @ jz=2 to F_NN @ jz=nz+1'
!  i.e. FF_NN @ jz=3 to FF_NN @ jz=nz+2'
! call mpi_sendrecv(FF_NN(1, 1, 2),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 8, &
!                   FF_NN(1, 1, nz + 1), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   8, &
!                   comm, status, ierr)
call mpi_sendrecv(FF_NN(1, 1, 3),      (nxt + 2)*(nynpy + 2), MPI_RPREC, down, 8, &
                  FF_NN(1, 1, nz + 2), (nxt + 2)*(nynpy + 2), MPI_RPREC, up,   8, &
                  comm, status, ierr)

if ( coordz == npz - 1 ) then
    ! FF_NN(:, :, nz + 1) = FF_NN(:, :, nz)
    FF_NN(:, :, nz + 2) = FF_NN(:, :, nz + 1)
end if
! end of witch craft
!        puts u_lag and and v_lag on cs nodes
u_temp = u_lag/real(cs_count, kind=rprec)
v_temp = v_lag/real(cs_count, kind=rprec)

!--not sure if u_lag, u_temp at 0 are needed yet
u_lag(:, :, 0) = BOGUS
v_lag(:, :, 0) = BOGUS

if ( coordz == 0 ) then
    u_lag(:, :, 1) = u_temp(:, :, 1)
    v_lag(:, :, 1) = v_temp(:, :, 1)

    u_lag(:, :, 1) = u_lag(:, :, 1) + ust(1)
    v_lag(:, :, 1) = v_lag(:, :, 1) + vst(1)
    
    jz_min = 2
else
    jz_min = 1
end if

if ( coordz == npz - 1 ) then

    u_lag(:, :, nz) = u_temp(:, :, nz-1)
    v_lag(:, :, nz) = v_temp(:, :, nz-1)

    u_lag(:, :, nz) = u_lag(:, :, nz-1) + ust(nz-1)
    v_lag(:, :, nz) = v_lag(:, :, nz-1) + vst(nz-1)

    jz_max = nz - 1
else
    jz_max = nz
end if


do jz = jz_min, jz_max
    u_lag(:, :, jz) = 0.5_rprec*(u_temp(:, :, jz) + u_temp(:, :, jz - 1))
    v_lag(:, :, jz) = 0.5_rprec*(v_temp(:, :, jz) + v_temp(:, :, jz - 1))

    u_lag(:, :, jz) = u_lag(:, :, jz) + 0.5_rprec*(ust(jz) + ust(jz - 1))
    v_lag(:, :, jz) = v_lag(:, :, jz) + 0.5_rprec*(vst(jz) + vst(jz - 1))
end do

w_lag = w_lag/real(cs_count, kind=rprec)

!--not sure if 0-level used
w_lag(:, :, 0) = BOGUS

if ( coordz == 0 ) then
    w_lag(:, :, 1) = 0.25_rprec*w_lag(:, :, 2)
end if

if ( coordz == npz - 1 ) then
    w_lag(:, :, nz) = 0.25_rprec*w_lag(:, :, nz-1)
end if

!     computes the 3-D inverse displacement arrays that describe
!     the location where the point was at the previous step
xp = -u_lag*dt*real(cs_count, kind=rprec)/dx !! is u already normalized
yp = -v_lag*dt*real(cs_count, kind=rprec)/dy
zp = -w_lag*dt*real(cs_count, kind=rprec)/dz

!--perhaps can remove 0-level altogether
xp(:, :, 0) = BOGUS
yp(:, :, 0) = BOGUS
zp(:, :, 0) = BOGUS

if ( coordz == 0 ) then
!        because first plane is on u,v,p nodes
!        this corrects for the fact that the first cell
!        in the z direction has height dz/2
!        it doubles the zp fraction if this fraction relates to the cell
!     beneath it
    do jx = 1, nxt
    do jy = 1, nynpy
        zp(jx, jy, 2) = zp(jx, jy, 2) + min(zp(jx, jy, 2), 0._rprec)
        zp(jx, jy, 1) = zp(jx, jy, 1) + max(zp(jx, jy, 1), 0._rprec)
    end do
    end do

    zp(:, :, 2) = min(1._rprec, zp(:, :, 2))
    zp(:, :, 1) = min(1._rprec, zp(:, :, 2))
    zp(:, :, 2) = max(-1._rprec, zp(:, :, 2))
    zp(:, :, 1) = max(-1._rprec, zp(:, :, 2))

end if

! if ( rank == 0 ) then
!     if (mod(jt, p_count) .eq. 0) print *, 'Lagrangian CFL condition= ', &
!           maxval(abs(xp(1:nxt, :, 1:nz)))
! end if

do jz = 1, nz
    jjz = jz + 1
    do jy = 1, nynpy
        jjy = jy + 1
        do jx = 1, nxt
            jjx = jx + 1
!           the are the values to add to the indices jx, jy, and jz
!           the are +1 or -1 depending on what cube should be used for interpolation
            addx = int(sign(1._rprec, xp(jx, jy, jz)))
            addy = int(sign(1._rprec, yp(jx, jy, jz)))
            addz = int(sign(1._rprec, zp(jx, jy, jz)))
            jxaddx = jjx + addx
            jyaddy = jjy + addy
            jzaddz = jjz + addz
!           computes the relative weights given to F_** in the cube depending on point location
            comp_x = abs(xp(jx, jy, jz))
            comp_y = abs(yp(jx, jy, jz))
            comp_z = abs(zp(jx, jy, jz))
            frac_x = 1._rprec - comp_x
            frac_y = 1._rprec - comp_y
            frac_z = 1._rprec - comp_z

!     computes interpolated F_LM

            F_LM(jx, jy, jz) = frac_x*frac_y* &
           (FF_LM(jjx, jjy, jjz)*frac_z + FF_LM(jjx, jjy, jzaddz)*comp_z) &
           + frac_x*comp_y* &
           (FF_LM(jjx, jyaddy, jjz)*frac_z + FF_LM(jjx, jyaddy, jzaddz)*comp_z) &
           + comp_x*frac_y* &
           (FF_LM(jxaddx, jjy, jjz)*frac_z + FF_LM(jxaddx, jjy, jzaddz)*comp_z) &
           + comp_x*comp_y* &
           (FF_LM(jxaddx, jyaddy, jjz)*frac_z + FF_LM(jxaddx, jyaddy, jzaddz)*comp_z)

!     computes interpolated F_MM
          F_MM(jx, jy, jz) = frac_x*frac_y* &
           (FF_MM(jjx, jjy, jjz)*frac_z + FF_MM(jjx, jjy, jzaddz)*comp_z) &
           + frac_x*comp_y* &
           (FF_MM(jjx, jyaddy, jjz)*frac_z + FF_MM(jjx, jyaddy, jzaddz)*comp_z) &
           + comp_x*frac_y* &
           (FF_MM(jxaddx, jjy, jjz)*frac_z + FF_MM(jxaddx, jjy, jzaddz)*comp_z) &
           + comp_x*comp_y* &
           (FF_MM(jxaddx, jyaddy, jjz)*frac_z + FF_MM(jxaddx, jyaddy, jzaddz)*comp_z)

!     computes interpolated F_QN
          F_QN(jx, jy, jz) = frac_x*frac_y* &
           (FF_QN(jjx, jjy, jjz)*frac_z + FF_QN(jjx, jjy, jzaddz)*comp_z) &
           + frac_x*comp_y* &
           (FF_QN(jjx, jyaddy, jjz)*frac_z + FF_QN(jjx, jyaddy, jzaddz)*comp_z) &
           + comp_x*frac_y* &
           (FF_QN(jxaddx, jjy, jjz)*frac_z + FF_QN(jxaddx, jjy, jzaddz)*comp_z) &
           + comp_x*comp_y* &
           (FF_QN(jxaddx, jyaddy, jjz)*frac_z + FF_QN(jxaddx, jyaddy, jzaddz)*comp_z)

!     computes interpolated F_NN
          F_NN(jx, jy, jz) = frac_x*frac_y* &
           (FF_NN(jjx, jjy, jjz)*frac_z + FF_NN(jjx, jjy, jzaddz)*comp_z) &
           + frac_x*comp_y* &
           (FF_NN(jjx, jyaddy, jjz)*frac_z + FF_NN(jjx, jyaddy, jzaddz)*comp_z) &
           + comp_x*frac_y* &
           (FF_NN(jxaddx, jjy, jjz)*frac_z + FF_NN(jxaddx, jjy, jzaddz)*comp_z) &
           + comp_x*comp_y* &
           (FF_NN(jxaddx, jyaddy, jjz)*frac_z + FF_NN(jxaddx, jyaddy, jzaddz)*comp_z)
        end do
    end do
end do

u_lag = 0._rprec
v_lag = 0._rprec
w_lag = 0._rprec

if (VERBOSE) write (*, *) 'finished interpolag_Sdep'
! ---
end subroutine interpolag_Sdep
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------