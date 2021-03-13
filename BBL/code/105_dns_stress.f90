!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!  
subroutine dns_stress(txx, txy, txz, tyy, tyz, tzz)
!--------------------------------------------------------------------!
! using the 'sgs' sign convention for stress, so there is a - sign                                             
!--------------------------------------------------------------------!  
use types, only: rprec
use param, only: ldx, nxt, nynpy, nz, z_i, u_scale, nu_molec, &
                 coordz, npz, BOGUS
use sim_param, only: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
implicit none
! ---  
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: txx, txy, txz, tyy, tyz, tzz

real(kind=rprec) :: S11, S12, S13, S22, S23, S33
real(kind=rprec) :: nu
integer :: jx, jy, jz, jz_min

! non-dimensional molecular viscosity
nu = nu_molec/(z_i*u_scale)

! uvp-nodes
do jz = 1, nz - 1
do jy = 1, nynpy
do jx = 1, nxt
    S11 = dudx(jx, jy, jz)
    S12 = 0.5_rprec*(dudy(jx, jy, jz) + dvdx(jx, jy, jz))
    S22 = dvdy(jx, jy, jz)
    S33 = dwdz(jx, jy, jz)
    txx(jx, jy, jz) = -2._rprec*nu*S11
    txy(jx, jy, jz) = -2._rprec*nu*S12
    tyy(jx, jy, jz) = -2._rprec*nu*S22
    tzz(jx, jy, jz) = -2._rprec*nu*S33
end do
end do
end do

!--if values not needed, set to bogus value (easier to catch errors)
if ( coordz == npz - 1 ) then
    ! top values of txx, txy, tyy, tzz not needed for stress free bc's
    !--if values not needed, set to bogus value (easier to catch errors)
    txx(:, :, nz) = BOGUS
    txy(:, :, nz) = BOGUS
    tyy(:, :, nz) = BOGUS
    tzz(:, :, nz) = BOGUS
end if

! w-nodes
if ( coordz == 0 ) then
    ! leave the wall level alone: taken care of with wall stress
    !--assume here that wall stress has already been set (for MPI)
    jz_min = 2
else
    jz_min = 1
end if

do jz = jz_min, nz - 1
do jy = 1, nynpy
do jx = 1, nxt
    S13 = 0.5_rprec*(dudz(jx, jy, jz) + dwdx(jx, jy, jz))
    S23 = 0.5_rprec*(dvdz(jx, jy, jz) + dwdy(jx, jy, jz))
    txz(jx, jy, jz) = -2._rprec*nu*S13
    tyz(jx, jy, jz) = -2._rprec*nu*S23
end do
end do
end do

if ( coordz < npz - 1 ) then
    !--nz here saves communication in MPI version: can only do this since
    !  dudz, dwdx, dvdz, dwdy are available at nz (not true w/ all components)
    do jy = 1, nynpy
    do jx = 1, nxt
        S13 = 0.5_rprec*(dudz(jx, jy, nz) + dwdx(jx, jy, nz))
        S23 = 0.5_rprec*(dvdz(jx, jy, nz) + dwdy(jx, jy, nz))
        txz(jx, jy, nz) = -2._rprec*nu*S13
        tyz(jx, jy, nz) = -2._rprec*nu*S23
    end do
    end do
end if
! ---
end subroutine dns_stress
