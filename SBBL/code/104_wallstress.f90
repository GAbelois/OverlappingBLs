!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!  
subroutine wallstress_dns_bottom()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only:rprec
use param
use sim_param, only:u, v, dudz, dvdz, txz, tyz
use stokes_drift, only: deg2rad
implicit none
! ---
integer :: jx, jy
! ---
if ( ocean_flag ) then
    if (flag_dynStress) then
        ! Add dynamical wind stress (Bicheng Chen 10/24/2016)
        open(unit=fid_ustar, file=fn_ustar, access='direct',recl=len_ustar)
        read(unit=fid_ustar, rec=nums+ttshift_dynStress+1) ustar_dyn, agl_stress
        close(fid_ustar)
        rad_stress = deg2rad(agl_stress)
        ustar_dyn = ustar_dyn / u_scale
        txz(1:nxt, :, 1) = ustar_dyn**2 * cos(rad_stress)
        tyz(1:nxt, :, 1) = ustar_dyn**2 * sin(rad_stress)
    else
        txz(1:nxt, :, 1) = u_star**2 * cos(rad_stress)   ! !DY Everything has been normalized by u_scale
        tyz(1:nxt, :, 1) = u_star**2 * sin(rad_stress)
    end if

    dudz(1:nxt, :, 1) = (dudz(1:nxt, :, 2)*cos(rad_stress) +  &
                        dvdz(1:nxt, :, 2)*sin(rad_stress)) * cos(rad_stress)
    dvdz(1:nxt, :, 1) = (dudz(1:nxt, :, 2)*cos(rad_stress) +  &
                        dvdz(1:nxt, :, 2)*sin(rad_stress)) * sin(rad_stress)
else 
    select case (lbc_mom)

    case ('wall')
        do jy = 1, nynpy
        do jx = 1, nxt
            txz(jx, jy, 1) = -nu_molec/(z_i*u_scale)*u(jx, jy, 1)/(0.5_rprec*dz)
            tyz(jx, jy, 1) = -nu_molec/(z_i*u_scale)*v(jx, jy, 1)/(0.5_rprec*dz)
            dudz(jx, jy, 1) = u(jx, jy, 1)/(0.5_rprec*dz)
            dvdz(jx, jy, 1) = v(jx, jy, 1)/(0.5_rprec*dz)
        end do
        end do
    case ('stress free')
        txz(:, :, 1)  = 0._rprec
        tyz(:, :, 1)  = 0._rprec
        dudz(:, :, 1) = 0._rprec
        dvdz(:, :, 1) = 0._rprec
    case default
        write (*, *) 'invalid lbc_mom'
        stop
    end select
end if
! ---
end subroutine wallstress_dns_bottom
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!  
subroutine wallstress_dns_top()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only:rprec
use param, only:nxt, nynpy, nz, nu_molec, z_i, u_scale, dz, lbc_mom, ubc_mom
use sim_param, only:u, v, dudz, dvdz, txz, tyz
implicit none
! ---
integer :: jx, jy
! ---
select case (ubc_mom)

case ('wall')
    do jy = 1, nynpy
    do jx = 1, nxt
        txz(jx, jy, nz) = nu_molec/(z_i*u_scale)*u(jx, jy, nz-1)/(0.5_rprec*dz)
        tyz(jx, jy, nz) = nu_molec/(z_i*u_scale)*v(jx, jy, nz-1)/(0.5_rprec*dz)
        dudz(jx, jy, nz) = -u(jx, jy, nz-1)/(0.5_rprec*dz)
        dvdz(jx, jy, nz) = -v(jx, jy, nz-1)/(0.5_rprec*dz)
    end do
    end do
case ('stress free')
        txz(:, :, nz)  = 0._rprec
        tyz(:, :, nz)  = 0._rprec
        dudz(:, :, nz) = 0._rprec
        dvdz(:, :, nz) = 0._rprec
case default
    write (*, *) 'invalid ubc_mom'
    stop
end select
! ---
end subroutine wallstress_dns_top
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine wallstress_bottom()
!--------------------------------------------------------------------!
! For use with staggered grid LES JDA, 23 Jan 96
! zo is nondimensionalized, zo1 not!
!--provides txz, tyz, dudz, dvdz at jz=1
!--------------------------------------------------------------------!  
use types, only:rprec
use param
use sim_param, only:u, v, dudz, dvdz, txz, tyz
use bottombc, only:zo, psi_m, phi_m, ustar_avg, d0, B
use test_filtermodule
use stokes_drift, only: deg2rad
implicit none
! ---  
integer :: jx, jy
real(kind=rprec), dimension(nxt,nynpy) :: u_avg, denom
real(kind=rprec), dimension(ldx,nynpy) :: u1, v1
real(kind=rprec) :: const
! ---

if (ocean_flag) then
    if (flag_dynStress) then
        ! Add dynamical wind stress (Bicheng Chen 10/24/2016)
        open(unit=fid_ustar, file=fn_ustar, access='direct',recl=len_ustar)
        read(unit=fid_ustar, rec=nums+ttshift_dynStress+1) ustar_dyn, agl_stress
        close(fid_ustar)
        rad_stress = deg2rad(agl_stress)
        ustar_dyn = ustar_dyn / u_scale
        txz(1:nxt, :, 1) = ustar_dyn**2 * cos(rad_stress)
        tyz(1:nxt, :, 1) = ustar_dyn**2 * sin(rad_stress)
    else
        txz(1:nxt, :, 1) = u_star**2 * cos(rad_stress)   ! !DY Everything has been normalized by u_scale
        tyz(1:nxt, :, 1) = u_star**2 * sin(rad_stress)
    end if

    dudz(1:nxt, :, 1) = (dudz(1:nxt, :, 2)*cos(rad_stress) +  &
                        dvdz(1:nxt, :, 2)*sin(rad_stress)) * cos(rad_stress)
    dvdz(1:nxt, :, 1) = (dudz(1:nxt, :, 2)*cos(rad_stress) +  &
                        dvdz(1:nxt, :, 2)*sin(rad_stress)) * sin(rad_stress)
else
    select case (lbc_mom)
    case ('wall')
        u1 = u(:,:,1)
        v1 = v(:,:,1)

        call test_filter(u1, G_test)
        call test_filter(v1, G_test)
        
        select case (walltype)
        case('smooth')
            denom = dlog(0.5_rprec*dz/(nu_molec/u_scale/z_i)) / vonk + B
        case('rough')
            denom = (dlog((0.5_rprec*dz-d0)/zo) - psi_m) / vonk
        case default
            write(*,*) 'invalid wall type'
            stop
        end select
        u_avg = sqrt(u1(1:nxt,1:nynpy)**2 + v1(1:nxt,1:nynpy)**2)
        ustar_avg = u_avg / denom

        do jy = 1, nynpy
        do jx = 1, nxt
            const = -(ustar_avg(jx,jy)**2) / u_avg(jx,jy)
            txz(jx,jy,1) = const * u1(jx,jy)
            tyz(jx,jy,1) = const * v1(jx,jy)
             !this is as in Moeng 84

            dudz(jx,jy,1) = ustar_avg(jx,jy) / ((0.5_rprec*dz-d0(jx,jy))*vonK)  &
                          * u(jx,jy,1) / u_avg(jx,jy) * phi_m(jx,jy)
            dvdz(jx,jy,1) = ustar_avg(jx,jy) / ((0.5_rprec*dz-d0(jx,jy))*vonK)  &
                          * v(jx,jy,1) / u_avg(jx,jy) * phi_m(jx,jy)
            ! dudz(jx,jy,1) = ustar_avg(jx,jy) / (0.5_rprec*dz*vonK)  &
            !               * u(jx,jy,1) / u_avg(jx,jy)
            ! dvdz(jx,jy,1) = ustar_avg(jx,jy) / (0.5_rprec*dz*vonK)  &
            !               * v(jx,jy,1) / u_avg(jx,jy)

            dudz(jx,jy,1) = merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
            dvdz(jx,jy,1) = merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
        end do
        end do
    case ('stress free')
        txz(:, :, 1)  = 0._rprec
        tyz(:, :, 1)  = 0._rprec
        dudz(:, :, 1) = 0._rprec
        dvdz(:, :, 1) = 0._rprec
    case default
        write (*, *) 'invalid lbc_mom'
        stop
    end select
    
end if
! ---
end subroutine wallstress_bottom
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------! 
subroutine wallstress_top()
!--------------------------------------------------------------------!
! For use with staggered grid LES JDA, 23 Jan 96
! zo is nondimensionalized, zo1 not!
!--provides txz, tyz, dudz, dvdz at jz=1
!--------------------------------------------------------------------!  
use types, only:rprec
use param
use sim_param, only:u, v, dudz, dvdz, txz, tyz
use bottombc, only:zo, psi_m, phi_m, d0, B
use test_filtermodule
use stokes_drift, only: deg2rad
implicit none
! ---  
integer :: jx, jy
real(kind=rprec), dimension(nxt,nynpy) :: u_avg, denom, tmp_avg
real(kind=rprec), dimension(ldx,nynpy) :: u1, v1
real(kind=rprec) :: const
! ---
select case (ubc_mom)

case ('wall')
    u1 = u(:,:,nz-1)
    v1 = v(:,:,nz-1)

    call test_filter(u1, G_test)
    call test_filter(v1, G_test)

    select case (walltype)
    case('smooth')
        denom = dlog(0.5_rprec*dz/(nu_molec/u_scale/z_i)) / vonk + B
    case('rough')
        denom = dlog((0.5_rprec*dz-d0)/zo) / vonk ! - psi_m of no use for upper boundary
    case default
        write(*,*) 'invalid wall type'
        stop
    end select
    u_avg = sqrt(u1(1:nxt,1:nynpy)**2 + v1(1:nxt,1:nynpy)**2)
    tmp_avg = u_avg / denom

    do jy = 1, nynpy
    do jx = 1, nxt
        const = -(tmp_avg(jx,jy)**2) / u_avg(jx,jy)
        txz(jx,jy,nz) = -const * u1(jx,jy)
        tyz(jx,jy,nz) = -const * v1(jx,jy)
        !this is as in Moeng 84

        dudz(jx,jy,nz) = -tmp_avg(jx,jy) / (0.5_rprec*dz*vonK)  &
                        * u(jx,jy,nz-1) / u_avg(jx,jy)
        dvdz(jx,jy,nz) = -tmp_avg(jx,jy) / (0.5_rprec*dz*vonK)  &
                        * v(jx,jy,nz-1) / u_avg(jx,jy)

        dudz(jx,jy,nz) = merge(0._rprec,dudz(jx,jy,nz),u(jx,jy,nz-1).eq.0._rprec)
        dvdz(jx,jy,nz) = merge(0._rprec,dvdz(jx,jy,nz),v(jx,jy,nz-1).eq.0._rprec)
    end do
    end do
case ('stress free')
    txz(:, :, nz)  = 0._rprec
    tyz(:, :, nz)  = 0._rprec
    dudz(:, :, nz) = 0._rprec
    dvdz(:, :, nz) = 0._rprec
case default
    write (*, *) 'invalid ubc_mom'
    stop
end select
! ---
end subroutine wallstress_top
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------! 