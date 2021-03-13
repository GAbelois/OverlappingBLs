module topbc
!--------------------------------------------------------------------!
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
!--------------------------------------------------------------------!     
use types, only: rprec
use param, only: nz, nzt, damping_method
implicit none
save
! ---
real(kind=rprec), dimension(:), allocatable :: sponge

contains

    subroutine setsponge()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!     
    use param
    implicit none
    ! ---
    real(kind=rprec) :: factor
    real(kind=rprec) :: z_d, cfrdmp, z_local
    integer :: k
    
    integer :: k_global
    real(rprec) :: sponge_top

    real(rprec) :: bv, g_hat, sigma_max, fac, sigma_z, ell
    character(80) :: fname 
    ! ---
    if (damping_method == 1) then
    ! sets relaxation term to vertical momentum equation in top quarter
    ! of domain, relaxation timescale on top is 50s with a factor of 5 if we
    ! had Nz=40 for the layers 40...31, Nieuwstadt et al. 1991, turbulent shear
    ! flows
        sponge = 0._rprec
        factor = 9._rprec / (nzt - 3*nzt/4 + 1)
        sponge_top = z_i / (50._rprec * u_scale)
        do k = 1, nz-1
            k_global = k + coordz * (nz - 1)
            if (k_global > 3*nzt/4 + 1) then
                sponge(k) = sponge_top * 5._rprec**((k_global-nzt) * factor)
            end if
        end do
    else if ( damping_method == 2 ) then
        z_d = 0.75_rprec * lz_tot
        cfrdmp = 3.9_rprec
        sponge = 0._rprec
        do k = 1, nz - 1
            if (zuvp(k) .ge. z_d .and. zuvp(k) .le. lz_tot) then
                sponge(k) = cfrdmp*0.5_rprec*(1._rprec - cos(pi*(zuvp(k) - z_d)/(lz_tot - z_d)))
            else
                sponge(k) = 0._rprec
            end if
        end do

        ! if ( coordz == npz - 1 ) then
        !     sponge(nz) = sponge(nz - 1)
        ! end if
    else if ( damping_method == 3 ) then
        sponge = 0._rprec
        g_hat = g*(z_i/(u_scale**2)) 
        if ( ocean_flag ) then
            fac = 0.22
            bv  = sqrt(g_hat*alpha_w*inv_strength/T_scale*z_i)
        else
            fac = 1.5
            bv = sqrt(g_hat*inv_strength/T_scale*z_i)
        end if

        sigma_z   = 0.8_rprec * lz_tot
        sigma_max = fac * bv / (8.0*atan(1.0))

        do k = 1, nz
            if (zuvp(k) .ge. sigma_z .and. zuvp(k) .le. lz_tot) then
                ell = (zuvp(k) - sigma_z) / (lz_tot - sigma_z)
                sponge(k) = ell*ell*sigma_max
            else
                sponge(k) = 0._rprec
            end if
        end do
    end if
    ! ---
    end subroutine setsponge
! ---    
end module topbc
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine add_damping_layer()
!--------------------------------------------------------------------!
!   Add damping terms to the momentum RHS                                  
!--------------------------------------------------------------------!     
use types, only: rprec
use param, only: nxt, nynpy, nz, ubc_mom, ubc, damping_method
use sim_param, only: u, v, w, RHSx, RHSy, RHSz
use topbc, only: sponge 
implicit none
! ---
integer :: jz
real(kind=rprec) :: u_mn, v_mn, w_mn
! --- 
if ( ubc_mom == 'stress free' ) then
    if (ubc == 1 .and. damping_method /= 1) then
        do jz = 1, nz - 1
            call calc_hor_avg3(u(1:nxt, 1:nynpy, jz), u_mn)
            call calc_hor_avg3(v(1:nxt, 1:nynpy, jz), v_mn)
            call calc_hor_avg3(w(1:nxt, 1:nynpy, jz), w_mn)

            RHSx(1:nxt, 1:nynpy, jz) = RHSx(1:nxt, 1:nynpy, jz) - 0.5*(sponge(jz) + sponge(jz + 1))* &
                                        (u(1:nxt, 1:nynpy, jz) - u_mn)
            RHSy(1:nxt, 1:nynpy, jz) = RHSy(1:nxt, 1:nynpy, jz) - 0.5*(sponge(jz) + sponge(jz + 1))* &
                                        (v(1:nxt, 1:nynpy, jz) - v_mn)
            RHSz(1:nxt, 1:nynpy, jz) = RHSz(1:nxt, 1:nynpy, jz) - sponge(jz) * (w(1:nxt, 1:nynpy, jz) - w_mn)
        end do
    else if (ubc == 1 .and. damping_method == 1) then
        do jz = 1, nz - 1
            RHSz(1:nxt, 1:nynpy, jz) = RHSz(1:nxt, 1:nynpy, jz) - sponge(jz) * w(1:nxt, 1:nynpy, jz)
        end do
    end if
end if
! ---
end subroutine add_damping_layer
