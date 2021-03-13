!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
module scalars_module2
!--------------------------------------------------------------------!
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
!--------------------------------------------------------------------!    
use types, only:rprec
use param
use bottombc !makes obukhov functions available
use sim_param, only:u, v, w, theta, q, pcon
use sgsmodule, only:Nu_t
use scalars_module, only:L, wstar, dTdz, dqdz, sgs_t3, sgs_q3,          &
                         DPCondz, sgs_PCon3, res_PCon3, Kc_t, Cs2Sc,    &
                         sgs_PCon1
use intermediate
implicit none
save
! ---
character(*), parameter :: fmt_5168 = "(1400(E14.5))"

contains

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
    subroutine flow_slice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To compute the budget of TKE, we need the following terms:
! <u'w'>*d/dz(<U>), <v'w'>*d/dz(<V>), (g/<T>)*<w'T'>, [0.5*d/dz(<w'e>) where
! e=<u'^2>+<v'^2>+<w'^2>], d/dz(<w'p'>), dissipation (Nu_t*S^2),d/dz(<u'\tau{13}>)+
! d/dz(v'\tau{23}), <\tau{13}>d/dz(<U>)+<\tau{23}>d/dz(<V>)
! Of the eight terms above, we can directly compute terms 1,2,3,8 from the variables
! calculated/outputted in flow_slice.f90 and theta_slice.f90
! So, the rest 4 terms will be computed/outputted here
! Similarly for temperature variance budget, we need
! <w'T'>*d/dz(<T>), d/dz(<w*T^2>) and we already have term 1 from
! theta_slice.f90. so we just need to compute term 2 here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use sim_param
    use param, only:path, p_count, c_count, jt
    use sgsmodule
    use intermediate 
    use stokes_drift, only: ust, vst
    implicit none
    ! ---
    real(rprec), dimension(nxt, nz-1) :: tu1, tv1, tw1, tp1, tdudt, tdvdt, tdwdt, tdudz, tdvdz
    real(rprec), dimension(nxt, nz-1) :: tu2, tv2, tw2, tp2, tuv, tuw, tvw, twp, tCs
    real(rprec), dimension(nxt, nz-1) :: ttxx, ttxy, ttxz, ttyy, ttyz, ttzz, tu3, tv3, tw3
    real(rprec), dimension(nxt, nz-1) :: tu2_fluc, tv2_fluc, tuv_fluc, twp_fluc
    real(rprec), dimension(nxt, nz-1) :: tu2w_fluc, tv2w_fluc, tdissip
    real(rprec), dimension(nxt, nz-1) :: tutxz_fluc, tvtyz_fluc, twtzz
    real(kind=rprec), dimension(0:nz-1) :: ubar_profile, vbar_profile, pbar_profile
    real(kind=rprec), dimension(ldx, nynpy, 0:nz) :: p_org
    real(rprec), dimension(nxt, nynpy, nz-1) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
    integer :: jz
    ! ---
    tu1 = 0._rprec; tv1 = 0._rprec; tw1 = 0._rprec; tp1 = 0._rprec
    tdudt = 0._rprec; tdvdt = 0._rprec; tdwdt = 0._rprec; tdudz = 0._rprec; tdvdz = 0._rprec
    tu2 = 0._rprec; tv2 = 0._rprec; tw2 = 0._rprec; tp2 = 0._rprec;
    tuv = 0._rprec; tuw = 0._rprec; tvw = 0._rprec; twp = 0._rprec; tCs = 0._rprec 
    ttxx = 0._rprec; ttxy = 0._rprec; ttxz = 0._rprec
    ttyy = 0._rprec; ttyz = 0._rprec; ttzz = 0._rprec
    tu3 = 0._rprec; tv3 = 0._rprec; tw3 = 0._rprec
    tu2_fluc = 0._rprec; tv2_fluc = 0._rprec; tuv_fluc = 0._rprec
    tu2w_fluc = 0._rprec; tv2w_fluc = 0._rprec; tdissip = 0._rprec
    tutxz_fluc = 0._rprec; tvtyz_fluc = 0._rprec; twtzz = 0._rprec
    ! ---
    do jz = 0, nz-1
        p_org(:, :,jz) = p(:, :,jz) - 0.5*(u(:, :, jz)*u(:, :, jz) +   &
                                           v(:, :, jz)*v(:, :, jz) +   &
                                           w(:, :, jz)*w(:, :, jz))    &
                                    -     (u(:, :, jz)*ust(jz)+v(:, :, jz)*vst(jz))   &
                                    - 0.5*(ust(jz)**2 + vst(jz)**2)
    end do
    ! ---
    if (coordz == 0) then
        arg1(1:nxt, 1:nynpy, 1) = 0._rprec
        arg2(1:nxt, 1:nynpy, 1) = 0._rprec
        arg1(1:nxt, 1:nynpy, 2:nz-1) = (u(1:nxt, 1:nynpy, 2:nz-1) + u(1:nxt, 1:nynpy, 1:nz-2)) / 2._rprec
        arg2(1:nxt, 1:nynpy, 2:nz-1) = (v(1:nxt, 1:nynpy, 2:nz-1) + v(1:nxt, 1:nynpy, 1:nz-2)) / 2._rprec
    else
        arg1(1:nxt, 1:nynpy, 1:nz-1) = (u(1:nxt, 1:nynpy, 1:nz-1) + u(1:nxt, 1:nynpy, 0:nz-2)) / 2._rprec
        arg2(1:nxt, 1:nynpy, 1:nz-1) = (v(1:nxt, 1:nynpy, 1:nz-1) + v(1:nxt, 1:nynpy, 0:nz-2)) / 2._rprec
    end if
    arg3(1:nxt, 1:nynpy, 1:nz-1) = (p_org(1:nxt, 1:nynpy, 1:nz-1) + p_org(1:nxt, 1:nynpy, 0:nz-2)) / 2._rprec
    ! ---
    call avg_y3d(u(1:nxt, 1:nynpy, 1:nz-1), tu1)
    call avg_y3d(v(1:nxt, 1:nynpy, 1:nz-1), tv1)
    call avg_y3d(w(1:nxt, 1:nynpy, 1:nz-1), tw1)
    call avg_y3d(p_org(1:nxt, 1:nynpy, 1:nz-1), tp1)
    call avg_y3d(dudt(1:nxt, 1:nynpy, 1:nz-1), tdudt)
    call avg_y3d(dvdt(1:nxt, 1:nynpy, 1:nz-1), tdvdt)
    call avg_y3d(dwdt(1:nxt, 1:nynpy, 1:nz-1), tdwdt)
    call avg_y3d(dudz(1:nxt, 1:nynpy, 1:nz-1), tdudz)
    call avg_y3d(dvdz(1:nxt, 1:nynpy, 1:nz-1), tdvdz)

    call avg_y3d(u(1:nxt, 1:nynpy, 1:nz-1)**2., tu2)
    call avg_y3d(v(1:nxt, 1:nynpy, 1:nz-1)**2., tv2)
    call avg_y3d(w(1:nxt, 1:nynpy, 1:nz-1)**2., tw2)
    call avg_y3d(p_org(1:nxt, 1:nynpy, 1:nz-1)**2., tp2)
    call avg_y3d(u(1:nxt, 1:nynpy, 1:nz-1)*v(1:nxt, 1:nynpy, 1:nz-1), tuv)
    call avg_y3d(arg1(1:nxt, 1:nynpy, 1:nz-1)*w(1:nxt, 1:nynpy, 1:nz-1), tuw)
    call avg_y3d(arg2(1:nxt, 1:nynpy, 1:nz-1)*w(1:nxt, 1:nynpy, 1:nz-1), tvw) 
    call avg_y3d(arg3(1:nxt, 1:nynpy, 1:nz-1)*w(1:nxt, 1:nynpy, 1:nz-1), twp)
    call avg_y3d(txx(1:nxt, 1:nynpy, 1:nz-1), ttxx)
    call avg_y3d(txy(1:nxt, 1:nynpy, 1:nz-1), ttxy)
    call avg_y3d(txz(1:nxt, 1:nynpy, 1:nz-1), ttxz)
    call avg_y3d(tyy(1:nxt, 1:nynpy, 1:nz-1), ttyy)
    call avg_y3d(tyz(1:nxt, 1:nynpy, 1:nz-1), ttyz)
    call avg_y3d(tzz(1:nxt, 1:nynpy, 1:nz-1), ttzz)

    call avg_y3d(u(1:nxt, 1:nynpy, 1:nz-1)**3., tu3)
    call avg_y3d(v(1:nxt, 1:nynpy, 1:nz-1)**3., tv3)
    call avg_y3d(w(1:nxt, 1:nynpy, 1:nz-1)**3., tw3)
    call avg_y3d(dissip(1:nxt, 1:nynpy, 1:nz-1), tdissip)
    call avg_y3d(sqrt(Cs_opt2(1:nxt, 1:nynpy, 1:nz-1)), tCs)
    ! ---
    call calc_hor_avg(u(:,:,:), ubar_profile)
    call calc_hor_avg(v(:,:,:), vbar_profile)
    ! call calc_hor_avg(p_org(:,:,:), pbar_profile)
    do jz = 1, nz-1
        arg4(1:nxt,1:nynpy,jz) = (u(1:nxt,1:nynpy,jz) - ubar_profile(jz))**2.
        arg5(1:nxt,1:nynpy,jz) = (v(1:nxt,1:nynpy,jz) - vbar_profile(jz))**2.
        arg6(1:nxt,1:nynpy,jz) = (u(1:nxt,1:nynpy,jz) - ubar_profile(jz))*(v(1:nxt,1:nynpy,jz) - vbar_profile(jz))
    end do
    call avg_y3d(arg4(1:nxt,1:nynpy,1:nz-1), tu2_fluc)
    call avg_y3d(arg5(1:nxt,1:nynpy,1:nz-1), tv2_fluc)
    call avg_y3d(arg6(1:nxt,1:nynpy,1:nz-1), tuv_fluc)

    if (coordz == 0) then
        arg7(1:nxt,1:nynpy,1) = 0._rprec
        arg8(1:nxt,1:nynpy,1) = 0._rprec
        do jz = 2, nz-1
            arg7(1:nxt,1:nynpy,jz) = (u(1:nxt,1:nynpy,jz)-ubar_profile(jz) +    &
                                      u(1:nxt,1:nynpy,jz-1)-ubar_profile(jz-1)) / 2._rprec
            arg8(1:nxt,1:nynpy,jz) = (v(1:nxt,1:nynpy,jz)-vbar_profile(jz) +    &
                                      v(1:nxt,1:nynpy,jz-1)-vbar_profile(jz-1)) / 2._rprec  
        end do
    else
        do jz = 1, nz-1
            arg7(1:nxt,1:nynpy,jz) = (u(1:nxt,1:nynpy,jz)-ubar_profile(jz) +    &
                                      u(1:nxt,1:nynpy,jz-1)-ubar_profile(jz-1)) / 2._rprec 
            arg8(1:nxt,1:nynpy,jz) = (v(1:nxt,1:nynpy,jz)-vbar_profile(jz) +    &
                                      v(1:nxt,1:nynpy,jz-1)-vbar_profile(jz-1)) / 2._rprec
        end do
    end if

    call avg_y3d(arg7(1:nxt,1:nynpy,1:nz-1)**2.*w(1:nxt, 1:nynpy, 1:nz-1), tu2w_fluc)
    call avg_y3d(arg8(1:nxt,1:nynpy,1:nz-1)**2.*w(1:nxt, 1:nynpy, 1:nz-1), tv2w_fluc)
    call avg_y3d(arg7(1:nxt,1:nynpy,1:nz-1)*txz(1:nxt, 1:nynpy, 1:nz-1), tutxz_fluc)
    call avg_y3d(arg8(1:nxt,1:nynpy,1:nz-1)*tyz(1:nxt, 1:nynpy, 1:nz-1), tvtyz_fluc)
    call avg_y3d(w(1:nxt, 1:nynpy, 1:nz-1)*tzz(1:nxt, 1:nynpy, 1:nz-1), twtzz)
    ! ---

    avg_u(:, :) = avg_u(:, :) + fr*tu1(:, :)
    avg_v(:, :) = avg_v(:, :) + fr*tv1(:, :)
    avg_w(:, :) = avg_w(:, :) + fr*tw1(:, :)
    avg_p(:, :) = avg_p(:, :) + fr*tp1(:, :)
    avg_dudt(:, :) = avg_dudt(:, :) + fr*tdudt(:, :)
    avg_dvdt(:, :) = avg_dvdt(:, :) + fr*tdvdt(:, :)
    avg_dwdt(:, :) = avg_dwdt(:, :) + fr*tdwdt(:, :)
    avg_dudz(:, :) = avg_dudz(:, :) + fr*tdudz(:, :)
    avg_dvdz(:, :) = avg_dvdz(:, :) + fr*tdvdz(:, :)

    avg_u2(:, :) = avg_u2(:, :) + fr*tu2(:, :)
    avg_v2(:, :) = avg_v2(:, :) + fr*tv2(:, :)
    avg_w2(:, :) = avg_w2(:, :) + fr*tw2(:, :)
    avg_p2(:, :) = avg_p2(:, :) + fr*tp2(:, :)
    avg_uv(:, :) = avg_uv(:, :) + fr*tuv(:, :)
    avg_uw(:, :) = avg_uw(:, :) + fr*tuw(:, :)
    avg_vw(:, :) = avg_vw(:, :) + fr*tvw(:, :)
    avg_wp(:, :) = avg_wp(:, :) + fr*twp(:, :)

    avg_txx(:, :) = avg_txx(:, :) + fr*ttxx(:, :)
    avg_txy(:, :) = avg_txy(:, :) + fr*ttxy(:, :) 
    avg_txz(:, :) = avg_txz(:, :) + fr*ttxz(:, :) 
    avg_tyy(:, :) = avg_tyy(:, :) + fr*ttyy(:, :) 
    avg_tyz(:, :) = avg_tyz(:, :) + fr*ttyz(:, :) 
    avg_tzz(:, :) = avg_tzz(:, :) + fr*ttzz(:, :)

    avg_Cs(:, :) = avg_Cs(:, :) + fr*tCs(:, :)
    
    avg_u3(:, :) = avg_u3(:, :) + fr*tu3(:, :)
    avg_v3(:, :) = avg_v3(:, :) + fr*tv3(:, :)
    avg_w3(:, :) = avg_w3(:, :) + fr*tw3(:, :)

    avg_u2_fluc(:, :) = avg_u2_fluc(:, :) + fr*tu2_fluc(:, :)
    avg_v2_fluc(:, :) = avg_v2_fluc(:, :) + fr*tv2_fluc(:, :)
    avg_uv_fluc(:, :) = avg_uv_fluc(:, :) + fr*tuv_fluc(:, :)

    avg_u2w_fluc(:, :) = avg_u2w_fluc(:, :) + fr*tu2w_fluc(:, :)
    avg_v2w_fluc(:, :) = avg_v2w_fluc(:, :) + fr*tv2w_fluc(:, :)
    avg_dissip(:, :) = avg_dissip(:, :) + fr*tdissip(:, :)

    avg_utxz_fluc(:, :) = avg_utxz_fluc(:, :) + fr*tutxz_fluc(:, :)
    avg_vtyz_fluc(:, :) = avg_vtyz_fluc(:, :) + fr*tvtyz_fluc(:, :)
    avg_wtzz(:, :) = avg_wtzz(:, :) + fr*twtzz(:, :)

    ! ---
    if (mod(jt,p_count)==0) then
        call write_avgslice(avg_u(:,:), 'u')
        call write_avgslice(avg_v(:,:), 'v')
        call write_avgslice(avg_w(:,:), 'w')
        call write_avgslice(avg_p(:,:), 'p')
        call write_avgslice(avg_dudt(:,:), 'dudt')
        call write_avgslice(avg_dvdt(:,:), 'dvdt')
        call write_avgslice(avg_dwdt(:,:), 'dwdt')
        call write_avgslice(avg_dudz(:,:), 'dudz')
        call write_avgslice(avg_dvdz(:,:), 'dvdz')

        call write_avgslice(avg_u2(:,:), 'u2')
        call write_avgslice(avg_v2(:,:), 'v2')
        call write_avgslice(avg_w2(:,:), 'w2')
        call write_avgslice(avg_p2(:,:), 'p2')
        call write_avgslice(avg_uv(:,:), 'uv')
        call write_avgslice(avg_uw(:,:), 'uw')
        call write_avgslice(avg_vw(:,:), 'vw')
        call write_avgslice(avg_wp(:,:), 'wp')

        call write_avgslice(avg_txx(:,:), 'txx')
        call write_avgslice(avg_txy(:,:), 'txy')
        call write_avgslice(avg_txz(:,:), 'txz')
        call write_avgslice(avg_tyy(:,:), 'tyy')
        call write_avgslice(avg_tyz(:,:), 'tyz')
        call write_avgslice(avg_tzz(:,:), 'tzz')

        call write_avgslice(avg_u3(:,:), 'u3')
        call write_avgslice(avg_v3(:,:), 'v3')
        call write_avgslice(avg_w3(:,:), 'w3')

        call write_avgslice(avg_u2_fluc(:,:), 'u2fluc')
        call write_avgslice(avg_v2_fluc(:,:), 'v2fluc')
        call write_avgslice(avg_uv_fluc(:,:), 'uvfluc')
        call write_avgslice(avg_u2w_fluc(:,:), 'u2wfluc')
        call write_avgslice(avg_v2w_fluc(:,:), 'v2wfluc')

        call write_avgslice(avg_utxz_fluc(:,:), 'utxz_fluc')
        call write_avgslice(avg_vtyz_fluc(:,:), 'vtyz_fluc')
        call write_avgslice(avg_wtzz(:,:), 'wtzz')
        
        call write_avgslice(avg_Cs(:,:), 'Cs')
        call write_avgslice(avg_dissip(:,:), 'dissip')

        avg_u = 0._rprec; avg_v = 0._rprec; avg_w = 0._rprec; avg_p = 0._rprec
        avg_dudt = 0._rprec; avg_dvdt = 0._rprec; avg_dwdt = 0._rprec
        avg_dudz = 0._rprec; avg_dvdz = 0._rprec
        avg_u2 = 0._rprec; avg_v2 = 0._rprec; avg_w2 = 0._rprec; avg_p2 = 0._rprec
        avg_uw = 0._rprec; avg_vw = 0._rprec; avg_uv = 0._rprec; avg_wp = 0._rprec
        avg_txx = 0._rprec; avg_txy = 0._rprec; avg_txz = 0._rprec
        avg_tyy = 0._rprec; avg_tyz = 0._rprec; avg_tzz = 0._rprec
        avg_u3 = 0._rprec; avg_v3 = 0._rprec; avg_w3 = 0._rprec
        avg_u2_fluc = 0._rprec; avg_v2_fluc = 0._rprec; avg_uv_fluc = 0._rprec
        avg_u2w_fluc = 0._rprec; avg_v2w_fluc = 0._rprec; avg_dissip = 0._rprec
        avg_utxz_fluc = 0._rprec; avg_vtyz_fluc = 0._rprec; avg_wtzz = 0._rprec; avg_Cs = 0._rprec
    end if
    ! ---
    end subroutine flow_slice
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!  
    subroutine theta_slice()
    !subroutine theta_slice(w,t,q,sgs_t3,sgs_q3,dTdz,beta_scal,Nu_t,jt)
    !c This is exactly the same like the subroutine avgslice with the
    !c only difference being that it averages the scalar variables
    !c to find the y-averaged instantaneous x-z slices of variables
    !c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
    use intermediate
    implicit none
    ! ---
    real(rprec), dimension(nxt, nz-1) :: ttheta1, ttheta2, ttheta2_fluc, tdTdz, tsgst, twt, tnu_t
    real(rprec), dimension(nxt, nynpy, nz-1) :: arg1, arg2
    real(rprec), dimension(0:nz-1) :: Tbar_profile
    integer :: jx, jz
    ! ---
    ttheta1 = 0._rprec; ttheta2 = 0._rprec; ttheta2_fluc = 0._rprec; tdTdz = 0._rprec
    tsgst = 0._rprec; twt = 0._rprec; tnu_t = 0._rprec

    arg1 = (theta(1:nxt, 1:nynpy, 1:nz-1) + theta(1:nxt, 1:nynpy, 0:nz-2)) / 2._rprec
    call avg_y3d(theta(1:nxt, 1:nynpy, 1:nz-1), ttheta1)
    call avg_y3d(theta(1:nxt, 1:nynpy, 1:nz-1)*theta(1:nxt, 1:nynpy, 1:nz-1), ttheta2)
    call avg_y3d(dTdz(1:nxt, 1:nynpy, 1:nz-1), tdTdz)
    call avg_y3d(sgs_t3(1:nxt, 1:nynpy, 1:nz-1), tsgst)
    call avg_y3d(w(1:nxt, 1:nynpy, 1:nz-1)*arg1(1:nxt, 1:nynpy, 1:nz-1), twt)
    call avg_y3d(Nu_t(1:nxt, 1:nynpy, 1:nz-1), tnu_t)

    ! ---
    call calc_hor_avg(theta(:,:,:), Tbar_profile)
    do jz = 1, nz-1
        arg2(1:nxt,1:nynpy,jz) = (theta(1:nxt,1:nynpy,jz) - Tbar_profile(jz))**2.
    end do
    call avg_y3d(arg2(1:nxt, 1:nynpy, 1:nz-1), ttheta2_fluc)
    ! ---

    avg_theta(:, :)    = avg_theta(:, :)    + fr*ttheta1(:, :)
    avg_theta2(:, :)   = avg_theta2(:, :)   + fr*ttheta2(:, :)
    avg_dTdz(:, :)     = avg_dTdz(:, :)     + fr*tdTdz(:, :)
    avg_sgsTheta(:, :) = avg_sgsTheta(:, :) + fr*tsgst(:, :)
    avg_wTheta(:, :)   = avg_wTheta(:, :)   + fr*twt(:, :)
    avg_Nut(:, :)      = avg_Nut(:, :)      + fr*tnu_t(:, :)
    avg_T2_fluc(:, :)  = avg_T2_fluc(:, :)  + fr*ttheta2_fluc(:, :)

    if (mod(jt,p_count)==0) then
        call write_avgslice(avg_theta(:,:), 'theta')
        call write_avgslice(avg_theta2(:,:), 'theta2')
        call write_avgslice(avg_T2_fluc(:,:), 'Tfluc')
        call write_avgslice(avg_dTdz(:,:), 'dTdz')
        call write_avgslice(avg_sgsTheta(:,:), 'sgs_t3')
        call write_avgslice(avg_wTheta(:,:), 'wt')
        call write_avgslice(avg_Nut(:,:), 'Nu_t')

        avg_theta = 0._rprec; avg_theta2 = 0._rprec; avg_dTdz = 0._rprec
        avg_sgsTheta = 0._rprec; avg_wTheta = 0._rprec; avg_Nut = 0._rprec
        avg_T2_fluc = 0._rprec
    end if
    ! ---
    end subroutine theta_slice
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
    subroutine pcon_slice()
    !c This is exactly the same like the subroutine avgslice with the
    !c only difference being that it averages the scalar variables
    !c to find the y-averaged instantaneous x-z slices of variables
    !c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
    !c It also outputs the average covariance between wt and wq
    !use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3
    !use output_slice,only: collocate_MPI_averages
    use intermediate
    implicit none
    ! ---
    integer:: ipcon
    real(rprec), dimension(nxt, nz-1) :: tPCon, tPCon2, tdPCondz, tsgs_PCon1, tuPCon, tsgs_PCon3, twPCon 
    real(rprec), dimension(nxt, nz-1) :: tKc_t, tCs2Sc
    real(rprec), dimension(nxt, nynpy, nz-1) :: arg1
    ! ---
    do ipcon = 1, npcon
        tPCon = 0._rprec; tPCon2 = 0._rprec; tdPCondz = 0._rprec
        tsgs_PCon1 = 0._rprec; tsgs_PCon3 = 0._rprec; tuPCon = 0._rprec; twPCon = 0._rprec 

        call avg_y3d(PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon), tPCon)
        call avg_y3d(PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon)*PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon), tPCon2)
        call avg_y3d(dPCondz(1:nxt, 1:nynpy, 1:nz-1, ipcon), tdPCondz)
        call avg_y3d(sgs_PCon1(1:nxt, 1:nynpy, 1:nz-1, ipcon), tsgs_PCon1)
        call avg_y3d(sgs_PCon3(1:nxt, 1:nynpy, 1:nz-1, ipcon), tsgs_PCon3)
        call avg_y3d(u(1:nxt, 1:nynpy, 1:nz-1)*PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon), tuPCon)

        if (PCon_scheme == 1) then
            arg1 = (PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon) + PCon(1:nxt, 1:nynpy, 0:nz-2, ipcon)) / 2._rprec
            call avg_y3d(w(1:nxt, 1:nynpy, 1:nz-1)*arg1(1:nxt, 1:nynpy, 1:nz-1), twPCon)
        else 
            call avg_y3d(res_PCon3(1:nxt, 1:nynpy, 1:nz-1, ipcon), twPCon)
        end if

        avg_PCon(:,:,ipcon)     = avg_PCon(:,:,ipcon)     + fr*tPCon(:,:)
        avg_PCon2(:,:,ipcon)    = avg_PCon2(:,:,ipcon)    + fr*tPCon2(:,:)
        avg_dPCondz(:,:,ipcon)  = avg_dPCondz(:,:,ipcon)  + fr*tdPCondz(:,:)
        avg_sgsPCon1(:,:,ipcon) = avg_sgsPCon1(:,:,ipcon) + fr*tsgs_PCon1(:,:)
        avg_sgsPCon3(:,:,ipcon) = avg_sgsPCon3(:,:,ipcon) + fr*tsgs_PCon3(:,:)
        avg_uPCon(:,:,ipcon)    = avg_uPCon(:,:,ipcon)    + fr*tuPCon(:,:)
        avg_wPCon(:,:,ipcon)    = avg_wPCon(:,:,ipcon)    + fr*twPCon(:,:)
    end do

    tKc_t = 0._rprec; tCs2Sc = 0._rprec
    call avg_y3d(Kc_t(1:nxt, 1:nynpy, 1:nz-1), tKc_t)
    call avg_y3d(Cs2Sc(1:nxt, 1:nynpy, 1:nz-1), tCs2Sc)
    avg_Kct(:,:)   = avg_Kct(:,:)   + fr*tKc_t(:,:)
    avg_Cs2Sc(:,:) = avg_Cs2Sc(:,:) + fr*tCs2Sc(:,:)
    ! ---
    if (mod(jt, p_count) == 0) then
        do ipcon = 1, npcon
            call write_avgslice(avg_PCon(:,:,ipcon), 'PCon')
            call write_avgslice(avg_PCon2(:,:,ipcon), 'PCon2')
            call write_avgslice(avg_dPCondz(:,:,ipcon), 'dPCondz')
            call write_avgslice(avg_sgsPCon1(:,:,ipcon), 'sgsPCon1')
            call write_avgslice(avg_sgsPCon3(:,:,ipcon), 'sgsPCon3')
            call write_avgslice(avg_uPCon(:,:,ipcon), 'uPCon')
            call write_avgslice(avg_wPCon(:,:,ipcon), 'wPCon')
        end do
        call write_avgslice(avg_Kct, 'Kct')
        call write_avgslice(avg_Cs2Sc, 'Cs2Sc')
        ! --- Reset quantities
        avg_PCon = 0._rprec; avg_PCon2 = 0._rprec; avg_dPCondz = 0._rprec
        avg_sgsPCon1 = 0._rprec; avg_sgsPCon3 = 0._rprec
        avg_uPCon = 0._rprec; avg_wPCon = 0._rprec; avg_Kct = 0._rprec; avg_Cs2Sc = 0._rprec
    end if
    ! ---
    end subroutine pcon_slice
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
    subroutine write_avgslice(arr2d, fn_str)
    !-----------------------------------------------------------------------
    ! For average_dim_num = 1, each line is the average in x direction.
    ! For average_dim_num = 2, each line is the vertical profile of average
    ! across the entire domain.
    ! No matter what value the average_dim_num value it is, the first data of
    ! each line is always the non-dimensional time.
    !
    ! Parameter
    ! ---------
    ! arr2d: 2D real array
    !   The quantity need to output.
    ! fn_str: string
    !   The string used in filename to identify the quantity arr2d.
    !
    ! Return
    !-----------------------------------------------------------------------
    implicit none
    real(kind=rprec), dimension(nxt, nz-1), intent(in) :: arr2d
    character(*), intent(in) :: fn_str
    ! ---
    integer, parameter :: fid = 9
    integer :: ind1, ind2
    character(len=256) :: fname
    real(kind=rprec), dimension(nxt, nzt-1) :: arr2d_glb
    ! ---
    fname = path//'result/aver_'//trim(fn_str)//'.out'
    arr2d_glb = 0.0
    ! ---
    if (average_dim_num == 1) then
        do ind1 = 1, nxt
            call collect_data_z(arr2d(ind1, 1:nz-1), arr2d_glb(ind1, 1:nzt-1))
        end do
    else if ( average_dim_num .eq. 2 ) then
        call collect_data_z(sum(arr2d(:,:), dim=1)/nxt, arr2d_glb(1, 1:nzt-1))
    end if
    
    if ( rank == 0 ) then
        open(fid, file=trim(fname), status="unknown", position="append")
        if (average_dim_num == 1) then
            do ind2 = 1, nzt-1
                write(fid, fmt_5168) nums*dt, (arr2d_glb(ind1, ind2), ind1=1, nxt)
            end do
        else if ( average_dim_num .eq. 2 ) then
            write(fid, fmt_5168) nums*dt, (arr2d_glb(1, ind1), ind1=1, nzt-1)
        end if
        close(fid)
    end if
    ! ---
    end subroutine write_avgslice
    !-----------------------------------------------------------------------
    !   
    !-----------------------------------------------------------------------
! ---
end module scalars_module2



