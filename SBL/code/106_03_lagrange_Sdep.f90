!--------------------------------------------------------------------!
! this is the w-node version
!--provides Cs_opt2 1:nz
!--MPI: required u,v on 0:nz, except bottom node 1:nz
!--------------------------------------------------------------------!
subroutine lagrange_Sdep(S11, S12, S13, S22, S23, S33)
!--------------------------------------------------------------------!
! standard dynamic model to calculate the Smagorinsky coefficient
! this is done layer-by-layer to save memory
! everything is done to be on uv-nodes
! -note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary. mem requirement ~ n^2, not n^3
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
!--------------------------------------------------------------------!
use types, only:rprec
use param
use sim_param, only:u, v, w
use sgsmodule, only:F_LM, F_MM, F_QN, F_NN, beta, Cs_opt2, opftime, Beta_avg, Betaclip_avg, Cs_Ssim
use test_filtermodule
use sgsmodule, only:Betaclip, Cs_opt2_2d, Cs_opt2_4d, &
     S, L11, L12, L13, L22, L23, L33, Q11, Q12, Q13, Q22, Q23, Q33, &
     M11, M12, M13, M22, M23, M33, N11, N12, N13, N22, N23, N33, LM, MM, QN, NN, &
     Tn, epsi, dumfac, LMvert, MMvert, QNvert, NNvert, S_bar, S11_bar, S12_bar, &
     S13_bar, S22_bar, S23_bar, S33_bar, S_S11_bar, S_S12_bar, S_S13_bar, &
     S_S22_bar, S_S23_bar, S_S33_bar, S_hat, S11_hat, S12_hat, S13_hat, S22_hat, &
     S23_hat, S33_hat, S_S11_hat, S_S12_hat, S_S13_hat, S_S22_hat, S_S23_hat, &
     S_S33_hat, u_bar, v_bar, w_bar, u_hat, v_hat, w_hat
implicit none
! ---  
real(kind=rprec), dimension(ldx, nynpy, nz), intent(in) :: S11, S12, S13, S22, S23, S33
! ---
integer :: jx, jy, jz
real(kind=rprec) :: tf1, tf2, tf1_2, tf2_2 ! Size if the second test filter
! integer :: counter1, counter2
! real(kind=rprec) :: fractus

real(kind=rprec) :: delta, const
real(kind=rprec) :: lagran_dt, opftdelta, powcoeff
real(kind=rprec), dimension(ldx, nynpy, 0:nz) :: u_temp, v_temp
real(kind=rprec) :: Cs_opt2_havg, Cs_opt2_2d_havg, Cs_opt2_4d_havg  

logical, save :: F_LM_MM_init = .false.
logical, save :: F_QN_NN_init = .false.
! ---
delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)

tf1 = 2._rprec
tf2 = 4._rprec
tf1_2 = tf1**2
tf2_2 = tf2**2

!TS Add the opftdelta
!TS FOR BUILDING (opftime=1*1.5)
opftdelta = opftime*delta
powcoeff = -1._rprec/8._rprec

!TS
u_temp = u
v_temp = v

call interpolag_Sdep()

!! (this is the time step used in the lagrangian computations)
lagran_dt = dt * real(cs_count, kind=rprec)
! fractus   = 1._rprec/real(nyt*nxt, kind=rprec)
const     = 2._rprec*(delta**2)
!        if (jt==cs_count) Beta=1.0
do jz = 1, nz
    ! using L_ij as temp storage here
    if ( coordz == 0 .and. jz == 1 ) then
        !!! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
        ! put on uvp-nodes
        u_bar(:, :) = u_temp(:, :, 1) ! first uv-node
        v_bar(:, :) = v_temp(:, :, 1) ! first uv-node
        w_bar(:, :) = 0.25_rprec*w(:, :, 2) ! first uv-node
    else if ( coordz == npz - 1 .and. jz == nz ) then
        u_bar(:, :) = u_temp(:, :, nz - 1) 
        v_bar(:, :) = v_temp(:, :, nz - 1) 
        w_bar(:, :) = 0.25_rprec*w(:, :, nz - 1)
    else ! w-nodes
        u_bar(:, :) = 0.5_rprec*(u_temp(:, :, jz) + u_temp(:, :, jz - 1))
        v_bar(:, :) = 0.5_rprec*(v_temp(:, :, jz) + v_temp(:, :, jz - 1))
        w_bar(:, :) = w(:, :, jz)
    end if
    u_hat = u_bar
    v_hat = v_bar
    w_hat = w_bar
    L11 = u_bar*u_bar
    L12 = u_bar*v_bar
    L13 = u_bar*w_bar
    L23 = v_bar*w_bar
    L22 = v_bar*v_bar
    L33 = w_bar*w_bar
    Q11 = u_bar*u_bar
    Q12 = u_bar*v_bar
    Q13 = u_bar*w_bar
    Q22 = v_bar*v_bar
    Q23 = v_bar*w_bar
    Q33 = w_bar*w_bar

    call test_filter(u_bar, G_test)
    call test_filter(v_bar, G_test)
    call test_filter(w_bar, G_test)

    call test_filter(L11, G_test) ! in-place filtering
    L11 = L11-u_bar*u_bar
    call test_filter(L12, G_test)
    L12 = L12-u_bar*v_bar
    call test_filter(L13, G_test)
    L13 = L13-u_bar*w_bar
    call test_filter(L22, G_test)
    L22 = L22-v_bar*v_bar
    call test_filter(L23, G_test)
    L23 = L23-v_bar*w_bar
    call test_filter(L33, G_test)
    L33 = L33-w_bar*w_bar

    call test_filter(u_hat, G_test_test)
    call test_filter(v_hat, G_test_test)
    call test_filter(w_hat, G_test_test)

    call test_filter(Q11, G_test_test)
    Q11 = Q11-u_hat*u_hat
    call test_filter(Q12, G_test_test)
    Q12 = Q12-u_hat*v_hat
    call test_filter(Q13, G_test_test)
    Q13 = Q13-u_hat*w_hat
    call test_filter(Q22, G_test_test)
    Q22 = Q22-v_hat*v_hat
    call test_filter(Q23, G_test_test)
    Q23 = Q23-v_hat*w_hat
    call test_filter(Q33, G_test_test)
    Q33 = Q33-w_hat*w_hat

    ! calculate |S|
    S(:, :) = sqrt(2._rprec*(S11(:, :, jz)**2 + S22(:, :, jz)**2 + S33(:, :, jz)**2 + &
                   2._rprec*(S12(:, :, jz)**2 + S13(:, :, jz)**2 + S23(:, :, jz)**2)))

    ! S_ij already on w-nodes
    S11_bar(:, :) = S11(:, :, jz)
    S12_bar(:, :) = S12(:, :, jz)
    S13_bar(:, :) = S13(:, :, jz)
    S22_bar(:, :) = S22(:, :, jz)
    S23_bar(:, :) = S23(:, :, jz)
    S33_bar(:, :) = S33(:, :, jz)

    S11_hat = S11_bar
    S12_hat = S12_bar
    S13_hat = S13_bar
    S22_hat = S22_bar
    S23_hat = S23_bar
    S33_hat = S33_bar

    call test_filter(S11_bar, G_test)
    call test_filter(S12_bar, G_test)
    call test_filter(S13_bar, G_test)
    call test_filter(S22_bar, G_test)
    call test_filter(S23_bar, G_test)
    call test_filter(S33_bar, G_test)

    call test_filter(S11_hat, G_test_test)
    call test_filter(S12_hat, G_test_test)
    call test_filter(S13_hat, G_test_test)
    call test_filter(S22_hat, G_test_test)
    call test_filter(S23_hat, G_test_test)
    call test_filter(S33_hat, G_test_test)

    S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
                 2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

    S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
                 2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

    S_S11_bar(:, :) = S(:, :)*S11(:, :, jz)
    S_S12_bar(:, :) = S(:, :)*S12(:, :, jz)
    S_S13_bar(:, :) = S(:, :)*S13(:, :, jz)
    S_S22_bar(:, :) = S(:, :)*S22(:, :, jz)
    S_S23_bar(:, :) = S(:, :)*S23(:, :, jz)
    S_S33_bar(:, :) = S(:, :)*S33(:, :, jz)

    S_S11_hat(:, :) = S_S11_bar(:, :)
    S_S12_hat(:, :) = S_S12_bar(:, :)
    S_S13_hat(:, :) = S_S13_bar(:, :)
    S_S22_hat(:, :) = S_S22_bar(:, :)
    S_S23_hat(:, :) = S_S23_bar(:, :)
    S_S33_hat(:, :) = S_S33_bar(:, :)

    call test_filter(S_S11_bar, G_test)
    call test_filter(S_S12_bar, G_test)
    call test_filter(S_S13_bar, G_test)
    call test_filter(S_S22_bar, G_test)
    call test_filter(S_S23_bar, G_test)
    call test_filter(S_S33_bar, G_test)

    call test_filter(S_S11_hat, G_test_test)
    call test_filter(S_S12_hat, G_test_test)
    call test_filter(S_S13_hat, G_test_test)
    call test_filter(S_S22_hat, G_test_test)
    call test_filter(S_S23_hat, G_test_test)
    call test_filter(S_S33_hat, G_test_test)

    M11 = const*(S_S11_bar - tf1_2*S_bar*S11_bar)
    M12 = const*(S_S12_bar - tf1_2*S_bar*S12_bar)
    M13 = const*(S_S13_bar - tf1_2*S_bar*S13_bar)
    M22 = const*(S_S22_bar - tf1_2*S_bar*S22_bar)
    M23 = const*(S_S23_bar - tf1_2*S_bar*S23_bar)
    M33 = const*(S_S33_bar - tf1_2*S_bar*S33_bar)

    N11 = const*(S_S11_hat - tf2_2*S_hat*S11_hat)
    N12 = const*(S_S12_hat - tf2_2*S_hat*S12_hat)
    N13 = const*(S_S13_hat - tf2_2*S_hat*S13_hat)
    N22 = const*(S_S22_hat - tf2_2*S_hat*S22_hat)
    N23 = const*(S_S23_hat - tf2_2*S_hat*S23_hat)
    N33 = const*(S_S33_hat - tf2_2*S_hat*S33_hat)

    LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
    MM = M11**2 + M22**2 + M33**2 + 2._rprec*(M12**2 + M13**2 + M23**2)
    QN = Q11*N11+Q22*N22+Q33*N33+2._rprec*(Q12*N12+Q13*N13+Q23*N23)
    NN = N11**2 + N22**2 + N33**2 + 2._rprec*(N12**2 + N13**2 + N23**2)

    if (inilag) then
        if ((.not. F_LM_MM_init) .and. (jt == cs_count .or. jt == DYN_init)) then
            !print *, 'F_MM and F_LM initialized'
            F_MM(:, :, jz) = MM
            F_LM(:, :, jz) = 0.03_rprec*MM
            F_MM(ldx - 1:ldx, :, jz) = 1._rprec
            F_LM(ldx - 1:ldx, :, jz) = 1._rprec

            if (jz == nz) F_LM_MM_init = .true.
        end if
    end if

    ! ! if (inflow) then      ! cyan
    ! !     Tn = merge(.1_rprec*const*S**2, MM, MM .le. .1_rprec*const*S**2)
    ! !     MM = Tn
    ! !     LM(jx_s - x_relax + 1:jx_s, 1:nynpy) = 0._rprec
    ! !     F_LM(jx_s - x_relax + 1:jx_s, 1:nynpy, jz) = 0._rprec
    ! !     Tn = merge(.1_rprec*const*S**2, NN, NN .le. .1_rprec*const*S**2)
    ! !     NN = Tn
    ! !     QN(jx_s - x_relax + 1:jx_s, 1:nynpy) = 0._rprec
    ! !     F_QN(jx_s - x_relax + 1:jx_s, 1:nynpy, jz) = 0._rprec
    ! ! end if

    ! if (inflow) then
    !     Tn = merge(.1_rprec*const*S**2, MM, MM .le. .1_rprec*const*S**2)
    !     MM = Tn
    !     LM(nxt-jx_relax+1:nxt, 1:nynpy) = 0._rprec
    !     F_LM(nxt-jx_relax+1:nxt, 1:nynpy, jz) = 0._rprec
    !     Tn = merge(.1_rprec*const*S**2, NN, NN .le. .1_rprec*const*S**2)
    !     NN = Tn
    !     QN(nxt-jx_relax+1:nxt, 1:nynpy) = 0._rprec
    !     F_QN(nxt-jx_relax+1:nxt, 1:nynpy, jz) = 0._rprec
    ! end if


    !TS FIX THE BUG IN COMPAQ COMPILER
    !TS ADD the following
    !--you mean workaround?
    Tn = max(F_LM(:, :, jz)*F_MM(:, :, jz), real(1.e-24, kind=rprec))
    Tn = opftdelta*(Tn**powcoeff)
    Tn(:, :) = max(real(1.e-24, kind=rprec), Tn(:, :)) 
    dumfac = lagran_dt/Tn
    epsi = dumfac/(1._rprec + dumfac)

    F_LM(:, :, jz) = (epsi*LM + (1._rprec - epsi)*F_LM(:, :, jz))
    F_MM(:, :, jz) = (epsi*MM + (1._rprec - epsi)*F_MM(:, :, jz))
    F_LM(:, :, jz) = max(real(1.e-24, kind=rprec), F_LM(:, :, jz)) 

    Cs_opt2_2d(:, :, jz) = F_LM(:, :, jz)/F_MM(:, :, jz)
!   Cs_opt2_2d(:,:,jz) = SUM(LM(:,:))/SUM(MM(:,:))
!--why set ldx-1:ldx to this?
    Cs_opt2_2d(ldx, :, jz) = real(1.e-24, kind=rprec)
    Cs_opt2_2d(ldx - 1, :, jz) = real(1.e-24, kind=rprec)
    Cs_opt2_2d(:, :, jz) = max(real(1.e-24, kind=rprec), Cs_opt2_2d(:, :, jz)) 

    if (inilag) then
        if ((.not. F_QN_NN_init) .and. (jt == cs_count .or. jt == DYN_init)) then
            !print *, 'F_NN and F_QN initialized'
            F_NN(:, :, jz) = NN
            F_QN(:, :, jz) = 0.03_rprec*NN
            F_NN(ldx - 1:ldx, :, jz) = 1._rprec
            F_QN(ldx - 1:ldx, :, jz) = 1._rprec

            if (jz == nz) F_QN_NN_init = .true.
        end if
    end if

    !TS FIX THE BUG IN COMPAQ COMPILER
    !TS ADD the following
    Tn = max(F_QN(:, :, jz)*F_NN(:, :, jz), real(1.e-24, kind=rprec)) 
    Tn = opftdelta*(Tn**powcoeff)
    Tn(:, :) = max(real(1.e-24, kind=rprec), Tn(:, :)) 
    dumfac = lagran_dt/Tn
    epsi = dumfac/(1._rprec + dumfac)

    F_QN(:, :, jz) = (epsi*QN + (1._rprec - epsi)*F_QN(:, :, jz))
    F_NN(:, :, jz) = (epsi*NN + (1._rprec - epsi)*F_NN(:, :, jz))
    F_QN(:, :, jz) = max(real(1.e-24, kind=rprec), F_QN(:, :, jz))

    Cs_opt2_4d(:, :, jz) = F_QN(:, :, jz)/F_NN(:, :, jz)
!   Cs_opt2_4d(:,:,jz) = SUM(QN(:,:))/SUM(NN(:,:))
    Cs_opt2_4d(ldx, :, jz) = real(1.e-24, kind=rprec)
    Cs_opt2_4d(ldx - 1, :, jz) = real(1.e-24, kind=rprec)
    Cs_opt2_4d(:, :, jz) = max(real(1.e-24, kind=rprec), Cs_opt2_4d(:, :, jz))

    Beta(:, :, jz) = &
       (Cs_opt2_4d(:, :, jz)/Cs_opt2_2d(:, :, jz))**(log(tf1)/(log(tf2) - log(tf1)))

    ! counter1 = 0
    ! counter2 = 0

    ! do jx = 1, nxt
    ! do jy = 1, nynpy
    !     if (Beta(jx, jy, jz) .le. 1._rprec/(tf1*tf2)) then
    !         counter1 = counter1 + 1
    !     end if
    ! end do
    ! end do

    !--MPI: this is not valid
    if ( coordz == npz - 1 .and. jz == nz ) then
        Beta(:, :, jz) = 1._rprec
    end if

    Betaclip(:, :, jz) = max(Beta(:, :, jz), 1._rprec/(tf1*tf2)) 
    Cs_opt2(:, :, jz) = Cs_opt2_2d(:, :, jz)/Betaclip(:, :, jz)
    Cs_opt2(ldx, :, jz) = real(1.e-24, kind=rprec)
    Cs_opt2(ldx - 1, :, jz) = real(1.e-24, kind=rprec)
    Cs_opt2(:, :, jz) = max(real(1.e-24, kind=rprec), Cs_opt2(:, :, jz)) 

    
    call calc_hor_avg3(Cs_opt2(1:nxt, 1:nynpy, jz),    Cs_opt2_havg)
    call calc_hor_avg3(Cs_opt2_2d(1:nxt, 1:nynpy, jz), Cs_opt2_2d_havg)
    call calc_hor_avg3(Cs_opt2_4d(1:nxt, 1:nynpy, jz), Cs_opt2_4d_havg)

    Beta_avg(jz)     = Cs_opt2_2d_havg/Cs_opt2_havg
    Betaclip_avg(jz) = Cs_opt2_4d_havg/Cs_opt2_2d_havg
    
    ! Beta_avg(jz) = sum(Cs_opt2_2d(1:nxt, 1:nynpy, jz))/sum(Cs_opt2(1:nxt, 1:nynpy, jz))
    ! Betaclip_avg(jz) = sum(Cs_opt2_4d(1:nxt, 1:nynpy, jz))/sum(Cs_opt2_2d(1:nxt, 1:nynpy, jz))

    if (calc_lag_Ssim) then !calc_lag_Ssim is specified in param.f90
    ! This basicallpy saves the Cs calculated using the scale-invariance assumption
    ! within the scale-dependent model to Cs_Ssim
        Cs_Ssim = Cs_opt2_2d
    end if

    if (mod(jt, 200) == 0) then
        call calc_hor_avg3(LM(1:nxt, 1:nynpy), LMvert(jz))
        call calc_hor_avg3(MM(1:nxt, 1:nynpy), MMvert(jz))
        call calc_hor_avg3(QN(1:nxt, 1:nynpy), QNvert(jz))
        call calc_hor_avg3(NN(1:nxt, 1:nynpy), NNvert(jz))
        ! LMvert(jz) = sum(sum(LM, DIM=1), DIM=1)/nynpy/nxt
        ! MMvert(jz) = sum(sum(MM, DIM=1), DIM=1)/nynpy/nxt
        ! QNvert(jz) = sum(sum(QN, DIM=1), DIM=1)/nynpy/nxt
        ! NNvert(jz) = sum(sum(NN, DIM=1), DIM=1)/nynpy/nxt
    end if
end do
! ---
end subroutine lagrange_Sdep
