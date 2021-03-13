!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
module sgsmodule
!-----------------------------------------------------------------------
!   Modified to fit namelist feature, add allocatable attribution
!  (Bicheng Chen 06/12/2016)
!-----------------------------------------------------------------------
use types, only: rprec
implicit none
!TS In genreal, ofttime is 1.5 (Meneveau et al., 1996)
real(kind=rprec), parameter :: opftime = 1.5_rprec
!## scaledep_dynamic.f90 and lagrange_Sdep.f90 and lagrange_Ssim.f90
real(rprec), dimension(:, :, :), allocatable :: visc
real(rprec), dimension(:, :, :), allocatable :: Cs_opt2_2d, Cs_opt2_4d
real(rprec), dimension(:, :), allocatable :: S, tempos
real(rprec), dimension(:, :), allocatable :: u_pr, v_pr, w_pr
real(rprec), dimension(:, :), allocatable :: w_nod
real(rprec), dimension(:, :), allocatable :: L11, L12, L13, L22, L23, L33
real(rprec), dimension(:, :), allocatable, target :: Q11, Q12, Q13, Q22, Q23, Q33
real(rprec), dimension(:, :), pointer :: M11, M12, M13, M22, M23, M33
real(rprec), dimension(:, :), allocatable :: N11, N12, N13, N22, N23, N33
real(rprec), dimension(:, :), allocatable :: LM, MM, QN, NN, Tn, epsi, dumfac
real(rprec), dimension(:, :), allocatable :: S_bar, S11_bar, S12_bar, &
   S13_bar, S22_bar, S23_bar, S33_bar, S_S11_bar, S_S12_bar, &
   S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
real(rprec), dimension(:, :), allocatable :: S_hat, S11_hat, S12_hat, &
   S13_hat, S22_hat, S23_hat, S33_hat, S_S11_hat, S_S12_hat, &
   S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat
real(rprec), dimension(:, :), allocatable :: u_bar, v_bar, w_bar
real(rprec), dimension(:, :), allocatable :: u_hat, v_hat, w_hat
real(rprec), dimension(:, :), allocatable :: fourbeta
real(rprec), dimension(:), allocatable :: beta_sd
real(rprec), dimension(:), allocatable :: LMvert, MMvert, QNvert, NNvert
!## interpolag_Sdep.f90 and interpolag_Ssim.f90
real(rprec), dimension(:, :, :), allocatable :: xp, yp, zp, u_temp, v_temp
real(rprec), dimension(:, :, :), allocatable :: FF_LM, FF_MM, FF_QN, FF_NN, Beta_t
real(rprec), dimension(:, :, :), allocatable :: F_LM, F_MM, F_QN, F_NN, Beta, Betaclip
real(rprec), dimension(:), allocatable :: Beta_avg, Betaclip_avg

!xxxx----- Added by Vij - 04/14/04--xxxx----------------
! Nu_t is needed for scalar sgs
! dissip is dissipation calculated in sgs_stag and outputted when needed
! For more details look into scalars_module.f90
real(rprec), dimension(:, :, :), allocatable :: Nu_t, dissip
!xxxx---- Vij change ends here --xxxx-------------
real(rprec), dimension(:, :, :), allocatable :: u_lag, v_lag, w_lag
!real(kind=rprec),dimension(ld,nyt,$lbz:nz)::u_lag1,v_lag1,w_lag1
integer ::jt_count
real(rprec), dimension(:, :, :), allocatable :: Cs_opt2, Cs_opt2_avg
real(rprec), dimension(:, :, :), allocatable :: Cs_Ssim

! Added by chamecki to calculate dynamic Schmidt number
! This is used to store |S| in sgs_stag.f90
real(rprec), dimension(:, :, :), allocatable :: magS

! Added by chamecki to save the value of epsilon_2delta from the lagrangian
! calculations for the velocity field to be used for the scalar field
real(rprec), dimension(:, :, :), allocatable :: epsilon_lag, epsilon_lag2

! These are the two vector contractions for the dynamic lagrangian scalar SGS calculations
real(rprec), dimension(:, :, :), allocatable :: F_KX, F_XX
real(rprec), dimension(:, :, :), allocatable :: F_KX2, F_XX2

! These are the interpolated previous locations for the lagrangian averaging
! They are calculated in intrepola_Sdep.f90 and stored to be used in dynamic_sc.f90
real(rprec), dimension(:, :, :), allocatable :: xlag, ylag, zlag

contains

    real(kind=rprec) function rtnewt(A, jz)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    use types, only: rprec
    integer, parameter :: jmax = 100
    real(kind=rprec) :: x1, x2, xacc
    integer :: j, jz
    real(kind=rprec) :: df, dx, f
    real(kind=rprec), dimension(0:5) :: A
    x1 = 0._rprec
    x2 = 15._rprec ! try to find the largest root first....hmm
    xacc = 0.001_rprec ! doesn't need to be that accurate...
    rtnewt = 0.5_rprec*(x1 + x2)
    do j = 1, jmax
      f = A(0) + rtnewt*(A(1) + rtnewt*(A(2) + rtnewt*(A(3) + rtnewt*(A(4) + rtnewt*A(5)))))
      df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) + &
          rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))
      dx = f/df
      rtnewt = rtnewt - dx
    !        if ((x1-rtnewt)*(rtnewt-x2) < 0.) STOP 'rtnewt out of bounds'
      if (abs(dx) < xacc) return
    end do
    rtnewt = 1._rprec ! if dont converge fast enough
    write (6, *) 'using beta=1 at jz= ', jz
    end function rtnewt
! ---
end module sgsmodule
