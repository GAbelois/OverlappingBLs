!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module intp
!-----------------------------------------------------------------------
! Interpolation from coarse simulation
!-----------------------------------------------------------------------
use types, only:rprec
implicit none
save
public
! ---
real(rprec), dimension(:, :, :), allocatable :: u_tmp, v_tmp, w_tmp
real(rprec), dimension(:, :, :), allocatable :: RHSx_tmp, RHSy_tmp, RHSz_tmp
real(rprec), dimension(:, :, :), allocatable :: Cs_opt2_tmp
real(rprec), dimension(:, :, :), allocatable :: F_LM_tmp, F_MM_tmp, F_QN_tmp, F_NN_tmp
real(rprec), dimension(:, :, :), allocatable :: theta_tmp, RHS_T_tmp
real(rprec), dimension(:, :), allocatable :: sgs_t3_tmp, psi_m_tmp
! ---
integer(kind=8) :: plan_xf_intp, plan_yf_intp, plan_x2b_intp, plan_y2b_intp
integer, dimension(3) :: gsize2, lsize2, start2   ! dimensions of global and local variables 
! ---
end module 
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 