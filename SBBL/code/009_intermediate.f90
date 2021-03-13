!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
module intermediate
!-----------------------------------------------------------------------
!   save the intermediate variables
!-----------------------------------------------------------------------
use types, only: rprec
!use param, only: nxt, nyt, nz

implicit none
save
public
!### Output Setup ###
 !## scalar_modules2.f90, io.f90
 ! The variable average_dim_select_flag generates the following values based
 ! >>on the value of average_dim_num in param.f90 :-
 ! >>a) average_dim_num = 2 :
 ! >>average_dim_select_flag = 0 ==> average the 3d array over x and y and
 ! >>output the z profile
 ! >>b) average_dim_num = 1 :
 ! >>average_dim_select_flag = 1 ==> average the 3d array over y and
 ! >>output the (x,z) profile
! --- scalar_module2.f90.f90, subroutine flow_slice
real(rprec), dimension(:, :), allocatable :: avg_u, avg_v, avg_w, avg_p
real(rprec), dimension(:, :), allocatable :: avg_dudt, avg_dvdt, avg_dwdt
real(rprec), dimension(:, :), allocatable :: avg_dudz, avg_dvdz

real(rprec), dimension(:, :), allocatable :: avg_u2, avg_v2, avg_w2, avg_p2
real(rprec), dimension(:, :), allocatable :: avg_uv, avg_uw, avg_vw
real(rprec), dimension(:, :), allocatable :: avg_up, avg_vp, avg_wp
real(rprec), dimension(:, :), allocatable :: avg_pdudx, avg_pdvdy, avg_pdwdz
real(rprec), dimension(:, :), allocatable :: avg_txx, avg_txy, avg_txz, avg_tyy, avg_tyz, avg_tzz

real(rprec), dimension(:, :), allocatable :: avg_u3, avg_v3, avg_w3
real(rprec), dimension(:, :), allocatable :: avg_u2v, avg_u2w, avg_uv2, avg_v2w, avg_uw2, avg_vw2, avg_uvw

real(rprec), dimension(:, :), allocatable :: avg_utxx, avg_utxy, avg_utxz
real(rprec), dimension(:, :), allocatable :: avg_vtxy, avg_vtyy, avg_vtyz
real(rprec), dimension(:, :), allocatable :: avg_wtxz, avg_wtyz, avg_wtzz, avg_dissip

real(rprec), dimension(:, :), allocatable :: avg_utxz_fluc, avg_vtyz_fluc

real(rprec), dimension(:, :), allocatable :: avg_Cs
! real(rprec), dimension(:, :), allocatable :: aCs_Ssim, abeta_sgs, abetaclip_sgs 

real(rprec), dimension(:, :), allocatable :: avg_u2_fluc, avg_v2_fluc, avg_T2_fluc
real(rprec), dimension(:, :), allocatable :: avg_uw_fluc, avg_vw_fluc, avg_uv_fluc 
real(rprec), dimension(:, :), allocatable :: avg_wp_fluc, avg_wT_fluc
real(rprec), dimension(:, :), allocatable :: avg_u2w_fluc, avg_v2w_fluc

! --- scalar_module2.f90, subroutine theta_slice and pcon_slice
real(rprec), dimension(:, :), allocatable :: avg_theta, avg_theta2, avg_dTdz
real(rprec), dimension(:, :), allocatable :: avg_sgsTheta, avg_wTheta, avg_Nut

real(rprec), dimension(:, :, :), allocatable :: avg_PCon, avg_PCon2, avg_dPCondz
real(rprec), dimension(:, :, :), allocatable :: avg_sgsPCon1, avg_sgsPCon3, avg_uPCon, avg_wPCon
real(rprec), dimension(:, :), allocatable :: avg_Kct, avg_Cs2Sc
! --- scalar_module2.f90, subroutine flow_slice


!### Flow Parameters ###
! --- convec.f90
real(rprec), dimension(:, :, :), allocatable :: cc_big
real(rprec), dimension(:, :, :), allocatable :: u_big, v_big, w_big
real(rprec), dimension(:, :, :), allocatable :: cross_x_big, cross_y_big, cross_z_big

real(rprec), dimension(:, :, :), allocatable :: cross_x, cross_y, cross_z
real(rprec), dimension(:, :, :), allocatable :: vort1, vort2, vort3
real(rprec), dimension(:, :, :), allocatable :: vort1_big, vort2_big, vort3_big
! --- scalar_module.f90
real(rprec), dimension(:, :, :), allocatable :: u_cell, v_cell
real(rprec), dimension(:, :, :, :), allocatable :: pu0_cell, pv0_cell
real(rprec), dimension(:, :, :, :), allocatable :: pu_cell, pv_cell
real(rprec), dimension(:, :), allocatable :: decay_cell
! --- press_stag_array.f90
real(rprec), dimension(:, :, :), allocatable, target :: rH_x, rH_y, rH_z
complex(rprec), dimension(:, :, :), allocatable :: H_x, H_y, H_z
real(rprec), dimension(:, :), allocatable :: rtopw, rbottomw
complex(rprec), dimension(:, :), allocatable :: topw, bottomw
 
!### Temperature Parameters ###
 !## dynamic_sc.f90, subroutine interpolag_scalar and
 !## >>interpolag_scalar_Sdep
real(rprec), dimension(:, :, :), allocatable :: FF_KX, FF_XX, FF_KX2, FF_XX2

real(rprec) :: fr
! ---
end module intermediate
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
