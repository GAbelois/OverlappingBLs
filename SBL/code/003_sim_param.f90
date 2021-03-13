!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
module sim_param
!--------------------------------------------------------------------!
!   purpose:                                          
!--------------------------------------------------------------------!  
use types, only:rprec
implicit none
save
public
! --- Flow Field 
real(rprec), dimension(:, :, :), allocatable :: u, v, w
real(rprec), dimension(:, :, :), allocatable :: ke, ke_temp, dragx, dragy, dragz
real(rprec), dimension(:, :, :), allocatable :: dudx, dudy, dudz,       &
                                                dvdx, dvdy, dvdz,       &
                                                dwdx, dwdy, dwdz,       &
                                                RHSx, RHSy, RHSz,       &
                                                RHSx_f, RHSy_f, RHSz_f
real(rprec), dimension(:, :, :), allocatable :: dudt, dvdt, dwdt
real(rprec), dimension(:, :, :), allocatable :: dkedx, dkedy, dkedz
real(rprec), dimension(:, :, :), allocatable :: txx, txy, txz, tyy, tyz, tzz
real(rprec), dimension(:, :, :), allocatable :: p
real(rprec), dimension(:, :, :), allocatable :: dpdx, dpdy, dpdz
real(rprec), dimension(:, :, :), allocatable :: theta, q
real(rprec), dimension(:, :, :), allocatable :: divtx, divty, divtz

! --- Concentration Field 
! --- Added for pollen - pollen concentration in grains/m3 Chamecki - 08/01/2006
real(rprec), dimension(:, :, :, :), allocatable :: pcon
! BC Added by Bicheng Chen for Chemical reaction with background species 05/21/2015
real(rprec), dimension(:, :, :, :), allocatable :: k_ttlCR ! total chemical reaction rate for concentraiton field k_ttlCR = sum(k_cR(i) * con_bg(i))
! ---
end module sim_param
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
