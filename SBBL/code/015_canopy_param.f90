!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module canopy
!-----------------------------------------------------------------------
! canopy-related parameters
!-----------------------------------------------------------------------
use types, only:rprec
implicit none
save
public
! --- nutrient uptake related parameters
real(rprec), parameter :: Vmax = 752._rprec/3600._rprec  ! maximum uptake rate, units: umoles m-2 s-1
real(rprec), parameter :: Km   = 10.2_rprec*1000._rprec  ! half-saturation constant, units: umoles/m-3, unit conversion from liter to m^3
! ---
real(kind=rprec), dimension(:,:), allocatable :: nutrient_inflow
! ---
end module 
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 