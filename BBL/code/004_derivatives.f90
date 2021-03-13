!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
module derivatives
!-----------------------------------------------------------------------
! Module to calculate the spatial derivatives
! 01/16/2019 - First created by Bicheng Chen (chabby@ucla.edu)
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: nxt, nyt, nz, lhx, lhy, ldx, nxnpy, nynpy
implicit none
save
! ---
contains
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine ddx(dfdx, fun)
    !!! Compute the spatial derivative in x-direction using FFT with
    !!! nxt data points
    use fft, only: plan_xf, plan_xb, eye, kx
    implicit none
    ! ---
    integer :: iy, iz
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(in) :: fun
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(inout) :: dfdx
    
    complex(rprec), dimension(lhx) :: arr
    real(rprec) :: const

    const = 1._rprec / nxt
    ! Loop through the horizontal slice
    do iz = 0, nz
    do iy = 1, nynpy
        call dfftw_execute_dft_r2c(plan_xf, const*fun(1:nxt, iy, iz), arr)
        arr = eye * kx * arr
        call dfftw_execute_dft_c2r(plan_xb, arr, dfdx(1:nxt, iy, iz))
    end do
    end do
    ! ---
    end subroutine ddx
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine ddy(dfdy, fun)
    !!! Compute the spatial derivative in y-direction using FFT with
    !!! nyt data points
    use fft, only: plan_yf, plan_yb, eye, ky
    implicit none
    ! ---
    integer :: ix, iz
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(in) :: fun
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(inout) :: dfdy
    ! ---
    real(rprec), dimension(nxnpy, nyt) :: tmp
    complex(rprec), dimension(lhy) :: arr
    real(rprec), dimension(nxnpy, nyt) :: dfdy_tmp 
    real(rprec) :: const

    const = 1._rprec / nyt
    ! ---
    ! do iz = 0, nz
    ! do ix = 1, nxt
    !     call collect_data_y(const*fun(ix, 1:nynpy, iz), tmp)
    !     call dfftw_execute_dft_r2c(plan_yf, tmp, arr)
    !     arr = eye * ky * arr
    !     call dfftw_execute_dft_c2r(plan_yb, arr, tmp)
    !     call scatter_y(tmp, dfdy(ix, 1:nynpy, iz))
    ! end do
    ! end do
    do iz = 0, nz
        call transpose_phy_y_to_x(fun(1:nxt, 1:nynpy, iz), tmp)
        do ix = 1, nxnpy
            call dfftw_execute_dft_r2c(plan_yf, const*tmp(ix,1:nyt), arr)
            arr = eye * ky * arr
            call dfftw_execute_dft_c2r(plan_yb, arr, dfdy_tmp(ix, 1:nyt))
        end do
        call transpose_phy_x_to_y(dfdy_tmp(1:nxnpy, 1:nyt), dfdy(1:nxt, 1:nynpy, iz))
    end do
    ! ---
    end subroutine ddy
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine ddz_uv(dfdz, fun)
    !!! Compute the spatial derivative in z-direction at uv-node
    !!! using finite difference method
    !!! 01/17/2019 - Modified to use array operation
    !!! >>by Bicheng Chen (chabby@ucla.edu)
    use param, only: dz, coordz, npz, BOGUS
    implicit none
    ! ---
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(in) :: fun
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(inout) :: dfdz

    real(rprec) :: const
    ! ---
    const = 1._rprec/dz

    dfdz(:, :, 0) = BOGUS
    if ( coordz > 0 ) then
        !--ghost node information is available here
        !--watch the z-dimensioning!
        !--if coordz == 0, dudz(1) will be set in wallstress
        dfdz(1:nxt, 1:nynpy, 1) = const*(fun(1:nxt, 1:nynpy, 1) - fun(1:nxt, 1:nynpy, 0))
    end if

    dfdz(1:nxt, 1:nynpy, 2:nz-1) = const*(fun(1:nxt, 1:nynpy, 2:nz-1) - fun(1:nxt, 1:nynpy, 1:nz-2))

    if ( coordz < npz - 1 ) then
        !--if coordz == npz-1, dudz(nz) will be set in wallstress
        dfdz(1:nxt, 1:nynpy, nz) = const*(fun(1:nxt, 1:nynpy, nz) - fun(1:nxt, 1:nynpy, nz-1))
    end if

    ! ---
    end subroutine ddz_uv
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine ddz_uv_con(dfdz, fun)
    !!! Compute the spatial derivative in z-direction at uv-node
    !!! using finite difference method
    !!! 01/17/2019 - Modified to use array operation
    !!! >>by Bicheng Chen (chabby@ucla.edu)
    use param, only: dz, coordz, npz, BOGUS
    implicit none
    ! ---
    real(rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(in) :: fun
    real(rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(inout) :: dfdz

    real(rprec) :: const
    ! ---
    const = 1._rprec/dz

    dfdz(:, :, 0) = BOGUS
    if ( coordz > 0 ) then
        !--ghost node information is available here
        !--watch the z-dimensioning!
        !--if coordz == 0, dudz(1) will be set in wallstress
        dfdz(1:nxt, 0:nynpy+1, 1) = const*(fun(1:nxt, 0:nynpy+1, 1) - fun(1:nxt, 0:nynpy+1, 0))
    end if

    dfdz(1:nxt, 0:nynpy+1, 2:nz-1) = const*(fun(1:nxt, 0:nynpy+1, 2:nz-1) - fun(1:nxt, 0:nynpy+1, 1:nz-2))

    if ( coordz < npz - 1 ) then
        !--if coordz == npz-1, dudz(nz) will be set in wallstress
        dfdz(1:nxt, 0:nynpy+1, nz) = const*(fun(1:nxt, 0:nynpy+1, nz) - fun(1:nxt, 0:nynpy+1, nz-1))
    end if

    ! ---
    end subroutine ddz_uv_con
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine ddz_w(dfdz, fun)
    !!! Compute the spatial derivative in z-direction at w-node using FDM
    !!! 01/17/2019 - Moved from ddz_uv.f90 by Bicheng Chen (chabby@ucla.edu)
    !!! 01/17/2019 - Modified to use array operation
    !!! >>by Bicheng Chen (chabby@ucla.edu)
    use param, only: dz, coordz, npz, BOGUS
    implicit none
    ! ---
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(in) :: fun
    real(rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: dfdz
    ! ---
    real(rprec) :: const

    const = 1._rprec/dz
    dfdz(1:nxt, 1:nynpy, 0:nz-1) = const*(fun(1:nxt, 1:nynpy, 1:nz) - fun(1:nxt, 1:nynpy, 0:nz-1))

    if ( coordz == 0 ) then
        !--bottom process cannot calculate dfdz(jz=0)
        dfdz(:, :, 0) = BOGUS
    end if
    
    if ( coordz == npz - 1 ) then
        dfdz(:, :, nz) = 0._rprec
    else
        dfdz(:, :, nz) = BOGUS
    end if

    ! ---
    end subroutine ddz_w
! ---
end module derivatives
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------